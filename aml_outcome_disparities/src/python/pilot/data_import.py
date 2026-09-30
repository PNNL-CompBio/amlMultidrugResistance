import gzip
import logging
from io import StringIO
from os import PathLike
from os.path import abspath, dirname, exists, join
from typing import Iterable

import mygene
import numpy as np
import pandas as pd
import requests
import rpy2.robjects as ro
import synapseclient as sc
from inmoose.pycombat import pycombat_seq
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale, StandardScaler
from synapseclient.entity import File
from synapseutils import syncFromSynapse

REPO_PATH = abspath(dirname(dirname(__file__)))
logging.getLogger("rpy2.rinterface").setLevel(logging.ERROR)


def syn_login(auth_path: PathLike | None = None) -> sc.Synapse:
    """
    Login to Synapse.

    Args:
        auth_path (PathLike | None, default: None): Path to authentication file.

    Returns:
        sc.Synapse: Logged-in Synapse client.
    """
    if auth_path is None:
        auth_path = join(REPO_PATH, "auth_token.txt")

    syn = sc.Synapse(debug=False, silent=True)
    with open(auth_path, "r") as f:
        auth_token = f.read()

    syn.login(authToken=auth_token)

    return syn


def call_mutations(
    syn: sc.Synapse | None = None
) -> pd.DataFrame:
    """
    Calls mutations from WES measurements.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
        new one if None.

    Returns:
        pd.DataFrame: Mutation calls for BeatAML 210 cohort.
    """
    if syn is None:
        syn = syn_login()

    # Load mutation data
    mutation_calls = pd.read_excel(syn.get("syn26427388").path)

    # # Load additional mutation calls
    calls_177 = pd.read_csv(
        syn.get("syn32533104").path,
        sep="\t"
    )
    calls_177.rename(
        columns={"labId": "original_id", "alt_ID": "seq_id"},
        inplace=True
    )
    mutation_calls = pd.concat([mutation_calls, calls_177])

    # Drop duplicate calls due to multiple calling methods
    mutation_calls = mutation_calls.loc[
        ~mutation_calls.loc[
            :,
            ["original_id", "hgvsp"]
        ].duplicated(keep="first")
    ]

    # Track all patients with calls
    # Those without any recurrent may be interesting to keep--as a WES, a lack
    # of calls should indicate a lack of mutations
    all_patients = pd.Index(mutation_calls.loc[:, "original_id"].unique())

    # Trim to nonsynonymous mutations
    mutation_calls = mutation_calls.loc[
        ~mutation_calls.loc[:, "hgvsp"].isna(),
        :
    ]

    # Trim to recurrent mutations
    recurrent = mutation_calls.loc[
        :,
        ["pos_start", "pos_end", "alt", "original_id", "symbol"]
    ]
    recurrent = recurrent.loc[~recurrent.duplicated(keep="first")]
    recurrent = recurrent.loc[
        recurrent.loc[:, ["pos_start", "pos_end", "alt"]].duplicated(
            keep=False
        ),
        :
    ]

    # Drop duplicate patient/gene pairs
    # Some patients have a gene with multiple mutations
    recurrent = recurrent.loc[
        ~recurrent.loc[
            :,
            ["original_id", "symbol"]
        ].duplicated(keep="first"),
        :
    ]

    # Pivot to wide format
    recurrent.loc[:, "value"] = "Mutant"
    recurrent = recurrent.pivot(
        index="original_id",
        columns="symbol",
        values="value"
    )

    # Adds patients in WES without recurrent mutations
    recurrent = pd.concat(
        [
            recurrent,
            pd.DataFrame(
                "WT",
                index=all_patients.difference(recurrent.index),
                columns=recurrent.columns
            )
        ]
    )

    return recurrent


def import_metabolites(
    syn: sc.Synapse | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads metabolite data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
        new one if None.

    Returns:
        pd.DataFrame: Metabolite data for BeatAML 210 cohort.
        pd.DataFrame: Metabolite data for Pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    # Import only RP
    # Seems that HILIC has some batch challenges not correctable via ComBat
    data = pd.read_csv(syn.get("syn71896311").path, index_col=0)

    # # Import and concatenate metabolite datasets
    # rp = pd.read_csv(
    #     syn.get("syn71896311").path,
    #     index_col=0
    # )
    # hilic = pd.read_csv(
    #     syn.get("syn71896297").path,
    #     index_col=0
    # )
    # data = pd.concat([rp, hilic])

    # Import Sample ID to universal Accession numbers
    ptrc_conversions = pd.read_excel(syn.get("syn25796769").path, index_col=0)
    pilot_conversions = pd.read_excel(syn.get("syn68835814").path, index_col=0)

    # Isolate Sample ID
    ptrc = data.loc[:, data.columns.str.startswith("BEAT_AML_PNL")]
    pilot = data.loc[:, data.columns.str.startswith("PTRC_exp26_Met")]
    ptrc.columns = ptrc.columns.str[13:16].astype(int)
    pilot.columns = pilot.columns.str[15:].astype(int)

    # Convert ID numbers to Accession
    ptrc_conversions = ptrc_conversions.loc[:, "labId"]
    ptrc_conversions.index = ptrc_conversions.index.str[11:].astype(int)
    pilot_conversions.set_index("Metabolomics ID", inplace=True)
    pilot_conversions = pilot_conversions.loc[:, "Accession/\nBarcode"]
    pilot_conversions.index = pilot_conversions.index.str[15:].astype(int)

    # Trims additional information from patient IDs
    pilot_conversions.loc[pilot_conversions.str.startswith("P")] = (
        pilot_conversions.loc[pilot_conversions.str.startswith("P")].str[:9]
    )

    ptrc.rename(columns=ptrc_conversions, inplace=True)
    pilot.rename(columns=pilot_conversions, inplace=True)

    # Drop samples without Accession numbers, relabel Bridging samples
    ptrc = ptrc.loc[:, ptrc.columns.intersection(ptrc_conversions)]
    pilot = pilot.loc[
        :,
        pilot.columns.intersection(pilot_conversions),
    ]
    pilot.loc[
        :, pilot.columns.intersection(ptrc_conversions.values)
    ].columns += "-Bridge"

    # Drop patients in both datasets--default to BeatAML
    pilot = pilot.drop(pilot.columns.intersection(ptrc.columns), axis=1)

    return ptrc.T, pilot.T


def import_lipids(
    syn: sc.Synapse | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads lipid data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.

    Returns:
        pd.DataFrame: Lipid data for BeatAML 210 cohort.
        pd.DataFrame: Lipid data for Pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    # Load lipid meta-data
    lipid_meta = pd.read_csv(
        syn.get("syn71896673").path
    )
    lipid_conversions = pd.Series(
        lipid_meta.loc[:, "Accession"].values,
        index=lipid_meta.loc[:, "SampleID"],
    )

    # Pull data, convert sample to cross-experiment Accession IDs
    data = pd.read_csv(syn.get("syn71896667").path, index_col=0)
    data.rename(columns=lipid_conversions, inplace=True)
    data = data.loc[:, ~data.columns.isna()]
    data = data.loc[:, ~data.columns.duplicated()]

    # Split out BeatAMl from Pilot cohorts, for consistency
    lipid_meta.set_index("Accession", inplace=True)
    lipid_meta = lipid_meta.loc[data.columns, :]
    lipid_meta = lipid_meta.loc[~lipid_meta.index.duplicated(), :]
    ptrc = data.loc[:, lipid_meta.loc[:, "study"] == "beataml"]
    pilot = data.loc[:, lipid_meta.loc[:, "study"] == "pilot"]

    return ptrc.T, pilot.T


def import_rna(
    syn: sc.Synapse | None = None,
    return_symbols: bool = True,
    batch_correct: bool = True,
    tpm: bool = True,
    normalize: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads and TPMs RNA data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        return_symbols (bool, default: True): Converts Ensembl numbers to
            UniProt symbols.
        batch_correct (bool, default: True): Corrects for batches via ComBat-seq
            through inmoose.
        tpm (bool, default: True): Converts counts to transcripts-per-million.
        normalize (bool, default: False): Z-scores genes, separating BM from
            PB samples.

    Returns:
         pd.DataFrame: RNA data for BeatAML 210 cohort.
         pd.DataFrame: RNA data for Pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    # Check if BioMart data is stored and, if not, downloads it
    if not exists(
        join(
            REPO_PATH,
            "data",
            "biomart_gene_lengths.csv.gz"
        )
    ):
        print("BioMart gene lengths not found. Downloading from BioMart...")
        # Load XML request
        with open(
            join(
                REPO_PATH,
                "data",
                "biomart_request.xml"
            ),
            "r"
        ) as xml_file:
            xml_req = xml_file.read()

        # Download from BioMart
        # This can take a few minutes!
        r = requests.get(
            "http://www.ensembl.org/biomart/martservice?query=" + xml_req
        )
        if not r.ok:
            print(f"BioMart request failed with status code {r.status_code}")

        # Save as csv.gz file
        with gzip.open(
            join(
                REPO_PATH,
                "data",
                "biomart_gene_lengths.csv.gz"
            ),
            "wb"
        ) as gz_file:
            gz_file.write(r.content)

        print("Download complete. Saving to python/data/")

    # Load gene lengths, symbols
    gene_lengths = pd.read_csv(
        join(REPO_PATH, "data", "biomart_gene_lengths.csv.gz")
    )
    gene_lengths.set_index("Gene stable ID", inplace=True, drop=True)
    gene_lengths.sort_values("UniProtKB Gene Name symbol", inplace=True)
    gene_lengths = gene_lengths.loc[
        ~gene_lengths.index.duplicated(keep="first")
    ]

    if batch_correct:
        # Setup Synapse file objects
        ptrc_file = File(
            join(REPO_PATH, "data", "RNAseq_beataml_batch_corrected.csv"),
            parentId="syn68820222",
            id="syn72665045"
        )
        pilot_file = File(
            join(REPO_PATH, "data", "RNAseq_pilot_batch_corrected.csv"),
            parentId="syn68820222",
            id="syn72665034"
        )
        try:
            # Syncs with Synapse files
            syncFromSynapse(
                syn,
                ptrc_file,
                path=join(REPO_PATH, "data")
            )
            syncFromSynapse(
                syn,
                pilot_file,
                path=join(REPO_PATH, "data")
            )

            # Loads cached results (if this has been run previously) as
            # ComBat-Seq can be time-intensive with race covariates included
            ptrc = pd.read_csv(ptrc_file.path, index_col=0)
            pilot = pd.read_csv(pilot_file.path, index_col=0)

        except FileNotFoundError:
            # Load in meta-data for larger cohort, including patients outside
            # pilot study
            meta = import_meta(syn, aux_meta=True)
            meta = meta.loc[
                ~meta.loc[:, "Race"].isin(["Declined", "Unknown"]), :
            ]

            # Load raw RNA data
            ptrc, pilot = import_rna(
                syn,
                return_symbols=False,
                batch_correct=False,
                tpm=False,
                normalize=False
            )

            # Only keep patients with meta-data
            ptrc = ptrc.loc[:, ptrc.columns.intersection(meta.index)]
            pilot = pilot.loc[:, pilot.columns.intersection(meta.index)]
            data = pd.concat([ptrc, pilot], axis=1)
            meta = meta.loc[data.columns, :]

            # Batch correct, including race as a covariate to preserve
            data = pycombat_seq(
                data,
                meta.loc[:, "Study"].values,
                covar_mod=(meta.loc[:, "Race"] == "Black").astype(int).values,
            )
            ptrc = data.loc[:, ptrc.columns]
            pilot = data.loc[:, pilot.columns]

            # Store locally
            ptrc.to_csv(ptrc_file.path)
            pilot.to_csv(pilot_file.path)

            # Sync to Synapse
            syn.store(ptrc_file)
            syn.store(pilot_file)
    else:
        # Load data from Synapse
        ptrc = pd.read_csv(syn.get("syn64126462").path, index_col=0, sep="\t")
        ptrc = ptrc.iloc[:, 3:].astype(float)
        pilot = pd.read_csv(syn.get("syn68820229").path, index_col=0, sep="\t")

        # Convert experiment IDs to universal Accession IDs
        pilot_conversion = pd.read_csv(syn.get("syn68820228").path)
        pilot_conversion = pd.Series(
            pilot_conversion.loc[:, "Accession"].values,
            index=pilot_conversion.loc[:, "RNA"],
        )
        pilot_conversion = pilot_conversion.loc[~pilot_conversion.index.isna()]
        pilot.rename(columns=pilot_conversion, inplace=True)

        ptrc_conversion = pd.read_excel(syn.get("syn64126463").path, index_col=0)
        ptrc_conversion.dropna(subset="dbgap_rnaseq_sample", inplace=True)
        ptrc_conversion = pd.Series(
            ptrc_conversion.loc[:, "labId"].values,
            index=ptrc_conversion.loc[:, "dbgap_rnaseq_sample"].values,
        )
        ptrc.rename(columns=ptrc_conversion, inplace=True)

    # Removes low coverage genes
    ptrc = ptrc.loc[ptrc.mean(axis=1) > 1, :]
    pilot = pilot.loc[pilot.mean(axis=1) > 1, :]

    # Trims to genes in both datasets
    shared_genes = ptrc.index.intersection(pilot.index)
    shared_genes = shared_genes.intersection(gene_lengths.index)
    pilot = pilot.loc[shared_genes, :]
    ptrc = ptrc.loc[shared_genes, :]
    gene_lengths = gene_lengths.loc[shared_genes]

    # Converts to transcripts-per-million (TPM)
    if tpm:
        ptrc = ptrc.T
        pilot = pilot.T

        ptrc = ptrc / (
            gene_lengths.loc[:, "Transcript length (including UTRs and CDS)"]
            / 1000
        )
        sums = ptrc.sum(axis=1) / 1e6
        ptrc = ptrc.T / sums

        pilot = pilot / (
            gene_lengths.loc[:, "Transcript length (including UTRs and CDS)"]
            / 1000
        )
        sums = pilot.sum(axis=1) / 1e6
        pilot = pilot.T / sums

    # Translate Ensembl IDs to gene symbols
    if return_symbols:
        gene_lengths = gene_lengths.loc[
            ~gene_lengths.loc[:, "UniProtKB Gene Name symbol"].isna(),
            :
        ]
        ptrc.rename(
            index=gene_lengths.loc[:, "UniProtKB Gene Name symbol"],
            inplace=True
        )
        pilot.rename(
            index=gene_lengths.loc[:, "UniProtKB Gene Name symbol"],
            inplace=True
        )
        ptrc = ptrc.loc[~ptrc.index.duplicated(), :]
        pilot = pilot.loc[~pilot.index.duplicated(), :]

    # Normalize genes across samples (z-score)
    if normalize:
        # Import meta-data to isolate BM and PB samples
        sample_types = import_meta(syn, aux_meta=True)
        sample_types = sample_types.loc[:, "Specimen_Type"]
        sample_types.dropna(inplace=True)

        # Concatenate datasets, trim to those with known sample types
        combined = pd.concat([ptrc, pilot], axis=1)
        sample_types = sample_types.loc[
            sample_types.index.intersection(combined.columns)
        ]
        combined = combined.loc[:, sample_types.index]

        # Z-score across samples, then separately for BM and PB samples
        combined.loc[:] = scale(combined, axis=0)
        combined.loc[:, sample_types == "BM"] = scale(
            combined.loc[:, sample_types == "BM"], axis=1
        )
        combined.loc[:, sample_types == "PB"] = scale(
            combined.loc[:, sample_types == "PB"], axis=1
        )

        pilot = combined.loc[:, pilot.columns]
        ptrc = combined.drop(columns=pilot.columns)

    return ptrc.T, pilot.T


def plex_correct_proteomics(
    ba_data: pd.DataFrame,
    pilot_data: pd.DataFrame,
    syn: sc.Synapse | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Batch-corrects proteomic data using PNNL's NA-compatible ComBat.

    Args:
        ba_data (pd.DataFrame): BeatAML proteomic data.
        pilot_data (pd.DataFrame): Pilot proteomic data.
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.

    Returns:
        pd.DataFrame: Batch-corrected BeatAML proteomic data.
        pd.DataFrame: Batch-corrected Pilot proteomic data.
    """
    # Z-score measurements, use StandardScaler to undo later
    scaler = StandardScaler()
    data = pd.concat([ba_data, pilot_data], axis=0)
    data.loc[:] = scaler.fit_transform(data)
    ba_data = data.loc[ba_data.index, :]
    pilot_data = data.loc[pilot_data.index, :]

    # Get plex IDs for each sample and convert to accession IDs for
    # harmonization with phosphoproteomic sample IDs
    ba_data_conversion, pilot_data_conversion = import_sample_conversion(syn)
    ba_data_plex = pd.read_csv(
        syn.get("syn26534982").path, index_col=0, sep="\t"
    )
    ba_data_plex.index = ba_data_plex.index.str[11:].astype(int)
    ba_data_plex.rename(index=ba_data_conversion, inplace=True)
    ba_data_plex = "B" + ba_data_plex.loc[:, "Plex"].astype(str)

    # Repeat plex indexing for pilot_data samples
    pilot_data_plex = pd.read_excel(
        syn.get("syn68835814").path, sheet_name="TMT"
    ).set_index("Accession")
    pilot_data_plex = "P" + pilot_data_plex.loc[:, "Plex"].astype(str)
    pilot_data_plex = pilot_data_plex.loc[
        ~pilot_data_plex.index.isin(["EMPTY", "PooledReference"])
    ]

    # Relabel bridging samples
    bridging = pilot_data_plex.loc[
        pilot_data_plex.index.intersection(ba_data_plex.index)
    ]
    bridging.loc[:] = bridging.index + "-Bridge"
    pilot_data_plex.rename(index=bridging, inplace=True)

    # Concatenate plex into one array
    plex = pd.concat([ba_data_plex, pilot_data_plex], axis=0)

    # Source R batch correction script
    r_source = ro.r["source"]
    r_source(
        join(
            REPO_PATH,
            "..",
            "r",
            "proteomics",
            "proteomic_batch_correction.R",
        )
    )
    bc_function = ro.globalenv["combat_proteomics"]

    # Convert DataFrames to R
    with localconverter(pandas2ri.converter):
        ba_data_r = ro.conversion.py2rpy(ba_data)
        pilot_data_r = ro.conversion.py2rpy(pilot_data)
        meta_r = ro.conversion.py2rpy(import_meta(syn, aux_meta=False))
        plex_r = ro.conversion.py2rpy(plex)

    result = bc_function(ba_data_r, pilot_data_r, meta_r, plex_r)
    with localconverter(pandas2ri.converter):
        ba_data, pilot_data = ro.conversion.rpy2py(result)

    # Rename R DataFrame columns
    # R appends "X" to the front of numeric columns and replaces "-"
    # with "."
    pilot_data.columns = pilot_data.columns.str.replace({".": "-", "X": ""})
    ba_data.columns = ba_data.columns.str.replace({".": "-", "X": ""})

    # Undo Z-score back to measurements
    ba_data = ba_data.T
    pilot_data = pilot_data.T
    ba_data.loc[:] = scaler.inverse_transform(ba_data)
    pilot_data.loc[:] = scaler.inverse_transform(pilot_data)
    ba_data = ba_data.T
    pilot_data = pilot_data.T

    return ba_data, pilot_data


def correct_ptm(
    ptrc: pd.DataFrame,
    pilot: pd.DataFrame,
    syn: sc.Synapse | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Corrects PTM via global measurements.

    Args:
        ptrc (pd.DataFrame): BeatAML dataset.
        pilot (pd.DataFrame): Pilot dataset.
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.

    Returns:
        pd.DataFrame: Corrected PTM data for BeatAML cohort.
        pd.DataFrame: Corrected PTM data for Pilot cohort.
    """
    ptrc, pilot = ptrc.T, pilot.T

    # Merge datasets, import global
    prot = pd.concat(import_global(syn)).T
    ptrc = ptrc.loc[:, ptrc.columns.intersection(prot.columns)]
    pilot = pilot.loc[:, pilot.columns.intersection(prot.columns)]
    merged = pd.concat([ptrc, pilot], axis=1)

    # Split parent genes from sites, trim to global measurement
    merged.index = merged.index.str.split("-", expand=True)
    merged = merged.loc[
        merged.index.get_level_values(0).isin(prot.index),
        merged.columns.intersection(prot.columns)
    ]
    prot = prot.loc[
        :,
        prot.columns
    ]

    # Initialize OLS
    ols = LinearRegression()

    # Iterate through genes, correct
    for p_site in merged.index:
        # Correct measurements where both phospho and global are present
        correct_index = np.logical_and(
            merged.loc[p_site, :].notna(),
            prot.loc[p_site[0], :].notna()
        )

        ols.fit(
            prot.loc[p_site[0], correct_index].values.reshape(-1, 1),
            merged.loc[p_site, correct_index].values.reshape(-1, 1)
        )
        merged.loc[p_site, correct_index] = merged.loc[
                                                p_site,
                                                correct_index
                                            ] - ols.predict(
            prot.loc[p_site[0], correct_index].values.reshape(-1, 1)
        ).flatten()
        merged.loc[p_site, ~correct_index] = np.nan

    # Restore index format and split out PTRC, Pilot
    merged.index = ["-".join(row) for row in merged.index]

    ptrc = merged.loc[:, ptrc.columns]
    pilot = merged.loc[:, pilot.columns]

    return ptrc, pilot


def import_phospho(
    syn: sc.Synapse | None = None,
    batch_correct: str | None = "study",
    global_correct: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads phosphoproteomic data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        batch_correct (str | None, default: "cohort"): Batch-correction method.
            Accepts 'cohort' (corrected at cohort-level), 'study' (corrected at
            study-level, post-cohort concatenation), or None.
        global_correct (bool, default=True): Correct measurements using global
            proteomic abundances.

    Returns:
        pd.DataFrame: Phosphoproteomic data for BeatAML 210 cohort.
        pd.DataFrame: Phosphoproteomic data for Pilot cohort.
    """
    if batch_correct not in ["cohort", "study", None]:
        raise ValueError(
            "'batch_correct' must be 'cohort', 'study', or None, "
            f"received {batch_correct}"
        )

    if syn is None:
        syn = syn_login()

    if global_correct:
        ptrc_file = File(
            join(
                REPO_PATH,
                "data",
                "phospho_beataml_batch_corrected_global_adjust.csv"
            ),
            parentId="syn25714185",
            id="syn75382181"
        )
        pilot_file = File(
            join(
                REPO_PATH,
                "data",
                "phospho_pilot_batch_corrected_global_adjust.csv"
            ),
            parentId="syn69075538",
            id="syn75382183"
        )
        try:
            # Syncs with Synapse files
            syncFromSynapse(
                syn,
                ptrc_file,
                path=join(REPO_PATH, "data")
            )
            syncFromSynapse(
                syn,
                pilot_file,
                path=join(REPO_PATH, "data")
            )

            # Loads cached results (if this has been run previously) as
            # ComBat-Seq can be time-intensive with race covariates included
            ptrc = pd.read_csv(ptrc_file.path, index_col=0)
            pilot = pd.read_csv(pilot_file.path, index_col=0)

        except (FileNotFoundError, ValueError):
            # Import raw phosphoproteomic measurements
            ptrc, pilot = import_phospho(
                syn,
                batch_correct="study",
                global_correct=False
            )

            # Correct via global
            ptrc, pilot = correct_ptm(ptrc, pilot, syn)

            # Store locally
            ptrc.to_csv(ptrc_file.path)
            pilot.to_csv(pilot_file.path)

            # Sync to Synapse
            syn.store(ptrc_file)
            syn.store(pilot_file)

    elif batch_correct == "study":
        # Setup Synapse file objects
        # This import is a little different from the other 'omic measurements
        # since this preprocessing was written in the present repo and steps
        # are kept here for documentation
        ptrc_file = File(
            join(REPO_PATH, "data", "phospho_beataml_batch_corrected.csv"),
            parentId="syn25714185",
            id="syn73754460"
        )
        pilot_file = File(
            join(REPO_PATH, "data", "phospho_pilot_batch_corrected.csv"),
            parentId="syn69075538",
            id="syn73754465"
        )
        try:
            # Syncs with Synapse files
            syncFromSynapse(
                syn,
                ptrc_file,
                path=join(REPO_PATH, "data")
            )
            syncFromSynapse(
                syn,
                pilot_file,
                path=join(REPO_PATH, "data")
            )

            # Loads cached results (if this has been run previously) as
            # ComBat-Seq can be time-intensive with race covariates included
            ptrc = pd.read_csv(ptrc_file.path, index_col=0)
            pilot = pd.read_csv(pilot_file.path, index_col=0)

        except FileNotFoundError:
            # Import raw phosphoproteomic measurements
            ptrc, pilot = import_phospho(syn, batch_correct=None)
            ptrc, pilot = plex_correct_proteomics(
                ptrc,
                pilot,
                syn=syn
            )

            # Store locally
            ptrc.to_csv(ptrc_file.path)
            pilot.to_csv(pilot_file.path)

            # Sync to Synapse
            syn.store(ptrc_file)
            syn.store(pilot_file)
    else:
        if batch_correct == "cohort":
            ptrc = pd.read_csv(
                syn.get("syn32528196").path, index_col=0, sep="\t"
            )
            pilot = pd.read_csv(
                syn.get("syn69075544").path, index_col=0, sep="\t"
            )
        else:
            ptrc = pd.read_csv(
                syn.get("syn25714936").path, index_col=0, sep="\t"
            )
            pilot = pd.read_csv(
                syn.get("syn69075545").path, index_col=0, sep="\t"
            )

            # Drop phosphosites with >50% missingness in either dataset
            ptrc = ptrc.loc[ptrc.isna().mean(axis=1) < 0.5, :]
            pilot = pilot.loc[pilot.isna().mean(axis=1) < 0.5, :]

        shared_phospho = ptrc.index.intersection(pilot.index)
        ptrc = ptrc.loc[shared_phospho, :]
        pilot = pilot.loc[shared_phospho, :]

        # Convert sample IDs to Accession IDs
        ptrc_conversion, pilot_conversion = import_sample_conversion(syn)
        ptrc.columns = ptrc.columns.astype(int)
        ptrc.rename(columns=ptrc_conversion, inplace=True)
        pilot.rename(columns=pilot_conversion, inplace=True)

        # Relabel bridging samples
        bridge_columns = pilot.columns.intersection(ptrc.columns)
        pilot.rename(
            columns=pd.Series(bridge_columns + "-Bridge", index=bridge_columns),
            inplace=True,
        )

    return ptrc.T, pilot.T


def import_global(
    syn: sc.Synapse | None = None,
    batch_correct: str | None = "study",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads global proteomic data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        batch_correct (str | None, default: "cohort"): Batch-correction method.
            Accepts 'cohort' (corrected at cohort-level), 'study' (corrected at
            study-level, post-cohort concatenation), or None.

    Returns:
        pd.DataFrame: Global data for BeatAML 210 cohort.
        pd.DataFrame: global data for Pilot cohort.
    """
    if batch_correct not in ["cohort", "study", None]:
        raise ValueError(
            "'batch_correct' must be 'cohort', 'study', or None, "
            f"received {batch_correct}"
        )

    if syn is None:
        syn = syn_login()

    if batch_correct == "study":
        # Setup Synapse file objects
        ptrc_file = File(
            join(REPO_PATH, "data", "global_beataml_batch_corrected.csv"),
            parentId="syn25714186",
            id="syn73755181"
        )
        pilot_file = File(
            join(REPO_PATH, "data", "global_pilot_batch_corrected.csv"),
            parentId="syn69075539",
            id="syn73755182"
        )
        try:
            # Syncs with Synapse files
            syncFromSynapse(
                syn,
                ptrc_file,
                path=join(REPO_PATH, "data")
            )
            syncFromSynapse(
                syn,
                pilot_file,
                path=join(REPO_PATH, "data")
            )

            # Loads cached results (if this has been run previously) as
            # ComBat-Seq can be time-intensive with race covariates included
            ptrc = pd.read_csv(ptrc_file.path, index_col=0)
            pilot = pd.read_csv(pilot_file.path, index_col=0)

        except FileNotFoundError:
            # Import raw global measurements
            ptrc, pilot = import_global(syn, batch_correct=None)
            ptrc, pilot = plex_correct_proteomics(
                ptrc,
                pilot,
                syn=syn
            )

            # Store locally
            ptrc.to_csv(ptrc_file.path)
            pilot.to_csv(pilot_file.path)

            # Sync to Synapse
            syn.store(ptrc_file)
            syn.store(pilot_file)
    else:
        if batch_correct == "cohort":
            ptrc = pd.read_csv(
                syn.get("syn25714248").path, index_col=0, sep="\t"
            )
            ptrc.columns = ptrc.columns.astype(int)
            pilot = pd.read_csv(
                syn.get("syn69075555").path, index_col=0, sep="\t"
            )
        else:
            # Raw measurements
            ptrc = pd.read_csv(
                syn.get("syn25714254").path, index_col=0, sep="\t"
            )
            ptrc.columns = ptrc.columns.astype(int)
            pilot = pd.read_csv(
                syn.get("syn69075554").path, index_col=0, sep="\t"
            )

            # Drop proteins with >50% missingness in either dataset
            ptrc = ptrc.loc[ptrc.isna().mean(axis=1) < 0.5, :]
            pilot = pilot.loc[pilot.isna().mean(axis=1) < 0.5, :]

        ptrc_conversion, pilot_conversion = import_sample_conversion(syn)
        pilot_conversion.loc[
            pilot_conversion.isin(ptrc_conversion)
        ] += "-Bridge"
        ptrc.rename(columns=ptrc_conversion, inplace=True)
        pilot.rename(columns=pilot_conversion, inplace=True)

    # Remove measurements missing across all of Pilot or BeatAML cohorts
    pilot = pilot.dropna(how="all", axis=0)
    ptrc = ptrc.dropna(how="all", axis=0)

    pilot = pilot.loc[pilot.index.intersection(ptrc.index), :]
    ptrc = ptrc.loc[pilot.index, :]

    return ptrc.T, pilot.T


def import_meta(
    syn: sc.Synapse | None = None,
    aux_meta: bool = False,
    use_wes: bool = True
) -> pd.DataFrame:
    """
    Loads merged meta-data from Synapse.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        aux_meta (bool, default: False): Load auxiliary meta-data for patients
            not included in manuscript analyses.
        use_wes (bool, default: False): Use WES for mutation calls, otherwise
            uses existing meta-data.

    Returns:
         pd.DataFrame: Updated meta-data across cohorts
    """
    if syn is None:
        syn = syn_login()

    meta = pd.read_csv(syn.get("syn69692583").path, index_col=0)

    # Rename bridging samples
    meta_index = meta.index.to_numpy()
    meta_index[
        np.logical_and(
            meta.loc[:, "source"] == "BeatAML",
            meta.loc[:, "Study"] == "pilotStudy",
        )
    ] += "-Bridge"
    meta.index = meta_index

    # Drop duplicates
    meta = meta.loc[~meta.index.duplicated(keep="first"), :]
    meta = meta.rename(columns={"source": "Source"})

    # Fill in missing/unknown with inferred race
    mapping = pd.read_excel(syn.get("syn64126463").path, index_col=0)
    mapping.reset_index(drop=True, inplace=True)
    ba_aux = pd.read_excel(
        syn.get("syn64126458").path, header=1, index_col=0, sheet_name="summary"
    )
    ba_aux.reset_index(inplace=True, drop=True)

    mapping_ids = mapping.loc[:, "dbgap_rnaseq_sample"].str[:-1]
    mapping_ids.loc[mapping_ids.isna()] = mapping.loc[
        :, "dbgap_dnaseq_sample"
    ].str[:-1]
    sample_conversions = pd.Series(
        mapping.loc[:, "labId"].values, index=mapping_ids
    )

    sample_names = ba_aux.loc[:, "dbgap_rnaseq_sample"].str[:-1]
    sample_names.loc[sample_names.isna()] = ba_aux.loc[
        :, "dbgap_dnaseq_sample"
    ].str[:-1]
    inferred_race = pd.Series(
        ba_aux.loc[:, "inferred_ethnicity"].values,
        index=sample_names.replace(sample_conversions),
    )
    inferred_race = inferred_race.loc[
        meta.loc[
            meta.loc[:, "Race"].isin([np.nan, "Unknown"]), :
        ].index.intersection(inferred_race.index)
    ]
    meta.loc[inferred_race.index, "Race"] = inferred_race

    if aux_meta:
        # Rename columns to match format of pilot cohort
        ba_aux.set_index("dbgap_rnaseq_sample", inplace=True, drop=True)
        ba_aux = ba_aux.loc[~ba_aux.index.isna(), :]
        ba_aux.rename(
            columns={
                "ageAtDiagnosis": "Age",
                "consensus_sex": "Sex",
                "reportedRace": "Race",
                "specimenType": "Specimen_Type",
            },
            inplace=True,
        )
        ba_aux = ba_aux.loc[:, ["Age", "Sex", "Race", "Specimen_Type"]]
        ba_aux.loc[:, ["Source", "Study"]] = "BeatAML"

        # Convert sample names to match BeatAML study sample IDs
        ptrc_conversion = pd.read_excel(
            syn.get("syn64126463").path, index_col=0
        )
        ptrc_conversion.dropna(subset="dbgap_rnaseq_sample", inplace=True)
        ptrc_conversion = pd.Series(
            ptrc_conversion.loc[:, "labId"].values,
            index=ptrc_conversion.loc[:, "dbgap_rnaseq_sample"].values,
        )
        ba_aux.rename(index=ptrc_conversion, inplace=True)

        # Concatenate auxiliary meta-data with existing meta-data, keep
        # auxiliary samples if duplicated in existing meta-data
        meta = pd.concat([ba_aux, meta])
        meta = meta.loc[~meta.index.duplicated(keep="first"), :]

        # Rename columns to match pilot cohort format
        pilot_aux = pd.read_csv(syn.get("syn73997552").path, index_col=0)
        pilot_aux.rename(
            columns={"Sex (1-male, 0-female)": "Sex", "Race_label": "Race"},
            inplace=True,
        )
        pilot_aux.loc[:, ["Source", "Study"]] = "pilotStudy"
        pilot_aux = pilot_aux.loc[
            :, pilot_aux.columns.intersection(meta.columns)
        ]

        # Concatenate pilot auxiliary meta-data with existing meta-data, keep
        # auxiliary samples if duplicated in existing meta-data
        meta = pd.concat([pilot_aux, meta])
        meta = meta.loc[~meta.index.duplicated(keep="first"), :]
        meta.loc[meta.index.str.contains("Bridge"), "Specimen_Type"] = meta.loc[
            meta.loc[meta.index.str.contains("Bridge")].index.str[:-7],
            "Specimen_Type",
        ].values

        # Standardize specimen types across meta data sets
        meta.replace(
            {
                "Specimen_Type": {
                    "Bone Marrow Aspirate": "BM",
                    "Peripheral Blood": "PB",
                    "Leukapheresis": "PB",
                    "Not Specified": "BM",
                }
            },
            inplace=True,
        )

        # Reorders columns, putting study/source meta-data to the left and
        # mutation data to the right
        col_order = list(meta.loc[:, :"Race"].columns) + ["Source", "Study"]
        col_order = col_order + sorted(
            list(meta.drop(col_order, axis=1).columns)
        )
        meta = meta.loc[:, col_order]

    # Standardizes Race column
    meta.replace(
        {"Race": {"Declined": "Unknown", np.nan: "Unknown"}}, inplace=True
    )

    if use_wes:
        # Gets WES calls
        mutations = call_mutations(syn)

        # Sets up WES calls for bridging patients
        bridge = meta.index[meta.index.str.endswith("-Bridge")]
        bridge = mutations.loc[
            mutations.index.intersection(bridge.str[:-7]),
            :
        ]
        bridge.index = bridge.index + "-Bridge"
        mutations = pd.concat([mutations, bridge])

        # Trims to mutations already in meta data
        preserved_mutations = mutations.columns.intersection(meta.columns)

        # Update mutation calls
        meta.loc[
            mutations.index.intersection(meta.index),
            preserved_mutations
        ] = mutations.loc[
            mutations.index.intersection(meta.index),
            preserved_mutations
        ]

        # Adds column to denote mutation call source
        meta.loc[:, "WES_call"] = False
        meta.loc[
            mutations.index.intersection(meta.index),
            "WES_call"
        ] = True

        # Overwrite non-Mutant cases as WT
        # This only works with WES cases since 'Not measured' is equivalent to
        # wild-type in WES measurements
        meta.loc[:, "ALT":] = meta.loc[:, "ALT":].replace(
            {
                np.nan: "WT",
                "Not measured": "WT"
            }
        )

    # Add FLT3-ITD to FLT3 mutation calls
    meta.loc[
        np.logical_or(
            meta.loc[:, "FLT3"] == "Mutant",
            meta.loc[:, "FLT3_ITD"] == "Mutant"
        ),
        "FLT3"
    ] = "Mutant"

    return meta


def import_acetyl(
    syn: sc.Synapse | None = None,
    global_correct: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads Acetylomics data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        global_correct (bool, default=True): Correct measurements using global
            proteomic abundances.

    Returns:
        pd.DataFrame: Acetylomics data.
    """
    if syn is None:
        syn = syn_login()

    if global_correct:
        ptrc_file = File(
            join(
                REPO_PATH,
                "data",
                "acetyl_beataml_global_adjust.csv"
            ),
            parentId="syn53484023",
            id="syn76033232"
        )
        pilot_file = File(
            join(
                REPO_PATH,
                "data",
                "acetyl_pilot_global_adjust.csv"
            ),
            parentId="syn69075540",
            id="syn76033235"
        )
        try:
            # Syncs with Synapse files
            syncFromSynapse(
                syn,
                ptrc_file,
                path=join(REPO_PATH, "data")
            )
            syncFromSynapse(
                syn,
                pilot_file,
                path=join(REPO_PATH, "data")
            )

            # Loads cached results (if this has been run previously) as
            # the amount of OLS being run can be time-intensive
            ptrc = pd.read_csv(ptrc_file.path, index_col=0)
            pilot = pd.read_csv(pilot_file.path, index_col=0)

        except (FileNotFoundError, ValueError):
            # Import raw phosphoproteomic measurements
            ptrc, pilot = import_acetyl(
                syn,
                global_correct=False
            )

            # Correct via global
            ptrc, pilot = correct_ptm(ptrc, pilot, syn)

            # Store locally
            ptrc.to_csv(ptrc_file.path)
            pilot.to_csv(pilot_file.path)

            # Sync to Synapse
            syn.store(ptrc_file)
            syn.store(pilot_file)

    else:
        ptrc = pd.read_csv(syn.get("syn53484994").path, index_col=0, sep="\t")
        pilot = pd.read_csv(syn.get("syn69075568").path, index_col=0, sep="\t")
        ptrc_conversion, pilot_conversion = import_sample_conversion(syn)

        ptrc.columns = [int(col.split("_")[-1]) for col in ptrc.columns]

        ptrc.rename(columns=ptrc_conversion, inplace=True)
        pilot.rename(columns=pilot_conversion, inplace=True)

        ptrc = ptrc.loc[~ptrc.index.str.contains("NULL"), :]
        pilot = pilot.loc[~pilot.index.str.contains("NULL"), :]

        ptrc = ptrc.loc[ptrc.index.intersection(pilot.index), :]
        pilot = pilot.loc[ptrc.index, :]

        bridge_columns = pilot.columns.intersection(ptrc.columns)
        pilot.rename(
            columns=pd.Series(bridge_columns + "-Bridge", index=bridge_columns),
            inplace=True,
        )

    return ptrc.T, pilot.T


def import_sample_conversion(
    syn: sc.Synapse | None = None,
) -> tuple[pd.Series, pd.Series]:
    """
    Loads conversions between sample and patient IDs.

    Args:
        syn (sc.Synapse): Logged-in synapse client.

    Returns:
        pd.Series: Series mapping PTRC sample to patient IDs.
        pd.Series: Series mapping Pilot study samples to patient IDs.
    """
    ptrc_conversion = pd.read_csv(
        syn.get("syn25807733").path, index_col=0, sep="\t"
    )
    pilot_conversion = pd.read_excel(
        syn.get("syn68835814").path, sheet_name="TMT"
    )

    ptrc_conversion = pd.Series(
        ptrc_conversion.loc[:, "Barcode.ID"].values,
        index=ptrc_conversion.loc[:, "SampleID.abbrev"],
    )
    pilot_conversion = pd.Series(
        pilot_conversion.loc[:, "Accession"].values,
        index=pilot_conversion.loc[:, "New Sample ID"],
    )

    # Remove PooledReference and EMPTY samples
    pilot_conversion = pilot_conversion.loc[
        ~pilot_conversion.isin(["PooledReference", "EMPTY"])
    ]

    return ptrc_conversion, pilot_conversion


def convert_gene_symbols(gene_list: Iterable) -> list[str]:
    """
    Converts ensembl genes to symbols.

    Args:
        gene_list (list[str]): ensembl gene symbols.

    Returns:
        list[str]: list of converted gene symbols.
    """
    mg = mygene.MyGeneInfo()
    result = mg.getgenes(gene_list, as_dataframe=True)
    result = result.loc[~result.index.duplicated(), "symbol"]
    result = result.loc[gene_list]

    return list(result)


def query_proteins(genes: Iterable[str]) -> pd.DataFrame:
    """
    Queries protein IDs for provided genes.

    Args:
        genes (Iterable[str]): Genes to lookup.

    Returns:
        pd.DataFrame: UniProt entries for each gene.
    """
    tsv = ""
    for gene in genes:
        response = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query=reviewed:true+AND+gene:{gene}+AND+organism_name:Human&"
            f"format=tsv&size=1"
        ).text
        response = response.split("\n")

        if len(tsv) == 0:
            tsv += response[0] + "\n"

        tsv += response[1] + "\n"

    uniprot_df = pd.read_csv(StringIO(tsv), sep="\t", index_col=0)

    return uniprot_df


def get_ksea_table() -> pd.DataFrame:
    """
    Loads phospho meta-data for KSEA.

    Args:
        None.

    Returns:
        pd.DataFrame: Phospho meta-data for KSEA.
    """
    syn = syn_login()
    ptrc, pilot = import_phospho(syn, pre_corrected=True)
    phospho_meta = pd.read_csv(
        join(REPO_PATH, "../data", "Concatenated_msgfplus_syn_plus_ascore.txt"),
        index_col=0,
        sep="\t",
    )
    conversions = pd.read_csv(
        syn.get("syn25714920").path, sep="\t", index_col=0
    )

    lookup = conversions.index.str.split("@", expand=True)
    lookup = pd.Series(
        lookup.get_level_values(1).values,
        index=lookup.get_level_values(0).values,
    )
    lookup = lookup.loc[lookup.isin(phospho_meta.loc[:, "Peptide"])]
    lookup = lookup.loc[~lookup.index.duplicated()]

    ksea_table = pd.DataFrame(
        list(pilot.columns.str.split("-", expand=True)),
        index=pilot.columns,
        columns=["Gene", "Residue.Both"],
    )

    ksea_table = ksea_table.loc[ksea_table.index.intersection(lookup.index), :]
    ksea_table.loc[:, "Peptide"] = lookup.loc[ksea_table.index]

    ksea_table.loc[:, "Residue.Both"] = (
        ksea_table.loc[:, "Residue.Both"]
        .str.replace(r"[a-z]+", ";", regex=True)
        .str[:-1]
    )
    ksea_table.loc[:, "Peptide"] = ksea_table.loc[:, "Peptide"].str.replace(
        r"[.*]+", "", regex=True
    )

    ksea_table.loc[:, "Protein"] = "NULL"
    ksea_table.loc[:, "p"] = "NULL"
    ksea_table.loc[:, "FC"] = None
    ksea_table = ksea_table.loc[
        :, ["Protein", "Gene", "Peptide", "Residue.Both", "p", "FC"]
    ]

    return ksea_table


def import_deconvolution(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads cell type deconvolution data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.

    Returns:
        pd.DataFrame: Cell-type deconvolution data.
    """
    if syn is None:
        syn = syn_login()

    data = pd.read_csv(syn.get("syn69907563").path, sep="\t")
    data = pd.pivot(
        data, columns="cell_type", index="Accession", values="value"
    )
    data = data.T / data.sum(axis=1)

    return data.T
