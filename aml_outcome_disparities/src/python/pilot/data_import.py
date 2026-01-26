import os.path
import pickle
import warnings
from io import StringIO
from os import PathLike
from os.path import abspath, dirname, join
from typing import Iterable

import mygene
import numpy as np
import pandas as pd
import requests
import synapseclient as sc
from inmoose.pycombat import pycombat_seq

# Suppress SettingWithCopyWarning
# Seems to be a false warning for some import operations
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
REPO_PATH = abspath(dirname(dirname(__file__)))


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

    syn = sc.Synapse()
    with open(auth_path, "r") as f:
        auth_token = f.read()

    syn.login(authToken=auth_token)

    return syn


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

    # TODO: Upload adjusted R code and meta-data file
    lipid_meta = pd.read_csv(
        join(REPO_PATH, "..", "r", "lipids", "f_data_lipid.csv")
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

    Returns:
         pd.DataFrame: RNA data for BeatAML 210 cohort.
         pd.DataFrame: RNA data for Pilot cohort.
    """
    if syn is None:
        syn = syn_login()

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
    pilot.rename(columns=pilot_conversion, inplace=True)

    ptrc_conversion = pd.read_excel(syn.get("syn64126463").path, index_col=0)
    ptrc_conversion.dropna(subset="dbgap_rnaseq_sample", inplace=True)
    ptrc_conversion = pd.Series(
        ptrc_conversion.loc[:, "labId"].values,
        index=ptrc_conversion.loc[:, "dbgap_rnaseq_sample"].values,
    )
    ptrc.rename(columns=ptrc_conversion, inplace=True)

    # Direct import from BioMart API -- needs attributes
    # "Transcript length (including UTRs and CDS)", and "Gene stable ID"
    # We do this even if were not TPMing measurements to ensure same genes are
    # considered across analyses
    # TODO: Setup direct XML-based import for reproducibility
    gene_lengths = pd.read_csv(
        join(REPO_PATH, "data", "biomart_gene_lengths.txt.gz")
    )
    gene_lengths.set_index("Gene stable ID", inplace=True, drop=True)
    gene_lengths.sort_values("UniProtKB Gene Name symbol", inplace=True)
    gene_lengths = gene_lengths.loc[
        ~gene_lengths.index.duplicated(keep="first")
    ]

    # Trims to genes in both datasets
    shared_genes = ptrc.index.intersection(pilot.index)
    shared_genes = shared_genes.intersection(gene_lengths.index)
    pilot = pilot.loc[shared_genes, :]
    ptrc = ptrc.loc[shared_genes, :]
    gene_lengths = gene_lengths.loc[shared_genes]

    # Correct raw counts via ComBat-Seq prior to normalization or TPM
    if batch_correct:
        # Loads cached results (if this has been run previously) as ComBat-Seq
        # can be time-intensive with race covariates included
        if os.path.exists(join(REPO_PATH, "data", "rna_combat_seq.pkl")):
            (ptrc, pilot) = pd.read_pickle(
                join(REPO_PATH, "data", "rna_combat_seq.pkl")
            )
        else:
            # Load in meta-data for larger cohort, including patients outside
            # pilot study
            meta = import_meta(syn, aux_meta=True)
            meta = meta.loc[
                ~meta.loc[:, "Race"].isin(["Declined", "Unknown"]), :
            ]

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

            # Store batch-corrected to .pkl for quicker loading
            with open(join(REPO_PATH, "data", "rna_combat_seq.pkl"), "wb") as f:
                pickle.dump((ptrc, pilot), f)

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
        ptrc.rename(
            index=gene_lengths.loc[:, "UniProtKB Gene Name symbol"],
            inplace=True,
        )
        pilot.rename(
            index=gene_lengths.loc[:, "UniProtKB Gene Name symbol"],
            inplace=True,
        )

    return ptrc.T, pilot.T


def import_phospho(
    syn: sc.Synapse | None = None, pre_corrected: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads phosphoproteomic data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        pre_corrected (bool, default: True): Load pre-corrected data; otherwise,
            load post-concatenation corrected data.

    Returns:
        pd.DataFrame: Phosphoproteomic data for BeatAML 210 cohort.
        pd.DataFrame: Phosphoproteomic data for Pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    if pre_corrected:
        ptrc = pd.read_csv(syn.get("syn32528196").path, index_col=0, sep="\t")
        pilot = pd.read_csv(syn.get("syn69075544").path, index_col=0, sep="\t")
        shared_phospho = ptrc.index.intersection(pilot.index)
        ptrc = ptrc.loc[shared_phospho, :]
        pilot = pilot.loc[shared_phospho, :]
    else:
        ptrc = pd.read_csv(
            join(REPO_PATH, "data", "ba_phospho_corrected.csv"), index_col=0
        )
        ptrc.columns = ptrc.columns.str[1:]
        pilot = pd.read_csv(
            join(REPO_PATH, "data", "pilot_phospho_corrected.csv"), index_col=0
        )

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


def import_meta(
    syn: sc.Synapse | None = None, aux_meta: bool = False
) -> pd.DataFrame:
    """
    Loads merged meta-data from Synapse.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        aux_meta (bool, default: False): Load auxiliary meta-data for patients
            not included in manuscript analyses.

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

    if aux_meta:
        # Import additional metadata from BeatAML study
        ba_aux = pd.read_excel(
            syn.get("syn64126458").path, index_col=0, sheet_name="summary"
        )
        ba_aux.columns = ba_aux.iloc[0, :]
        ba_aux = ba_aux.iloc[1:, :]

        # Rename columns to match format of pilot cohort
        ba_aux.set_index("dbgap_rnaseq_sample", inplace=True, drop=True)
        ba_aux = ba_aux.loc[~ba_aux.index.isna(), :]
        ba_aux.rename(
            columns={
                "ageAtDiagnosis": "Age",
                "consensus_sex": "Sex",
                "reportedRace": "Race",
            },
            inplace=True,
        )
        ba_aux = ba_aux.loc[:, ["Age", "Sex", "Race"]]
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
        meta = pd.concat([meta, ba_aux])
        meta = meta.loc[~meta.index.duplicated(keep="first"), :]

        # Rename columns to match pilot cohort format
        pilot_aux = pd.read_excel(syn.get("syn62750469").path, index_col=0)
        pilot_aux.drop(pilot_aux.index.intersection(meta.index), inplace=True)
        pilot_aux.rename(
            columns={"Sex (1-male, 0-female)": "Sex", "Race_label": "Race"},
            inplace=True,
        )
        pilot_aux.loc[:, ["Source", "Study"]] = "pilotStudy"
        pilot_aux = pilot_aux.loc[
            :, pilot_aux.columns.intersection(meta.columns)
        ]
        meta = pd.concat([meta, pilot_aux])

    return meta


def import_acetyl(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads Acetylomics data.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.

    Returns:
        pd.DataFrame: Acetylomics data.
    """
    if syn is None:
        syn = syn_login()

    acetyl = pd.read_csv(syn.get("syn69075568").path, index_col=0, sep="\t")
    _, pilot_conversion = import_sample_conversion(syn)

    acetyl.rename(columns=pilot_conversion, inplace=True)
    acetyl = acetyl.loc[~acetyl.index.str.contains("NULL"), :]

    return acetyl.T


def import_global(
    syn: sc.Synapse | None = None, pre_corrected: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads global proteomic measurements.

    Args:
        syn (sc.Synapse | None, default: None): Logged-in Synapse object; loads
            new one if None.
        pre_corrected (bool, default: True): Load pre-corrected data; otherwise,
            load post-concatenation corrected data.

    Returns:
        pd.DataFrame: Global proteomics for 210 cohort.
        pd.DataFrame: Global proteomics for pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    if pre_corrected:
        ptrc = pd.read_csv(syn.get("syn25714248").path, index_col=0, sep="\t")
        ptrc.columns = ptrc.columns.astype(int)
        pilot = pd.read_csv(syn.get("syn69075555").path, index_col=0, sep="\t")
    else:
        ptrc = pd.read_csv(
            join(REPO_PATH, "data", "ba_global_corrected.csv"), index_col=0
        )
        ptrc.columns = ptrc.columns.str[1:].astype(int)
        pilot = pd.read_csv(
            join(REPO_PATH, "data", "pilot_global_corrected.csv"), index_col=0
        )

    ptrc_conversion, pilot_conversion = import_sample_conversion(syn)
    pilot_conversion.loc[pilot_conversion.isin(ptrc_conversion)] += "-Bridge"
    ptrc.rename(columns=ptrc_conversion, inplace=True)
    pilot.rename(columns=pilot_conversion, inplace=True)

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
