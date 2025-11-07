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

REPO_PATH = abspath(dirname(dirname(__file__)))


def syn_login(auth_path: PathLike | None = None) -> sc.Synapse:
    """
    Login to Synapse.

    Args:
        auth_path (PathLike | None): Path to authentication file.

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


def batch_correct(
    data: pd.DataFrame, batches: Iterable, race: Iterable | None = None
) -> pd.DataFrame:
    """
    Corrects for batches via ComBat-seq through inmoose.

    Args:
        data (pd.DataFrame): Data to batch-correct.
        batches (Iterable): Batches for samples in data.
        race (Iterable): Race covariate.

    Returns:
        pd.DataFrame: Corrected data.
    """
    if race is None:
        data = pycombat_seq(data, batches)
    else:
        data = pycombat_seq(data, batches, covar_mod=race)

    return data.T


def import_rna(syn: sc.Synapse | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads and TPMs RNA data.

    Args:
        syn (sc.Synapse): Logged-in Synapse object; loads new one if None.

    Returns:
         pd.DataFrame: TPM RNA data.
    """
    if syn is None:
        syn = syn_login()

    ptrc = pd.read_csv(syn.get("syn64126462").path, index_col=0, sep="\t")
    pilot = pd.read_csv(syn.get("syn68820229").path, index_col=0, sep="\t")
    gene_lengths = pd.read_csv(join(REPO_PATH, "data", "read_lengths.txt"), index_col=0)
    gene_lengths = gene_lengths.loc[
        ~gene_lengths.loc[:, "Ensembl Canonical"].isna(),
        "Transcript length (including UTRs and CDS)",
    ].squeeze()

    # Trims to genes in both datasets
    shared_genes = ptrc.index.intersection(pilot.index)
    shared_genes = shared_genes.intersection(gene_lengths.index)
    pilot = pilot.loc[shared_genes, :].T
    ptrc = ptrc.loc[shared_genes, :].T
    gene_lengths = gene_lengths.loc[shared_genes]

    # TPM data
    ptrc /= gene_lengths / 1000
    sums = ptrc.sum(axis=1) / 1e6
    ptrc = ptrc.T / sums

    pilot /= gene_lengths / 1000
    sums = pilot.sum(axis=1) / 1e6
    pilot = pilot.T / sums

    return ptrc.T, pilot.T


def import_phospho(
    syn: sc.Synapse | None = None, corrected: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads phosphoproteomic data.

    Args:
        syn (sc.Synapse): Logged-in Synapse object; loads new one if None.
        corrected (bool, default: True): Load corrected data.

    Returns:
        pd.DataFrame: Phosphoproteomic data across studies.
    """
    if syn is None:
        syn = syn_login()

    if corrected:
        ptrc = pd.read_csv(syn.get("syn32528196").path, index_col=0, sep="\t")
        pilot = pd.read_csv(
            join(REPO_PATH, "data", "ptrc_ex26_crosstab_phospho_siteid_corrected.txt"),
            index_col=0,
            sep="\t",
        )
    else:
        ptrc = pd.read_csv(syn.get("syn25714936").path, index_col=0, sep="\t")
        pilot = pd.read_csv(
            join(REPO_PATH, "data", "ptrc_ex26_crosstab_phospho_siteid_original.txt"),
            index_col=0,
            sep="\t",
        )

    ptrc_conversion, pilot_conversion = import_sample_conversion(syn)

    ptrc.columns = ptrc.columns.astype(int)
    ptrc.rename(columns=ptrc_conversion, inplace=True)
    pilot.rename(columns=pilot_conversion, inplace=True)

    shared_phospho = ptrc.index.intersection(pilot.index)
    ptrc = ptrc.loc[shared_phospho, :]
    pilot = pilot.loc[shared_phospho, :]

    bridge_columns = pilot.columns.intersection(ptrc.columns)
    pilot.rename(
        columns=pd.Series(bridge_columns + "-Bridge", index=bridge_columns),
        inplace=True,
    )

    return ptrc.T, pilot.T


def import_meta(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads merged meta-data from Synapse.

    Args:
        syn (sc.Synapse): Logged-in Synapse object; loads new one if None.

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
            meta.loc[:, "source"] == "BeatAML", meta.loc[:, "study"] == "Pilot"
        )
    ] += "-Bridge"
    meta.index = meta_index

    # Drop duplicates
    meta = meta.loc[~meta.duplicated(), :]

    # Collate with patients that have measurements
    ptrc, pilot = import_phospho(syn)
    phospho = pd.concat([ptrc, pilot])
    meta = meta.reindex(phospho.index)

    meta = meta.rename(columns={"source": "Source"})

    return meta


def import_acetyl(syn: sc.Synapse | None = None) -> pd.DataFrame:
    """
    Loads Acetylomics data.

    Args:
        syn (sc.Synapse): Logged-in Synapse client.

    Returns:
        pd.DataFrame: Acetylomics data.
    """
    if syn is None:
        syn = syn_login()

    acetyl = pd.read_csv(syn.get("syn69075568").path, index_col=0, sep="\t")
    _, pilot_conversion = import_sample_conversion(syn)
    acetyl.rename(columns=pilot_conversion, inplace=True)

    return acetyl


def import_global(syn: sc.Synapse | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads global proteomic measurements.

    Args:
        syn (sc.Synapse): Logged-in synapse client.

    Returns:
        pd.DataFrame: Global proteomics for 210 cohort.
        pd.DataFrame: Global proteomics for pilot cohort.
    """
    if syn is None:
        syn = syn_login()

    ptrc = pd.read_csv(syn.get("syn25714248").path, index_col=0, sep="\t")
    ptrc.columns = ptrc.columns.astype(int)
    pilot = pd.read_csv(syn.get("syn69075555").path, index_col=0, sep="\t")

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
    ptrc_conversion = pd.read_csv(syn.get("syn25807733").path, index_col=0, sep="\t")
    pilot_conversion = pd.read_excel(syn.get("syn68835814").path, sheet_name="TMT")

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
    ptrc, pilot = import_phospho(syn, corrected=True)
    phospho_meta = pd.read_csv(
        join(REPO_PATH, "data", "Concatenated_msgfplus_syn_plus_ascore.txt"),
        index_col=0,
        sep="\t",
    )
    conversions = pd.read_csv(syn.get("syn25714920").path, sep="\t", index_col=0)

    lookup = conversions.index.str.split("@", expand=True)
    lookup = pd.Series(
        lookup.get_level_values(1).values, index=lookup.get_level_values(0).values
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
        syn (sc.Synapse): Logged-in Synapse client. Creates new one if None
            is provided.

    Returns:
        pd.DataFrame: Cell-type deconvolution data.
    """
    if syn is None:
        syn = syn_login()

    data = pd.read_csv(syn.get("syn69907563").path, sep="\t")
    data = pd.pivot(data, columns="cell_type", index="Accession", values="value")
    data = data.T / data.sum(axis=1)

    return data.T
