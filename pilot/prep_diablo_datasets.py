import pandas as pd

from pilot.data_import import (import_global, import_meta, import_phospho,
                               syn_login)


def main():
    syn = syn_login()
    ptrc_phospho, pilot_phospho = import_phospho(syn, corrected=True)
    ptrc_global, pilot_global = import_global(syn)
    phospho = pd.concat([ptrc_phospho, pilot_phospho])
    global_prot = pd.concat([ptrc_global, pilot_global])

    meta = import_meta(syn)
    meta = meta.loc[meta.loc[:, "Race"].isin(["Black", "White"]), "Race"]

    phospho = phospho.loc[meta.index, :]
    global_prot = global_prot.loc[meta.index, :]

    phospho.to_csv("data/phospho.txt.gz", index_label=False)
    global_prot.to_csv("data/global.txt.gz", index_label=False)
    meta.to_csv("data/race_label.txt.gz", index_label=False)


if __name__ == "__main__":
    main()
