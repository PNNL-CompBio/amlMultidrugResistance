from magine.data.experimental_data import ExperimentalData
from pathlib import Path
_file = Path(__file__).parent.joinpath('data_ptrc2_for_magine.csv').__str__()
exp_data = ExperimentalData(_file)