from pysb import *

Model()

Monomer('S')
Monomer('ER')
Monomer('LR')

Initial(S(), Parameter('S_0', 1000))

Observable('S_tot', S())
Observable('ER_tot', ER())
Observable('LR_tot', LR())
Observable('Tot_cells', LR() + S() + ER())

Parameter('k_S_div_0', 0.1)
Parameter('k_S_die_0', 0.2)

Rule('S_div', S() >> S() + S(), k_S_div_0)
Rule('S_die', S() >> None, k_S_die_0)

Parameter('k_ER_div_0', 0.05)
Parameter('k_ER_die', 0.1)

Rule('ER_div', ER() >> ER() + ER(), k_ER_div_0)
Rule('ER_die', ER() >> None, k_ER_die)

Parameter('k_LR_div', 0.1)
Parameter('k_LR_die', 0.01)
Rule('LR_div', LR() >> LR() + LR(), k_LR_div)
Rule('LR_die', LR() >> None, k_LR_die)

Parameter('k_S_diff_ER', 0.001)
Parameter('k_ER_diff_LR', 0.001)
Rule('S_diff_ER', S() >> ER(), k_S_diff_ER)
Rule('ER_diff_LR', ER() >> LR(), k_ER_diff_LR)
