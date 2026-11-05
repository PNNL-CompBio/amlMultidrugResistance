"""
    define a common color palette to use across the different figures
"""


from enum import StrEnum


class COLORS(StrEnum):
    BLACK_WT = "#3B31FC"
    BLACK_NPM1 = "#1D99B8"
    BLACK_FLT3ITD = "#0CA186"
    BLACK_DOUBLE = "#0E3162"
    WHITE_WT = "#E3B505"
    WHITE_NPM1 = "#D47F00"
    WHITE_FLT3ITD = "#E8487B"
    WHITE_DOUBLE = "#751323"



# map race/mutation labels to colors
# key: (race, mutation)
# value: color
RACE_MUT_TO_COLOR = {
    ("Black", "NPM1-/FLT3-ITD-"): COLORS.BLACK_WT,
    ("Black", "NPM1+/FLT3-ITD-"): COLORS.BLACK_NPM1,
    ("Black", "NPM1-/FLT3-ITD+"): COLORS.BLACK_FLT3ITD,
    ("Black", "NPM1+/FLT3-ITD+"): COLORS.BLACK_DOUBLE,
    ("White", "NPM1-/FLT3-ITD-"): COLORS.WHITE_WT,
    ("White", "NPM1+/FLT3-ITD-"): COLORS.WHITE_NPM1,
    ("White", "NPM1-/FLT3-ITD+"): COLORS.WHITE_FLT3ITD,
    ("White", "NPM1+/FLT3-ITD+"): COLORS.WHITE_DOUBLE,
}
