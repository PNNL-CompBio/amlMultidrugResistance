#!/bin/zsh

# NPM1 vs. WT, black patients

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/B_NPM1-WT_Lipidomics.csv \
        _cache/diablo/B_NPM1-WT_Metabolomics.csv \
        _cache/diablo/B_NPM1-WT_Phosphoproteomics.csv \
        _cache/diablo/B_NPM1-WT_Proteomics.csv  \
    --block-names \
        Lipidomics \
        Metabolomics \
        Phosphoproteomics \
        Proteomics  \
    --design-csv _cache/diablo/B_NPM1-WT_null-design.csv  \
    --target-csv _cache/diablo/B_NPM1-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/B_NPM1-WT/null/  \
    --n-components 3  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/B_NPM1-WT_null-keepX.csv  \
| tee _cache/diablo/B_NPM1-WT_final-null.log

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/B_NPM1-WT_Lipidomics.csv \
        _cache/diablo/B_NPM1-WT_Metabolomics.csv \
        _cache/diablo/B_NPM1-WT_Phosphoproteomics.csv \
        _cache/diablo/B_NPM1-WT_Proteomics.csv  \
    --block-names \
        Lipidomics \
        Metabolomics \
        Phosphoproteomics \
        Proteomics  \
    --design-csv _cache/diablo/B_NPM1-WT_plsda-design.csv  \
    --target-csv _cache/diablo/B_NPM1-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/B_NPM1-WT/plsda/  \
    --n-components 2  \
    --distance-metric max.dist  \
    --keepx-csv _cache/diablo/B_NPM1-WT_plsda-keepX.csv  \
| tee _cache/diablo/B_NPM1-WT_final-plsda.log

# NPM1 vs. WT, white patients

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/W_NPM1-WT_Transcriptomics.csv \
        _cache/diablo/W_NPM1-WT_Phosphoproteomics.csv \
        _cache/diablo/W_NPM1-WT_Proteomics.csv  \
    --block-names \
        Transcriptomics\
        Phosphoproteomics \
        Proteomics  \
    --design-csv _cache/diablo/W_NPM1-WT_null-design.csv  \
    --target-csv _cache/diablo/W_NPM1-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/W_NPM1-WT/null/  \
    --n-components 5  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/W_NPM1-WT_null-keepX.csv  \
| tee _cache/diablo/W_NPM1-WT_final-null.log

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/W_NPM1-WT_Phosphoproteomics.csv \
        _cache/diablo/W_NPM1-WT_Transcriptomics.csv \
        _cache/diablo/W_NPM1-WT_Proteomics.csv  \
    --block-names \
        Phosphoproteomics \
        Transcriptomics\
        Proteomics  \
    --design-csv _cache/diablo/W_NPM1-WT_plsda-design.csv  \
    --target-csv _cache/diablo/W_NPM1-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/W_NPM1-WT/plsda/  \
    --n-components 7  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/W_NPM1-WT_plsda-keepX.csv  \
| tee _cache/diablo/W_NPM1-WT_final-plsda.log

# NRAS vs. WT, black patients

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/B_NRAS-WT_Metabolomics.csv \
        _cache/diablo/B_NRAS-WT_Phosphoproteomics.csv \
        _cache/diablo/B_NRAS-WT_Proteomics.csv  \
        _cache/diablo/B_NRAS-WT_Acetylomics.csv \
        _cache/diablo/B_NRAS-WT_Lipidomics.csv \
    --block-names \
        Metabolomics \
        Phosphoproteomics \
        Proteomics  \
        Acetylomics \
        Lipidomics \
    --design-csv _cache/diablo/B_NRAS-WT_null-design.csv  \
    --target-csv _cache/diablo/B_NRAS-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/B_NRAS-WT/null/  \
    --n-components 3  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/B_NRAS-WT_null-keepX.csv  \
| tee _cache/diablo/B_NRAS-WT_final-null.log

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/B_NRAS-WT_Acetylomics.csv \
        _cache/diablo/B_NRAS-WT_Metabolomics.csv \
        _cache/diablo/B_NRAS-WT_Phosphoproteomics.csv \
        _cache/diablo/B_NRAS-WT_Proteomics.csv  \
        _cache/diablo/B_NRAS-WT_Lipidomics.csv \
    --block-names \
        Acetylomics \
        Metabolomics \
        Phosphoproteomics \
        Proteomics  \
        Lipidomics \
    --design-csv _cache/diablo/B_NRAS-WT_plsda-design.csv  \
    --target-csv _cache/diablo/B_NRAS-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/B_NRAS-WT/plsda/  \
    --n-components 2  \
    --distance-metric max.dist  \
    --keepx-csv _cache/diablo/B_NRAS-WT_plsda-keepX.csv  \
| tee _cache/diablo/B_NRAS-WT_final-plsda.log

# NRAS vs. WT, white patients

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/W_NRAS-WT_Transcriptomics.csv \
        _cache/diablo/W_NRAS-WT_Phosphoproteomics.csv \
        _cache/diablo/W_NRAS-WT_Proteomics.csv  \
    --block-names \
        Transcriptomics\
        Phosphoproteomics \
        Proteomics  \
    --design-csv _cache/diablo/W_NRAS-WT_null-design.csv  \
    --target-csv _cache/diablo/W_NRAS-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/W_NRAS-WT/null/  \
    --n-components 2  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/W_NRAS-WT_null-keepX.csv  \
| tee _cache/diablo/W_NRAS-WT_final-null.log

./MomDiablo.R final \
    --block-csvs \
        _cache/diablo/W_NRAS-WT_Phosphoproteomics.csv \
        _cache/diablo/W_NRAS-WT_Transcriptomics.csv \
        _cache/diablo/W_NRAS-WT_Proteomics.csv  \
    --block-names \
        Phosphoproteomics \
        Transcriptomics\
        Proteomics  \
    --design-csv _cache/diablo/W_NRAS-WT_plsda-design.csv  \
    --target-csv _cache/diablo/W_NRAS-WT_Targets.csv  \
    --target-csv-column Target \
    --plot-dir _figures/diablo/W_NRAS-WT/plsda/  \
    --n-components 4  \
    --distance-metric mahalanobis.dist  \
    --keepx-csv _cache/diablo/W_NRAS-WT_plsda-keepX.csv  \
| tee _cache/diablo/W_NRAS-WT_final-plsda.log
