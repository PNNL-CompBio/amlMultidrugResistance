library(Biobase)
library(limma)

meta <- read.csv("data/patient_meta.csv", row.names = 1)
meta = meta[meta$Race %in% c("White", "Black"),]
meta[!meta$FLT3_ITD %in% c("WT", "Mutant"), "FLT3_ITD"] = NA

FLT3_ITD_status <- paste(meta$Race, meta$FLT3_ITD, sep=".")
FLT3_ITD_status <- factor(
  FLT3_ITD_status, 
  levels=c(
    "Black.Mutant", "Black.WT",
    "White.Mutant", "White.WT"
  )
)
design <- model.matrix(
  ~0 + FLT3_ITD_status,
  data = meta
)
colnames(design) <- levels(FLT3_ITD_status)

phospho <- read.csv("data/phospho.csv", row.names = 1)
phospho <- phospho[row.names(design),]
phospho <- t(as.matrix(phospho))
p_eset <- ExpressionSet(phospho)

global <- read.csv("data/global.csv", row.names = 1)
global <- global[row.names(design),]
global <- t(as.matrix(global))
g_eset <- ExpressionSet(global)

cont.matrix <- makeContrasts(
  Black = Black.Mutant - Black.WT,
  White = White.Mutant - White.WT, 
  Diff = (Black.Mutant - Black.WT) - (White.Mutant - White.WT),
  levels = design
)

p_fit <- lmFit(p_eset, design)
g_fit <- lmFit(g_eset, design)

p_contrasts <- contrasts.fit(p_fit, cont.matrix)
p_contrasts <- eBayes(p_contrasts)
p_results <- topTable(p_contrasts, number = Inf, adjust.method = "BH")
p_tests <- decideTests(p_contrasts)

g_contrasts <- contrasts.fit(g_fit, cont.matrix)
g_contrasts <- eBayes(g_contrasts)
g_results <- topTable(g_contrasts, number = Inf, adjust.method = "BH")
g_tests <- decideTests(g_contrasts)

both_over <- intersect(
  row.names(g_tests[g_tests[,"White"] == 1,]), 
  row.names(g_tests[g_tests[,"Black"] == 1,])
)
both_under <- intersect(
  row.names(g_tests[g_tests[,"White"] == -1,]), 
  row.names(g_tests[g_tests[,"Black"] == -1,])
)

black_over <- intersect(
  row.names(g_tests[g_tests[,"White"] != 1,]), 
  row.names(g_tests[g_tests[,"Black"] == 1,])
)
black_under <- intersect(
  row.names(g_tests[g_tests[,"White"] != -1,]), 
  row.names(g_tests[g_tests[,"Black"] == -1,])
)
white_over <- intersect(
  row.names(g_tests[g_tests[,"White"] == 1,]), 
  row.names(g_tests[g_tests[,"Black"] != 1,])
)
white_under <- intersect(
  row.names(g_tests[g_tests[,"White"] == -1,]), 
  row.names(g_tests[g_tests[,"Black"] != -1,])
)
