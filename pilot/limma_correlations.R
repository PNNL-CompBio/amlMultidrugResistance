library(Biobase)
library(limma)

meta <- read.csv("data/patient_meta.csv", row.names = 1)
design <- model.matrix(
  ~0 + Race + Study + Age,
  data = meta
)
colnames(design)[
  which(colnames(design) == "RacePacific Islander")
] = "RacePacificIslander"

phospho <- read.csv("data/phospho.csv", row.names = 1)
phospho <- phospho[row.names(design),]
phospho <- t(as.matrix(phospho))
p_eset <- ExpressionSet(phospho)

global <- read.csv("data/global.csv", row.names = 1)
global <- global[row.names(design),]
global <- t(as.matrix(global))
g_eset <- ExpressionSet(global)

cont.matrix <- makeContrasts(RaceBlack - RaceWhite, levels = design)
p_fit <- lmFit(p_eset, design)
g_fit <- lmFit(g_eset, design)

p_contrasts <- contrasts.fit(p_fit, cont.matrix)
p_contrasts <- eBayes(p_contrasts)
p_results <- topTable(p_contrasts, number = Inf, adjust.method = "BH")

g_contrasts <- contrasts.fit(g_fit, cont.matrix)
g_contrasts <- eBayes(g_contrasts)
g_results <- topTable(g_contrasts, number = Inf, adjust.method = "BH")

print(dim(p_results[p_results$adj.P.Val < 0.05 & p_results$logFC < 0,]))
print(dim(p_results[p_results$adj.P.Val < 0.05 & p_results$logFC > 0,]))
print(dim(g_results[g_results$adj.P.Val < 0.05 & g_results$logFC < 0,]))
print(dim(g_results[g_results$adj.P.Val < 0.05 & g_results$logFC > 0,]))