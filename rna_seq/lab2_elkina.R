library(DESeq2)
library(biomaRt)

sample_name <- 'E-GEOD-78936.sdrf.txt'
counts_name <- 'E-GEOD-78936-raw-counts.tsv'

# считаем и обработаем данные

smpls <- read.table(sample_name, sep = "\t", header = T)
counts <- read.table(counts_name, sep = '\t', header = T, row.names = 1)[, -1]

samples <- smpls

# возьмем только важные нам столбцы
smpls <- unique(smpls[, c("Scan.Name", "FactorValue..organism.part.", "FactorValue..disease.")])
colnames(smpls) <- c("Scan", "Brain_area", "Disease")
rownames(smpls) <- smpls$Scan
smpls$Scan <- NULL

# поскольку строки в smpls уже отсортированы, отсортируем столбцы в counts
# чтобы номер столбца был равен номеру строки в smpls

counts <- counts[, order(names(counts))]

# использование DESeq2
# тестируем простую линейную модель тестируемая часть мозга + заболевание
# то есть будем тестировать влияние на заболевание фактора части мозга
des_obj <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = smpls,
                                   design = ~Brain_area + Disease)

# построим для наглядности данные в pca (наглядность особо ничего не даст, жаль)
# blind = TRUE because of 
# comparing samples in an manner unbiased by prior information on samples
vstdata <- vst(des_obj, blind = T)
plotPCA(vstdata, intgroup = c("Brain_area", "Disease"))

# проведем сам анализ: он делается следующим образом
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

DES <- DESeq(des_obj)

# все возможные имена
resultsNames(DES)
# выберем биполярное расстройство и шизофрению
res <- results(DES, name = "Disease_normal_vs_bipolar.disorder", alpha = 0.05)
res <- res[complete.cases(res), ]
background <- res
res <- res[res[, "padj"] < 0.05, ]
res <- res[order(res$padj), ]

library(EnsDb.Hsapiens.v79)

gene_names <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                     keys = rownames(res),
                                     keytype = "GENEID", columns = c("GENEID","SYMBOL"))

##########################################################################################

gene_entrez <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                keys = rownames(res),
                                keytype = "GENEID", columns = c("GENEID", "ENTREZID","SYMBOL"))
                                