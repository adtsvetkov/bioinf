library(ArrayExpress)
library(oligo)
library(limma)
library(AnnotationDbi)
library(reactome.db)
setwd("C:/microarrayslab")

# ПУНКТ ПЕРВЫЙ - ЗАГРУЗКА ДАННЫХ

# загрузка данных. Сохраним их в файл с датой

# geod15765 <- getAE("E-GEOD-15765", type = "full")
# save(geod15765, file=geod_file)
load("geod15765.RData")
aeset <- ae2bioc(mageFiles = geod15765) # A-AFFY-37

# смотрим, какие факторы поставлены в соответствие образцам
# это те факторы, которые должны присутствовать в анализе
colnames1 <- colnames(phenoData(aeset))
fac <- colnames1[grep("Factor", colnames1)]
fac

# выводим данные с этими факторами
groups <- phenoData(aeset)[, fac]
pData(groups)

# итак единственный фактор отвечает за тип злокачественного новообразования

# ПУНКТ ВТОРОЙ - ПРОВЕРКА КАЧЕСТВА ДАННЫХ

# строим модель данных
dataPLM <- fitProbeLevelModel(aeset)
par(mar=c(8, 4.1, 4.1, 2.1))
NUSE(dataPLM, las = 2)
nuse <- NUSE(dataPLM, type="stats")

# по графику NUSE видим выбросы:
# GSM395733.CEL - 22
# GSM395727.CEL - 31
# GSM395736.CEL - 41
# GSM395714.CEL - 45
# GSM395731.CEL - 49
# GSM395732.CEL - 55
# GSM395737.CEL - 57
# GSM395734.CEL - 70
# GSM395721.CEL - 74
# GSM395724.CEL - 89

RLE(dataPLM, las = 2)
rle <- RLE(dataPLM, type = "stats")

# теперь для этих же данных сделаем MA-плоты
# уберем те файлы, в которых есть выбросы

extremes <- c(22, 31, 41, 45, 49, 55, 57, 70, 74, 89)

new_aeset <- aeset[, -extremes]

# нарисуем MA-плоты для нормированных и ненормированных данных

normed_new_aeset <- rma(new_aeset)

arrays <- c(16, 32, 48, 64, 80)
for (a in arrays)
{
  png(paste("Original, №", a, ".png", sep=""), width = 720)
  MAplot(new_aeset, which = a, main = "Original data")
  dev.off()
  png(paste("Normed, №", a, ".png", sep=""), width = 720)
  MAplot(normed_new_aeset, which = a, main = "Normalized")
  dev.off()
}

new_groups <- groups[-extremes, ]

# ПУНКТ ТРЕТИЙ - ДИФФЕРЕНЦИАЛЬНАЯ ЭКСПРЕССИЯ

# создаем матрицу дизайна

factorname <- factor(new_groups$Factor.Value..TISSUE.)
designmat <- model.matrix(~0 + factorname)


# cholangiocarcinoma -> cc
# combined hepatocellular carcinoma and cholangiocarcinoma -> hc_cc
# hepatocellular carcinoma -> hc

colnames(designmat) <- c("cc", "hc_cc", "hc")
designmat

print(is.fullrank(designmat))

contrmat <- makeContrasts(hc_cc - cc, hc_cc - hc, levels = designmat)
contrmat

print(is.fullrank(contrmat))

fitted <- lmFit(normed_new_aeset, designmat)
contrfit <- contrasts.fit(fitted, contrmat)
new_fitted <- eBayes(contrfit)

png("volcano1.png", height = 720)
volcanoplot(new_fitted, coef = "hc_cc - cc", highlight = 15)
dev.off()
png("volcano2.png", height = 720)
volcanoplot(new_fitted, coef = "hc_cc - hc", highlight = 15)
dev.off()

head(new_fitted$coefficients)

# по наибольшему p-value, использованному в статье,
# отрежем данные 

# для hc_cc - cc

p_thresh_1 = 0.05

res1 <- topTable(new_fitted, coef = "hc_cc - cc", number = nrow(normed_new_aeset))
filtered_res1 <- res1[res1$adj.P.Val < p_thresh_1, ]

hc_cc_cc <- c(nrow(filtered_res1), nrow(filtered_res1[filtered_res1$logFC <= 0, ]), 
              nrow(filtered_res1[filtered_res1$logFC > 0, ]))

p_thresh_2 = 0.001

res2 <- topTable(new_fitted, coef = "hc_cc - hc", number = nrow(normed_new_aeset))
filtered_res2 <- res2[res2$adj.P.Val < p_thresh_2, ]

hc_cc_hc <- c(nrow(filtered_res2), nrow(filtered_res2[filtered_res2$logFC <= 0, ]), 
              nrow(filtered_res2[filtered_res2$logFC > 0, ]))

resframe <- rbind(hc_cc_cc, hc_cc_hc)
colnames(resframe) <- c("total", "<=0", ">0")
resframe

# нарисуем тепловые карты найденных диф. экспрессированных генов

res1_rownames <- rownames(filtered_res1)
res2_rownames <- rownames(filtered_res2)

hmap1 <- exprs(normed_new_aeset[res1_rownames, ])
colnames(hmap1) <- unlist(pData(new_groups)[1])
colnames(hmap1)[colnames(hmap1) == 
                  "combined hepatocellular carcinoma and cholangiocarcinoma"] <-"combined"
png("heatmap1.png", width = 1080, height = 1080)
heatmap(hmap1, margins = c(10, 10))
dev.off()

hmap2 <- exprs(normed_new_aeset[res2_rownames, ])
colnames(hmap2) <- unlist(pData(new_groups)[1])
colnames(hmap2)[colnames(hmap2) == 
                  "combined hepatocellular carcinoma and cholangiocarcinoma"] <-"combined"
png("heatmap2.png", width = 1080, height = 1080)
heatmap(hmap2, margins = c(10, 10))
dev.off()

# ПУНКТ №4 - АНАЛИЗ ОБОГАЩЕННОСТИ

# сохраним данные для анализа в давиде

# надо сохранить rownames(filtered_res1) rownames(filtered_res2)

write(rownames(res1), file = 'background1.txt')
write(rownames(filtered_res1), file = "genes1.txt")
write(rownames(res2), file = 'background2.txt')
write(rownames(filtered_res2), file = "genes2.txt")

# загружаем файлы давида

res1_david <- read.table("ALL_res1.txt", sep = "\t", header = T)
res1_david$Genes <- NULL

# hc_cc - cc
# Для этого в статье указаны development/differentiation- or metastasis/adhesion-related functions

# Найдем результаты, связанные с этим

match1 <- c("develop", "differ", "metastas", "adh")

matched_res1 <- res1_david[grep(paste(match1, collapse = "|"), res1_david$Term), ]

library(AnnotationDbi)
library(reactome.db)

react_annotations1 <- select(reactome.db, 
                             gsub(":.*","", 
                                  res1_david[res1_david$Category == 
                                               "REACTOME_PATHWAY", ]$Term),
       keytype = "PATHID",
       c("PATHNAME"))
react_annotations1 <- react_annotations1[complete.cases(react_annotations1), ]
react_annotations1$PATHNAME <- tolower(react_annotations1$PATHNAME)

matched_res1_react <- react_annotations1[grep(paste(match1, collapse = "|"), 
                                              react_annotations1$PATHNAME),]
matched_res1_react

# ничего не подошло

# аналогично поступаем с вторым файлом

res2_david <- read.table("ALL_res2.txt", sep = "\t", header = T)
res2_david$Genes <- NULL

# hc_cc - hc
# Для этого в статье указаны metabolism- and immune-related functions

# Найдем результаты, связанные с этим

match2 <- c("metabolism", "immun")

matched_res2 <- res2_david[grep(paste(match2, collapse = "|"), res2_david$Term), ]

library(AnnotationDbi)
library(reactome.db)

react_annotations2 <- select(reactome.db, 
                             gsub(":.*","", 
                                  res2_david[res2_david$Category == 
                                               "REACTOME_PATHWAY", ]$Term),
                            keytype = "PATHID",
                            c("PATHNAME"))
react_annotations2 <- react_annotations2[complete.cases(react_annotations2), ]
react_annotations2$PATHNAME <- tolower(react_annotations2$PATHNAME)

matched_res2_react <- react_annotations2[grep(paste(match2, collapse = "|"), react_annotations2$PATHNAME),]
matched_res2_react

res2_david[res2_david$Term == paste(matched_res2_react$PATHID, matched_res2_react$PATHID, sep=":"), ]

# рассмотрим отдельно ген TP53 из статьи

TP53 <- read.csv("TP53.csv", header = T)
TP53_GO <- TP53$Accession

res1_david[substring(res1_david$Term, 1, 10) %in% TP53_GO, ]
res2_david[substring(res2_david$Term, 1, 10) %in% TP53_GO, ]