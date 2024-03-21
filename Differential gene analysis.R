library(stringr)
library(edgeR)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(cols4all)
library(GSEABase)
library(enrichplot)
library(stats)
library(forcats)
library(ggstance)

rawcount <- read.table("data.txt",row.names = 1,  sep = "\t", header = T)
colnames(rawcount)


dim(rawcount)

group <- read.table("data/group.txt",
                    header = T,sep = "\t",quote = "\"")

group <- group[match(colnames(rawcount),group$sample), c("sample","group")]
group

group_list <- group$group
group_list

keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)

filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

express_cpm <- log2(cpm(filter_count)+ 1)
express_cpm[1:6,1:6]


save(filter_count,express_cpm,group,
     file = "data/Step01-airwayData.Rdata")



lname <- load(file = "data/Step01-airwayData.Rdata")
lname

exprSet <- filter_count
dim(exprSet)
exprSet[1:4,1:4]

group_list <- group[match(colnames(filter_count),group$sample),2]
group_list

comp <- unlist(strsplit("day1_VS_day6",split = "_VS_"))
group_list <- factor(group_list)
group_list <- relevel(group_list,ref = "day6")

group_list
table(group_list)

design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(exprSet)
colnames(design) <- levels(factor(group_list))
design

DEG <- DGEList(counts=exprSet, 
               group=factor(group_list))

DEG <- calcNormFactors(DEG)

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)

lrt <- glmLRT(fit, contrast=c(1,-1)) 

DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

fc_cutoff <- 2
pvalue <- 0.05

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),
                    which(DEG_edgeR$PValue<pvalue))
loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)

id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENSEMBL", 
                  toType = "SYMBOL", 
                  OrgDb = org.Hs.eg.db)
head(id2symbol)

DEG_edgeR <- cbind(GeneID=rownames(DEG_edgeR),DEG_edgeR)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all.y=T)
head(DEG_edgeR_symbol)

DEG_edgeR_symbol_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")

write.table(DEG_edgeR_symbol,"results/DEG_edgeR_all.xls", row.names = F, sep="\t",quote = F)
write.table(DEG_edgeR_symbol_Sig,"results/DEG_edgeR_Sig.xls", row.names = F, sep="\t",quote = F)
save(DEG_edgeR_symbol,DEG_edgeR_symbol_Sig,file = "data/Step03-edgeR_nrDEG.Rdata")


load(file = "data/Step03-edgeR_nrDEG.Rdata")
ls()

DEG <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",2]
DEG <- as.character(na.omit(DEG))
head(DEG)

DEG_up <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated=="up",2]
DEG_up <- as.character(na.omit(DEG_up))
head(DEG_up)
DEG_down <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated=="down",2]
DEG_down <- as.character(na.omit(DEG_down))
head(DEG_down)


ego_BP_up <- enrichGO(gene=DEG_up, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 0.05,qvalueCutoff= 0.05)

ego_BP_down <- enrichGO(gene=DEG_down, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 0.05,qvalueCutoff= 0.05)


ego_BP_up_dt <- as.data.frame(ego_BP_up)
ego_BP_down_dt <- as.data.frame(ego_BP_down)
p <- barplot(ego_BP_up )
p

p111 <- barplot(ego_BP_down )
p111


down_ego_BP <- ego_BP_down_dt[ego_BP_down_dt$p.adjust<0.05,]

up_ego_BP<-ego_BP_up_dt[ego_BP_up_dt$p.adjust<0.05,]

down_ego_BP$change<- "down"
up_ego_BP$change<- "up"
dt <-rbind(up_ego_BP,down_ego_BP) 
dt$'-log10p.adjust' <- ifelse(dt$change == 'up', -log10(dt$p.adjust), -(-log10(dt$p.adjust)))
dt <- dt[order(-dt$`-log10p.adjust`), ]
level_up <- dt$Description[dt$change=='up']
level_down <- dt$Description[dt$change=='down']
level <- c(level_up,level_down)

dt$Description <- factor(dt$Description,
                         levels = rev(level))

dt$abs_log10p_adjust <- abs(dt$`-log10p.adjust`)


dt <- dt[order(dt$abs_log10p_adjust, decreasing = TRUE), ]


upregulated <- dt[dt$`-log10p.adjust` > 0, ]
top10_upregulated <- upregulated

downregulated <- dt[dt$`-log10p.adjust` < 0, ]
top10_downregulated <- downregulated[,]

combined <- rbind(top10_upregulated, top10_downregulated)

combined$regulation <- ifelse(combined$`-log10p.adjust` > 0, "Upregulated", "Downregulated")

p <- ggplot(combined,aes(x = `-log10p.adjust`,y = Description,fill = `-log10p.adjust`)) + #数据映射
  geom_col() + 
  theme_classic()
p


p1<-p + scale_fill_continuous_c4a_div('sunset', mid = 0) 
p1


p2 <- ggplot(combined,aes(x = `-log10p.adjust`,y = Description)) +
  geom_col(width = 0.1, 
           fill = 'black') +
  theme_classic()
p2
p3 <- p2 + geom_point(aes(size = abs(`-log10p.adjust`)),
                      color = 'black') +
  scale_size_continuous(range = c(3, 10)) 
p4 <- ggplot(combined,aes(x = `-log10p.adjust`,y = Description)) +
  geom_col(aes(fill = `-log10p.adjust`), width = 0.1) +
  geom_point(aes(size = Count,
                 color = `-log10p.adjust`)) +
  scale_size_continuous(range = c(2, 7)) +
  scale_color_continuous_c4a_div('div_bu_wh_rd', mid = 0, reverse = T) +
  scale_fill_continuous_c4a_div('div_bu_wh_rd', mid = 0, reverse = T) +
  theme_classic()
p4

p5 <- p4 + scale_x_continuous(breaks = seq(-5,5,by = 2.5),
                              labels = abs(seq(-5,5,by = 2.5))) + 
  ylab('') +
  theme(
    axis.text = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)  )
p5


load("data/Step03-edgeR_nrDEG.Rdata")


DEG <- DEG_edgeR_symbol

DEG <- DEG[!is.na(DEG$SYMBOL),]

DEG <- DEG[!duplicated(DEG$SYMBOL),]
geneList <- DEG$logFC
names(geneList) <- DEG$SYMBOL
geneList <- sort(geneList,decreasing = T)

geneset1 <- read.gmt("data/MsigDB/h.all.v2023.1.Hs.symbols.gmt")
geneset <- rbind(geneset1)

egmt <- GSEA(geneList, TERM2GENE=geneset, pvalueCutoff = 0.05, verbose = T)
x <- data.frame(egmt@result)

x <- x[x$pvalue<0.05 & x$p.adjust<0.25 & abs(x$NES)>1,]

y <- arrange(x, desc(abs(NES))) %>% group_by(sign(NES))

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue)) + 
  geom_barh(stat='identity') + 
  scale_fill_continuous(low='red', high='blue') + 
  theme_bw() + ylab(NULL) + 
  theme(axis.text.y = element_text(size=8, face="bold"))