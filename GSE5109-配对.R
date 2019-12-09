rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
gse = "GSE5109"
eSet <- getGEO("GSE5109", 
               destdir = '.', 
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
exp = log2(exp+1)
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation


identical(rownames(pd),colnames(exp))

pd = pd[match(colnames(exp),rownames(pd)),]

group_list = ifelse(str_detect(pd$title,"Pre"),"pre","post")
group_list=factor(group_list,levels = c("pre","post"),ordered = T)
pairinfo = factor(c(1,2,1,3,2,3))
exp = log2(exp+1)

#PCA
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)#画主成分分析图需要加载这两个包
  library(factoextra) 
  # pca的统一操作走起
  dat.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave(paste0(gse,"PCA.png"))
}

exp[1:4,1:4]
boxplot(exp[1,]~group_list) 

#差异分析，用limma包来做
#需要表达矩阵和group_list，其他都不要动
{
  library(limma)
  design=model.matrix(~group_list+pairinfo)
  fit=lmFit(exp,design)
  fit=eBayes(fit)
  
  #差异基因排名
  deg=topTable(fit,coef=2,number = Inf)
  head(deg)
  
  #为deg数据框添加几列
  #1.加probe_id列，把行名变成一列
  library(dplyr)
  deg <- mutate(deg,probe_id=rownames(deg))
  #tibble::rownames_to_column()
  head(deg)
}


#2.加symbol列，火山图要用
#id转换，查找芯片平台对应的包
eSet[[1]]@annotation
#http://www.bio-info-trainee.com/1399.html
#hgu133plus2
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
{
  
  #merge
  deg <- inner_join(deg,ids,by="probe_id")
  deg <- deg[!duplicated(deg$symbol),]
  head(deg)
  #3.加change列：上调或下调，火山图要用
  
  logFC_t=0 #不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
  change=ifelse(deg$P.Value>0.05,'stable', 
                ifelse( deg$logFC >logFC_t,'up', 
                        ifelse( deg$logFC < -logFC_t,'down','stable') )
  )
  deg <- mutate(deg,change)
  head(deg)
  table(deg$change)
  deg <- mutate(deg,v = -log10(P.Value))
  
  #4.加ENTREZID列，后面富集分析要用
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  s2e <- bitr(unique(deg$symbol), fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Hs.eg.db)
  head(s2e)
  head(deg)
  deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
  
  head(deg)
  
  library(dplyr)
  dat <- mutate(deg,v=-log10(P.Value))
  head(dat)
}

{
  p <- ggplot(data = deg, 
              aes(x = logFC, 
                  y = v)) +
    geom_point(alpha=0.4, size=3.5, 
               aes(color=change)) +
    ylab("-log10(Pvalue)")+
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    theme_bw()
  for_label <- deg %>% 
    filter(abs(logFC) >0.7& P.Value< 0.05)
  p +
    geom_point(size = 3, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = symbol),
      data = for_label,
      color="black"
    )
  ggsave(paste0(gse,"volcano.png"))
}


x=deg$logFC 
names(x)=deg$probe_id 
cg=c(names(head(sort(x),30)),
     names(tail(sort(x),30)))
#热图实现了配对画图，pre在前，post在后
library(pheatmap)
test = data.frame(gsm = colnames(exp),group_list,pairinfo)
test
col = (arrange(test,pairinfo,group_list))$gsm
od = match(col,colnames(exp))
n=exp[cg,od]

annotation_col=data.frame(group= as.character(group_list)[od],
                          pair = as.character(pairinfo)[od])
rownames(annotation_col)=colnames(n) 

pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4)
         ) 
dev.off()
#保存

png(file = "heatmap.png")
pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4)) 
dev.off()

# 配对样本的箱线图

exp2 = exp[match(deg$probe_id,rownames(exp)),]
dim(exp2)
rownames(exp2) = deg$symbol
dat <- data.frame(pairinfo=pairinfo,group=group_list,t(exp2))
library(ggplot2)
colnames(dat)[3]
ggplot(dat, aes(group,STK26,fill=group)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of",colnames(dat)[3]))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(paste0(colnames(dat)[3],"box.png"))
### 文章报道的基因
c("GRB14","GPD1","GDF8") %in% rownames(exp2)

pn2 = exp2[c("GRB14","GPD1"),]

ggplot(dat, aes(group,GRB14,fill=group)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of","GRB14"))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(paste0("GRB14","box.png"))

ggplot(dat, aes(group,GPD1,fill=group)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of","GPD1"))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(paste0("GPD1","box.png"))

# 差异基因表格
test = deg[deg$P.Value<0.05,]
nrow(test)
range(test$logFC)
table(test$change)
genes = deg$symbol[deg$change != "stable"]

#用arrange间接实现了配对画图，pre在前，post在后
library(pheatmap)
test = data.frame(gsm = colnames(exp),group_list,pairinfo)
test

library(dplyr)
col = (arrange(test,pairinfo,group_list))$gsm
od = match(col,colnames(exp))

n=exp2[genes,od]

annotation_col=data.frame(group= as.character(group_list)[od],
                          pair = as.character(pairinfo)[od])
rownames(annotation_col)=colnames(n) 
png("26genes.png")
pheatmap(n,show_colnames =F,
         #show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4)) 
dev.off()

c("GRB14","GPD1","GDF8") %in% rownames(n)
cgs = rownames(n)
save(cgs,file = "cgs.Rdata")


#用arrange间接实现了配对画图，pre在前，post在后
library(pheatmap)
test = data.frame(gsm = colnames(exp),group_list,pairinfo)
test

library(dplyr)
col = (arrange(test,pairinfo,group_list))$gsm
od = match(col,colnames(exp))
load("cgin.Rdata")
n=exp2[cgin,od]

annotation_col=data.frame(group= as.character(group_list)[od],
                          pair = as.character(pairinfo)[od])
rownames(annotation_col)=colnames(n) 
png("474genes.png")
pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4)) 
dev.off()
