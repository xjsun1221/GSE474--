#1.getdata----
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
library(AnnoProbe)
gse = "GSE5109"

if(require(AnnoProbe)){
  if(!file.exists(paste0(gse,"_eSet.Rdata"))) geoChina(gse)
  load(paste0(gse,"_eSet.Rdata"))
  eSet <- gset
  rm(gset)
}else{
  eSet <- getGEO(gse, 
                 destdir = '.', 
                 getGPL = F)
}

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
exp = log2(exp+1)
#(2)提取临床信息
pd <- pData(eSet[[1]])
#(3)提取芯片平台编号
gpl <- eSet[[1]]@annotation

p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#2.group----
group_list = ifelse(str_detect(pd$title,"Pre"),"pre","post")
group_list=factor(group_list,levels = c("pre","post"),ordered = T)
pairinfo = factor(c(1,2,1,3,2,3))

#3.PCA----
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)#画主成分分析图需要加载这两个包
  library(factoextra) 
  # pca的统一操作走起
  dat.pca <- PCA(dat, graph = FALSE)
  pca_plot <- fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  print(pca_plot)
}
ggsave(plot = pca_plot,filename = paste0(gse,"PCA.png"))
#4.差异分析----
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

#5.芯片注释----
eSet[[1]]@annotation
#http://www.bio-info-trainee.com/1399.html
#hgu133plus2
if(F){
  if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db")
  ids <- toTable(hgu133plus2SYMBOL)
  head(ids)
}else if(F){
  getGEO(gpl)
  ids = data.table::fread(paste0(gpl,".soft"),header = T,skip = "ID",data.table = F)
  ids = ids[c("ID","Gene Symbol"),]
  colnames(ids) = c("probe_id")
}else if(T){
  ids = idmap(gpl,type = "bioc")
}
#6.加列----
{
  #merge
  deg <- inner_join(deg,ids,by="probe_id")
  deg <- deg[!duplicated(deg$symbol),]
  head(deg)
  #3.加change列：上调或下调，火山图要用
  logFC_t=mean(deg$logFC)+2*sd(deg$logFC)
  #logFC_t=2
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
}
#7.火山图----
library(dplyr)
dat <- mutate(deg,v=-log10(P.Value))
head(dat)
{
  for_label <- deg%>% filter(abs(logFC)>logFC_t & P.Value < 0.05) %>%head(10)
  p <- ggplot(data = dat, 
              aes(x = logFC, 
                  y = v)) +
    geom_point(alpha=0.4, size=3.5, 
               aes(color=change)) +
    ylab("-log10(Pvalue)")+
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    theme_bw()
volcano_plot <- p +
    geom_point(size = 3, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = symbol),
      data = for_label,
      color="black"
    )
}
ggsave(plot = volcano_plot,filename = paste0(gse,"volcano.png"))
#8.配对样本的差异基因热图----
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
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =F,
                                   show_rownames = F,
                                   scale = "row",
                                   cluster_cols = F, 
                                   annotation_col=annotation_col,
                                   gaps_col = c(2,4)
))
dev.off()
#保存

png(file = paste0(gse,"heatmap.png"))
pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4)) 
dev.off()

#9.配对样本的箱线图----
expm = exp[match(deg$probe_id,rownames(exp)),]
rownames(expm) = deg$symbol
#配对样本箱线图批量绘制,画10张玩玩
n = 10
library(ggplot2)
x = deg$symbol[order(head(deg$logFC,n))]
dat <- data.frame(pairinfo=pairinfo,group=group_list,t(expm[x,]))
pl = list()
for(i in 1:n){
  pl[[i]] = ggplot(dat, aes_string("group",colnames(dat)[i+2],fill="group")) +
    geom_boxplot() +
    geom_point(size=2, alpha=0.5) +
    geom_line(aes(group=pairinfo), colour="black", linetype="11") +
    xlab("") +
    ylab(paste("Expression of",colnames(dat)[i+2]))+
    theme_classic()+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#00AFBB", "#E7B800"))
}
#拼图
library(patchwork)
pb_top10= wrap_plots(pl,nrow = 2)+plot_annotation(tag_levels = 'A')
pb_top10
ggsave(plot = pb_top10,filename = paste0(gse,"box.png"),width = 15,height = 6)

(pca_plot + volcano_plot +heatmap_plot)/ pb_top10 + plot_annotation(tag_levels = "A")

table(deg$change)

