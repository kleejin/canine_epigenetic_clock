library(ggplot2)
library(wesanderson)

##################################################################
## Load results
##################################################################
load("glm.loo.atac_alpha0.5_20230331.rda")
load("glm.loo.meth_alpha0.5_20230331.rda")
load("glm.loo.comb_alpha0.5_20230331.rda")

tmp = unlist(lapply(glm.loo.meth, function(x) x$coef))
tmp2 = sort(table(tmp))
clockfeat.meth = tmp2[tmp2 > 1]
metafeat.meth = tmp2[!grepl("_", names(tmp2))]

tmp = unlist(lapply(glm.loo.atac, function(x) x$coef))
tmp2 = sort(table(tmp))
clockfeat.atac = tmp2[tmp2 > 1]
metafeat.atac = tmp2[!grepl("_", names(tmp2))]

tmp = unlist(lapply(glm.loo.comb, function(x) x$coef))
tmp2 = sort(table(tmp))
clockfeat.comb = tmp2[tmp2 > 1]
metafeat.comb = tmp2[!grepl("_", names(tmp2))]

## Summarize into DF
tmp1 = data.frame(metafeat.atac)
tmp1$type = "ATAC"
tmp2 = data.frame(metafeat.meth)
tmp2$type = "DNAm"
tmp3 = data.frame(metafeat.comb)
tmp3$type = "Combined"
meta.summary.df = rbind(tmp1, tmp2, tmp3)
names(meta.summary.df) = c("Cell type", "Freq", "Data type")

meta.summary.df$`Cell type` = as.character(meta.summary.df$`Cell type`)
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD62L..CD8.T.cells"] = "CD62L+ CD8 T cells"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD62L..CD8.T.cells.1"] = "CD62L- CD8 T cells"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD94..T.cells.1"] = "CD94- T cells"

meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD62L..DP.T.cells.1"] = "CD62L- DP T cells"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD62L..DN.T.cells.1"] = "CD62L- DN T cells"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "weight.cat"] = "Weight category"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "non.T.non.B"] = "non T non B"
meta.summary.df$`Cell type`[meta.summary.df$`Cell type` == "CD4.T.cells"] = "CD4 T cells"
meta.summary.df$`Data type` = factor(meta.summary.df$`Data type`, levels = c("ATAC", "DNAm", "Combined"))
meta.summary.df$`Cell type` = factor(meta.summary.df$`Cell type`)


## Add zeros for plotting
for(i in unique(meta.summary.df$`Cell type`)){
  
  print(i)
  
  tmp.df = data.frame(`Cell type` = i, Freq = 0, `Data type` = c("ATAC", "DNAm", "Combined"))
  names(tmp.df) = c("Cell type", "Freq", "Data type")
  tmp.summary = meta.summary.df[meta.summary.df$`Cell type` == i,]
  tmp.keep = setdiff(tmp.df$`Data type`, tmp.summary$`Data type`)
  tmp.df = tmp.df[tmp.df$`Data type` %in% tmp.keep,]
  
  if(nrow(tmp.df) >= 1){
    meta.summary.df = rbind(tmp.df, meta.summary.df)
  } else{
    next
  }
  
}

pal <- wes_palette("FantasticFox1", 3, type = "continuous")
names(pal) = c("ATAC", "DNAm", "Combined")

meta.summary.df$`Data type` = factor(meta.summary.df$`Data type`, levels = c("Combined", "DNAm", "ATAC"))
g = meta.summary.df %>%
  ggplot(aes(x = Freq, y = `Cell type`, fill = `Data type`)) +
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_manual(values = pal)+
  # facet_wrap(~`Data type`, nrow = 1)+
  xlab("Frequency selected") + ylab("Meta feature")+
  theme_bw()

ggsave(g, filename = "meta_feature_barplot.pdf", width = 6, height = 3)
