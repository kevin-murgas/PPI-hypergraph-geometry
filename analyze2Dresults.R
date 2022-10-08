# analysis of 2D ppi network model curvature results

library(tidyverse)
library(Matrix)
library(ggpubr)
select <- dplyr::select

# specify data folder and experiment subfolder
setwd("Analysis")
#setwd(your_data_directory)
datadir = getwd()


### FIGURE 3A,B: GSE75748_3 hESC cell type
expdir = "GSE75748_3"

# load in csv results from matlab curvature algorithm
adj = read_csv(file.path(datadir, expdir, "adj.csv"))
#expr = read_csv(file.path(datadir, expdir, "expr.csv")) %>% select(-Row)
piDist = read_csv(file.path(datadir, expdir, "piDist.csv")) %>% select(-Row)
nodeEntropy = read_csv(file.path(datadir, expdir, "nodeEntropy.csv")) %>% select(-Row)
curv1D_in = read_csv(file.path(datadir, expdir, "1Dcurv_in.csv")) %>% select(-Row)
curv1D_out = read_csv(file.path(datadir, expdir, "1Dcurv_out.csv")) %>% select(-Row)
curv2D_in = read_csv(file.path(datadir, expdir, "2Dcurv_in.csv")) %>% select(-Row)
curv2D_out = read_csv(file.path(datadir, expdir, "2Dcurv_out.csv")) %>% select(-Row)

# get sample (column) and gene (row) names
sampName = colnames(piDist)
geneName = adj$Row
adj = select(adj,-Row) %>% as.matrix %>% Matrix(sparse=T)
deg_out = rowSums(adj)
deg_in = colSums(adj)

# compute summary measures (global averages)
curv1D = curv1D_in/deg_in - curv1D_out/deg_out
curv2D = curv2D_in/deg_in - curv2D_out/deg_out
Fglob1 = colSums( piDist * curv1D_in/deg_out )
Fglob2 = colSums( piDist * curv2D_in/deg_out )
Sglob  = colSums( piDist * nodeEntropy )

# maximum possible entropy
# power method to determine first eigenvalue of adjacency matrix
b0 = deg_out*0+1
b1 = b0 %>% as.matrix
adj = as.matrix(adj)
for (i in 1:10) {
  b1 = adj %*% b1  # run this line 10 times or so to iterate power method until stable
}
lam = sum( (adj %*% b1) * b1) / sum(b1 * b1) # compute first eigenvalue
#a = eigen(adj, only.values = T); lam = a$values[1] # SLOW
Smax = log(lam) # from power method

# labels for 75748_3
sampGroup = str_extract(sampName,"[^_]+")
sampGroup[sampGroup %in% c("H1","H9")] = "hESC"
sampGroup[sampGroup %in% c("DEC")] = "DEP"
sampGroup = factor(sampGroup, levels= c("hESC","NPC","DEP","TB","HFF","EC"))
# violins
plot_df = data.frame(i = sampName, g = sampGroup,
                     s=Sglob/Smax, f1 = Fglob1, f2 = Fglob2) %>% arrange(as.numeric(g))
p1 = plot_df %>% ggplot(aes(x=g,y=s)) + geom_violin(aes(fill=g)) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.99) +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p2 = plot_df %>% ggplot(aes(x=g,y=f1)) + geom_violin(aes(fill=g)) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.99) +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p3 = plot_df %>% ggplot(aes(x=g,y=f2)) + geom_violin(aes(fill=g)) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.99) +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p4 = plot_df %>% ggplot(aes(x=f1,y=s,color=g)) + geom_point() +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p5 = plot_df %>% ggplot(aes(x=f2,y=s,color=g)) + geom_point() +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p6 = plot_df %>% ggplot(aes(x=f1,y=f2,color=g)) + geom_point() +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "none")
p_all = ggarrange(p1,p4,p2,p5,p3,p6,nrow=3, ncol=2)
p_all


#### FIGURE 3C: GSE75748_4 differentiation time course
expdir = "GSE75748_4"

# load in csv results from matlab curvature algorithm
adj = read_csv(file.path(datadir, expdir, "adj.csv"))
#expr = read_csv(file.path(datadir, expdir, "expr.csv")) %>% select(-Row)
piDist = read_csv(file.path(datadir, expdir, "piDist.csv")) %>% select(-Row)
nodeEntropy = read_csv(file.path(datadir, expdir, "nodeEntropy.csv")) %>% select(-Row)
curv1D_in = read_csv(file.path(datadir, expdir, "1Dcurv_in.csv")) %>% select(-Row)
curv1D_out = read_csv(file.path(datadir, expdir, "1Dcurv_out.csv")) %>% select(-Row)
curv2D_in = read_csv(file.path(datadir, expdir, "2Dcurv_in.csv")) %>% select(-Row)
curv2D_out = read_csv(file.path(datadir, expdir, "2Dcurv_out.csv")) %>% select(-Row)

# get sample (column) and gene (row) names
sampName = colnames(piDist)
geneName = adj$Row
adj = select(adj,-Row) %>% as.matrix %>% Matrix(sparse=T)
deg_out = rowSums(adj)
deg_in = colSums(adj)

# compute summary measures (global averages)
curv1D = curv1D_in/deg_in - curv1D_out/deg_out
curv2D = curv2D_in/deg_in - curv2D_out/deg_out
Fglob1 = colSums( piDist * curv1D_in/deg_out )
Fglob2 = colSums( piDist * curv2D_in/deg_out )
Sglob  = colSums( piDist * nodeEntropy )

# maximum possible entropy
# power method to determine first eigenvalue of adjacency matrix
b0 = deg_out*0+1
b1 = b0 %>% as.matrix
adj = as.matrix(adj)
for (i in 1:10) {
  b1 = adj %*% b1  # run this line 10 times or so to iterate power method until stable
}
lam = sum( (adj %*% b1) * b1) / sum(b1 * b1) # compute first eigenvalue
#a = eigen(adj, only.values = T); lam = a$values[1] # SLOW
Smax = log(lam) # from power method

# labels for 75748_4
sampGroup = str_extract(sampName,"(?<=_)[^_]+(?=h)")
sampType = str_extract(sampName,"[^_]+")

# plot means and errorbars
plot_df = data.frame(i = sampName, g = as.numeric(sampGroup),
                     s=Sglob/Smax, f1 = Fglob1, f2 = Fglob2)
mean_df = plot_df %>% group_by(g) %>%
  summarize(s_m = mean(s), f1_m = mean(f1), f2_m = mean(f2),
            s_s = sd(s), f1_s = sd(f1), f2_s = sd(f2))

p1 = ggplot() + geom_violin(data=plot_df, aes(x=g,y=s,fill=as.factor(g))) +
  geom_point(data=mean_df, aes(x=g,y=s_m),shape=21,fill="white", size=3) +
  geom_line(data=mean_df, aes(x=g,y=s_m),linetype="dotted",size=0.5) +
  theme_classic() + theme(legend.position = "none", aspect.ratio = 1) +
  scale_x_continuous(name="Time (h)", breaks = c(0,12,24,36,72,96))
  #geom_errorbar(data=mean_df, aes(x=g,y=s_m,ymin=s_m-s_s,ymax=s_m+s_s), width=5, size=1)
p2 = ggplot() + geom_violin(data=plot_df, aes(x=g,y=f1,fill=as.factor(g))) +
  geom_point(data=mean_df, aes(x=g,y=f1_m),shape=21,fill="white", size=3) +
  geom_line(data=mean_df, aes(x=g,y=f1_m),linetype="dotted",size=0.5) +
  theme_classic() + theme(legend.position = "none", aspect.ratio = 1) +
  scale_x_continuous(name="Time (h)", breaks = c(0,12,24,36,72,96))
p3 = ggplot() + geom_violin(data=plot_df, aes(x=g,y=f2,fill=as.factor(g))) +
  geom_point(data=mean_df, aes(x=g,y=f2_m),shape=21,fill="white", size=3) +
  geom_line(data=mean_df, aes(x=g,y=f2_m),linetype="dotted",size=0.5) +
  theme_classic() + theme(legend.position = "none", aspect.ratio = 1) +
  scale_x_continuous(name="Time (h)", breaks = c(0,12,24,36,72,96))
p_all = ggarrange(p1,p2,p3,nrow=3,ncol=1)
p_all


### FIGURE 4A,B: GSE72056 Melanoma results
expdir = "GSE72056"

# load in csv results from matlab curvature algorithm
adj = read_csv(file.path(datadir, expdir, "adj.csv"))
expr = read_csv(file.path(datadir, expdir, "expr.csv")) %>% select(-Row)
piDist = read_csv(file.path(datadir, expdir, "piDist.csv")) %>% select(-Row)
nodeEntropy = read_csv(file.path(datadir, expdir, "nodeEntropy.csv")) %>% select(-Row)
curv1D_in = read_csv(file.path(datadir, expdir, "1Dcurv_in.csv")) %>% select(-Row)
curv1D_out = read_csv(file.path(datadir, expdir, "1Dcurv_out.csv")) %>% select(-Row)
curv2D_in = read_csv(file.path(datadir, expdir, "2Dcurv_in.csv")) %>% select(-Row)
curv2D_out = read_csv(file.path(datadir, expdir, "2Dcurv_out.csv")) %>% select(-Row)

# get sample (column) and gene (row) names
sampName = colnames(piDist)
geneName = adj$Row
adj = select(adj,-Row) %>% as.matrix %>% Matrix(sparse=T)
deg_out = rowSums(adj)
deg_in = colSums(adj)

# compute summary measures (global averages)
curv1D = curv1D_in/deg_in - curv1D_out/deg_out
curv2D = curv2D_in/deg_in - curv2D_out/deg_out
Fglob1 = colSums( piDist * curv1D_in/deg_out )
Fglob2 = colSums( piDist * curv2D_in/deg_out )
Sglob  = colSums( piDist * nodeEntropy )

# maximum possible entropy
# power method to determine first eigenvalue of adjacency matrix
b0 = deg_out*0+1
b1 = b0 %>% as.matrix
adj = as.matrix(adj)
for (i in 1:10) {
  b1 = adj %*% b1  # run this line 10 times or so to iterate power method until stable
}
lam = sum( (adj %*% b1) * b1) / sum(b1 * b1) # compute first eigenvalue
#a = eigen(adj, only.values = T); lam = a$values[1] # SLOW
Smax = log(lam) # from power method

sampGroup = str_extract(sampName,"(?<=x)[^_]+") # patient label
sampType = str_extract(sampName,"(?<=_)[^_]+") # cell type: 1=normal, 2=tumor
samp_df = data.frame(i = sampName, g = sampGroup, t = sampType) # number of samples in each group/type
n_normal = samp_df %>% filter(t==1) %>% group_by(g) %>% summarize(nn = n())
n_tumor = samp_df %>% filter(t==2) %>% group_by(g) %>% summarize(nt = n())
# (72056) keep only groups with minimum 5 normal and tumor cells
# keep: 53 59 60 71 79 80 81 82 84 88 89 94, exclude: 58 65 67 72 74 75 78
keep_groups = full_join(n_normal,n_tumor) %>% filter(nn>5 & nt>5) %>% pull(g)
plot_df = data.frame(i = sampName, g = sampGroup, t= sampType,
                     s=Sglob/Smax, f1 = Fglob1, f2 = Fglob2) %>% filter( g %in% keep_groups )

# patient-split boxplots (Supplementary Figure S1)
plot_df %>% ggplot(aes(x=g,y=s)) + geom_boxplot(aes(fill=t), outlier.size = 0) + theme_classic() + theme(aspect.ratio = 1/4)
plot_df %>% ggplot(aes(x=g,y=f1)) + geom_boxplot(aes(fill=t), outlier.size = 0) + theme_classic() + theme(aspect.ratio = 1/4)
plot_df %>% ggplot(aes(x=g,y=f2)) + geom_boxplot(aes(fill=t), outlier.size = 0) + theme_classic() + theme(aspect.ratio = 1/4)
# all cells violins
plot_df %>% ggplot(aes(x=t, y=s)) + geom_violin(aes(fill=t)) +
  geom_boxplot( width=0.1, outlier.size = 0) + theme_classic() + theme(aspect.ratio = 2)
plot_df %>% ggplot(aes(x=t, y=f1)) + geom_violin(aes(fill=t)) +
  geom_boxplot( width=0.1, outlier.size = 0) + theme_classic() + theme(aspect.ratio = 2)
plot_df %>% ggplot(aes(x=t, y=f2)) + geom_violin(aes(fill=t)) +
  geom_boxplot( width=0.1, outlier.size = 0) + theme_classic() + theme(aspect.ratio = 2)

# stats: compare normal vs tumor for each measure
shapiro.test(plot_df %>% filter(t==1) %>% pull(s)) # data is not normally distributed
wilcox.test(plot_df %>% filter(t==1) %>% pull(s), plot_df %>% filter(t==2) %>% pull(s))
wilcox.test(plot_df %>% filter(t==1) %>% pull(f1), plot_df %>% filter(t==2) %>% pull(f1))
wilcox.test(plot_df %>% filter(t==1) %>% pull(f2), plot_df %>% filter(t==2) %>% pull(f2))

# Supplementary Figure S2: compute ROC AUC for each measure
# considering normal = negative, tumor = positive, to compute TPR and FPR
library(pROC)
ROC_all_s = roc(plot_df$t ~ plot_df$s)
ROC_all_f1 = roc(plot_df$t ~ plot_df$f1)
ROC_all_f2 = roc(plot_df$t ~ plot_df$f2)

ROC_all = plot(ROC_all_s, col="blue", print.auc=T, print.auc.y = .4, print.auc.x = .3)
ROC_all = plot(ROC_all_f1, col="red", print.auc=T, print.auc.y = .3, print.auc.x = .3, add=T)
ROC_all = plot(ROC_all_f2, col="green", print.auc=T, print.auc.y = .5, print.auc.x = .3, add=T)
text(.3,c(.3,.4,.5)-0.024,c("FR_1D","SR","FR_2D"), col = c("red","blue","green"), pos=2)

# compare AUC
roc.test(response = plot_df$t, predictor1 = plot_df$f2, predictor2 = plot_df$s)
roc.test(response = plot_df$t, predictor1 = plot_df$f2, predictor2 = plot_df$f1)
roc.test(response = plot_df$t, predictor1 = plot_df$f1, predictor2 = plot_df$s)

# differential expression and curvature for each gene
ind_N = which(sampType == 1); ind_T = which(sampType == 2)
diff_df = data.frame(gene = geneName,
                     ex_tum = rowMeans(expr[,ind_T]), ex_norm = rowMeans(expr[,ind_N]),
                     m1_tum = rowMeans(curv1D[,ind_T]), m1_norm = rowMeans(curv1D[,ind_N]),
                     m2_tum = rowMeans(curv2D[,ind_T]), m2_norm = rowMeans(curv2D[,ind_N])) %>%
  mutate(diffEx = ex_tum-ex_norm, diff1D = m1_tum - m1_norm, diff2D = m2_tum - m2_norm,
         pvalE = 0, pval1 = 0, pval2=0)
# test for significance with nonparametric wilcoxon t-test
diff_df$pvalE = expr %>% apply(1,function(x) wilcox.test(x[ind_T],x[ind_N])$p.value)
diff_df$pval1 = curv1D %>% apply(1,function(x) wilcox.test(x[ind_T],x[ind_N])$p.value)
diff_df$pval2 = curv2D %>% apply(1,function(x) wilcox.test(x[ind_T],x[ind_N])$p.value)
diff_df$padjE = p.adjust(diff_df$pvalE, "BH")
diff_df$padj1 = p.adjust(diff_df$pval1, "BH")
diff_df$padj2 = p.adjust(diff_df$pval2, "BH")
diff_df = diff_df %>% mutate(ratEx = diffEx / ex_norm,
                             rat1D = diff1D / m1_norm,
                             rat2D = diff2D / m2_norm)

# Figure 5: 
diff_df %>% arrange(padj2,-abs(diff2D)) %>% filter(abs(rat2D)>2)
indE = which(diff_df$padjE < 0.05 & abs(diff_df$ratEx)>2)
ind1 = which(diff_df$padj1 < 0.05 & abs(diff_df$rat1D)>2)
ind2 = which(diff_df$padj2 < 0.05 & abs(diff_df$rat2D)>2)
set = list(ind1,indE,  ind2)
names(set) <- c("FR_1D","DE","FR_2D")
library(RVenn)
ggvenn(Venn(set), slice = 1:3) + theme_void()

# Table 2/3: Increased/Decreased PPI Curvature Pathways in Melanoma
# the significant genes are fed into Reactome online pathway analysis tool
diff_df = diff_df %>% mutate(sig_up_ex = (padjE<0.05 & diffEx>0 & abs(ratEx)>2),
                             sig_dn_ex = (padjE<0.05 & diffEx<0 & abs(ratEx)>2),
                             sig_up_1D = (padj1<0.05 & diff1D>0 & abs(rat1D)>2),
                             sig_dn_1D = (padj1<0.05 & diff1D<0 & abs(rat1D)>2),
                             sig_up_2D = (padj2<0.05 & diff2D>0 & abs(rat2D)>2),
                             sig_dn_2D = (padj2<0.05 & diff2D<0 & abs(rat2D)>2))
write_csv(diff_df, "gse72056_siggenes.csv")


### FIGURE 4C: GSE81861 colorectal cancer
# GSE81861 colorectal cancer heterogeneity from 11 patient tumors and matching normal mucosa
expdir = "GSE81861"

# load in csv results from matlab curvature algorithm
adj = read_csv(file.path(datadir, expdir, "adj.csv"))
#expr = read_csv(file.path(datadir, expdir, "expr.csv")) %>% select(-Row)
piDist = read_csv(file.path(datadir, expdir, "piDist.csv")) %>% select(-Row)
#nodeEntropy = read_csv(file.path(datadir, expdir, "nodeEntropy.csv")) %>% select(-Row)
curv1D_in = read_csv(file.path(datadir, expdir, "1Dcurv_in.csv")) %>% select(-Row)
curv1D_out = read_csv(file.path(datadir, expdir, "1Dcurv_out.csv")) %>% select(-Row)
curv2D_in = read_csv(file.path(datadir, expdir, "2Dcurv_in.csv")) %>% select(-Row)
curv2D_out = read_csv(file.path(datadir, expdir, "2Dcurv_out.csv")) %>% select(-Row)

# get sample (column) and gene (row) names
sampName = colnames(piDist)
geneName = adj$Row
adj = select(adj,-Row) %>% as.matrix %>% Matrix(sparse=T)
deg_out = rowSums(adj)
deg_in = colSums(adj)

# compute summary measures (global averages)
curv1D = curv1D_in/deg_in - curv1D_out/deg_out
curv2D = curv2D_in/deg_in - curv2D_out/deg_out
Fglob1 = colSums( piDist * curv1D_in/deg_out )
Fglob2 = colSums( piDist * curv2D_in/deg_out )

# labels for 81861
sampType = str_extract(sampName,"(?<=_)[^_]+") # 7 cell types + NA (8)
sampGroup = str_extract(sampName,"[^_]+$") # normal or tumor

samp_df = data.frame(t = sampType, g = sampGroup)
samp_df %>% group_by(g,t) %>% summarize(n=n())
samp_df %>% group_by(t) %>% summarize(nn=sum(g=="normal"), nt=sum(g=="tumor"))

keep_types = samp_df %>% group_by(t) %>%
  summarize(nn=sum(g=="normal"), nt=sum(g=="tumor")) %>%
  filter(nn>5 & nt>5) %>% pull(t)

# filter to keep only samples with n_normal > 5 and n_tumor > 5
plot_df = data.frame(i = sampName, g = sampGroup, t= sampType,
                     f1 = Fglob1, f2 = Fglob2)
plot_df = plot_df %>% dplyr::filter( t %in% keep_types ) %>%
  mutate(t2 = if_else(t=="Epithelial", "epi", "not epi"))

p2 = plot_df %>% ggplot() + geom_boxplot(aes(x=t2, fill=g, y=f1)) +
  theme_classic() + theme(aspect.ratio = 1/3, legend.position = "none") +
  scale_y_continuous(breaks=seq(-1800,-600,400))
p3 = plot_df %>% ggplot() + geom_boxplot(aes(x=t2, fill=g, y=f2)) +
  theme_classic() + theme(aspect.ratio = 1/3, legend.position = "none")
p_all = ggarrange(p2,p3,ncol=1,nrow=2)
p_all

# statistical test comparisons
kruskal.test(f2 ~ t2, data=plot_df %>% filter(g=="normal"))  # n.s.
kruskal.test(f2 ~ t2, data=plot_df %>% filter(g=="tumor"))   # p<0.001 ***
kruskal.test(f2 ~ g, data=plot_df %>% filter(t2=="epi"))     # p<0.001 ***
kruskal.test(f2 ~ g, data=plot_df %>% filter(t2=="not epi")) # n.s.
x = plot_df$f2
g = paste(plot_df$t2, plot_df$g, sep="_")
pairwise.wilcox.test(x,g,p.adjust.method="BH")

# ANOVA on cell type (in tumor cells is significant)
lm.f2 <- lm(f2 ~ t, data = plot_df %>% filter(g=="tumor")) # one-way for cell type in tumor only
lm.f2 <- lm(f2 ~ t2 * g , data = plot_df) # two-way for cell type x normal/tumor
av.f2 <- aov(lm.f2)
summary(av.f2)
tukey.test <- TukeyHSD(av.f2)
tukey.test
plot(tukey.test)


##### Fig 4D: GSE130019 results
# Ewing cell line (A673/TR/shEF) treated with DOX, time-course
expdir = "GSE130019"

# load in csv results from matlab curvature algorithm
adj = read_csv(file.path(datadir, expdir, "adj.csv"))
#expr = read_csv(file.path(datadir, expdir, "expr.csv")) %>% select(-Row)
piDist = read_csv(file.path(datadir, expdir, "piDist.csv")) %>% select(-Row)
#nodeEntropy = read_csv(file.path(datadir, expdir, "nodeEntropy.csv")) %>% select(-Row)
curv1D_in = read_csv(file.path(datadir, expdir, "1Dcurv_in.csv")) %>% select(-Row)
curv1D_out = read_csv(file.path(datadir, expdir, "1Dcurv_out.csv")) %>% select(-Row)
curv2D_in = read_csv(file.path(datadir, expdir, "2Dcurv_in.csv")) %>% select(-Row)
curv2D_out = read_csv(file.path(datadir, expdir, "2Dcurv_out.csv")) %>% select(-Row)

# get sample (column) and gene (row) names
sampName = colnames(piDist)
geneName = adj$Row
adj = select(adj,-Row) %>% as.matrix %>% Matrix(sparse=T)
deg_out = rowSums(adj)
deg_in = colSums(adj)

# compute summary measures (global averages)
curv1D = curv1D_in/deg_in - curv1D_out/deg_out
curv2D = curv2D_in/deg_in - curv2D_out/deg_out
Fglob1 = colSums( piDist * curv1D_in/deg_out )
Fglob2 = colSums( piDist * curv2D_in/deg_out )

sampLabels = read_delim(file.path(datadir, expdir, "GSE130019_sampLabels.txt"))
sampType = sampLabels[1,-c(1,2)] %>% as.character()
sampType = factor(sampType, levels = unique(sampType))

plot_df = data.frame(i = sampName, t= sampType,
                     f1 = Fglob1, f2 = Fglob2) %>%
  mutate(timept = str_extract(t,"(?<=\\+)[^+]+$")) %>%
  mutate(timept2 = if_else(is.na(timept),"0",timept) %>% as.numeric)
mean_df = plot_df %>% group_by(timept2) %>%
  summarize(f1_m = mean(f1), f2_m = mean(f2),
            f1_s = sd(f1), f2_s = sd(f2))

p1 = ggplot() + geom_violin(data=plot_df, aes(x=timept2, y=f1, fill=as.factor(timept2) )) +
  geom_point(data=mean_df, aes(x=timept2,y=f1_m),shape=21,fill="white", size=3) +
  geom_line(data=mean_df, aes(x=timept2,y=f1_m),linetype="dotted",size=0.5) +
  theme_classic() + theme(legend.position = "none", aspect.ratio = 1/3) +
  scale_x_continuous(name="Time (d)", breaks = c(0,2,3,4,7,10,15,22))
p2 = ggplot() + geom_violin(data=plot_df, aes(x=timept2, y=f2, fill=as.factor(timept2) )) +
  geom_point(data=mean_df, aes(x=timept2,y=f2_m),shape=21,fill="white", size=3) +
  geom_line(data=mean_df, aes(x=timept2,y=f2_m),linetype="dotted",size=0.5) +
  theme_classic() + theme(legend.position = "none", aspect.ratio = 1/3) +
  scale_x_continuous(name="Time (d)", breaks = c(0,2,3,4,7,10,15,22))
p_all = ggarrange(p1,p2,ncol=1,nrow=2)
p_all

# test with ANOVA
lm.f2 = lm(f2 ~ timept2, data = plot_df %>% filter(!str_detect(t,"_Xeno")))
av.f2 = aov(lm.f2)
summary(lm.f2)
summary(av.f2)
lm.f1 = lm(f1 ~ timept2, data = plot_df %>% filter(!str_detect(t,"_Xeno")))
summary(lm.f1)
