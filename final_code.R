eig.val <- get_eigenvalue(PCA_Overall)
eig.val

Conv <- slice(prot_cov_2, 22:31)

controlconv <- slice(prot_cov_2, 1:31)
sevmild <- slice(prot_cov_2, 32:68)

View(Conv)

prot_cov_2$group_col = NULL

View(Conv)

PCA_1 <- PCA(severe2, graph = FALSE)
PCA_2 <- PCA(mild, graph = FALSE)
PCA_3 <- PCA(HC, graph = FALSE)
PCA_4 <- PCA(Conv, graph = FALSE)
PCA_A <- PCA(controlconv, graph = FALSE)
PCA_B <- PCA(sevmild, graph = FALSE)

PCA_Severe_contribution <- PCA(severe3, graph = FALSE)
PCA_HCConv <- PCA(convhc, graph = FALSE)
PCA_Mildsev <- PCA(mildsev, graph = FALSE)

severe3 <- t(severe2)

fviz_contrib(PCA_Severe_contribution, choice = "ind", axes = 1:2)

fviz_pca_var(PCA_1, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

sev <- fviz_eig(PCA_1, addlabels = TRUE, ylim = c(0, 50))
ggpubr::ggpar(sev,
              title = "Principal Component Analysis",
              subtitle = "Severe Eigenvalues",
              ggtheme = theme_excel(), palette = "gdocs_pal"
)
fviz_eig(PCA_2, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_3, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_4, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_A, addlabels = TRUE, ylim = c(0, 50))
ggpubr::ggpar(fviz_eig(PCA_A, addlabels = TRUE, ylim = c(0, 50))
              ,
              title = "Principal Component Analysis",
              subtitle = "Mild/Severe Cos2",
              ggtheme = theme_excel(), palette = "economist"
)
fviz_eig(PCA_B, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_Overall, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_HCConv, addlabels = TRUE, ylim = c(0, 50))
fviz_eig(PCA_Mildsev, addlabels = TRUE, ylim = c(0, 50))

#library defining----

library(factoextra)
library(FactoMineR)
library(tidyverse)
library(BSDA)
library(corrplot)
library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(enrichR)
library(ggthemes)

#PCA Formatting --------------------------------------------------------------------------------------------------------------------------


convhc <- slice(prot_cov_2, 1:31)
mildsev <- slice(prot_cov_2, 32:69)
------------------------------------
#Misc Data cleaning----
Severe_var <- get_pca_var(PCA_1)
Severe_var
fviz_pca_var(PCA_1,col.var="black")


library("corrplot")
corrplot(Severe_var$contrib, is.corr=FALSE)

#variables contributions to PC-----
fviz_contrib(PCA_Overall, choice="var", axes = 1, top = 10)
fviz_contrib(PCA_Overall, choice="var", axes = 2, top = 10)
sevcontrib <- fviz_contrib(PCA_A, choice="var", axes = 1:3, top = 10)
ggpubr::ggpar(sevcontrib,
              title = "Principal Component Analysis",
              subtitle = "Convalescent/Control Cos2",
              ggtheme = theme_excel(), palette = "economist"
)
fviz_contrib(PCA_Overall, choice="var", axes = 4, top = 10)
#PCA3d experimentation----
library(pca3d)
data(PCA_1)
pca_sev <- prcomp(severe3[,-1], scale.=TRUE)
gr <- factor(severe2[,1])    
summary(gr)

pca3d(pca_sev, group=gr, bg = "white", axes.color = "black", show.labels=TRUE, show.shapes = FALSE)
snapshotPCA3d(file="severe_plot.png")
?pca3d
#Var defining-------------------------------
mild_var <- get_pca_var(PCA_2)
HC_Var <- get_pca_var(PCA_3)
Conv_var <- get_pca_var(PCA_4)
main_var <- get_pca_var(PCA_Overall)

View(severe2)

#Cos2 Graphs-----------------------------------------------------------------------------------------------------------------------------

corrplot(var$cos2, is.corr=FALSE)
sev_cos2 <- fviz_cos2(PCA_1, choice = "var", axes = 1:3, top = 10) #Cos2 graph for severe

ggpubr::ggpar(sev_cos2,
              title = "Principal Component Analysis",
              subtitle = "Severe Cos2",
              ggtheme = theme_excel(), palette = "economist"
)

mild_cos2 <- fviz_cos2(PCA_2, choice = "var", axes = 1:3) #Cos2 graph for mild
hc_cos2 <- fviz_cos2(PCA_3, choice = "var", axes = 1:3) #Cos2 graph for HC
conv_cos2 <- fviz_cos2(PCA_4, choice = "var", axes = 1:3) #Cos2 graph for convalescent
fviz_cos2(PCA_Overall, choice = "var", axes = 1:3, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_Overall, choice = "var", axes = 1:2, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_Overall, choice = "var", top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_HCConv, choice = "var", axes = 1:3, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_HCConv, choice = "var", axes = 1:2, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_HCConv, choice = "var", top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_Mildsev, choice = "var", axes = 1:3, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_Mildsev, choice = "var", axes = 1:2, top = 10) #Cos2 graph for convalescent
fviz_cos2(PCA_Mildsev, choice = "var", top = 10) #Cos2 graph for convalescent

fviz_cos2(PCA_1, choice = "var", axes = 1:3, top = 30) #Cos2 graph for severe
ggpubr::ggpar(fviz_cos2(PCA_B, choice = "var", axes = 1:3, top = 10)
,
              title = "Principal Component Analysis",
              subtitle = "Mild/Severe Cos2",
              ggtheme = theme_excel(), palette = "economist"
)
fviz_cos2(PCA_B, choice = "var", axes = 1:3, top = 10) #Cos2 graph for severe


#Biplot Graphs----------------------------------------------------------------------------------------------------------------------------

fviz_pca_biplot(PCA_1, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(PCA_2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(PCA_3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(PCA_4, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_pca_biplot(PCA_Overall, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

#Correlation Matrices---------------------------------------------------------------------------------------------------
corrplot(Severe_var$cos2, is.corr=FALSE)
corrplot(mild_var$cos2, is.corr=FALSE)
corrplot(HC_Var$cos2, is.corr=FALSE)
corrplot(Conv_var$cos2, is.corr=FALSE)
corrplot(main_var$cos2, is.corr=FALSE)

#PC_1 - Severe
#PC_2 - Mild
#PC_3 - HC
#Pc_4 - Convalescent


#Correlation Circle Graphs--------------------------------------------------------------------------------------------------------------------
fviz_pca_var(PCA_1, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

fviz_pca_var(PCA_2, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
)

fviz_pca_var(PCA_3, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

fviz_pca_var(PCA_4, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

fviz_pca_var(PCA_Overall, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)


#Violin Plots + T tests-----------------------------------------------------------------------------------------------------------------------------------
view(prot_cov_2_melted)

GZMH <- slice(prot_cov_2_melted, 2109:2176)

GZMHg <- ggplot(GZMH, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                         position=position_dodge(1)) + labs(title="gzmh",x="Protein", y = "Fold increase")
S_GZMH <- slice(GZMH, 58:68)
M_GZMH <- slice(GZMH, 32:57)
mean(S_GZMH$value)
mean(M_GZMH$value)
sd(S_GZMH$value)
sd(M_GZMH$value)

tsum.test(
  mean.x = 5.423583,
  s.x = 0.5584862,
  n.x = 11,
  mean.y = 5.570023,
  s.y = 0.8143145,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

GZMHg

CXCL10 <- slice(prot_cov_2_melted, 3741:3808)

CXCL10g <- ggplot(CXCL10, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                         position=position_dodge(1)) + labs(title="CXCL10",x="Protein", y = "Fold increase")

CXCL10g


CD28 <- slice(prot_cov_2_melted, 2789:2856)

CD28g <- ggplot(CD28, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                         position=position_dodge(1)) + labs(title="CD28",x="Protein", y = "Fold increase")

CD28g

summary(CD28t)
View(CD28)
CD28t <- t(CD28)

CD28$Protein = NULL

S_CD28 <- slice(CD28, 58:68)
M_CD28 <- slice(CD28, 32:57)

length(S_CD28)
length(M_CD28)

sd(M_CD28$value)

tsum.test(
  mean.x = 1.0587,
  s.x = 0.4517157,
  n.x = 11,
  mean.y = 0.9121,
  s.y = 0.3995486,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)



CD83 <- slice(prot_cov_2_melted, 6053:6120)

CD83g <- ggplot(CD83, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                         position=position_dodge(1)) + labs(title="CD83",x="Protein", y = "Fold increase")

CD83g

fviz_pca_biplot(PCA_Overall, 
                col.ind = prot_cov_2$group_col, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species") 


IL72 <- slice(prot_cov_2_melted, 613:680)

write.csv(IL72, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\IL72.csv")

IL7g <- ggplot(IL72, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                         position=position_dodge(1)) + labs(title="IL7",x="Protein", y = "Fold increase")

IL7g

S_IL7 <- slice(IL7, 58:68)
M_IL7 <- slice(IL7, 32:57)

mean(S_IL7$value)
mean(M_IL7$value)
sd(S_IL7$value)
sd(M_IL7$value)

tsum.test(
  mean.x = 3.717664,
  s.x = 0.7378655,
  n.x = 11,
  mean.y = 3.204504,
  s.y = 0.6693909,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)



IL10 <- slice(prot_cov_2_melted, 3877:3944)


IL10g <- ggplot(IL10, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="IL10",x="Protein", y = "Fold increase")

IL10g


CD83g

install.packages("BSDA")
library(BSDA)


CCL4 <- slice(prot_cov_2_melted, 2993:3060)

S_CCL4 <- slice(CCL4, 58:68)
M_CCL4 <- slice(CCL4, 32:57)

mean(S_CCL4$value)
mean(M_CCL4$value)
sd(S_CCL4$value)
sd(M_CCL4$value)

tsum.test(
  mean.x = 6.70577,
  s.x = 0.930827,
  n.x = 11,
  mean.y = 6.339841,
  s.y = 0.5397448,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

MCP_3_b <- slice(prot_cov_2_melted, 205:272)
S_MCP3 <- slice(MCP_3, 58:68)
M_MCP3 <- slice(MCP_3, 32:57)

write.csv(MCP_3_b, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\MCP_3_b.csv")

LTB1_2g <- ggplot(LTB1_2, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + labs(title="L TGF Beta",x="Protein", y = "Fold increase")
MCP_3g
mean(S_MCP3$value)
mean(M_MCP3$value)
sd(S_MCP3$value)
sd(M_MCP3$value)

tsum.test(
  mean.x = 4.961365,
  s.x = 1.468397,
  n.x = 11,
  mean.y = 3.371262,
  s.y = 1.246916,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)



LTB1_2 <- slice(prot_cov_2_melted, 2245:2312)
write.csv(MCP_3_b, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\MCP_3_c.csv")

S_LTB1 <- slice(LTB1, 58:68)
M_LTB1 <- slice(LTB1, 32:57)
CC_LTB1 <- slice(LTB1, 1:31)
MS_LTB1 <- slice(LTB1, 32:68)


LTB13g <- ggplot(LTB13, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + labs(title="L TGF Beta",x="Protein", y = "Fold increase")

IL73g <- ggplot(IL73, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + labs(title="IL7",x="Protein", y = "Fold increase")

MCP_3_cg <- ggplot(MCP_3_c_csv, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + labs(title="MCP 3",x="Protein", y = "Fold increase")

mean(S_LTB1$value)
sd(S_LTB1$value)
mean(M_LTB1$value)
sd(M_LTB1$value)


MCP_3_c_csv$...1 = NULL

mean(MS_LTB1$value)
sd(MS_LTB1$value)
mean(CC_LTB1$value)
sd(CC_LTB1$value)

tsum.test(
  mean.x = 8.783336,
  s.x = 0.58567,
  n.x = 37,
  mean.y = 8.001794,
  s.y = 0.6124713,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

LTB1g <- ggplot(LTB1, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="LTB1",x="Protein", y = "Fold increase")
LTB1g

PDGF <- slice(prot_cov_2_melted, 2585:2652)

S_PDGF <- slice(PDGF, 58:68)
M_PDGF <- slice(PDGF, 32:57)

mean(S_PDGF$value)
sd(S_PDGF$value)
mean(M_PDGF$value)
sd(M_PDGF$value)


tsum.test(
  mean.x = 9.744619,
  s.x = 0.5664372,
  n.x = 11,
  mean.y = 9.690443,
  s.y = 0.4475549,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

PDCD1 <- slice(prot_cov_2_melted, 2653:2720)

S_PDCD1 <- slice(PDCD1, 58:68)
M_PDCD1 <- slice(PDCD1, 32:57)

mean(S_PDCD1$value)
sd(S_PDCD1$value)
mean(M_PDCD1$value)
sd(M_PDCD1$value)

tsum.test(
  mean.x = 3.753248,
  s.x = 0.7745835,
  n.x = 11,
  mean.y = 3.636043,
  s.y = 0.4660622,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CD5 <- slice(prot_cov_2_melted, 4081:4148)

S_CD5 <- slice(CD5, 58:68)
M_CD5 <- slice(CD5, 32:57)

mean(S_CD5$value)
sd(S_CD5$value)
mean(M_CD5$value)
sd(M_CD5$value)

tsum.test(
  mean.x = 4.383215,
  s.x = 0.4436853,
  n.x = 11,
  mean.y = 4.491174,
  s.y = 0.3051814,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

TNFRSF21 <- slice(prot_cov_2_melted, 4489:4556)

S_TNFRSF21 <- slice(TNFRSF21, 58:68)
M_TNFRSF21 <- slice(TNFRSF21, 32:57)

mean(S_TNFRSF21$value)
sd(S_TNFRSF21$value)
mean(M_TNFRSF21$value)
sd(M_TNFRSF21$value)

tsum.test(
  mean.x = 6.827106,
  s.x = 0.2594822,
  n.x = 11,
  mean.y = 6.822098,
  s.y = 0.3106433,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

S_CXCL10 <- slice(CXCL10, 58:68)
M_CXCL10 <- slice(CXCL10, 32:57)

mean(S_CXCL10$value)
sd(S_CXCL10$value)
mean(M_CXCL10$value)
sd(M_CXCL10$value)

tsum.test(
  mean.x = 11.47284,
  s.x = 1.41107,
  n.x = 11,
  mean.y = 11.24849,
  s.y = 1.421795,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

S_MCP1 <- slice(MCP1, 58:68)
M_MCP1 <- slice(MCP1, 32:57)

mean(S_MCP1$value)
sd(S_MCP1$value)
mean(M_MCP1$value)
sd(M_MCP1$value)

tsum.test(
  mean.x = 11.37601,
  s.x = 1.004056,
  n.x = 11,
  mean.y = 10.92102,
  s.y = 0.7916269,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

MCP2 <- slice(prot_cov_2_melted, 2925:2992)

S_MCP2 <- slice(MCP2, 58:68)
M_MCP2 <- slice(MCP2, 32:57)

mean(S_MCP2$value)
sd(S_MCP2$value)
mean(M_MCP2$value)
sd(M_MCP2$value)

tsum.test(
  mean.x = 8.240035,
  s.x = 0.9965753,
  n.x = 11,
  mean.y = 8.077885,
  s.y = 1.096826,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

IL6 <- slice(prot_cov_2_melted, 749:816)
S_IL6 <- slice(IL6, 58:68)
M_IL6 <- slice(IL6, 32:57)

mean(S_IL6$value)
sd(S_IL6$value)
mean(M_IL6$value)
sd(M_IL6$value)

tsum.test(
  mean.x = 4.772535,
  s.x = 1.550149,
  n.x = 11,
  mean.y = 3.614998,
  s.y = 1.35072,
  n.y = 26,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)


PSBg <- ggplot(PSB, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="PDGF Subunit B",x="Protein", y = "Fold increase")

write_csv(APT1, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\APT1.csv")
write_csv(PSB, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\PSB.csv")



#Diseased vs Nondiseased testing----------------------

LTB1 <- slice(prot_cov_2_melted, 2245:2312)

CC_LTB1 <- slice(LTB1, 1:31)
MS_LTB1 <- slice(LTB1, 32:68)




mean(MS_LTB1$value)
sd(MS_LTB1$value)
mean(CC_LTB1$value)
sd(CC_LTB1$value)

tsum.test(
  mean.x = 8.783336,
  s.x = 0.58567,
  n.x = 37,
  mean.y = 8.001794,
  s.y = 0.6124713,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

EGF <- slice(prot_cov_2_melted, 477:544)


CC_EGF <- slice(EGF, 1:31)
MS_EGF <- slice(EGF, 32:68)


mean(MS_EGF$value)
sd(MS_EGF$value)
mean(CC_EGF$value)
sd(CC_EGF$value)

tsum.test(
  mean.x = 8.533139,
  s.x = 1.112316,
  n.x = 37,
  mean.y = 7.164375,
  s.y = 1.042188,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

APT1 <- slice(prot_cov_2_melted, 545:612)


CC_APT1 <- slice(APT1, 1:31)
MS_APT1 <- slice(APT1, 32:68)


mean(MS_APT1$value)
sd(MS_APT1$value)
mean(CC_APT1$value)
sd(CC_APT1$value)

tsum.test(
  mean.x = 8.80007,
  s.x = 0.7822165,
  n.x = 37,
  mean.y = 7.921524,
  s.y = 1.126028,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

PSB <- slice(prot_cov_2_melted, 2585:2652)


CC_PSB <- slice(PSB, 1:31)
MS_PSB <- slice(PSB, 32:68)


mean(MS_PSB$value)
sd(MS_PSB$value)
mean(CC_PSB$value)
sd(CC_PSB$value)

tsum.test(
  mean.x = 9.706549,
  s.x = 0.4783894,
  n.x = 37,
  mean.y = 8.848897,
  s.y = 0.8987264,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CC_TNFRSF21 <- slice(TNFRSF21, 1:31)
MS_TNFRSF21 <- slice(TNFRSF21, 32:68)

mean(MS_TNFRSF21$value)
sd(MS_TNFRSF21$value)
mean(CC_TNFRSF21$value)
sd(CC_TNFRSF21$value)

tsum.test(
  mean.x = 6.823587,
  s.x = 0.2927829,
  n.x = 37,
  mean.y = 6.987327,
  s.y = 0.2305873,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

TNFRSF21g <-  ggplot(TNFRSF21, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="TNFRSF21",x="Protein", y = "Fold increase")


CC_IL7 <- slice(IL7, 1:31)
MS_IL7 <- slice(IL7, 32:68)

mean(MS_IL7$value)
sd(MS_IL7$value)
mean(CC_IL7$value)
sd(CC_IL7$value)

tsum.test(
  mean.x = 3.357065,
  s.x = 0.7203789,
  n.x = 37,
  mean.y = 2.435057,
  s.y = 0.53701,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)


CC_CD83 <- slice(CD83, 1:31)
MS_CD83 <- slice(CD83, 32:68)

mean(MS_CD83$value)
sd(MS_CD83$value)
mean(CC_CD83$value)
sd(CC_CD83$value)

tsum.test(
  mean.x = 2.473592,
  s.x = 0.352235,
  n.x = 37,
  mean.y = 2.288846,
  s.y = 0.2655307,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

TNF9 <- slice(prot_cov_2_melted, 69:136)


CC_TNF9 <- slice(TNF9, 1:31)
MS_TNF9 <- slice(TNF9, 32:68)

mean(MS_TNF9$value)
sd(MS_TNF9$value)
mean(CC_TNF9$value)
sd(CC_TNF9$value)

tsum.test(
  mean.x = 5.323564,
  s.x = 0.4909854,
  n.x = 37,
  mean.y = 5.286691,
  s.y = 0.3638851,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CC_MCP2 <- slice(MCP2, 1:31)
MS_MCP2 <- slice(MCP2, 32:68)

mean(MS_MCP2$value)
sd(MS_MCP2$value)
mean(CC_MCP2$value)
sd(CC_MCP2$value)

tsum.test(
  mean.x = 8.126092,
  s.x = 1.056863,
  n.x = 37,
  mean.y = 6.906903,
  s.y = 0.5526626,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CXCL11 <- slice(prot_cov_2_melted, 1021:1088)

CC_CXCL11 <- slice(CXCL11, 1:31)
MS_CXCL11 <- slice(CXCL11, 32:68)

mean(MS_CXCL11$value)
sd(MS_CXCL11$value)
mean(CC_CXCL11$value)
sd(CC_CXCL11$value)

tsum.test(
  mean.x = 9.384652,
  s.x = 1.206476,
  n.x = 37,
  mean.y = 7.228065,
  s.y = 0.8281565,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

VEGFA <- slice(prot_cov_2_melted, 5441:5508)

CC_VEGFA <- slice(VEGFA, 1:31)
MS_VEGFA <- slice(VEGFA, 32:68)

mean(MS_VEGFA$value)
sd(MS_VEGFA$value)
mean(CC_VEGFA$value)
sd(CC_VEGFA$value)

tsum.test(
  mean.x = 7.97876,
  s.x = 0.5094833,
  n.x = 37,
  mean.y = 7.167766,
  s.y = 0.4220001,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CX3CL1 <- slice(prot_cov_2_melted, 3673:3740)

CC_CX3CL1 <- slice(CX3CL1, 1:31)
MS_CX3CL1 <- slice(CX3CL1, 32:68)

mean(MS_CX3CL1$value)
sd(MS_CX3CL1$value)
mean(CC_CX3CL1$value)
sd(CC_CX3CL1$value)

tsum.test(
  mean.x = 6.300394,
  s.x = 0.5522002,
  n.x = 37,
  mean.y = 5.723699,
  s.y = 0.3438834,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CX3CL1 <- slice(prot_cov_2_melted, 3673:3740)

CC_CD5 <- slice(CD5, 1:31)
MS_CD5 <- slice(CD5, 32:68)

mean(MS_CD5$value)
sd(MS_CD5$value)
mean(CC_CD5$value)
sd(CC_CD5$value)

tsum.test(
  mean.x = 4.459078,
  s.x = 0.3490881,
  n.x = 37,
  mean.y = 4.504135,
  s.y = 0.3523925,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CC_MCP1 <- slice(MCP1, 1:31)
MS_MCP1 <- slice(MCP1, 32:68)

mean(MS_MCP1$value)
sd(MS_MCP1$value)
mean(CC_MCP1$value)
sd(CC_MCP1$value)

tsum.test(
  mean.x = 11.05629,
  s.x = 0.8715924,
  n.x = 37,
  mean.y = 10.14185,
  s.y = 0.4797713,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CC_MCP3 <- slice(MCP_3, 1:31)
MS_MCP3 <- slice(MCP_3, 32:68)

mean(MS_MCP3$value)
sd(MS_MCP3$value)
mean(CC_MCP3$value)
sd(CC_MCP3$value)

tsum.test(
  mean.x = 3.843995,
  s.x = 1.490487,
  n.x = 37,
  mean.y = 0.8827113,
  s.y = 0.4622532,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

MCP2g <- TNFRSF21g <-  ggplot(MCP2, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="MCP2g",x="Protein", y = "Fold increase")

MCP1g <- TNFRSF21g <-  ggplot(MCP1, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="MCP1g",x="Protein", y = "Fold increase")


PDL1 <- slice(prot_cov_2_melted, 3197:3264)

PDL1g <- TNFRSF21g <-  ggplot(PDL1, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="PDL1",x="Protein", y = "Fold increase")

CC_PDL1 <- slice(PDL1, 1:31)
MS_PDL1 <- slice(PDL1, 32:68)

mean(MS_PDL1$value)
sd(MS_PDL1$value)
mean(CC_PDL1$value)
sd(CC_PDL1$value)

tsum.test(
  mean.x = 5.260294,
  s.x = 0.5121187,
  n.x = 37,
  mean.y = 4.141884,
  s.y = 0.3007477,
  n.y = 31,
  alternative = "two.sided",
  mu = 0,
  var.equal = FALSE,
  conf.level = 0.95
)

CD40 <- slice(prot_cov_2_melted, 273:340)

CD40g <- ggplot(CD40, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="CD40",x="Protein", y = "Fold increase")

GAL9 <- slice(prot_cov_2_melted, 1837:1904)

GAL9g <- ggplot(GAL9, aes(x=Protein, y=value, fill=group_col)) +
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1)) + labs(title="GAL9",x="Protein", y = "Fold increase")

#Graphs_2------
LTB1g
write.csv(LTB1, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\LTB1.csv")
EGFg
write.csv(EGF, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\EGF.csv")
IL7g1
write.csv(IL7, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\IL7.csv")
CD83g
write.csv(CD83, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\CD83.csv")
TNF9g
write.csv(TNF9, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\TNF9.csv")
MCP_3g
write.csv(MCP_3, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\MCP_3.csv")
PDL1g
write.csv(PDL1, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\PDL1.csv")
CD40g
write.csv(CD40, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\CD40.csv")
GAL9g
write.csv(GAL9, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\GAL9.csv")
IL6g1
write.csv(IL6, "C:\\Users\\netra\\Documents\\~\\Desktop\\r-novice-inflammation\\IL6.csv")
?reorder

