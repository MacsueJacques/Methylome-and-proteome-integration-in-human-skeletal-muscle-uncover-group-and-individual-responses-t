#Analyses for paper on DNA methylation and proteomics using repeted testing#
####Load data sets and perform random effects mixed models####
library(tidyverse)

#Because of variability when imputing make sure to use the same data always for analyses
dir=c("D:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Proteomics and methylation changes at the individual level/New Analyses 06092021/Data used in the analysis")
setwd(dir)
data_final=read_csv("MetaData_norm_and_batch_3.csv")
data_final=as.data.frame(data_final)
rownames(data_final)=data_final$...1
data_final=dplyr::select(data_final, SI_CON_1:SI_12WP_17)
data_final=t(data_final)

library(readxl)
Pheno_Second_intervention=read_excel("Gene SMART SG100-SG195.xlsx",sheet="Second intervention",col_names=TRUE,
                                     col_types=c("text","text","date","date","date","date","date","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                                                 "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"),skip=1)

#Subset variables of interest
library(tidyverse)
SI_subset=Pheno_Second_intervention%>%
  dplyr::select(ID, "New ID", Age_PRE:Age_12WP, Kcals_inter, Wpeak_rel_PRE, Wpeak_rel_4WP, Wpeak_rel_8WP, Wpeak_rel_12WP, 
                LT_poly_rel_PRE, LT_poly_rel_4WP, LT_poly_rel_8WP, LT_poly_rel_12WP,
                VO2max_rel_PRE, VO2max_rel_4WP, VO2max_rel_8WP, VO2max_rel_12WP, 
                "Type I fibre %_PRE":"SOD3_12WP")

#change data to a longer format
SI_subset_long <- pivot_longer(data = SI_subset,
                               cols = -c("ID","New ID", Kcals_inter),
                               names_to = c(".value","Timepoint"), #indicates that component of the name
                               #defines the name of the column containing the cell values, overriding values_to.
                               names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",
                               values_drop_na = 'FALSE')


#add timepoint as a number so models will work further on
SI_subset_long=SI_subset_long%>%
  mutate(time=rep(c(1,2,3,4), times=20))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))


#Transform data so you can fit loop
data_irs_log3=data_final[9:76,]

#Load experiment design
exp_design<-read.table("exp_design_SI.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

exp_design=exp_design[c(-33,-73) ,1:5] #Remove reference channels and bad samples

exp_design$condition<-factor(exp_design$condition, 
                             levels = c("SI_CON", "SI_PRE", "SI_4WP", "SI_8WP", "SI_12WP"))

exp_des3=drop_na(exp_design)
IDs=exp_des3[c(9:76),2]
IDs=c("SI3_PRE","SI4_PRE","SI5_PRE","SI6_PRE","SI8_PRE","SI9_PRE","SI10_PRE","SI11_PRE","SI12_PRE","SI14_PRE","SI16_PRE",
      "SI14_PRE","SI18_PRE","SI19_PRE","SI20_PRE","SI20_PRE","SI13_4WP","SI4_4WP","SI4_4WP","SI5_4WP","SI6_4WP","SI9_4WP",
      "SI10_4WP","SI11_4WP","SI12_4WP","SI3_4WP","SI14_4WP","SI15_4WP","SI16_4WP","SI17_4WP","SI18_4WP","SI19_4WP","SI19_4WP",
      "SI20_4WP","SI10_8WP","SI11_8WP","SI12_8WP","SI13_8WP","SI14_8WP","SI15_8WP","SI16_8WP","SI17_8WP","SI18_8WP","SI19_8WP",
      "SI20_8WP","SI13_8WP","SI13_8WP","SI4_8WP","SI5_8WP","SI6_8WP","SI8_8WP","SI9_8WP","SI3_12WP","SI10_12WP","SI11_12WP",
      "SI12_12WP","SI14_12WP","SI15_12WP","SI17_12WP","SI18_12WP","SI18_12WP","SI19_12WP","SI20_12WP","SI4_12WP","SI6_12WP",
      "SI8_12WP","SI9_12WP","SI13_12WP")

Timepoint=exp_des3[c(9:76),3]
data_irs_log3=data.frame(data_irs_log3, Timepoint, IDs)
data_irs_log3=as_tibble(data_irs_log3)
#remove duplicates
data_irs_log3=data_irs_log3[c(-12,-16,-18,-32,-46,-47,-61),]

SI_subset_long2=unite(SI_subset_long, "IDs", "New ID",Timepoint, sep = "_")

data_final=inner_join(data_irs_log3, SI_subset_long2, by="IDs")

cols=data_final[,1:2365]

data_final

data_final=rename(data_final, Type_I_perc=`Type I fibre %`)

#Do imputations for missing VO2 values
library(mice)

#Use SI_zscore data subset
md.pattern(data_final)

#For mice to work, we need to ensure the variables we want to impute are
#in the right format (numeric, factor, etc.)
data_final

#examine the number of missing cases in a table format and proportion of missing
missing.indicator=data.frame(is.na(data_final))
propMissing=apply(missing.indicator, 2, mean)

#To impute the missing data, use the mice() function
pred_mat=quickpred(data_final, mincor = 0.25) #quickly selects data that has min corr of 0.25
imputed <- mice(data = data_final, #dataset containing missing values
                m =  5, #number of imputations you want (default = 5)
                #method = NULL. It tells the algorithm HOW it will impute the data, depending on data type (continuous variable, factor with 1 or more levels)
                print=F,
                pred=pred_mat
)

#Have a look at the summary of imputations
imputed

#To access imputed values:
imputed$imp
plot(imputed)

data_final=complete(imputed, 3)
rownames(data_final)=data_final$IDs

data_final[,c("IDs", "VO2max_rel")]

VO2max_PRE=c(45.43, 45.45,69.61, 46.95,46.97,58.24,48.90,50.74,65.57,58.53,37.32,58.62,39.85,
             63.39,59.99,45.45,69.60,46.95,58.24,48.90,50.74,65.57,45.43,58.53,64.26,37.32,
             37.30,58.62,39.85,63.39,48.90,50.74,65.57,59.99,58.54,64.26,37.32,37.30,58.62,
             39.84,63.39,45.45,69.60,46.95,46.97,58.24,45.43,48.90,50.74,65.57,58.53,64.26,58.62,
             39.85,63.39,45.45,46.95,46.97,58.24,59.99)
data_final=cbind(data_final,VO2max_PRE)

#load lmertest package
library(lmerTest)

#Create a matrix to save results
Proteins=c(colnames(data_final[,1:2365]))

ResultsVO2max=matrix(NA, nrow = 2365, ncol=3)
colnames(ResultsVO2max)=c("Estimate", "SE", "p-value")
rownames(ResultsVO2max)=Proteins

ResultsTime=matrix(NA, nrow = 2365, ncol=3)
colnames(ResultsTime)=c("Estimate", "SE", "p-value")
rownames(ResultsTime)=Proteins

ResultsRandom=matrix(NA, nrow = 2365, ncol=3)
colnames(ResultsRandom)=c("Residual Variance", "Residual SE", "p-value")
rownames(ResultsRandom)=Proteins

cols=c(colnames(data_final[,1:2365]))
#data_final=data_final[-53,]

#loop by cols 
for (c in cols){
  model=lmer(data_final[,c]~data_final[,'VO2max_PRE']+
               data_final[,'time']+
               data_final[,'Age']+
               (1+data_final[,'time']|data_final[,'ID']))
  
  ResultsVO2max[c,"Estimate"]=summary(model)$coefficients[2,"Estimate"]
  ResultsTime[c,"Estimate"]=summary(model)$coefficients[3,"Estimate"]
   
  ResultsVO2max[c,"SE"]=summary(model)$coefficients[2,"Std. Error"]
  ResultsTime[c,"SE"]=summary(model)$coefficients[3,"Std. Error"]
  
  ResultsVO2max[c,"p-value"]=summary(model)$coefficients[2,"Pr(>|t|)"]
  ResultsTime[c,"p-value"]=summary(model)$coefficients[3,"Pr(>|t|)"]
   
  ResultsRandom[c,"Residual Variance"]=as.data.frame(VarCorr(model))[4,4]
  ResultsRandom[c,"Residual SE"]=as.data.frame(VarCorr(model))[4,5]
  ResultsRandom[c,"p-value"]=as.data.frame(rand(model))[2,6]
  
  print(summary(model))
  
}

#create vectors for p-values to adjust
p_VO2max=c(ResultsVO2max[,3])
p_Time=c(ResultsTime[,3])
p_Random=c(ResultsRandom[,3])

p_adj_VO2max=p.adjust(p_VO2max, method = "fdr")
p_adj_Time=p.adjust(p_Time, method = "fdr")
p_adj_Random=p.adjust(p_Random, method = "fdr")

#calculate and add p-values
ResultsTimeSig=cbind(ResultsTime,p_adj_Time)

ResultsRandomSig=cbind(ResultsRandom,p_adj_Random)

ResultsVO2maxSig=cbind(ResultsVO2max,p_adj_VO2max)

setwd=dir
#write.csv(ResultsRandomSig, "ResultsRandom_log_Proteomics2109.csv")
#write.csv(ResultsTimeSig, "ResultsTime_log_Proteomics2109.csv")
#write.csv(ResultsVO2maxSig, "ResultsVO2max_log_Proteomics2109.csv")

ResultsRandomSig=read.csv("ResultsRandom_log_Proteomics2109.csv")
ResultsTimeSig=read.csv("ResultsTime_log_Proteomics2109.csv")
ResultsVO2maxSig=read.csv("ResultsVO2max_log_Proteomics2109.csv")


setwd("D:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses")
beta <- read.delim("Combined_First_Second_B.txt")
M <- read.delim("Combined_First_Second_M.txt")
library(readxl)

#Subset variables of interest
library(tidyverse)
library(limma)

#Subset M matrix for analyses
M_subset=M%>%
  dplyr::select(contains("SI"), -"SI2_8WP", -contains("CON"))
M_subset=dplyr::rename(M_subset, SI18_4W=SI18_4W1)

####Run loop with random effects to look at trainability####
library("lmerTest")
library("tidyverse")
library("broom")

#Transform data so you can fit loop
M_subset_transposed=t(M_subset)

CpGs=c(colnames(M_subset_transposed))

sapply(unlist(M_subset_transposed), as.numeric)

M_subset_transposed=as_tibble(M_subset_transposed)

Timepoint=c(SI_subset_long$Timepoint)

SI_subset_long2=SI_subset_long%>%
  unite(Names, "New ID", Timepoint)

SI_subset_long2=rename(SI_subset_long2, Type_I_perc=`Type I fibre %`)

#Do imputations for missing VO2 values
library(mice)

#Use SI_zscore data subset
md.pattern(SI_subset_long2)

#For mice to work, we need to ensure the variables we want to impute are
#in the right format (numeric, factor, etc.)
SI_subset_long2

#examine the number of missing cases in a table format and proportion of missing
missing.indicator=data.frame(is.na(SI_subset_long2))
propMissing=apply(missing.indicator, 2, mean)

#To impute the missing data, use the mice() function
pred_mat=quickpred(SI_subset_long2, mincor = 0.25) #quickly selects data that has min corr of 0.25
imputed <- mice(data = SI_subset_long2, #dataset containing missing values
                m =  5, #number of imputations you want (default = 5)
                #method = NULL. It tells the algorithm HOW it will impute the data, depending on data type (continuous variable, factor with 1 or more levels)
                print=F,
                pred=pred_mat
)

#Have a look at the summary of imputations
imputed

#To access imputed values:
imputed$imp
plot(imputed)

SI_subset_long2=complete(imputed, 3)


colnames(M_subset)
Names=c("SI20_PRE", "SI17_4WP", "SI14_8WP", "SI18_12WP", "SI19_4WP","SI8_8WP","SI15_PRE", "SI15_12WP","SI17_PRE", "SI11_4WP", "SI9_8WP",  
        "SI18_4WP", "SI12_8WP", "SI15_4WP", "SI17_12WP", "SI14_12WP", "SI4_8WP",  "SI19_PRE","SI19_12WP", "SI16_PRE",  "SI6_8WP", "SI18_PRE", 
        "SI13_4WP", "SI16_8WP", "SI14_PRE", "SI11_12WP", "SI12_4WP", "SI13_12WP", "SI20_8WP", "SI11_PRE", "SI14_4WP", "SI13_PRE", "SI12_PRE", 
        "SI17_8WP", "SI10_12WP", "SI9_PRE", "SI13_8WP", "SI12_12WP", "SI10_PRE","SI9_12WP",  "SI8_12WP", "SI3_PRE",  "SI10_4WP",  "SI8_4WP",   
        "SI4_PRE", "SI20_4WP", "SI6_PRE", "SI4_12WP", "SI9_4WP",  "SI18_8WP", "SI4_4WP",  "SI19_8WP", "SI3_12WP", "SI5_PRE", "SI11_8WP", 
        "SI3_4WP",  "SI8_PRE","SI5_4WP", "SI15_8WP",  "SI6_4WP",  "SI6_12WP",  "SI5_8WP", "SI16_4WP",  "SI20_12WP", "SI10_8WP")
M_subset_transposed=cbind(M_subset_transposed, Names)
SI_subset_long2=unite(SI_subset_long2, Names, sep="_", New_ID,Timepoint)
SI_subset_long2=filter(SI_subset_long2, Names %in% Names)

Pheno_M700k=left_join(SI_subset_long2, M_subset_transposed, by="Names")

VO2max_PRE_Meth=c(37.30,37.30,37.30,37.30,45.45,45.45,45.45,45.45,58.62,58.62,58.62,58.62,
                  65.57,65.57,65.57,65.57, 59.99,59.99,59.99,59.99,58.53,58.53,58.53,58.53,
                  46.97,46.97,46.97,46.97,46.52,46.52,46.52,46.52,50.74,50.74,50.74,50.74,
                  69.60,69.60,69.60,69.60,46.95,46.95,46.95,46.95,33.89,33.89,33.89,33.89,
                  45.43,45.43,45.43,45.43,37.32,37.32,37.32,37.32,64.26,64.26,64.26,64.26,
                  58.24,58.24,58.24,58.24,63.39,63.39,63.39,63.39,48.90,48.90,48.90,48.90,
                  42.45,42.45,42.45,42.45,39.84,39.84,39.84,39.84)
Pheno_M700k=cbind(Pheno_M700k, VO2max_PRE_Meth)

####Create a matrix to save results####
CpGs=c(rownames(M))
ResultsTime_Meth=matrix(NA, nrow = 721543, ncol=3)
colnames(ResultsTime_Meth)=c("Estimate", "SE", "p-value")
rownames(ResultsTime_Meth)=CpGs

ResultsVO2max_Meth=matrix(NA, nrow = 721543, ncol=3)
colnames(ResultsVO2max_Meth)=c("Estimate", "SE", "p-value")
rownames(ResultsVO2max_Meth)=CpGs

ResultsRandom_Meth=matrix(NA, nrow = 721543, ncol=3)
colnames(ResultsRandom_Meth)=c("Residual Variance", "Residual SE", "p-value")
rownames(ResultsRandom_Meth)=CpGs

colnames(Pheno_M700k)=c("ID", "Names","Kcals_inter","Age","VO2max_rel", "Type_I_perc", "time", 
                        rownames(M), "VO2max_PRE_Meth")
cols=c(colnames(Pheno_M700k[,8:721550]))

Pheno_M700k=drop_na(Pheno_M700k)

#loop by cols 
for (c in cols){
  model=lmer(Pheno_M700k[,c]~Pheno_M700k[,'VO2max_PRE_Meth']+
              Pheno_M700k[,'time']+
               Pheno_M700k[,'Age']+
               (1+Pheno_M700k[,'time']|Pheno_M700k[,'ID']))
  
  ResultsTime_Meth[c,"Estimate"]=summary(model)$coefficients[3,"Estimate"]
  ResultsVO2max_Meth[c,"Estimate"]=summary(model)$coefficients[2,"Estimate"]
  
  ResultsTime_Meth[c,"SE"]=summary(model)$coefficients[3,"Std. Error"]
  ResultsVO2max_Meth[c,"SE"]=summary(model)$coefficients[2,"Std. Error"]
  
  ResultsTime_Meth[c,"p-value"]=summary(model)$coefficients[3,"Pr(>|t|)"]
  ResultsVO2max_Meth[c,"p-value"]=summary(model)$coefficients[2,"Pr(>|t|)"]
  
  ResultsRandom_Meth[c,"Residual Variance"]=as.data.frame(VarCorr(model))[4,4]
  ResultsRandom_Meth[c,"Residual SE"]=as.data.frame(VarCorr(model))[4,5]
  ResultsRandom_Meth[c,"p-value"]=as.data.frame(rand(model))[2,6]
  
  print(summary(model))
  
}

#create vectors for p-values to adjust
p_VO2max_Meth=c(ResultsVO2max_Meth[,3])
p_Time_Meth=c(ResultsTime_Meth[,3])
p_Random_Meth=c(ResultsRandom_Meth[,3])

p_adj_VO2max_Meth=p.adjust(p_VO2max_Meth, method = "fdr")
p_adj_Time_Meth=p.adjust(p_Time_Meth, method = "fdr")
p_adj_Random_Meth=p.adjust(p_Random_Meth, method = "fdr")

#calculate and add p-values
ResultsTimeSig_Meth=cbind(ResultsTime_Meth,p_adj_Time_Meth)

ResultsRandomSig_Meth=cbind(ResultsRandom_Meth,p_adj_Random_Meth)

ResultsVO2maxSig_Meth=cbind(ResultsVO2max_Meth,p_adj_VO2max_Meth)

#write.csv(ResultsVO2maxSig_Meth, "ResultsVO2max_Meth2109.csv")
#write.csv(ResultsRandomSig_Meth, "ResultsRandom_Meth2109.csv")
#write.csv(ResultsTimeSig_Meth, "ResultsTime_Meth2109.csv")

ResultsVO2maxSig_Meth=read.csv("ResultsVO2max_Meth2109.csv")
ResultsRandomSig_Meth=read.csv("ResultsRandom_Meth2109.csv")
ResultsTimeSig_Meth=read.csv("ResultsTime_Meth2109.csv")

#Filter results based on significance
ResultsTimeSig_Meth=as.data.frame(ResultsTimeSig_Meth)
ResultsTimeSig_Meth2=ResultsTimeSig_Meth%>%
  filter(p_adj_Time_Meth<0.05)

hist_meth=ggplot(ResultsTimeSig_Meth, aes(x=`p.value`))+
  geom_histogram(fill="#CCCCCC",color="#3399FF")+
  theme_minimal()

tiff('histogram methylation.tiff', width =6, height = 4, units = 'in', res=600)
hist_meth
dev.off()

#Plot PCA for residuals
#Get CpGs from model
ResultsTimeSig_Meth=ResultsTimeSig_Meth[order(ResultsTimeSig_Meth$`p.value`),][1:10000,]

CpGs=rownames(ResultsTimeSig_Meth)

####Run a for loop without time as fixed effects to get residuals####
library("lmerTest")
rownames(M_subset)=rownames(M)
M_subset_transposed=t(M_subset)

M_S_t_700k=M_subset_transposed[,c(CpGs)]
Names=c(rownames(M_S_t_700k))

sapply(unlist(M_S_t_700k), as.numeric)

M_S_t_700k=as_tibble(M_S_t_700k)

Timepoint=c(SI_subset_long$Timepoint)

SI_subset_long2=SI_subset_long%>%
  unite(Names, "New ID", Timepoint)

Names=c("SI20_PRE", "SI17_4WP", "SI14_8WP", "SI18_12WP", "SI19_4WP","SI8_8WP","SI15_PRE", "SI15_12WP","SI17_PRE", "SI11_4WP", "SI9_8WP",  
        "SI18_4WP", "SI12_8WP", "SI15_4WP", "SI17_12WP", "SI14_12WP", "SI4_8WP",  "SI19_PRE","SI19_12WP", "SI16_PRE",  "SI6_8WP", "SI18_PRE", 
        "SI13_4WP", "SI16_8WP", "SI14_PRE", "SI11_12WP", "SI12_4WP", "SI13_12WP", "SI20_8WP", "SI11_PRE", "SI14_4WP", "SI13_PRE", "SI12_PRE", 
        "SI17_8WP", "SI10_12WP", "SI9_PRE", "SI13_8WP", "SI12_12WP", "SI10_PRE","SI9_12WP",  "SI8_12WP", "SI3_PRE",  "SI10_4WP",  "SI8_4WP",   
        "SI4_PRE", "SI20_4WP", "SI6_PRE", "SI4_12WP", "SI9_4WP",  "SI18_8WP", "SI4_4WP",  "SI19_8WP", "SI3_12WP", "SI5_PRE", "SI11_8WP", 
        "SI3_4WP",  "SI8_PRE","SI5_4WP", "SI15_8WP",  "SI6_4WP",  "SI6_12WP",  "SI5_8WP", "SI16_4WP",  "SI20_12WP", "SI10_8WP")
M_S_t_700k=cbind(M_S_t_700k, Names)

Pheno_M700k=full_join(SI_subset_long2, M_S_t_700k, by="Names")
Pheno_M700k=Pheno_M700k[,-6]
Pheno_M700k=drop_na(Pheno_M700k)

####Create a matrix to save results####
ResultsResiduals=matrix(NA, nrow = 61, ncol=10000)
colnames(ResultsResiduals)=CpGs
rownames(ResultsResiduals)=Pheno_M700k$Names

cols=c(colnames(Pheno_M700k[,7:10006]))

#loop by cols 
for (c in cols){
  model=lmer(pull(Pheno_M700k[,c])~
               pull(Pheno_M700k[,'VO2max_rel'])+
               pull(Pheno_M700k[,'Age'])+
               (1+pull(Pheno_M700k[,'time'])|pull(Pheno_M700k[,'ID'])))
  
  ResultsResiduals[,c]=data.frame(resid(model))[,1]
  
  print(summary(model))
  
}

#write.csv(ResultsResiduals, "ResultsResiduals22092021.csv")
Residuals=read.csv("ResultsResiduals22092021.csv")

#Run PCA analyses on residuals
library(FactoMineR)
library(factoextra)
PC=PCA(ResultsResiduals)
print(PC)

#extract eigenvalue 
#amount of variation retained by each principal component
#usually larger for the first PCs
eig.val<- get_eigenvalue(PC)
fviz_eig(PC, addlabels = TRUE, ylim=c(0,50))

#to extract results for variables from PCA output
#the components can be used as followed
#var$coord: coordinates of variables to create a scatter plot
#var$cos2: represents the quality of representation for variables on the factor map.
#var$contrib: contains the contributions (%) of the variables to the principal components
var=get_pca_var(PC)
fviz_pca_var(PC, col.var = "black")

#to get individuals results
ind<-get_pca_ind(PC)
coordinates<-data.frame(ind$coord)

coordinates=coordinates%>%
  mutate(time=c(1,2,3,1,2,3,4,1,
                2,3,4,1,2,3,4,1, 
                2,3,4,1,2,3,4,1, 
                2,3,4,1,2,3,4,1,
                2,3,4,1,2,4,1,2,   
                3,1,2,3,4,1,2,3,
                4,1,2,3,4,1,2,3,
                4,1,2,3,4))

#get mean values for each timepoint
Mean_PRE=coordinates%>%
  filter(time==1)
Mean_PRE=mean(Mean_PRE$Dim.1)

Mean_4WP=coordinates%>%
  filter(time==2)
Mean_4WP=mean(Mean_4WP$Dim.1)

Mean_8WP=coordinates%>%
  filter(time==3)
Mean_8WP=mean(Mean_8WP$Dim.1)

Mean_12WP=coordinates%>%
  filter(time==4)
Mean_12WP=mean(Mean_12WP$Dim.1)

Mean_PRE2=coordinates%>%
  filter(time==1)
Mean_PRE2=mean(Mean_PRE2$Dim.2)

Mean_4WP2=coordinates%>%
  filter(time==2)
Mean_4WP2=mean(Mean_4WP2$Dim.2)

Mean_8WP2=coordinates%>%
  filter(time==3)
Mean_8WP2=mean(Mean_8WP2$Dim.2)

Mean_12WP2=coordinates%>%
  filter(time==4)
Mean_12WP2=mean(Mean_12WP2$Dim.2)

#Run linear model to see if changes in methylome are significant overtime
model=lm(Dim.1~time, data=coordinates)
summary(model)


fviz_pca_ind(PC, col.ind = "cos2",
             gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

#Color by groups
Timepoint=coordinates$time
Timepoint=ifelse(Timepoint==1, "PRE", Timepoint)
Timepoint=ifelse(Timepoint==2, "4WP", Timepoint)
Timepoint=ifelse(Timepoint==3, "8WP", Timepoint)
Timepoint=ifelse(Timepoint==4, "12WP", Timepoint)
fviz_pca_ind(PC, geom.ind="point", col.ind = Timepoint,
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#660099"),
             addEllipses = TRUE,
             legend.title = "Timepoint")


Timepoint=factor(Timepoint, levels = c("PRE", "4WP", "8WP", "12WP"))

plot<-ggplot(coordinates, aes(x=Dim.1, y=Dim.2))+
  theme(panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(size=0.5, linetype = "solid", colour = "grey"))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(color=factor(Timepoint, levels = c("PRE", "4WP", "8WP", "12WP"))))+
  annotate("text", x=c(48.31,1.63,-27.25,-27.88), y=c(5.69,10.02,-10.04,-7.19), 
           label=c("Mean Pre", "Mean 4WP", "Mean 8WP", "Mean 12WP"), 
           color=c("#FC4E07", "#669900", "#00AFBB","#660099"))+
  labs(color="Timepoint")


plot

tiff('PCA methylation overtime.tiff', width =6, height = 4, units = 'in', res=600)
plot
dev.off()

####Plot significant proteins overtime as a heatmap####
ResultsTimeSig=as.data.frame(ResultsTimeSig)
ResultsTimeSig2=filter(ResultsTimeSig, p_adj_Time<=0.05)
ResultsTimeSig3=filter(ResultsTimeSig, abs(Estimate)>=0.5)

ResultsTimeSig4=ResultsTimeSig2$X

data_heatTime=select(data_irs_log3, c(ResultsTimeSig4))
rownames(data_heatTime)=data_irs_log3$IDs
data_heatTime=data.matrix(data_heatTime)
data_heatTime_t=t(data_heatTime)

means=rowMeans(data_heatTime_t, na.rm = TRUE)

df=data.matrix(data_heatTime_t-means)

library(ComplexHeatmap)
library(circlize)

#Separate matrix by timepoints
ptn_df1=data.matrix(df[,grep("PRE", colnames(df))])
ptn_df2=data.matrix(df[,grep("4WP", colnames(df))])
ptn_df3=data.matrix(df[,grep("8WP", colnames(df))])
ptn_df4=data.matrix(df[,grep("12WP", colnames(df))])

ht_list=Heatmap(ptn_df1, name="centered data",km=2, show_column_names= FALSE, show_row_names = FALSE,
                show_column_dend=FALSE, column_title = "PRE")+
  Heatmap(ptn_df2,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "4WP",show_heatmap_legend = FALSE,show_row_names = FALSE)+
  Heatmap(ptn_df3,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "8WP",show_heatmap_legend = FALSE,show_row_names = FALSE)+
  Heatmap(ptn_df4,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "12WP",show_heatmap_legend = FALSE,show_row_names = FALSE)
draw(ht_list, row_title="Significant proteins overtime",
     annotation_legend_side = "left", heatmap_legend_side = "left")


tiff('Heatmap Proteins overtime.tiff', width =6, height = 8, units = 'in', res=600)
draw(ht_list, row_title="Significant proteins overtime",
     annotation_legend_side = "left", heatmap_legend_side = "left")
dev.off()

####Kmeans clustering ####
library(reshape2)
set.seed(1235334)

#how many number of clusters 
library(factoextra)

#First have to scale data
data_final_m=data.matrix(data_final[,1:2365])
protein_scaled<-scale(data_final_m, scale = TRUE)
protein_scaled=data.frame(protein_scaled)
protein_scaled=rename(protein_scaled, AMPD3=AMPD3.x)
protein_scaled=rename(protein_scaled, USP2=USP2.x)

ptn_scaled=protein_scaled%>%
  select(ResultsTimeSig2$X)%>%
  mutate(Timepoint=data_final$Timepoint, ID=data_final$ID)%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))%>%
  arrange(Timepoint)%>%
  pivot_longer(!c("Timepoint","ID"),
               names_to = 'proteins', #indicates that component of the name
               #defines the name of the column containing the cell values, overriding values_to.
               values_to="intensity",
               values_drop_na = 'FALSE')%>%
  pivot_wider(id_cols = proteins, names_from=Timepoint, values_from=intensity, values_fn=list(intensity=mean))

#is.na(data_final_wide) %>% sum()

qc_cluster<-fviz_nbclust(ptn_scaled[,2:5], kmeans, method = "wss", k.max = 10)
print(qc_cluster)

#protein_kmeans<-kmeans(protein_scaled[1:3],centers = 4) 
protein_kmeans<-kmeans(ptn_scaled[,2:5],centers = 6) 

protein_cluster<-data.frame(ptn_scaled, 
                            cluster=factor(protein_kmeans$cluster)) 


colnames(protein_cluster)<-gsub("_protein", "", colnames(protein_cluster))

#protein_cluster<-
protein_ribbon<-melt(protein_cluster, variable.name = "stage", value.name = "intensity") %>% 
  group_by(stage, cluster) %>% 
  summarise(clust_min=min(intensity), 
            clust_max=max(intensity),
            clust_mean=mean(intensity))

protein_long<-melt(protein_cluster, variable.name = "stage", value.name = "intensity")

protein_cluster_merged<-merge(protein_long,protein_ribbon, 
                              by=c("stage","cluster"))

levels(protein_cluster_merged$stage)<-c("PRE","4WP","8WP","12WP" )

protein_kmeans$size ## 50 42 82 69 112 106

protein_cluster_merged$cluster_size<-protein_cluster_merged$cluster

levels(protein_cluster_merged$cluster_size)<-protein_kmeans$size

#write.csv(protein_cluster, "kmeans protein clusters.csv")

plot<-ggplot(protein_cluster_merged, aes(x=stage, group=1, color=cluster_size)) +
  geom_ribbon(aes(ymin=clust_min, ymax= clust_max, fill= cluster_size), alpha=0.4 ) +
  geom_line(aes(y=clust_mean))+
  facet_wrap(~cluster_size)+
  labs(title= "K-means Clustering", 
       x="Timepoint", y="Scaled Expression")+
  theme_bw()+
  # scale_fill_brewer(palette = "Dark2")+
  # scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text = element_text(size=10, face = "bold"))+ # facetwrap font and bold
  theme(legend.position="none")
print(plot)

tiff('Kmeans Proteins overtime.tiff', width =8, height = 6, units = 'in', res=600)
print(plot)
dev.off()

####Investigate trainability####
####Select all proteins that change overtime and present trainability####
ResultsRandomWpeak=read.csv("ResultsRandom_log_ProteomicsWpeak.csv")
ResultsRandomLT=read.csv("ResultsRandom_log_ProteomicsLT.csv")
ResultsRandomVO2max=read.csv("ResultsRandom_log_ProteomicsVO2max.csv")
ResultsRandomTime=read.csv("ResultsRandom_log_ProteomicsTime.csv")

ResultsRandomWpeak=filter(ResultsRandomWpeak, p_adj_Random<=0.05)
ResultsRandomLT=filter(ResultsRandomLT, p_adj_Random<=0.05)
ResultsRandomVO2max=filter(ResultsRandomVO2max, p_adj_Random<=0.05)
ResultsRandomTime=filter(ResultsRandomTime, p_adj_Random<=0.05)

ResultsRandomCombined=rbind(ResultsRandomWpeak, ResultsRandomLT,ResultsRandomVO2max,ResultsRandomTime)
ResultsRandomCombined=ResultsRandomCombined[order(ResultsRandomCombined$p_adj_Random),]
ResultsRandomCombined=ResultsRandomCombined[!duplicated(ResultsRandomCombined$X),]

rownames(ResultsTimeSig)=ResultsTimeSig$X
ResultsGT=ResultsTimeSig[ResultsRandomCombined$X,]
ResultsGT=filter(ResultsGT, p_adj_Time<0.05)
ResultsGT2=filter(ResultsGT, abs(Estimate)>=0.5)

data_final=dplyr::rename(data_final, AMPD3=AMPD3.x, USP2=USP2.x, SOD3=SOD3.x)

PTN_Random=dplyr::select(data_final, rownames(ResultsGT2),Timepoint:time)

PTN_Random<- PTN_Random%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

PTN_Random2=pivot_longer(data = PTN_Random,
                         cols = -c(Timepoint:time),
                         names_to = "Proteins") #indicates that component of the name

library(ggplot2)
ptn2=ggplot(data = PTN_Random2,
            mapping = aes(x = Timepoint,
                          y = value))+
  geom_point(aes(fill=ID, group=seq_along(time)),col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "Intensity")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+
  facet_wrap(~Proteins)

ptn2

####run over-representation analyses for enrichment for each k-mean of proteins using clusterprofiler####
library(mitch)
library(qusage)
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- read.gmt("ReactomePathways.gmt")

library(clusterProfiler)
library(enrichplot)

#separate the clusters for the analysis
de_k1Table = filter(protein_cluster, cluster==1)
de_k1Table= de_k1Table$proteins

de_k2Table = filter(protein_cluster, cluster==2)
de_k2Table= de_k2Table$proteins

de_k3Table = filter(protein_cluster, cluster==3)
de_k3Table= de_k3Table$proteins

de_k4Table = filter(protein_cluster, cluster==4)
de_k4Table= de_k4Table$proteins

de_k5Table = filter(protein_cluster, cluster==5)
de_k5Table= de_k5Table$proteins

de_k6Table = filter(protein_cluster, cluster==6)
de_k6Table= de_k6Table$proteins

de_bg=unique(protein_cluster$proteins)

de_trainability=unique(ResultsRandomCombined$X)

o_k1_res <- as.data.frame(enricher(gene = de_k1Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k1 <- rownames(subset(o_k1_res, p.adjust < 0.05))

o_k2_res <- as.data.frame(enricher(gene = de_k2Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k2 <- rownames(subset(o_k2_res, p.adjust < 0.05))

o_k3_res <- as.data.frame(enricher(gene = de_k3Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k3 <- rownames(subset(o_k3_res, p.adjust < 0.05))

o_k4_res <- as.data.frame(enricher(gene = de_k4Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k4 <- rownames(subset(o_k4_res, p.adjust < 0.05))
  
o_k5_res <- as.data.frame(enricher(gene = de_k5Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k5 <- rownames(subset(o_k5_res, p.adjust < 0.05))

o_k6_res <- as.data.frame(enricher(gene = de_k6Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 ))
o_k6 <- rownames(subset(o_k6_res, p.adjust < 0.05))

o_com_res <- as.data.frame(enricher(gene = de_bg, 
                                    universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                    pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))
o_com <- rownames(subset(o_com_res, p.adjust < 0.05))


o_tranability_res <- as.data.frame(enricher(gene = de_trainability, 
                                    universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                    pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))
o_trainability <- rownames(subset(o_tranability_res, p.adjust < 0.05))


o_k1_res2 <-enricher(gene = de_k1Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_k2_res2 <-enricher(gene = de_k2Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                                   pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_k3_res2 <-enricher(gene = de_k3Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_k4_res2 <-enricher(gene = de_k4Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_k5_res2 <-enricher(gene = de_k5Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_k6_res2 <-enricher(gene = de_k6Table, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_com_res2 <-enricher(gene = de_bg, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                     pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

o_trainability_res2 <-enricher(gene = de_trainability, universe = colnames(data_final),  maxGSSize = 5000, TERM2GENE = genesets, 
                      pAdjustMethod="fdr", pvalueCutoff = 1, qvalueCutoff = 1 )

length(o_k1)
length(o_k2)
length(o_k3)
length(o_k4)
length(o_k5)
length(o_k6)

length(o_k1) + length(o_k2)+ length(o_k3)+ length(o_k4)+ length(o_k5)+ length(o_k6)
length(o_com)

dotplot(o_k1_res2, showCategory=20)
dotplot(o_k2_res2, showCategory=20)
dotplot(o_k3_res2, showCategory=20)
dotplot(o_k4_res2, showCategory=18)
dotplot(o_k5_res2, showCategory=6)
dotplot(o_k6_res2, showCategory=20)
dotplot(o_com_res2, showCategory=15)
dotplot(o_trainability_res2, showCategory=10)

tiff('dotplot ORA enrichment all de proteins.tiff', width =9, height = 11, units = 'in', res=600)
dotplot(o_com_res2, showCategory=14)
dev.off()

tiff('dotplot ORA enrichment k4.tiff', width =9, height = 13, units = 'in', res=600)
dotplot(o_k4_res2, showCategory=18)
dev.off()

tiff('dotplot ORA enrichment k5.tiff', width =9, height = 6, units = 'in', res=600)
dotplot(o_k5_res2, showCategory=6)
dev.off()

library("topGO")
plotGOgraph(o_k4_res2)

####Plot significant individual Proteins overtime####
PTN_Random<- PTN_Random%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

ptn=colnames(PTN_Random[1:7])

library(ggplot2)

for (pn in ptn) {
  p1=ggplot(data = PTN_Random,
            mapping = aes_string(x = "Timepoint",
                                 y = pn))+
    geom_point(aes(fill=ID), col="black",
               pch = 21,
               size = 3, show.legend=FALSE)+
    theme_classic(base_size = 14)+
    labs(y = "Intensity")+
    geom_line(aes(group = ID,
                  col = ID),
              show.legend=FALSE)+ggtitle(pn)
  ggsave(p1, file=paste0("plot_",pn,"_log.png"))
  
}



#Plot significant proteins overtime as a heatmap
ResultsRandomCombined2=filter(ResultsRandomCombined, p_adj_Random<=0.05)
ResultsRandomSig3=ResultsRandomCombined$X

data_heatTime=dplyr::select(data_irs_log3, c(ResultsRandomSig3))
rownames(data_heatTime)=data_irs_log3$IDs
data_heatTime=data.matrix(data_heatTime)
data_heatTime_t=t(data_heatTime)

means=rowMeans(data_heatTime_t, na.rm = TRUE)

df=data.matrix(data_heatTime_t-means)

library(ComplexHeatmap)
library(circlize)

#Separate matrix by timepoints
ptn_df1=data.matrix(df[,grep("PRE", colnames(df))])
ptn_df2=data.matrix(df[,grep("4WP", colnames(df))])
ptn_df3=data.matrix(df[,grep("8WP", colnames(df))])
ptn_df4=data.matrix(df[,grep("12WP", colnames(df))])

ht_list=Heatmap(ptn_df1, name="centered data",km=4, show_column_names= FALSE, show_row_names = FALSE,
                show_column_dend=FALSE, column_title = "PRE")+
  Heatmap(ptn_df2,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "4WP",show_heatmap_legend = FALSE,show_row_names = FALSE)+
  Heatmap(ptn_df3,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "8WP",show_heatmap_legend = FALSE,show_row_names = FALSE)+
  Heatmap(ptn_df4,  show_column_names= FALSE,
          show_column_dend=FALSE, column_title = "12WP",show_heatmap_legend = FALSE,show_row_names = TRUE)
draw(ht_list, row_title="Significant proteins overtime",
     annotation_legend_side = "left", heatmap_legend_side = "left")


tiff('Heatmap Proteins Trainability.tiff', width =10, height = 16, units = 'in', res=600)
draw(ht_list, row_title="Trainability proteins overtime",
     annotation_legend_side = "left", heatmap_legend_side = "left")
dev.off()

#Individual CpGs
annotation=read.delim('Annotation_new.txt')
rownames(annotation)=annotation[,1]

annotation_subset=filter(annotation, genesUniq %in% rownames(ResultsGT2))

M_subset2=M_subset[annotation_subset$probeID,]
rownames(M_subset2)=rownames(annotation_subset)
M_subset2=as.data.frame(M_subset2)
M_subset2=drop_na(M_subset2)
cpgs=rownames(M_subset2)

####Run loop with random effects to look at trainability####
library("lmerTest")
library("tidyverse")
library("broom")


#Transform data so you can fit loop
M_subset_transposed=t(M_subset2)
M_subset_transposed=as.data.frame(M_subset_transposed)
M_subset_transposed=select(M_subset_transposed, -contains("NA"))

sapply(unlist(M_subset_transposed), as.numeric)

CpGs=c(colnames(M_subset_transposed))

Timepoint2=c(SI_subset_long$Timepoint)
SI_subset_long=cbind(SI_subset_long, Timepoint2)

SI_subset_long2=SI_subset_long%>%
  unite(Names, "New ID", Timepoint)

Names=c("SI20_PRE", "SI17_4WP", "SI14_8WP", "SI18_12WP", "SI19_4WP","SI8_8WP","SI15_PRE", "SI15_12WP","SI17_PRE", "SI11_4WP", "SI9_8WP",  
        "SI18_4WP", "SI12_8WP", "SI15_4WP", "SI17_12WP", "SI14_12WP", "SI4_8WP",  "SI19_PRE","SI19_12WP", "SI16_PRE",  "SI6_8WP", "SI18_PRE", 
        "SI13_4WP", "SI16_8WP", "SI14_PRE", "SI11_12WP", "SI12_4WP", "SI13_12WP", "SI20_8WP", "SI11_PRE", "SI14_4WP", "SI13_PRE", "SI12_PRE", 
        "SI17_8WP", "SI10_12WP", "SI9_PRE", "SI13_8WP", "SI12_12WP", "SI10_PRE","SI9_12WP",  "SI8_12WP", "SI3_PRE",  "SI10_4WP",  "SI8_4WP",   
        "SI4_PRE", "SI20_4WP", "SI6_PRE", "SI4_12WP", "SI9_4WP",  "SI18_8WP", "SI4_4WP",  "SI19_8WP", "SI3_12WP", "SI5_PRE", "SI11_8WP", 
        "SI3_4WP",  "SI8_PRE","SI5_4WP", "SI15_8WP",  "SI6_4WP",  "SI6_12WP",  "SI5_8WP", "SI16_4WP",  "SI20_12WP", "SI10_8WP")
M_subset_transposed=cbind(M_subset_transposed, Names)

Pheno_M700k=full_join(SI_subset_long2, M_subset_transposed, by="Names")

####Create a matrix to save results####
ResultsTime=matrix(NA, nrow = 188, ncol=3)
colnames(ResultsTime)=c("Estimate", "SE", "p-value")
rownames(ResultsTime)=cpgs

ResultsRandom=matrix(NA, nrow = 188, ncol=3)
colnames(ResultsRandom)=c("Residual Variance", "Residual SE", "p-value")
rownames(ResultsRandom)=CpGs

cols=c(colnames(Pheno_M700k[,9:196]))

Pheno_M700k=as_tibble(Pheno_M700k)
Pheno_M700k=Pheno_M700k[,-c(3,5,6)]
Pheno_M700k=drop_na(Pheno_M700k)

#loop by cols 
for (c in cols){
  model=lmer(pull(Pheno_M700k[,c])~
               pull(Pheno_M700k[,'time'])+
               pull(Pheno_M700k[,'Age'])+
               (1+pull(Pheno_M700k[,'time'])|pull(Pheno_M700k[,'ID'])))
  
  ResultsTime[c,"Estimate"]=summary(model)$coefficients[2,"Estimate"]
  ResultsTime[c,"SE"]=summary(model)$coefficients[2,"Std. Error"]
  ResultsTime[c,"p-value"]=summary(model)$coefficients[2,"Pr(>|t|)"]
  
  ResultsRandom[c,"Residual Variance"]=as.data.frame(VarCorr(model))[4,4]
  ResultsRandom[c,"Residual SE"]=as.data.frame(VarCorr(model))[4,5]
  ResultsRandom[c,"p-value"]=as.data.frame(rand(model))[2,6]
  
  print(summary(model))
  
}

p_Time=c(ResultsTime[,3])
p_RandomTime=c(ResultsRandom[,3])

#calculate and add p-values
p_adj_Time=p.adjust(p_Time, method = "fdr")

p_adj_RandomTime=p.adjust(p_RandomTime, method = "fdr")

ResultsTime=cbind(ResultsTime, p_adj_Time)
ResultsRandom=cbind(ResultsRandom, p_adj_RandomTime)

setwd("D:/PhD related/Proteomics project/Data for analyses/Methylation for the 8 trainability proteins")
write.csv(ResultsTime, "Results time 8 proteins Group.csv")
write.csv(ResultsRandom, "Results random 8 proteins individual.csv")

####Validation for PCRs####
#Subset variables of interest
library(tidyverse)
SI_subset=Pheno_Second_intervention%>%
  dplyr::select(ID, "New ID", Age_PRE:Age_12WP, Kcals_inter, Wpeak_rel_PRE, Wpeak_rel_4WP, Wpeak_rel_8WP, Wpeak_rel_12WP, 
                LT_poly_rel_PRE, LT_poly_rel_4WP, LT_poly_rel_8WP, LT_poly_rel_12WP,
                VO2max_rel_PRE, VO2max_rel_4WP, VO2max_rel_8WP, VO2max_rel_12WP, 
                "Type_I_PRE":"SOD3_12WP")%>%
  dplyr::rename("New_ID"="New ID")

#change data to a longer format
SI_subset_long <- pivot_longer(data = SI_subset,
                               cols = -c("ID","New_ID", Kcals_inter),
                               names_to = c(".value","Timepoint"), #indicates that component of the name
                               #defines the name of the column containing the cell values, overriding values_to.
                               names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",
                               values_drop_na = 'FALSE')


#add timepoint as a number so models will work further on
SI_subset_long=SI_subset_long%>%
  mutate(time=rep(c(1,2,3,4), times=20))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))


#Do imputations for missing VO2 values
library(mice)

#Use SI_zscore data subset
md.pattern(SI_subset)

#For mice to work, we need to ensure the variables we want to impute are
#in the right format (numeric, factor, etc.)
SI_subset

#examine the number of missing cases in a table format and proportion of missing
missing.indicator=data.frame(is.na(SI_subset))
propMissing=apply(missing.indicator, 2, mean)

#To impute the missing data, use the mice() function
pred_mat=quickpred(SI_subset, mincor = 0.25) #quickly selects data that has min corr of 0.25
imputed <- mice(data = SI_subset, #dataset containing missing values
                m =  5, #number of imputations you want (default = 5)
                #method = NULL. It tells the algorithm HOW it will impute the data, depending on data type (continuous variable, factor with 1 or more levels)
                print=F,
                pred=pred_mat
)

#Have a look at the summary of imputations
imputed

#To access imputed values:
imputed$imp
plot(imputed)

SI_subset1=complete(imputed, 3)
SI_subset2=complete(imputed,4)
SI_subset3=complete(imputed, 5)

#pivot longer for the model
SI_subset_long1 <- pivot_longer(data = SI_subset1,
                               cols = -c("ID","New_ID", Kcals_inter),
                               names_to = c(".value","Timepoint"), #indicates that component of the name
                               #defines the name of the column containing the cell values, overriding values_to.
                               names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",

                                                              values_drop_na = 'FALSE')
#add timepoint as a number so models will work further on
SI_subset_long1=SI_subset_long1%>%
  mutate(time=rep(c(1,2,3,4), times=20))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))

#pivot longer for the model
SI_subset_long2 <- pivot_longer(data = SI_subset2,
                               cols = -c("ID","New_ID", Kcals_inter),
                               names_to = c(".value","Timepoint"), #indicates that component of the name
                               #defines the name of the column containing the cell values, overriding values_to.
                               names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",
                               
                               values_drop_na = 'FALSE')
#add timepoint as a number so models will work further on
SI_subset_long2=SI_subset_long2%>%
  mutate(time=rep(c(1,2,3,4), times=20))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))


#pivot longer for the model
SI_subset_long3 <- pivot_longer(data = SI_subset3,
                                cols = -c("ID","New_ID", Kcals_inter),
                                names_to = c(".value","Timepoint"), #indicates that component of the name
                                #defines the name of the column containing the cell values, overriding values_to.
                                names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",
                                
                                values_drop_na = 'FALSE')
#add timepoint as a number so models will work further on
SI_subset_long3=SI_subset_long3%>%
  mutate(time=rep(c(1,2,3,4), times=20))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))

library(lmerTest)
AMPD3_model=lmer(AMPD3~time+Age+(1+time|ID), SI_subset_long)
summary(AMPD3_model)
rand(AMPD3_model)

USP2_model=lmer(USP2~time+Age+(1+time|ID), SI_subset_long)
summary(USP2_model)
rand(USP2_model)

SOD3_model=lmer(SOD3~time+Age+(1+time|ID), SI_subset_long)
summary(SOD3_model)
rand(SOD3_model)


SI_subset_long<- SI_subset_long%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

ggplot(data = SI_subset_long,
       mapping = aes_string(x = "Timepoint",
                            y = "AMPD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("AMPD3")


ggplot(data = SI_subset_long,
       mapping = aes_string(x = "Timepoint",
                            y = "USP2"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("USP2")

ggplot(data = SI_subset_long,
       mapping = aes_string(x = "Timepoint",
                            y = "SOD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("SOD3")


AMPD3_model1=lmer(AMPD3~time+Age+(1+time|ID), SI_subset_long1)
summary(AMPD3_model1)
rand(AMPD3_model1)

USP2_model1=lmer(USP2~time+Age+(1+time|ID), SI_subset_long1)
summary(USP2_model1)
rand(USP2_model1)

SOD3_model1=lmer(SOD3~time+Age+(1+time|ID), SI_subset_long1)
summary(SOD3_model1)
rand(SOD3_model1)

SI_subset_long1<- SI_subset_long1%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

ggplot(data = SI_subset_long1,
       mapping = aes_string(x = "Timepoint",
                            y = "AMPD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("AMPD3")


ggplot(data = SI_subset_long1,
       mapping = aes_string(x = "Timepoint",
                            y = "USP2"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("USP2")

ggplot(data = SI_subset_long1,
       mapping = aes_string(x = "Timepoint",
                            y = "SOD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("SOD3")


AMPD3_model2=lmer(AMPD3~time+Age+(1+time|ID), SI_subset_long2)
summary(AMPD3_model2)
rand(AMPD3_model2)

USP2_model2=lmer(USP2~time+Age+(1+time|ID), SI_subset_long2)
summary(USP2_model2)
rand(USP2_model2)

SOD3_model2=lmer(SOD3~time+Age+(1+time|ID), SI_subset_long2)
summary(SOD3_model2)
rand(SOD3_model2)

SI_subset_long2<- SI_subset_long2%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

ggplot(data = SI_subset_long2,
       mapping = aes_string(x = "Timepoint",
                            y = "AMPD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("AMPD3")


ggplot(data = SI_subset_long2,
       mapping = aes_string(x = "Timepoint",
                            y = "USP2"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("USP2")

ggplot(data = SI_subset_long2,
       mapping = aes_string(x = "Timepoint",
                            y = "SOD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("SOD3")



AMPD3_model3=lmer(AMPD3~time+Age+(1+time|ID), SI_subset_long3)
summary(AMPD3_model3)
rand(AMPD3_model3)

USP2_model3=lmer(USP2~time+Age+(1+time|ID), SI_subset_long3)
summary(USP2_model3)
rand(USP2_model3)

SOD3_model3=lmer(SOD3~time+Age+(1+time|ID), SI_subset_long3)
summary(SOD3_model3)
rand(SOD3_model3)

SI_subset_long3<- SI_subset_long3%>%
  mutate(Timepoint = factor(Timepoint,
                            levels = c("PRE","4WP","8WP","12WP")))

ggplot(data = SI_subset_long3,
       mapping = aes_string(x = "Timepoint",
                            y = "AMPD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("AMPD3")


ggplot(data = SI_subset_long3,
       mapping = aes_string(x = "Timepoint",
                            y = "USP2"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("USP2")

ggplot(data = SI_subset_long3,
       mapping = aes_string(x = "Timepoint",
                            y = "SOD3"))+
  geom_point(aes(fill=ID), col="black",
             pch = 21,
             size = 3, show.legend=FALSE)+
  theme_classic(base_size = 14)+
  labs(y = "logFC")+
  geom_line(aes(group = ID,
                col = ID),
            show.legend=FALSE)+ggtitle("SOD3")

#Transform imputed data to long format 
#Calculate MHI in the imputed data
long1=complete(imputed, action = 'long', include = TRUE)

#pivot longer for the model
long2 <- pivot_longer(data = long1,
                                cols = -c(".imp",".id","ID","New_ID", Kcals_inter),
                                names_to = c(".value","Timepoint"), #indicates that component of the name
                                #defines the name of the column containing the cell values, overriding values_to.
                                names_pattern = "(.+)_(PRE|4WP|8WP|12WP)$",
                                
                                values_drop_na = 'FALSE')
#add timepoint as a number so models will work further on
long2=long2%>%
  mutate(time=rep(c(1,2,3,4), times=120))%>%
  mutate(".id"=rep(c(1:80), times=6))%>%
  mutate(Timepoint=factor(Timepoint, levels=c("PRE", "4WP", "8WP", "12WP")))

imputed=as.mids(long2)

#convert into a list of datasets
datlist2=miceadds::mids2datlist(imputed)

#Build a predictive model
library(lme4)
AMPD3model=with(datlist2, lme4::lmer(AMPD3~time+Age+(1+time|ID)))
AMPD3mod=miceadds::lmer_pool(AMPD3model)
AMPD3mo=summary(AMPD3mod)

USP2model=with(datlist2, lme4::lmer(USP2~time+Age+(1+time|ID)))
USP2mod=miceadds::lmer_pool(USP2model)
USP2mo=summary(USP2mod)

SOD3model=with(datlist2, lme4::lmer(SOD3~time+Age+(1+time|ID)))
SOD3mod=miceadds::lmer_pool(SOD3model)
SOD3mo=summary(SOD3mod)

#Proteins and transcriptome correlation
USP2_Cor_prot_trans=ggscatter(data=data_final, x="USP2.x", y="USP2.y",
                       add = "reg.line",  # Add regressin line
                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                       conf.int = TRUE, # Add confidence interval
                       cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                       cor.coeff.args = list(method = "spearman", label.sep = "\n"),
                       #cor.coef.coord = c(-800, 400),
                       xlab="USP2 Protein",
                       ylab="USP2 Transcript",
                       font.x=12,
                       font.y=12,
                       #xlim=c(-900,1000),
                       font.xtickslab=12,
                       font.ytickslab=12)

AMPD3_Cor_prot_trans=ggscatter(data=data_final, x="AMPD3.x", y="AMPD3.y",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE, # Add confidence interval
                              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                              cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                              #cor.coef.coord = c(-800, 400),
                              xlab="USP2 Protein",
                              ylab="USP2 Transcript",
                              font.x=12,
                              font.y=12,
                              #xlim=c(-900,1000),
                              font.xtickslab=12,
                              font.ytickslab=12)

SOD3_Cor_prot_trans=ggscatter(data=data_final, x="SOD3.x", y="SOD3.y",
                               add = "reg.line",  # Add regressin line
                               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                               conf.int = TRUE, # Add confidence interval
                               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                               cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                               #cor.coef.coord = c(-800, 400),
                               xlab="USP2 Protein",
                               ylab="USP2 Transcript",
                               font.x=12,
                               font.y=12,
                               #xlim=c(-900,1000),
                               font.xtickslab=12,
                               font.ytickslab=12)


####Integration with mixed-omics####
#first integrate group results by time
library(mixOmics)
Residuals_CpGs=colnames(Residuals[2:10001])
M_subset3=data.matrix(M_subset)
Residuals2=M_subset3[c(rownames(M_subset3)%in%Residuals_CpGs),data_irs_log3$IDs]

Prot_mixomis=data_irs_log3[,1:2365]
rownames(Prot_mixomis)=data_irs_log3$IDs

X=t(Residuals2)
Y=Prot_mixomis

dim(X)
dim(Y)

#Generate repeated measure table
repeat.indiv <- rownames(X)
repeat.indiv <- c(1,   2,  3,  4,   5, 6, 7, 8, 9,  10,  11,  12,
                  13, 14,  15,  2,   3,  4,   6,  7, 8,  9, 1, 10, 
                  16,  11, 17,  12,  13, 14, 7, 8,  9, 15, 10, 16, 
                  11,  17,  12, 13, 14, 2,   3,  4,   5,  6,  1, 7,
                  8,9, 10, 16, 17, 12,13, 14, 2,  4,  5,  6,
                  15)
design <- data.frame(sample = repeat.indiv) # load this into a dataframe

summary(as.factor(repeat.indiv)) # 16 rats, 3-4 measurements each

BPPARAM <- BiocParallel::SnowParam(workers = max(parallel::detectCores()-1, 2))
spls.muscle <- spls(X, Y, ncomp = 5, mode = 'regression')
perf.spls.muscle <- perf(spls.muscle, validation = 'Mfold', folds=10, nrepeat = 5, BPPARAM = BPPARAM)
plot(perf.spls.muscle, criterion='Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(50, 500, 50))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(seq(20, 300, 10))

tune.spls.muscle <- tune.spls(X, Y, ncomp = 4,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor', BPPARAM = BPPARAM) 
plot(tune.spls.muscle)         # use the correlation measure for tuning

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.muscle$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.muscle$choice.keepY

optimal.ncomp <-  length(optimal.keepX)

tunespls=tune.splslevel(X, Y, multilevel = design, ncomp = 1, test.keepX=c(50,50,300,50), test.keepY=c(300,30,280,110))

spls.muscle.multilevel <- spls(X, Y, # generate a tuned sPLS model
                              multilevel = design,
                              ncomp = optimal.ncomp,
                              keepX = optimal.keepX, 
                              keepY = optimal.keepY,
                              mode = 'canonical')

#Plots
# project onto averaged components, use time as colour, use dosage as symbol
setwd("D:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Proteomics and methylation changes at the individual level")
tiff('mixomics PCA individual.tiff', width =9, height = 8, units = 'in', res=600)
plotIndiv(spls.muscle.multilevel, rep.space = "XY", 
          group = c("PRE", "PRE","PRE","PRE", "PRE","PRE","PRE","PRE","PRE","PRE","PRE","PRE","PRE","PRE","4WP","4WP","4WP", "4WP", "4WP", "4WP","4WP","4WP","4WP", 
                    "4WP","4WP","4WP","4WP","4WP","4WP","4WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP","8WP", "8WP", "8WP", "8WP", "8WP", 
                    "12WP","12WP","12WP","12WP","12WP", "12WP", "12WP", "12WP","12WP", "12WP", "12WP","12WP","12WP","12WP","12WP"), 
          col.per.group = color.mixo(1:4), pch=16,
          legend = TRUE, legend.title = 'Timepoint')
dev.off()


tiff('mixomics correlation circle plot.tiff', width =8, height = 8, units = 'in', res=600)
plotVar(spls.muscle.multilevel, 
        var.names = TRUE, 
        cex = c(2,2))
dev.off()

heat=cim(spls.muscle.multilevel)

tiff('mixomics CIM plot.tiff', width =8, height = 8, units = 'in', res=600)
cim(spls.muscle.multilevel)
dev.off()

color.edge <- color.GreenRed(50) 
network(spls.muscle.multilevel, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS muscle Network Plot')

#Investigate correlations
Negative_Correlation=c("cg13415390","cg05130518","cg01346152","cg11337289","cg01906055","cg03699601","cg04864807","cg00171421",
                       "cg17274072","cg20512501","cg18345115","cg18995817","cg09943102","cg01894498","cg15763020","cg03967627",
                       "cg13957413","cg05103987","cg06631784","cg02852238","cg22662233","cg09414426","cg14353148","cg11518301",
                       "cg20941258","cg16691177","cg08517040","cg15709989","cg06296503","cg01639524","cg06085432","cg24707573",
                       "cg17666539","cg17885791","cg26394775","cg18022883","cg23726566","cg20324207","cg16867973","cg26843227",
                       "cg05554000","cg01975453","cg12513066","cg03957481","cg25882366")

Positive_Correlation=c("cg19977501","cg13106484","cg01330596","cg09830360","cg23701513","cg04978197","cg19990438","cg07775917",
                       "cg27236246","cg12964546","cg04499119","cg17936901","cg01454349","cg15146200","cg08492753","cg17411709",
                       "cg11099756","cg08689708","cg09005159","cg02707927","cg24691042","cg15005210","cg24120162","cg26581261",
                       "cg06799664","cg19249188","cg21834754","cg22846423","cg01362115","cg03498995","cg25386282","cg16233136",
                       "cg23840275","cg26537607","cg26537640","cg07581999","cg17161130","cg03286609","cg26692016","cg15492900",
                       "cg06716980","cg09005221","cg10666248","cg22323421","cg14314529")

Proteins_Correlation=c("GYG1","L2HGDH","MRPS16","MPC2","MRPL16","TMEM109","DHODH","HSD17B12","RTN4IP1","GCDH","MRPL40",
                       "ALDH1B1","NDUFAF2","SPRYD4","NAMPT","OXCT1","SF3B2","ACY1","BPNT1","SLC16A1","GADD45GIP1","GLUD1",
                       "TRAP1","FAHD1","MRPS27","WARS","MT.ND1","AK2","ATP5IF1","SUPV3L1","MRPS23","SDHD","LETM1","RAB5A",
                       "t.","SCO1","HSD17B10","MRPL37","GPT","ACAD8","NME2","ETFA","IBA57","CISD1","ADSSL1","PRXL2A", 
                       "GOT1","COQ3","SAMM50","ECHS1","GYS1","NDUFB11","SLC25A4","VDAC2","SLC25A11","NDUFB5","UQCRC1","ATP5PD",
                       "CKMT2","NDUFS4","COX4I1","LRPPRC","SLC25A3","HSPA9","NDUFB3","DLAT","NDUFB10","NDUFS2","NDUFA12",
                       "UQCRC2","PDHB","UQCRFS1","AK3","NDUFB6","CPT1B","ADHFE1","NDUFB7","NDUFS6","IMMT","ATP5F1B","ATP5F1A",
                       "GATD3B","MT.CO2","CRAT","NDUFS1","NDUFA7","HADHB","PRDX3","COX7A1","NDUFA4","ALDH5A1","DLST","SUCLA2",
                       "ACAA2","PDHX","CYC1","AIFM1","PHB","NDUFS8","IDH3B","SUCLG1","NIPSNAP3B","PERM1","LDHD","COX7A2","MP68",
                       "ETFDH","ATP5MG","TIMM44","HINT2","PDPR","ACADS","C1QBP","VWA8","SLIRP","GSTK1","NDUFV3","PDK2","TMEM65",
                       "IDH3G","MCCC2","HSD17B8","ATP5ME","MRPS36","UQCRQ","HADH","PTGES2","DECR1","SDHB","NDUFB8","PHB2",
                       "ACOT1","NDUFA10","UQCRB","VDAC3","OPA1","NDUFS5","NDUFA8","NIPSNAP2","ATP5PB","MECR","TST","GPD1L",
                       "ALDH6A1","TMLHE","ACSF2","TFAM","VDAC1","COQ5","ECI1","PPIF","CHCHD3","HSDL2","MDH2","MAOB","LACTB",
                       "ECHDC3","SLC25A12","DLD","PDHA1","HADHA","PCCB","ACOT13","HSPD1","OGDH","CS","ACAT1","MT.ATP6","ECHDC1",
                       "HIBADH","FABP3","MMUT","LDHB","PKM","CKB","BCS1L","ACO1","TXN2","NDUFAB1","NDUFA6", "ATP5MD","ATP5F1D",
                       "MT.ND5","ECH1","COX5B","FXN","ATP5PO","SUCLG2","IARS2","ACADM","NDUFC2","ATP5F1C","NDUFV2","ACADVL",
                       "NNT","ACSL1","CYCS","MICOS13","NDUFA9","ACO2","COX5A","GRHPR","TUFM","FH","SDHA","NDUFS3","GOT2","ETFB",
                       "MCEE","NDUFA13","NDUFV1","HSPE1","NDUFS7","SOD2","AFG3L2","IDH3A")

Positive_Correlation=annotation[c(rownames(annotation)%in%Positive_Correlation),]

Negative_Correlation=annotation[c(rownames(annotation)%in%Negative_Correlation),]


