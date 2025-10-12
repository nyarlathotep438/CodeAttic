#Set_Work_Dictionary####
setwd("C:/Users/96328/Documents/assessment_data") 

#Image_Save&Load####
save.image("image_assessment_1.RData")
load("image_assessment_1.RData")

#Library####
library(ggplot2)
library(ggrepel)

#Function####
get_PCA_plot = function(pca_data, group, group_title)
{
  ggp = ggplot(pca_data, aes(x=PC1, y= PC2, colour = group)) +
    geom_point() +
    labs(x= x_axis_label, y= y_axis_label, colour = group_title) +
    theme_bw()
  return(ggp)
}

save_plot= function(plot, plot.path, plot.height, plot.width)
{
  tiff(plot.path, width = plot.width, height = plot.height)
  print(plot)
  dev.off()
}

#Input_Data####
em = read.table("em.csv", header = TRUE, row.names = 1, sep="\t")
ss = read.table("ss.csv", header = TRUE, row.names = 1, sep="\t")
ss_per_patient = read.table("ss_per_patient.csv", header = TRUE, row.names = 1, sep="\t")

#Data clean####

## Basic clinic information
#Define Discrete Data
ss_per_patient$sample_group = as.factor(ss_per_patient$sample_group)
ss_per_patient$censored = as.factor(ss_per_patient$censored)
ss_per_patient$batch = as.factor(ss_per_patient$batch)

#Divide into two group:GB_1(metastasis) & GB_2(non metastasis)
ss_GB1 = subset(ss_per_patient,sample_group == "GB_1")
ss_GB2 = subset(ss_per_patient,sample_group == "GB_2")

#Summary Table
summary_of_ss = summary(ss_per_patient)
summary_of_ss_GB1 = summary(ss_GB1)
summary_of_ss_GB2 = summary(ss_GB2)

#Standard Error
age_sd = round(sd(ss_per_patient$age),2)
bmi_sd = round(sd(ss_per_patient$bmi),2)
days_survived_sd = round(sd(ss_per_patient$days_survived),2)
slide_area_sd = round(sd(ss_per_patient$slide_area),2)

age_GB1_sd = round(sd(ss_GB1$age),2)
bmi_GB1_sd = round(sd(ss_GB1$bmi),2)
days_survived_GB1_sd = round(sd(ss_GB1$days_survived),2)
slide_area_GB1_sd = round(sd(ss_GB1$slide_area),2)

age_GB2_sd = round(sd(ss_GB2$age),2)
bmi_GB2_sd = round(sd(ss_GB2$bmi),2)
days_survived_GB2_sd = round(sd(ss_GB2$days_survived),2)
slide_area_GB2_sd = round(sd(ss_GB2$slide_area),2)

#P-value
shapiro.test(ss_GB1$age)
shapiro.test(ss_GB1$bmi)
shapiro.test(ss_GB1$days_survived)
shapiro.test(ss_GB1$slide_area)

shapiro.test(ss_GB2$age)
shapiro.test(ss_GB2$bmi)
shapiro.test(ss_GB2$days_survived)
shapiro.test(ss_GB2$slide_area)

t_test_age = t.test(ss_GB1$age,ss_GB2$age)
t_test_slide_area = t.test(ss_GB1$slide_area,ss_GB2$slide_area)
t_test_days_survived = t.test(ss_GB1$days_survived,ss_GB2$days_survived)
t_test_BMI = t.test(ss_GB1$bmi,ss_GB2$bmi)

batch_Chisq_Test = chisq.test(ss_per_patient$batch,ss_per_patient$sample_group)
censored_Chisq_Test = chisq.test(ss_per_patient$censored,ss_per_patient$sample_group)

##PCA plot
# scale data
em_scaled = na.omit(data.frame(t(scale(t(em)))))

#PCA
xx = prcomp(t(em_scaled))
pca_coordinates = data.frame(xx$x)

#get % variation
vars = apply(xx$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100
x_axis_label = paste0("PC1 (" ,prop_x, "%)")
y_axis_label = paste0("PC2 (" ,prop_y, "%)")

#Continuous data grouping
age_discrete = cut(ss$age,2, labels = c("low", "high")) 
bmi_discrete = cut(ss$bmi,2, labels = c("low", "high"))
days_survived_discrete = cut(ss$days_survived,2, labels = c("low", "high"))
slide_area_discrete = cut(ss$slide_area,2, labels = c("low", "high"))

ss$age_discrete = age_discrete
ss$bmi_discrete = bmi_discrete
ss$days_survived_discrete = days_survived_discrete
ss$slide_area_discrete = slide_area_discrete

#Plot
ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$sample_group, group_title = "Sample Group")
save_plot(ggp,"ggp_sample_group.tiff",549,537)

ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$batch, group_title = "BATCH")
save_plot(ggp,"ggp_batch.tiff",549,537)

ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$age_discrete, group_title = "AGE")
save_plot(ggp,"ggp_age.tiff",549,537)

ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$bmi_discrete, group_title = "BMI")
save_plot(ggp,"ggp_bmi.tiff",549,537)

ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$days_survived_discrete, group_title = "Survived Days")
save_plot(ggp,"ggp_days_survived.tiff",549,537)

ggp = get_PCA_plot(pca_data = pca_coordinates,group = ss$slide_area_discrete, group_title = "Slide Size")
save_plot(ggp,"ggp_slide_area.tiff",549,537)
