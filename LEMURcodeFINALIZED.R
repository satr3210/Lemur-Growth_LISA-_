####lemur analysis -- updated version with good subset

#load needed packages 
library(dplyr)
library(ggplot2)
library(bootstrap)

#read in data 
getwd() #determine working directory 
setwd("/Users/sandratredinnick/Downloads/doi_10") #set working directory 
data_dlcweight= read.csv("DataRecord_3_DLC_Weight_File_05Feb2019.csv")
data_dlclist=read.csv("DataRecord_2_DLC_Animal_List_05Feb2019.csv")
data_dlctable=read.csv("DataRecord_1b_DLC_LH_Table_Analysis_05Feb2019.csv")


##filtering data set 
# just want lemurs
LCATdata = data_dlcweight %>% filter(Taxon == "LCAT")

#want to see which lemurs we have weight when they are 2 weeks or younger
# and 5 years or older

# first determine lemurs with 20+ observations 
lemurnames = table(LCATdata$Name) # creating table with just lemur names
# from the table lemurnames, select the name of the lemur that occurs more than 20 times
# finds rows of 'Name' in LCATdata that matches the > 20 names
# result is a logical vector that can be used to subset the data frame LCATdata
# subset data frame by a logical vector
# 'TRUE' rows are kept, 'FALSE' rows are removed.
# assign the result to df2
# link where code was found to do this: 
# https://stackoverflow.com/questions/20204257/subset-data-frame-based-on-number-of-rows-per-group
df2 = subset(LCATdata, Name %in% names(lemurnames[lemurnames > 20]))
length(unique(df2$Name)) #100 lemurs with 20+ observations 


# see how many lemurs have a birth weight (saying birth weight is first weight when they are 2 weeks old or younger)
df3 = df2 %>% filter(AgeAtWt_wk <= 2)
length(unique(df3$Name)) # 56 with weight recorded when they were 2 weeks old or younger

# see how many lemurs have weight recorded when they are 5 years or older 
df4 = df2 %>% filter(AgeAtWt_y >= 5)
length(unique(df4$Name)) # 70 with weight recorded after age 5

# find which lemurs have both age <= 2 weeks and age >= 5 years
youngerlemurnames = unique(df3$Name)
olderlemurnames = unique(df4$Name)
# Reduce takes a function f of two arguments and a list or vector x which is to be ‘reduced’ using f
# in this case, the function f is intersect and x is the list of lemur names that meet specific age criteria
# it finds the intersection (overlapping) names between the two lists and returns a a vector 
# of the names 
samenames = Reduce(intersect, list(youngerlemurnames, olderlemurnames)) # gives vector of 26 names
samenames

# filter df2 dataset with the lemurs from samenames 
LCATsub = df2 %>% filter(Name == "Ivy" | 
                           Name == "SAPPHO" |
                           Name == "Ginger" | 
                           Name == "Jones" | 
                           Name == "PHOTIUS" |
                           Name == "Alena" |
                           Name == "Herodotus" |
                           Name == "Dorieus" |
                           Name == "Stewart" |
                           Name == "Alexander" |
                           Name == "Cebes" |
                           Name == "Gretl" |
                           Name == "Artemesia" |
                           Name == "Justine" |
                           Name == "Fern" |
                           Name == "Nikos" |
                           Name == "Nicaea" |
                           Name == "Lilah" |
                           Name == "Sierra Mist" |
                           Name == "Onyx" |
                           Name == "Tellus" |
                           Name == "Persephone" |
                           Name == "Hibernia" |
                           Name == "Berisades" |
                           Name == "Sophia" |
                           Name == "Brigitta")

LCATsub= LCATsub %>% arrange(Name)
  
  

#sanity check to make sure there are 26 individual lemurs are in data set
length(unique(LCATsub$Name)) # 26!

# double check that there are at least 20 observations for each lemur 
observationcount = LCATsub %>% count(Name) # all have more than 20! 

#looking at summary statistics of age at weight in years for subset datset 
summary(LCATsub$AgeAtWt_y) # mean is 5.593!

#statistics from the whole data set
summary(LCATdata$AgeAtWt_y) # mean of whole LCAT data set is 8 years old 

#doesnt appear to be horrible variation in subset and whole lemur dataset!


#plot subset -LCATsub
qplot(AgeAtWt_y, Weight_g, colour = Name, data = LCATsub,
      xlab = "Age in Year", ylab = "Weight in Grams",
      main = "Weight vs Age for Subset Dataset of Lemurs")




########lets do some modified universal growth curve modeling

##need birth weight - youngest age recorded (<= 2 weeks)
##need average adult weight 
##need a parameters, tau parameters, and r parameters 

#to do universal fitting we will fit each individual lemur in the subset seperately (by subject)
#and then average the coefficients in each subject based modified curve to get one general (represenative)
#modifid growth curve for the lemur population

#now that we have the subset datset of lemurs, lets do modeling!
###store needed parameters for modified universal modeling 
a_BMR = 45.1
a_RMR = 651.4

############first model SAPPHO!!!!

#filter subset for just SAPPHO
dfSAPPHO=LCATsub%>%filter(Name=="SAPPHO")

#birth weight for SHAPPO
birthweightSAPPHO=min(dfSAPPHO$Weight_g)

#average adult weight 
#get df when SAPPHO is an adult
adultdfsapp= dfSAPPHO%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightsapp= mean(adultdfsapp$Weight_g)

#contruct tau and r columns for sappho
dfSAPPHO$r_sapp = (dfSAPPHO$Weight_g/avgadweightsapp)^0.25
dfSAPPHO$tauBMR_sapp=((a_BMR*dfSAPPHO$AgeAtWt_y)/(4*avgadweightsapp^0.25))-log(1-(birthweightSAPPHO/avgadweightsapp)^0.25)
dfSAPPHO$tauRMR_sapp=((a_RMR*dfSAPPHO$AgeAtWt_y)/(4*avgadweightsapp^0.25))-log(1-(birthweightSAPPHO/avgadweightsapp)^0.25)

#modeling sappho
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestsapp_BMR=lm(r_sapp~log(tauBMR_sapp), data=dfSAPPHO)
summary(betaestsapp_BMR)
betasapp_BMR_sappho=0.12222 

modelsapp_BMR= function(tau){1-exp(-betasapp_BMR_sappho*tau)}
plot(dfSAPPHO$tauBMR_sapp, dfSAPPHO$r_sapp, xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight", main="Modified Universal Curve For SAPPHO") 
curve(modelsapp_BMR,0, 13, add=TRUE) #LSE didnt do well, so estiate by eye
#estimate betas by eye 
betasapp_BMR2=1.0
modelsapp_BMR2= function(tau){1-exp(-betasapp_BMR2*tau)}
plot(dfSAPPHO$tauBMR_sapp, dfSAPPHO$r_sapp, xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight", main="Modified Universal Curve For SAPPHO") 
curve(modelsapp_BMR2,0, 13, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestsapp_RMR=lm(r_sapp~log(tauRMR_sapp), data=dfSAPPHO)
summary(betaestsapp_RMR)
betasapp_RMR=0.098210 
###the right way to plot it! using a_RMR
modelsapp_RMR= function(tau){1-exp(-betasapp_RMR*tau)}
plot(dfSAPPHO$tauRMR_sapp, dfSAPPHO$r_sapp,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For SAPPHO") 
curve(modelsapp_RMR,0, 178, add=TRUE)



############second model Jones!!!!

#filter subset for just Jones
dfjones=LCATsub%>%filter(Name=="Jones")

#birth weight for SHAPPO
birthweightjones=min(dfjones$Weight_g)

#average adult weight 
#get df when SAPPHO is an adult
adultdfjones= dfjones%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightjones= mean(adultdfjones$Weight_g)

#contruct tau and r columns for sappho
dfjones$r_jones = (dfjones$Weight_g/avgadweightjones)^0.25
dfjones$tauBMR_jones=((a_BMR*dfjones$AgeAtWt_y)/(4*avgadweightjones^0.25))-log(1-(birthweightjones/avgadweightjones)^0.25)
dfjones$tauRMR_jones=((a_RMR*dfjones$AgeAtWt_y)/(4*avgadweightjones^0.25))-log(1-(birthweightjones/avgadweightjones)^0.25)

#modeling jones
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestjones_BMR=lm(r_jones~log(tauBMR_jones), data=dfjones)
summary(betaestjones_BMR)
betajones_BMR=0.168522 
###the right way to plot it! using a_BMR
modeljones_BMR= function(tau){1-exp(-betajones_BMR*tau)}
plot(dfjones$tauBMR_jones, dfjones$r_jones,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Jones") 
curve(modeljones_BMR,0, 14, add=TRUE) 
#estimate betas by eye 
betajones_BMR2=.9
modeljones_BMR2= function(tau){1-exp(-betajones_BMR2*tau)}
plot(dfjones$tauBMR_jones, dfjones$r_jones,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Jones") 
curve(modeljones_BMR2,0, 14, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestjones_RMR=lm(r_jones~log(tauRMR_jones), data=dfjones)
summary(betaestjones_RMR)
betajones_RMR=0.116682
###the right way to plot it! using a_RMR
modeljones_RMR= function(tau){1-exp(-betajones_RMR*tau)}
plot(dfjones$tauRMR_jones, dfjones$r_jones,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Jones") 
curve(modeljones_RMR,0, 175, add=TRUE)




############third model Justine!!!!

#filter subset for just Justine
dfjustine=LCATsub%>%filter(Name=="Justine")

#birth weight for Justine
birthweightjustine=min(dfjustine$Weight_g)

#average adult weight 
#get df when justine is an adult
adultdfjustine= dfjustine%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightjustine= mean(adultdfjustine$Weight_g)

#contruct tau and r columns for justine
dfjustine$r_justine = (dfjustine$Weight_g/avgadweightjustine)^0.25
dfjustine$tauBMR_justine=((a_BMR*dfjustine$AgeAtWt_y)/(4*avgadweightjustine^0.25))-log(1-(birthweightjustine/avgadweightjustine)^0.25)
dfjustine$tauRMR_justine=((a_RMR*dfjustine$AgeAtWt_y)/(4*avgadweightjustine^0.25))-log(1-(birthweightjustine/avgadweightjustine)^0.25)

#modeling justine
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestjustine_BMR=lm(r_justine~log(tauBMR_justine), data=dfjustine)
summary(betaestjustine_BMR)
betajustine_BMR=0.068688
###the right way to plot it! using a_BMR
modeljustine_BMR= function(tau){1-exp(-betajustine_BMR*tau)}
plot(dfjustine$tauBMR_justine, dfjustine$r_justine,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Justine") 
curve(modeljones_BMR,0, 25, add=TRUE) 
#estimate betas by eye 
betajustine_BMR2=.6
modeljustine_BMR2= function(tau){1-exp(-betajustine_BMR2*tau)}
plot(dfjustine$tauBMR_justine, dfjustine$r_justine,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Justine") 
curve(modeljustine_BMR2,0, 25, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestjustine_RMR=lm(r_justine~log(tauRMR_justine), data=dfjustine)
summary(betaestjustine_RMR)
betajustine_RMR=0.062486
###the right way to plot it! using a_RMR
modeljustine_RMR= function(tau){1-exp(-betajustine_RMR*tau)}
plot(dfjustine$tauRMR_justine, dfjustine$r_justine,
     xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Justine") 
curve(modeljustine_RMR,0, 340, add=TRUE)




############forth model Lilah!!!!

#filter subset for just lilah
dflilah=LCATsub%>%filter(Name=="Lilah")

#birth weight for lilah
birthweightlilah=min(dflilah$Weight_g)

#average adult weight 
#get df when lilah is an adult
adultdflilah= dflilah%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightlilah= mean(adultdflilah$Weight_g)

#contruct tau and r columns for lilah
dflilah$r_lilah = (dflilah$Weight_g/avgadweightlilah)^0.25
dflilah$tauBMR_lilah=((a_BMR*dflilah$AgeAtWt_y)/(4*avgadweightlilah^0.25))-log(1-(birthweightlilah/avgadweightlilah)^0.25)
dflilah$tauRMR_lilah=((a_RMR*dflilah$AgeAtWt_y)/(4*avgadweightlilah^0.25))-log(1-(birthweightlilah/avgadweightlilah)^0.25)

#modeling lilah
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestlilah_BMR=lm(r_lilah~log(tauBMR_lilah), data=dflilah)
summary(betaestlilah_BMR)
betalilah_BMR=0.085059
###the right way to plot it! using a_BMR
modellilah_BMR= function(tau){1-exp(-betalilah_BMR*tau)}
plot(dflilah$tauBMR_lilah, dflilah$r_lilah,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Lilah") 
curve(modellilah_BMR,0, 25, add=TRUE) 
#estimate betas by eye 
betalilah_BMR2=.7
modellilah_BMR2= function(tau){1-exp(-betalilah_BMR2*tau)}
plot(dflilah$tauBMR_lilah, dflilah$r_lilah,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Lilah") 
curve(modellilah_BMR2,0, 25, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestlilah_RMR=lm(r_lilah~log(tauRMR_lilah), data=dflilah)
summary(betaestlilah_RMR)
betalilah_RMR=0.073023
###the right way to plot it! using a_RMR
modellilah_RMR= function(tau){1-exp(-betalilah_RMR*tau)}
plot(dflilah$tauRMR_lilah, dflilah$r_lilah,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Lilah") 
curve(modellilah_RMR,0, 333, add=TRUE)







############five model Nicaea!!!!

#filter subset for just Nicaea
dfnicaea=LCATsub%>%filter(Name=="Nicaea")

#birth weight for nicaea
birthweightnicaea=min(dfnicaea$Weight_g)

#average adult weight 
#get df when nicaea is an adult
adultdfnicaea= dfnicaea%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightnicaea= mean(adultdfnicaea$Weight_g)

#contruct tau and r columns for nicaea
dfnicaea$r_nicaea = (dfnicaea$Weight_g/avgadweightnicaea)^0.25
dfnicaea$tauBMR_nicaea=((a_BMR*dfnicaea$AgeAtWt_y)/(4*avgadweightnicaea^0.25))-log(1-(birthweightnicaea/avgadweightnicaea)^0.25)
dfnicaea$tauRMR_nicaea=((a_RMR*dfnicaea$AgeAtWt_y)/(4*avgadweightnicaea^0.25))-log(1-(birthweightnicaea/avgadweightnicaea)^0.25)

#modeling nicaea
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestnicaea_BMR=lm(r_nicaea~log(tauBMR_nicaea), data=dfnicaea)
summary(betaestnicaea_BMR)
betanicaea_BMR=0.163964
###the right way to plot it! using a_BMR
modelnicaea_BMR= function(tau){1-exp(-betanicaea_BMR*tau)}
plot(dfnicaea$tauBMR_nicaea, dfnicaea$r_nicaea,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nicaea") 
curve(modelnicaea_BMR,0, 15, add=TRUE) 
#estimate betas by eye 
betanicaea_BMR2=.7
modelnicaea_BMR2= function(tau){1-exp(-betanicaea_BMR2*tau)}
plot(dfnicaea$tauBMR_nicaea, dfnicaea$r_nicaea,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nicaea") 
curve(modelnicaea_BMR2,0, 15, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestnicaea_RMR=lm(r_nicaea~log(tauRMR_nicaea), data=dfnicaea)
summary(betaestnicaea_RMR)
betanicaea_RMR=0.117327
###the right way to plot it! using a_RMR
modelnicaea_RMR= function(tau){1-exp(-betanicaea_RMR*tau)}
plot(dfnicaea$tauRMR_nicaea, dfnicaea$r_nicaea,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nicaea") 
curve(modelnicaea_RMR,0, 185, add=TRUE)






############sixth model Nikos!!!!

#filter subset for just Nikos
dfnikos=LCATsub%>%filter(Name=="Nikos")

#birth weight for nicaea
birthweightnikos=min(dfnikos$Weight_g)

#average adult weight 
#get df when nicaea is an adult
adultdfnikos= dfnikos%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightnikos= mean(adultdfnikos$Weight_g)

#contruct tau and r columns for nikos
dfnikos$r_nikos = (dfnikos$Weight_g/avgadweightnikos)^0.25
dfnikos$tauBMR_nikos=((a_BMR*dfnikos$AgeAtWt_y)/(4*avgadweightnikos^0.25))-log(1-(birthweightnikos/avgadweightnikos)^0.25)
dfnikos$tauRMR_nikos=((a_RMR*dfnikos$AgeAtWt_y)/(4*avgadweightnikos^0.25))-log(1-(birthweightnikos/avgadweightnikos)^0.25)

#modeling nikos
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestnikos_BMR=lm(r_nikos~log(tauBMR_nikos), data=dfnikos)
summary(betaestnikos_BMR)
betanikos_BMR=0.17898
###the right way to plot it! using a_BMR
modelnikos_BMR= function(tau){1-exp(-betanikos_BMR*tau)}
plot(dfnikos$tauBMR_nikos, dfnikos$r_nikos,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nikos") 
curve(modelnikos_BMR,0, 11, add=TRUE) 
#estimate betas by eye 
betanikos_BMR2=.98
modelnikos_BMR2= function(tau){1-exp(-betanikos_BMR2*tau)}
plot(dfnikos$tauBMR_nikos, dfnikos$r_nikos,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nikos") 
curve(modelnikos_BMR2,0, 11, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestnikos_RMR=lm(r_nikos~log(tauRMR_nikos), data=dfnikos)
summary(betaestnikos_RMR)
betanikos_RMR=0.122424
###the right way to plot it! using a_RMR
modelnikos_RMR= function(tau){1-exp(-betanikos_RMR*tau)}
plot(dfnikos$tauRMR_nikos, dfnikos$r_nikos,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Nikos") 
curve(modelnikos_RMR,0, 137, add=TRUE)





############seventh model Onyx!!!!

#filter subset for just Onyx
dfonyx=LCATsub%>%filter(Name=="Onyx")

#birth weight for onyx
birthweightonyx=min(dfonyx$Weight_g)

#average adult weight 
#get df when onyx is an adult
adultdfonyx= dfonyx%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightonyx= mean(adultdfonyx$Weight_g)

#contruct tau and r columns for onyx
dfonyx$r_onyx = (dfonyx$Weight_g/avgadweightonyx)^0.25
dfonyx$tauBMR_onyx=((a_BMR*dfonyx$AgeAtWt_y)/(4*avgadweightonyx^0.25))-log(1-(birthweightonyx/avgadweightonyx)^0.25)
dfonyx$tauRMR_onyx=((a_RMR*dfonyx$AgeAtWt_y)/(4*avgadweightonyx^0.25))-log(1-(birthweightonyx/avgadweightonyx)^0.25)

#modeling onyx
#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestonyx_BMR=lm(r_onyx~log(tauBMR_onyx), data=dfonyx)
summary(betaestonyx_BMR)
betaonyx_BMR=0.19727
###the right way to plot it! using a_BMR
modelonyx_BMR= function(tau){1-exp(-betaonyx_BMR*tau)}
plot(dfonyx$tauBMR_onyx, dfonyx$r_onyx,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Onyx") 
curve(modelonyx_BMR,0, 11, add=TRUE) 
#estimate betas by eye 
betaonyx_BMR2=.98
modelonyx_BMR2= function(tau){1-exp(-betaonyx_BMR2*tau)}
plot(dfonyx$tauBMR_onyx, dfonyx$r_onyx,xlab="Dimensionless Time using BMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Onyx") 
curve(modelonyx_BMR2,0, 11, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestonyx_RMR=lm(r_onyx~log(tauRMR_onyx), data=dfonyx)
summary(betaestonyx_RMR)
betaonyx_RMR=0.129018
###the right way to plot it! using a_RMR
modelonyx_RMR= function(tau){1-exp(-betaonyx_RMR*tau)}
plot(dfonyx$tauRMR_onyx, dfonyx$r_onyx,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Onyx") 
curve(modelonyx_RMR,0, 130, add=TRUE)





############eigith model Persephone!!!!

#filter subset for just Persephone
dfpersephone=LCATsub%>%filter(Name=="Persephone")

#birth weight for persephone
birthweightpersephone=min(dfpersephone$Weight_g)

#average adult weight 
#get df when persephone is an adult
adultdfpersephone= dfpersephone%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightpersephone= mean(adultdfpersephone$Weight_g)

#contruct tau and r columns for persephone
dfpersephone$r_persephone = (dfpersephone$Weight_g/avgadweightpersephone)^0.25
dfpersephone$tauBMR_persephone=((a_BMR*dfpersephone$AgeAtWt_y)/(4*avgadweightpersephone^0.25))-log(1-(birthweightpersephone/avgadweightpersephone)^0.25)
dfpersephone$tauRMR_persephone=((a_RMR*dfpersephone$AgeAtWt_y)/(4*avgadweightpersephone^0.25))-log(1-(birthweightpersephone/avgadweightpersephone)^0.25)

#modeling persephone
#model = (r=1-e^-beta*tau)

####BMR case persephone
#estimate betas for BMR
betaestpersephone_BMR=lm(r_persephone~log(tauBMR_persephone), data=dfpersephone)
summary(betaestpersephone_BMR)
betapersephone_BMR=0.09320
###the right way to plot it! using a_BMR
modelpersephone_BMR= function(tau){1-exp(-betapersephone_BMR*tau)}
plot(dfpersephone$tauBMR_persephone, dfpersephone$r_persephone,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For Persephone") 
curve(modelpersephone_BMR,0, 25, add=TRUE) 
#estimate betas by eye 
betapersephone_BMR2=.98
  modelpersephone_BMR2= function(tau){1-exp(-betapersephone_BMR2*tau)}
  plot(dfpersephone$tauBMR_persephone, dfpersephone$r_persephone,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Persephone") 
  curve(modelpersephone_BMR2,0, 25, add=TRUE) #LSE didnt do well, so estiate by eye
  
  
  
  ###RMR case 
  #estimate betas for RMR
  betaestpersephone_RMR=lm(r_persephone~log(tauRMR_persephone), data=dfpersephone)
  summary(betaestpersephone_RMR)
  betapersephone_RMR=0.078228
  ###the right way to plot it! using a_RMR
  modelpersephone_RMR= function(tau){1-exp(-betapersephone_RMR*tau)}
  plot(dfpersephone$tauRMR_persephone, dfpersephone$r_persephone,
       xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Persephone") 
  curve(modelpersephone_RMR,0, 325, add=TRUE)
  
  
  
  
  
  
  ############ninth model PHOTIUS!!!!
  
  #filter subset for just PHOTIUS
  dfPHOTIUS=LCATsub%>%filter(Name=="PHOTIUS")
  
  #birth weight for PHOTIUS
  birthweightPHOTIUS=min(dfPHOTIUS$Weight_g)
  
  #average adult weight 
  #get df when PHOTIUS is an adult
  adultdfPHOTIUS= dfPHOTIUS%>%filter(Age_Category=="adult")
  #average the adult weights
  avgadweightPHOTIUS= mean(adultdfPHOTIUS$Weight_g)
  
  #contruct tau and r columns for PHOTIUS
  dfPHOTIUS$r_PHOTIUS = (dfPHOTIUS$Weight_g/avgadweightPHOTIUS)^0.25
  dfPHOTIUS$tauBMR_PHOTIUS=((a_BMR*dfPHOTIUS$AgeAtWt_y)/(4*avgadweightPHOTIUS^0.25))-log(1-(birthweightPHOTIUS/avgadweightPHOTIUS)^0.25)
  dfPHOTIUS$tauRMR_PHOTIUS=((a_RMR*dfPHOTIUS$AgeAtWt_y)/(4*avgadweightPHOTIUS^0.25))-log(1-(birthweightPHOTIUS/avgadweightPHOTIUS)^0.25)
  
  #modeling PHOTIUS
  #model = (r=1-e^-beta*tau)
  
  ####BMR case PHOTIUS
  #estimate betas for BMR
  betaestPHOTIUS_BMR=lm(r_PHOTIUS~log(tauBMR_PHOTIUS), data=dfPHOTIUS)
  summary(betaestPHOTIUS_BMR)
  betaPHOTIUS_BMR=0.12588
  ###the right way to plot it! using a_BMR
  modelPHOTIUS_BMR= function(tau){1-exp(-betaPHOTIUS_BMR*tau)}
  plot(dfPHOTIUS$tauBMR_PHOTIUS, dfPHOTIUS$r_PHOTIUS,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For PHOTIUS") 
  curve(modelPHOTIUS_BMR,0, 13, add=TRUE) 
  #estimate betas by eye 
  betaPHOTIUS_BMR2=1.01
  modelPHOTIUS_BMR2= function(tau){1-exp(-betaPHOTIUS_BMR2*tau)}
  plot(dfPHOTIUS$tauBMR_PHOTIUS, dfPHOTIUS$r_PHOTIUS,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For PHOTIUS") 
  curve(modelPHOTIUS_BMR2,0, 13, add=TRUE) #LSE didnt do well, so estiate by eye
  
  
  
  ###RMR case 
  #estimate betas for RMR
  betaestPHOTIUS_RMR=lm(r_PHOTIUS~log(tauRMR_PHOTIUS), data=dfPHOTIUS)
  summary(betaestPHOTIUS_RMR)
  betaPHOTIUS_RMR=0.097900
  ###the right way to plot it! using a_RMR
  modelPHOTIUS_RMR= function(tau){1-exp(-betaPHOTIUS_RMR*tau)}
  plot(dfPHOTIUS$tauRMR_PHOTIUS, dfPHOTIUS$r_PHOTIUS,
       xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For PHOTIUS") 
  curve(modelPHOTIUS_RMR,0, 170, add=TRUE)
  
  
  
  
  
  
  ############tenth  model Sierra Mist!!!!
  
  #filter subset for just Sierra Mist
  dfSierraMist=LCATsub%>%filter(Name=="Sierra Mist")
  
  #birth weight for PHOTIUS
  birthweightSierraMist=min(dfSierraMist$Weight_g)
  
  #average adult weight 
  #get df when SierraMist is an adult
  adultdfSierraMist= dfSierraMist%>%filter(Age_Category=="adult")
  #average the adult weights
  avgadweightSierraMist= mean(adultdfSierraMist$Weight_g)
  
  #contruct tau and r columns for SierraMist
  dfSierraMist$r_SierraMist = (dfSierraMist$Weight_g/avgadweightSierraMist)^0.25
  dfSierraMist$tauBMR_SierraMist=((a_BMR*dfSierraMist$AgeAtWt_y)/(4*avgadweightSierraMist^0.25))-log(1-(birthweightSierraMist/avgadweightSierraMist)^0.25)
  dfSierraMist$tauRMR_SierraMist=((a_RMR*dfSierraMist$AgeAtWt_y)/(4*avgadweightSierraMist^0.25))-log(1-(birthweightSierraMist/avgadweightSierraMist)^0.25)
  
  #modeling SierraMist
  #model = (r=1-e^-beta*tau)
  
  ####BMR case SierraMist
  #estimate betas for BMR
  betaestSierraMist_BMR=lm(r_SierraMist~log(tauBMR_SierraMist), data=dfSierraMist)
  summary(betaestSierraMist_BMR)
  betaSierraMist_BMR=0.131487
  ###the right way to plot it! using a_BMR
  modelSierraMist_BMR= function(tau){1-exp(-betaSierraMist_BMR*tau)}
  plot(dfSierraMist$tauBMR_SierraMist, dfSierraMist$r_SierraMist,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Sierra Mist") 
  curve(modelSierraMist_BMR,0, 20, add=TRUE) 
  #estimate betas by eye 
  betaSierraMist_BMR2=1.01
  modelSierraMist_BMR2= function(tau){1-exp(-betaSierraMist_BMR2*tau)}
  plot(dfSierraMist$tauBMR_SierraMist, dfSierraMist$r_SierraMist,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Sierra Mist") 
  curve(modelSierraMist_BMR2,0, 20, add=TRUE) #LSE didnt do well, so estiate by eye
  
  
  
  ###RMR case 
  #estimate betas for RMR
  betaestSierraMist_RMR=lm(r_SierraMist~log(tauRMR_SierraMist), data=dfSierraMist)
  summary(betaestSierraMist_RMR)
  betaSierraMist_RMR=0.100012
  ###the right way to plot it! using a_RMR
  modelSierraMist_RMR= function(tau){1-exp(-betaSierraMist_RMR*tau)}
  plot(dfSierraMist$tauRMR_SierraMist, dfSierraMist$r_SierraMist,
       xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Sierra Mist") 
  curve(modelSierraMist_RMR,0, 255, add=TRUE)
  
  
  
  
  
  ############eleventh  model Sophia!!!!
  
  #filter subset for just Sophia
  dfSophia=LCATsub%>%filter(Name=="Sophia")
  
  #birth weight for Sophia
  birthweightSophia=min(dfSophia$Weight_g)
  
  #average adult weight 
  #get df when Sophia is an adult
  adultdfSophia= dfSophia%>%filter(Age_Category=="adult")
  #average the adult weights
  avgadweightSophia= mean(adultdfSophia$Weight_g)
  
  #contruct tau and r columns for Sophia
  dfSophia$r_Sophia = (dfSophia$Weight_g/avgadweightSophia)^0.25
  dfSophia$tauBMR_Sophia=((a_BMR*dfSophia$AgeAtWt_y)/(4*avgadweightSophia^0.25))-log(1-(birthweightSophia/avgadweightSophia)^0.25)
  dfSophia$tauRMR_Sophia=((a_RMR*dfSophia$AgeAtWt_y)/(4*avgadweightSophia^0.25))-log(1-(birthweightSophia/avgadweightSophia)^0.25)
  
  #modeling Sophia
  #model = (r=1-e^-beta*tau)
  
  ####BMR case Sophia
  #estimate betas for BMR
  betaestSophia_BMR=lm(r_Sophia~log(tauBMR_Sophia), data=dfSophia)
  summary(betaestSophia_BMR)
  betaSophia_BMR=0.074083
  ###the right way to plot it! using a_BMR
  modelSophia_BMR= function(tau){1-exp(-betaSophia_BMR*tau)}
  plot(dfSophia$tauBMR_Sophia, dfSophia$r_Sophia,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Sophia") 
  curve(modelSophia_BMR,0, 25, add=TRUE) 
  #estimate betas by eye 
  betaSophia_BMR2=0.99
  modelSophia_BMR2= function(tau){1-exp(-betaSophia_BMR2*tau)}
  plot(dfSophia$tauBMR_Sophia, dfSophia$r_Sophia,
       xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
       main="Modified Universal Curve For Sophia") 
  curve(modelSierraMist_BMR2,0, 25, add=TRUE) #LSE didnt do well, so estiate by eye
  
  
  
  ###RMR case 
  #estimate betas for RMR
  betaestSophia_RMR=lm(r_Sophia~log(tauRMR_Sophia), data=dfSophia)
  summary(betaestSophia_RMR)
  betaSophia_RMR=0.067290
  ###the right way to plot it! using a_RMR
modelSophia_RMR= function(tau){1-exp(-betaSophia_RMR*tau)}
plot(dfSophia$tauRMR_Sophia, dfSophia$r_Sophia,xlab="Dimensionless Time using RMR",
     ylab="Dimensionless Weight",main="Modified Universal Curve For Sophia") 
curve(modelSophia_RMR,0, 345, add=TRUE)
  



############twelfth model Stewart!!!!

#filter subset for just Stewart
dfStewart=LCATsub%>%filter(Name=="Stewart")

#birth weight for Stewart
birthweightStewart=min(dfStewart$Weight_g)

#average adult weight 
#get df when Stewart is an adult
adultdfStewart= dfStewart%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightStewart= mean(adultdfStewart$Weight_g)

#contruct tau and r columns for Stewart
dfStewart$r_Stewart = (dfStewart$Weight_g/avgadweightStewart)^0.25
dfStewart$tauBMR_Stewart=((a_BMR*dfStewart$AgeAtWt_y)/(4*avgadweightStewart^0.25))-log(1-(birthweightStewart/avgadweightStewart)^0.25)
dfStewart$tauRMR_Stewart=((a_RMR*dfStewart$AgeAtWt_y)/(4*avgadweightStewart^0.25))-log(1-(birthweightStewart/avgadweightStewart)^0.25)

#modeling Stewart
#model = (r=1-e^-beta*tau)

####BMR case Stewart
#estimate betas for BMR
betaestStewart_BMR=lm(r_Stewart~log(tauBMR_Stewart), data=dfStewart)
summary(betaestStewart_BMR)
betaStewart_BMR=0.163813
###the right way to plot it! using a_BMR
modelStewart_BMR= function(tau){1-exp(-betaStewart_BMR*tau)}
plot(dfStewart$tauBMR_Stewart, dfStewart$r_Stewart,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Stewart") 
curve(modelStewart_BMR,0, 14, add=TRUE) 
#estimate betas by eye 
betaStewart_BMR2=1.01
modelStewart_BMR2= function(tau){1-exp(-betaStewart_BMR2*tau)}
plot(dfStewart$tauBMR_Stewart, dfStewart$r_Stewart,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Stewart") 
curve(modelStewart_BMR2,0, 14, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestStewart_RMR=lm(r_Stewart~log(tauRMR_Stewart), data=dfStewart)
summary(betaestStewart_RMR)
betaStewart_RMR=0.113112
###the right way to plot it! using a_RMR
modelStewart_RMR= function(tau){1-exp(-betaStewart_RMR*tau)}
plot(dfStewart$tauRMR_Stewart, dfStewart$r_Stewart,
     xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Stewart") 
curve(modelStewart_RMR,0, 345, add=TRUE)





############thirfteenth model Tellus!!!!

#filter subset for just Tellus
dfTellus=LCATsub%>%filter(Name=="Tellus")

#birth weight for Tellus
birthweightTellus=min(dfTellus$Weight_g)

#average adult weight 
#get df when Tellus is an adult
adultdfTellus= dfTellus%>%filter(Age_Category=="adult")
#average the adult weights
avgadweightTellus= mean(adultdfTellus$Weight_g)

#contruct tau and r columns for Tellus
dfTellus$r_Tellus = (dfTellus$Weight_g/avgadweightTellus)^0.25
dfTellus$tauBMR_Tellus=((a_BMR*dfTellus$AgeAtWt_y)/(4*avgadweightTellus^0.25))-log(1-(birthweightTellus/avgadweightTellus)^0.25)
dfTellus$tauRMR_Tellus=((a_RMR*dfTellus$AgeAtWt_y)/(4*avgadweightTellus^0.25))-log(1-(birthweightTellus/avgadweightTellus)^0.25)

#modeling Tellus
#model = (r=1-e^-beta*tau)

####BMR case Tellus
#estimate betas for BMR
betaestTellus_BMR=lm(r_Tellus~log(tauBMR_Tellus), data=dfTellus)
summary(betaestTellus_BMR)
betaTellus_BMR=0.116262
###the right way to plot it! using a_BMR
modelTellus_BMR= function(tau){1-exp(-betaTellus_BMR*tau)}
plot(dfTellus$tauBMR_Tellus, dfTellus$r_Tellus,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Tellus") 
curve(modelTellus_BMR,0, 25, add=TRUE) 
#estimate betas by eye 
betaTellus_BMR2=1.01
modelTellus_BMR2= function(tau){1-exp(-betaTellus_BMR2*tau)}
plot(dfTellus$tauBMR_Tellus, dfTellus$r_Tellus,
     xlab="Dimensionless Time using BMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Tellus") 
curve(modelTellus_BMR2,0, 25, add=TRUE) #LSE didnt do well, so estiate by eye



###RMR case 
#estimate betas for RMR
betaestTellus_RMR=lm(r_Tellus~log(tauRMR_Tellus), data=dfTellus)
summary(betaestTellus_RMR)
betaTellus_RMR=0.094058
###the right way to plot it! using a_RMR
modelTellus_RMR= function(tau){1-exp(-betaTellus_RMR*tau)}
plot(dfTellus$tauRMR_Tellus, dfTellus$r_Tellus,
     xlab="Dimensionless Time using RMR",ylab="Dimensionless Weight",
     main="Modified Universal Curve For Tellus") 
curve(modelTellus_RMR,0, 345, add=TRUE)


# Alena!

#filter subset for just Alena
dfAlena = LCATsub %>% filter(Name =="Alena")

#birth weight 
birthweightAlena = min(dfAlena$Weight_g)

#average adult weight 
#get df when Alena is an adult
adultdfAlena = dfAlena %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightAlena = mean(adultdfAlena$Weight_g)

#contruct tau and r columns 
dfAlena$r = (dfAlena$Weight_g/avgadweightAlena)^0.25
dfAlena$tauBMR = ((a_BMR*dfAlena$AgeAtWt_y)/(4*avgadweightAlena^0.25))-log(1-(birthweightAlena/avgadweightAlena)^0.25)
dfAlena$tauRMR=  ((a_RMR*dfAlena$AgeAtWt_y)/(4*avgadweightAlena^0.25))-log(1-(birthweightAlena/avgadweightAlena)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestAlena_BMR = lm(r ~ log(tauBMR), data = dfAlena)
summary(betaestAlena_BMR)
betaAlena_BMR = 0.106645 

###the right way to plot it! using a_BMR
modelAlena_BMR = function(tau){1-exp(-betaAlena_BMR*tau)}
plot(dfAlena$tauBMR, dfAlena$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alena") 
curve(modelAlena_BMR, 0, 25, add = TRUE) #didnt do well, so estimate by eye

#estimate beta by eye 
betaAlena_BMR2 = 1.3
modelAlena_BMR2 = function(tau){1-exp(-betaAlena_BMR2*tau)}
plot(dfAlena$tauBMR, dfAlena$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alena") 
curve(modelAlena_BMR2, 0, 25, add = TRUE)



###RMR case 
#estimate betas for RMR
betaestAlena_RMR = lm(r~log(tauRMR), data = dfAlena)
summary(betaestAlena_RMR)
betaAlena_RMR = 0.086339 
###the right way to plot it! using a_RMR
modelAlena_RMR= function(tau){1-exp(-betaAlena_RMR*tau)}
plot(dfAlena$tauRMR, dfAlena$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alena") 
curve(modelAlena_RMR, 0, 327, add=TRUE)



# Alexander!

#filter subset for just Alexander
dfAlex = LCATsub %>% filter(Name == "Alexander")

#birth weight 
birthweightAlex = min(dfAlex$Weight_g)

#average adult weight 
#get df when Alex is an adult
adultdfAlex= dfAlex %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightAlex = mean(adultdfAlex$Weight_g)

#contruct tau and r columns 
dfAlex$r = (dfAlex$Weight_g/avgadweightAlex)^0.25
dfAlex$tauBMR = ((a_BMR*dfAlex$AgeAtWt_y)/(4*avgadweightAlex^0.25))-log(1-(birthweightAlex/avgadweightAlex)^0.25)
dfAlex$tauRMR=  ((a_RMR*dfAlex$AgeAtWt_y)/(4*avgadweightAlex^0.25))-log(1-(birthweightAlex/avgadweightAlex)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestAlex_BMR = lm(r ~ log(tauBMR), data = dfAlex)
summary(betaestAlex_BMR)
betaAlex_BMR = 0.145650  

###the right way to plot it! using a_BMR
modelAlex_BMR = function(tau){1-exp(-betaAlex_BMR*tau)}

plot(dfAlex$tauBMR, dfAlex$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alexander") 
curve(modelAlex_BMR, 0, 16, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaAlex_BMR2 = 1.2
modelAlex_BMR2 = function(tau){1-exp(-betaAlex_BMR2*tau)}
plot(dfAlex$tauBMR, dfAlex$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alexander") 
curve(modelAlex_BMR2,0, 16, add = TRUE)



###RMR case 
#estimate betas for RMR
betaestAlex_RMR = lm(r~log(tauRMR), data = dfAlex)
summary(betaestAlex_RMR)
betaAlex_RMR = 0.106267
###the right way to plot it! using a_RMR
modelAlex_RMR= function(tau){1-exp(-betaAlex_RMR*tau)}
plot(dfAlex$tauRMR, dfAlex$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alexander") 
curve(modelAlex_RMR,0, 200, add=TRUE)



# Artemesia!

#filter subset for just Artemesia
dfArt = LCATsub %>% filter(Name == "Artemesia")

#birth weight 
birthweightArt = min(dfArt$Weight_g)

#average adult weight 
#get df when Artemesia is an adult
adultdfArt= dfArt %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightArt = mean(adultdfArt$Weight_g)

#contruct tau and r columns 
dfArt$r = (dfArt$Weight_g/avgadweightArt)^0.25
dfArt$tauBMR = ((a_BMR*dfArt$AgeAtWt_y)/(4*avgadweightArt^0.25))-log(1-(birthweightArt/avgadweightArt)^0.25)
dfArt$tauRMR=  ((a_RMR*dfArt$AgeAtWt_y)/(4*avgadweightArt^0.25))-log(1-(birthweightArt/avgadweightArt)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestArt_BMR = lm(r ~ log(tauBMR), data = dfArt)
summary(betaestArt_BMR)
betaArt_BMR = 0.17313    

###the right way to plot it! using a_BMR
modelArt_BMR = function(tau){1-exp(-betaArt_BMR*tau)}

plot(dfArt$tauBMR, dfArt$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Artemesia") 
curve(modelArt_BMR, 0, 12, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaArt_BMR2 = 1.2
modelArt_BMR2 = function(tau){1-exp(-betaArt_BMR2*tau)}
plot(dfArt$tauBMR, dfArt$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Artemesia") 
curve(modelArt_BMR2,0, 12, add = TRUE)



###RMR case 
#estimate betas for RMR
betaestArt_RMR = lm(r~log(tauRMR), data = dfArt)
summary(betaestArt_RMR)
betaArt_RMR = 0.118137 
###the right way to plot it! using a_RMR
modelArt_RMR= function(tau){1-exp(-betaArt_RMR*tau)}
plot(dfArt$tauRMR, dfArt$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Artemesia") 
curve(modelArt_RMR, 0, 130, add=TRUE)




# Berisades!

#filter subset for just Berisades
dfBeri = LCATsub %>% filter(Name == "Berisades")

#birth weight 
birthweightBeri = min(dfBeri$Weight_g)

#average adult weight 
#get df when Berisades is an adult
adultdfBeri = dfBeri %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightBeri = mean(adultdfBeri$Weight_g)

#contruct tau and r columns 
dfBeri$r = (dfBeri$Weight_g/avgadweightBeri)^0.25
dfBeri$tauBMR = ((a_BMR*dfBeri$AgeAtWt_y)/(4*avgadweightBeri^0.25))-log(1-(birthweightBeri/avgadweightBeri)^0.25)
dfBeri$tauRMR=  ((a_RMR*dfBeri$AgeAtWt_y)/(4*avgadweightBeri^0.25))-log(1-(birthweightBeri/avgadweightBeri)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestBeri_BMR = lm(r ~ log(tauBMR), data = dfBeri)
#names(summary(betaestBeri_BMR))
betaBeri_BMR = summary(betaestBeri_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelBeri_BMR = function(tau){1-exp(-betaBeri_BMR*tau)}

plot(dfBeri$tauBMR, dfBeri$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Berisades") 
curve(modelBeri_BMR, 0, 27, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaBeri_BMR2 = 1
modelBeri_BMR2 = function(tau){1-exp(-betaBeri_BMR2*tau)}
plot(dfBeri$tauBMR, dfBeri$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Berisades") 
curve(modelBeri_BMR2,0, 27, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestBeri_RMR = lm(r~log(tauRMR), data = dfBeri)
summary(betaestBeri_RMR)
betaBeri_RMR = summary(betaestBeri_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelBeri_RMR= function(tau){1-exp(-betaBeri_RMR*tau)}
plot(dfBeri$tauRMR, dfBeri$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Berisades") 
curve(modelBeri_RMR, 0, 350, add=TRUE)




# Brigitta! 

#filter subset for just Brigitta
dfBrig = LCATsub %>% filter(Name == "Brigitta")

#birth weight 
birthweightBrig = min(dfBrig$Weight_g)

#average adult weight 
#get df when Brigitta is an adult
adultdfBrig = dfBrig %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightBrig = mean(adultdfBrig$Weight_g)

#contruct tau and r columns 
dfBrig$r = (dfBrig$Weight_g/avgadweightBrig)^0.25
dfBrig$tauBMR = ((a_BMR*dfBrig$AgeAtWt_y)/(4*avgadweightBrig^0.25))-log(1-(birthweightBrig/avgadweightBrig)^0.25)
dfBrig$tauRMR=  ((a_RMR*dfBrig$AgeAtWt_y)/(4*avgadweightBrig^0.25))-log(1-(birthweightBrig/avgadweightBrig)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestBrig_BMR = lm(r ~ log(tauBMR), data = dfBrig)
#names(summary(betaestBeri_BMR))
betaBrig_BMR = summary(betaestBrig_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelBrig_BMR = function(tau){1-exp(-betaBrig_BMR*tau)}

plot(dfBrig$tauBMR, dfBrig$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Brigitta") 
curve(modelBrig_BMR, 0, 15, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaBrig_BMR2 = 1.1
modelBrig_BMR2 = function(tau){1-exp(-betaBrig_BMR2*tau)}
plot(dfBrig$tauBMR, dfBrig$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Brigitta") 
curve(modelBrig_BMR2,0, 25, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestBrig_RMR = lm(r~log(tauRMR), data = dfBrig)
summary(betaestBrig_RMR)
betaBrig_RMR = summary(betaestBrig_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelBrig_RMR= function(tau){1-exp(-betaBrig_RMR*tau)}
plot(dfBrig$tauRMR, dfBrig$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Brigitta") 
curve(modelBrig_RMR, 0, 185, add=TRUE)




# Cebes!

#filter subset for just Cebes
dfCebes = LCATsub %>% filter(Name == "Cebes")

#birth weight 
birthweightCebes = min(dfCebes$Weight_g)

#average adult weight 
#get df when Cebes is an adult
adultdfCebes = dfCebes %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightCebes = mean(adultdfCebes$Weight_g)

#contruct tau and r columns 
dfCebes$r = (dfCebes$Weight_g/avgadweightCebes)^0.25
dfCebes$tauBMR = ((a_BMR*dfCebes$AgeAtWt_y)/(4*avgadweightCebes^0.25))-log(1-(birthweightCebes/avgadweightCebes)^0.25)
dfCebes$tauRMR=  ((a_RMR*dfCebes$AgeAtWt_y)/(4*avgadweightCebes^0.25))-log(1-(birthweightCebes/avgadweightCebes)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestCebes_BMR = lm(r ~ log(tauBMR), data = dfCebes)
#names(summary(betaestBeri_BMR))
betaCebes_BMR = summary(betaestCebes_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelCebes_BMR = function(tau){1-exp(-betaCebes_BMR*tau)}

plot(dfCebes$tauBMR, dfCebes$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Cebes") 
curve(modelCebes_BMR, 0, 12, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaCebes_BMR2 = 1.1
modelCebes_BMR2 = function(tau){1-exp(-betaCebes_BMR2*tau)}
plot(dfCebes$tauBMR, dfCebes$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Cebes") 
curve(modelCebes_BMR2,0, 12, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestCebes_RMR = lm(r~log(tauRMR), data = dfCebes)
summary(betaestCebes_RMR)
betaCebes_RMR = summary(betaestCebes_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelCebes_RMR= function(tau){1-exp(-betaCebes_RMR*tau)}
plot(dfCebes$tauRMR, dfCebes$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Cebes") 
curve(modelCebes_RMR, 0, 140, add = TRUE)




# Dorieus!

#filter subset for just Dorieus
dfDor = LCATsub %>% filter(Name == "Dorieus")

#birth weight 
birthweightDor = min(dfDor$Weight_g)

#average adult weight 
#get df when Dorieus is an adult
adultdfDor = dfDor %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightDor = mean(adultdfDor$Weight_g)

#contruct tau and r columns 
dfDor$r = (dfDor$Weight_g/avgadweightDor)^0.25
dfDor$tauBMR = ((a_BMR*dfDor$AgeAtWt_y)/(4*avgadweightDor^0.25))-log(1-(birthweightDor/avgadweightDor)^0.25)
dfDor$tauRMR=  ((a_RMR*dfDor$AgeAtWt_y)/(4*avgadweightDor^0.25))-log(1-(birthweightDor/avgadweightDor)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestDor_BMR = lm(r ~ log(tauBMR), data = dfDor)
betaDor_BMR = summary(betaestDor_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelDor_BMR = function(tau){1-exp(-betaDor_BMR*tau)}

plot(dfDor$tauBMR, dfDor$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Dorieus") 
curve(modelDor_BMR, 0, 33, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaDor_BMR2 = 0.9
modelDor_BMR2 = function(tau){1-exp(-betaDor_BMR2*tau)}
plot(dfDor$tauBMR, dfDor$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Dorieus") 
curve(modelDor_BMR2,0, 33, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestDor_RMR = lm(r~log(tauRMR), data = dfDor)
summary(betaestDor_RMR)
betaDor_RMR = summary(betaestDor_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelDor_RMR= function(tau){1-exp(-betaDor_RMR*tau)}
plot(dfDor$tauRMR, dfDor$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Dorieus") 
curve(modelDor_RMR, 0, 440, add = TRUE)




# Fern!

#filter subset for just Fern
dfFern = LCATsub %>% filter(Name == "Fern")

#birth weight 
birthweightFern = min(dfFern$Weight_g)

#average adult weight 
#get df when Fern is an adult
adultdfFern = dfFern %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightFern = mean(adultdfFern$Weight_g)

#contruct tau and r columns 
dfFern$r = (dfFern$Weight_g/avgadweightFern)^0.25
dfFern$tauBMR = ((a_BMR*dfFern$AgeAtWt_y)/(4*avgadweightFern^0.25))-log(1-(birthweightFern/avgadweightFern)^0.25)
dfFern$tauRMR=  ((a_RMR*dfFern$AgeAtWt_y)/(4*avgadweightFern^0.25))-log(1-(birthweightFern/avgadweightFern)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestFern_BMR = lm(r ~ log(tauBMR), data = dfFern)
betaFern_BMR = summary(betaestFern_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelFern_BMR = function(tau){1-exp(-betaFern_BMR*tau)}

plot(dfFern$tauBMR, dfFern$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Fern") 
curve(modelFern_BMR, 0, 13, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaFern_BMR2 = 1.1
modelFern_BMR2 = function(tau){1-exp(-betaFern_BMR2*tau)}
plot(dfFern$tauBMR, dfFern$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Fern") 
curve(modelFern_BMR2,0, 13, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestFern_RMR = lm(r~log(tauRMR), data = dfFern)
summary(betaestFern_RMR)
betaFern_RMR = summary(betaestFern_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelFern_RMR= function(tau){1-exp(-betaFern_RMR*tau)}
plot(dfFern$tauRMR, dfFern$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Fern") 
curve(modelFern_RMR, 0, 147, add = TRUE)



# Ginger!

#filter subset for just Ginger
dfGinger = LCATsub %>% filter(Name == "Ginger")

#birth weight 
birthweightGinger = min(dfGinger$Weight_g)

#average adult weight 
#get df when Ginger is an adult
adultdfGinger = dfGinger %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightGinger = mean(adultdfGinger$Weight_g)

#contruct tau and r columns 
dfGinger$r = (dfGinger$Weight_g/avgadweightGinger)^0.25
dfGinger$tauBMR = ((a_BMR*dfGinger$AgeAtWt_y)/(4*avgadweightGinger^0.25))-log(1-(birthweightGinger/avgadweightGinger)^0.25)
dfGinger$tauRMR=  ((a_RMR*dfGinger$AgeAtWt_y)/(4*avgadweightGinger^0.25))-log(1-(birthweightGinger/avgadweightGinger)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestGinger_BMR = lm(r ~ log(tauBMR), data = dfGinger)
betaGinger_BMR = summary(betaestGinger_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelGinger_BMR = function(tau){1-exp(-betaGinger_BMR*tau)}

plot(dfGinger$tauBMR, dfGinger$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ginger") 
curve(modelGinger_BMR, 0, 15, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaGinger_BMR2 = 1.1
modelGinger_BMR2 = function(tau){1-exp(-betaGinger_BMR2*tau)}
plot(dfGinger$tauBMR, dfGinger$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ginger") 
curve(modelGinger_BMR2,0, 15, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestGinger_RMR = lm(r~log(tauRMR), data = dfGinger)
summary(betaestGinger_RMR)
betaGinger_RMR = summary(betaestGinger_RMR)$coefficients[2]
###the right way to plot it! using a_RMR
modelGinger_RMR= function(tau){1-exp(-betaGinger_RMR*tau)}
plot(dfGinger$tauRMR, dfGinger$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ginger") 
curve(modelGinger_RMR, 0, 177, add = TRUE)



# Gretl!

#filter subset for just Gretl
dfGretl = LCATsub %>% filter(Name == "Gretl")

#birth weight 
birthweightGretl = min(dfGretl$Weight_g)

#average adult weight 
#get df when Gretl is an adult
adultdfGretl = dfGretl %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightGretl = mean(adultdfGretl$Weight_g)

#contruct tau and r columns 
dfGretl$r = (dfGretl$Weight_g/avgadweightGretl)^0.25
dfGretl$tauBMR = ((a_BMR*dfGretl$AgeAtWt_y)/(4*avgadweightGretl^0.25))-log(1-(birthweightGretl/avgadweightGretl)^0.25)
dfGretl$tauRMR=  ((a_RMR*dfGretl$AgeAtWt_y)/(4*avgadweightGretl^0.25))-log(1-(birthweightGretl/avgadweightGretl)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestGretl_BMR = lm(r ~ log(tauBMR), data = dfGretl)
betaGretl_BMR = summary(betaestGretl_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelGretl_BMR = function(tau){1-exp(-betaGretl_BMR*tau)}

plot(dfGretl$tauBMR, dfGretl$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Gretl") 
curve(modelGretl_BMR, 0, 15, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaGretl_BMR2 = 1.2
modelGretl_BMR2 = function(tau){1-exp(-betaGretl_BMR2*tau)}
plot(dfGretl$tauBMR, dfGretl$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Gretl") 
curve(modelGretl_BMR2,0, 15, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestGretl_RMR = lm(r~log(tauRMR), data = dfGretl)
betaGretl_RMR = summary(betaestGretl_RMR)$coefficients[2]
betaGretl_RMR
###the right way to plot it! using a_RMR
modelGretl_RMR= function(tau){1-exp(-betaGretl_RMR*tau)}
plot(dfGretl$tauRMR, dfGretl$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Gretl") 
curve(modelGretl_RMR, 0, 160, add = TRUE)



# Herodotus!

#filter subset for just Herodotus
dfHero = LCATsub %>% filter(Name == "Herodotus")

#birth weight 
birthweightHero = min(dfHero$Weight_g)

#average adult weight 
#get df when Herodotus is an adult
adultdfHero = dfHero %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightHero = mean(adultdfHero$Weight_g)

#contruct tau and r columns 
dfHero$r = (dfHero$Weight_g/avgadweightHero)^0.25
dfHero$tauBMR = ((a_BMR*dfHero$AgeAtWt_y)/(4*avgadweightHero^0.25))-log(1-(birthweightHero/avgadweightHero)^0.25)
dfHero$tauRMR=  ((a_RMR*dfHero$AgeAtWt_y)/(4*avgadweightHero^0.25))-log(1-(birthweightHero/avgadweightHero)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestHero_BMR = lm(r ~ log(tauBMR), data = dfHero)
betaHero_BMR = summary(betaestHero_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelHero_BMR = function(tau){1-exp(-betaHero_BMR*tau)}

plot(dfHero$tauBMR, dfHero$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Herodotus") 
curve(modelHero_BMR, 0, 15, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaHero_BMR2 = 1.2
modelHero_BMR2 = function(tau){1-exp(-betaHero_BMR2*tau)}
plot(dfHero$tauBMR, dfHero$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Herodotus") 
curve(modelHero_BMR2,0, 15, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestHero_RMR = lm(r~log(tauRMR), data = dfHero)
betaHero_RMR = summary(betaestHero_RMR)$coefficients[2]
betaHero_RMR
###the right way to plot it! using a_RMR
modelHero_RMR= function(tau){1-exp(-betaHero_RMR*tau)}
plot(dfHero$tauRMR, dfHero$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Herodotus") 
curve(modelHero_RMR, 0, 140, add = TRUE)



# Hibernia!

#filter subset for just Hibernia
dfHib = LCATsub %>% filter(Name == "Hibernia")

#birth weight 
birthweightHib = min(dfHib$Weight_g)

#average adult weight 
#get df when Hibernia is an adult
adultdfHib = dfHib %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightHib = mean(adultdfHib$Weight_g)

#contruct tau and r columns 
dfHib$r = (dfHib$Weight_g/avgadweightHib)^0.25
dfHib$tauBMR = ((a_BMR*dfHib$AgeAtWt_y)/(4*avgadweightHib^0.25))-log(1-(birthweightHib/avgadweightHib)^0.25)
dfHib$tauRMR=  ((a_RMR*dfHib$AgeAtWt_y)/(4*avgadweightHib^0.25))-log(1-(birthweightHib/avgadweightHib)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestHib_BMR = lm(r ~ log(tauBMR), data = dfHib)
betaHib_BMR = summary(betaestHib_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelHib_BMR = function(tau){1-exp(-betaHib_BMR*tau)}

plot(dfHib$tauBMR, dfHib$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Hibernia") 
curve(modelHib_BMR, 0, 17, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaHib_BMR2 = 1.1
modelHib_BMR2 = function(tau){1-exp(-betaHib_BMR2*tau)}
plot(dfHib$tauBMR, dfHib$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Hibernia") 
curve(modelHib_BMR2,0, 17, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestHib_RMR = lm(r~log(tauRMR), data = dfHib)
betaHib_RMR = summary(betaestHib_RMR)$coefficients[2]
betaHib_RMR
###the right way to plot it! using a_RMR
modelHib_RMR= function(tau){1-exp(-betaHib_RMR*tau)}
plot(dfHib$tauRMR, dfHib$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Hibernia") 
curve(modelHib_RMR, 0, 229, add = TRUE)



# Ivy!

#filter subset for just Ivy
dfIvy = LCATsub %>% filter(Name == "Ivy")

#birth weight 
birthweightIvy = min(dfIvy$Weight_g)

#average adult weight 
#get df when Ivy is an adult
adultdfIvy = dfIvy %>% filter(Age_Category == "adult")
#average the adult weights
avgadweightIvy = mean(adultdfIvy$Weight_g)

#contruct tau and r columns 
dfIvy$r = (dfIvy$Weight_g/avgadweightIvy)^0.25
dfIvy$tauBMR = ((a_BMR*dfIvy$AgeAtWt_y)/(4*avgadweightIvy^0.25))-log(1-(birthweightIvy/avgadweightIvy)^0.25)
dfIvy$tauRMR=  ((a_RMR*dfIvy$AgeAtWt_y)/(4*avgadweightIvy^0.25))-log(1-(birthweightIvy/avgadweightIvy)^0.25)

#model = (r=1-e^-beta*tau)

####BMR case 
#estimate betas for BMR
betaestIvy_BMR = lm(r ~ log(tauBMR), data = dfIvy)
betaIvy_BMR = summary(betaestIvy_BMR)$coefficients[2]   

###the right way to plot it! using a_BMR
modelIvy_BMR = function(tau){1-exp(-betaIvy_BMR*tau)}

plot(dfIvy$tauBMR, dfIvy$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ivy") 
curve(modelIvy_BMR, 0, 15, add = TRUE) #didnt do well, so estimate by eye

#estimate betas by eye 
betaIvy_BMR2 = 1
modelIvy_BMR2 = function(tau){1-exp(-betaIvy_BMR2*tau)}
plot(dfIvy$tauBMR, dfIvy$r, xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ivy") 
curve(modelIvy_BMR2,0, 15, add = TRUE)


###RMR case 
#estimate betas for RMR
betaestIvy_RMR = lm(r~log(tauRMR), data = dfIvy)
betaIvy_RMR = summary(betaestIvy_RMR)$coefficients[2]
betaIvy_RMR
###the right way to plot it! using a_RMR
modelIvy_RMR= function(tau){1-exp(-betaIvy_RMR*tau)}
plot(dfIvy$tauRMR, dfIvy$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ivy") 
curve(modelIvy_RMR, 0, 172, add = TRUE)




#######now that we have universal modified growth curve for all of the 26 lemurs
######we need to create the averaged representative growth curve 

#first do RMR case 
#need to average all of the RMR beta coefficients 

#make array of all the Beta essimates using RMR variable in a equation
dfRMRbetas= c(betaAlena_RMR,betaAlex_RMR,
              betaArt_RMR,betaBeri_RMR,betaBrig_RMR,betaCebes_RMR,betaDor_RMR,
              betaFern_RMR,betaGinger_RMR,betaGretl_RMR,betaHero_RMR,betaHib_RMR,
              betaIvy_RMR, betajones_RMR,betajustine_RMR,betalilah_RMR,betanicaea_RMR,
              betanikos_RMR,betaonyx_RMR,betapersephone_RMR,betaPHOTIUS_RMR,betasapp_RMR,betaSierraMist_RMR,
              betaSophia_RMR,betaStewart_RMR,betaTellus_RMR)


#averaging beta coefficents 
avgRMR= mean(dfRMRbetas)

#creating cumulative/represenative modified universal model and plotting it
repmodel_RMR= function(tau){1-exp(-avgRMR*tau)}
###need df of the constructed tau_RMR and r values including all 26 lemurs for plotting purposes
vectauRMR= c(dfAlena$tauRMR,dfAlex$tauRMR, dfArt$tauRMR,dfBeri$tauRMR,
             dfBrig$tauRMR, dfCebes$tauRMR, dfDor$tauRMR, dfFern$tauRMR, dfGinger$tauRMR,
             dfGretl$tauRMR, dfHero$tauRMR, dfHib$tauRMR, dfIvy$tauRMR, 
             dfjones$tauRMR_jones, dfjustine$tauRMR_justine, dflilah$tauRMR_lilah,
             dfnicaea$tauRMR_nicaea,dfnikos$tauRMR_nikos,dfonyx$tauRMR_onyx,
             dfpersephone$tauRMR_persephone,dfPHOTIUS$tauRMR_PHOTIUS, 
             dfSAPPHO$tauRMR_sapp, dfSierraMist$tauRMR_SierraMist, 
             dfSophia$tauRMR_Sophia,dfStewart$tauRMR_Stewart, dfTellus$tauRMR_Tellus)


vecR=c(dfAlena$r,dfAlex$r, dfArt$r,dfBeri$r,
       dfBrig$r, dfCebes$r, dfDor$r, dfFern$r, dfGinger$r,
       dfGretl$r, dfHero$r, dfHib$r, dfIvy$r, dfjones$r_jones, dfjustine$r_justine, dflilah$r_lilah,
       dfnicaea$r_nicaea,dfnikos$r_nikos,dfonyx$r_onyx,dfpersephone$r_persephone,
       dfPHOTIUS$r_PHOTIUS, dfSAPPHO$r_sapp, dfSierraMist$r_SierraMist, dfSophia$r_Sophia,dfStewart$r_Stewart, 
       dfTellus$r_Tellus)


#make vector a df to plot 
dftauRMR = as.data.frame(vectauRMR)
dfR=as.data.frame(vecR)

plot(dftauRMR$vectauRMR, dfR$vecR, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Representative Modified Universal Growth Curve for Lemur Subset") 
curve(repmodel_RMR, 0, 437, add = TRUE,col='red')


#add columns 
LCATsub$tauRMRrep= dftauRMR$vectauRMR
LCATsub$Rrep= dfR$vecR

#multiple ways to plot 
LCATsub %>%
  ggplot(aes(x=tauRMRrep, y=Rrep, colour=Name)) +
  geom_line() + #geom_line connects the data points for each lemur together 
  geom_point() + stat_function(fun=repmodel_RMR) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using RMR") + ylab("Dimensionless Weight")


LCATsub %>%
  ggplot(aes(x=tauRMRrep, y=Rrep, colour=Name))  +
  geom_point() + stat_function(fun=repmodel_RMR, aes(colour = "The Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using RMR") + ylab("Dimensionless Weight")

LCATsub %>%
  ggplot(aes(x=tauRMRrep, y=Rrep))  +
  geom_point() + stat_function(fun=repmodel_RMR, aes(colour = "Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using RMR") + ylab("Dimensionless Weight") + labs(color='Legend') 


#now do BMR case using properly estimated betas 
#need to average all of the BMR beta coefficients 

#make array of all the good Beta estimates using RMR variable in a parameter equation
dfBMRbetas= c(betaAlena_BMR2,betaAlex_BMR2,
              betaArt_BMR2,betaBeri_BMR2,betaBrig_BMR2,betaCebes_BMR2,betaDor_BMR2,
              betaFern_BMR2,betaGinger_BMR2,betaGretl_BMR2,betaHero_BMR2,betaHib_BMR2,
              betaIvy_BMR2, betajones_BMR2,betajustine_BMR2,betalilah_BMR2,betanicaea_BMR2,
              betanikos_BMR2,betaonyx_BMR2,betapersephone_BMR2,betaPHOTIUS_BMR2,betasapp_BMR2, betaSierraMist_BMR2,
              betaSophia_BMR2,betaStewart_BMR2,betaTellus_BMR2)


#averaging beta coefficents 
avgBMR= mean(dfBMRbetas)

#creating cumulative/represenative modified universal model and plotting it
#curve uses the good betas 
repmodel_BMR= function(tau){1-exp(-avgBMR*tau)}
###need df of the constructed tau_BMR and r values including all 26 lemurs for plotting purposes
vectauBMR= c(dfAlena$tauBMR,dfAlex$tauBMR, dfArt$tauBMR,dfBeri$tauBMR,
             dfBrig$tauBMR, dfCebes$tauBMR, dfDor$tauBMR, dfFern$tauBMR, dfGinger$tauBMR,
             dfGretl$tauBMR, dfHero$tauBMR, dfHib$tauBMR, dfIvy$tauBMR, dfjones$tauBMR_jones, dfjustine$tauBMR_justine, dflilah$tauBMR_lilah,
             dfnicaea$tauBMR_nicaea,dfnikos$tauBMR_nikos,dfonyx$tauBMR_onyx,dfpersephone$tauBMR_persephone,
             dfPHOTIUS$tauBMR_PHOTIUS, dfSAPPHO$tauBMR_sapp,dfSierraMist$tauBMR_SierraMist, dfSophia$tauBMR_Sophia,dfStewart$tauBMR_Stewart, 
             dfTellus$tauBMR_Tellus)

###note that the way i stacked all of the values are in alphabeltical order to easily combine dfs for plotting 

#make vector a df to plot 
dftauBMR = as.data.frame(vectauBMR)

plot(dftauBMR$vectauBMR, dfR$vecR, 
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Representative Modified Universal Growth Curve for Lemur Subset") 
curve(repmodel_BMR, 0, 31, add = TRUE,col='red')



#add columns to LCATsub dataframe for plotting purposes
LCATsub$tauBMRrep= dftauBMR$vectauBMR

LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep, colour=Name)) +
  geom_line() + #geom_line connects the data points for each lemur together 
  geom_point() + stat_function(fun=repmodel_BMR) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight")

#multiple ways to plot 
LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep, colour=Name))  +
  geom_point() + stat_function(fun=repmodel_BMR, aes(colour = "The Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight")

LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep))  +
  geom_point() + stat_function(fun=repmodel_BMR, aes(colour = "Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve for Lemur Subset") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight") + labs(color='Legend') 


#now do BMR case using the poorly fitting estimated betas 
#need to average all of the BMR beta coefficients 

#make array of all the Beta essimates using RMR variable in a equation
dfBMRbadbetas= c(betaAlena_BMR,betaAlex_BMR,
              betaArt_BMR,betaBeri_BMR,betaBrig_BMR,betaCebes_BMR,betaDor_BMR,
              betaFern_BMR,betaGinger_BMR,betaGretl_BMR,betaHero_BMR,betaHib_BMR,
              betaIvy_BMR, betajones_BMR,betajustine_BMR,betalilah_BMR,betanicaea_BMR,
              betanikos_BMR,betaonyx_BMR,betapersephone_BMR,betaPHOTIUS_BMR,
              betasapp_BMR_sappho, betaSierraMist_BMR,
              betaSophia_BMR,betaStewart_BMR,betaTellus_BMR)


#averaging beta coefficents 
avgBMRbad= mean(dfBMRbadbetas)

#creating cumulative/represenative modified universal model and plotting it
repmodel_BMRbad= function(tau){1-exp(-avgBMRbad*tau)}


###note that the way i stacked all of the values are in alphabeltical order to easily combine dfs for plotting 


plot(dftauBMR$vectauBMR, dfR$vecR, 
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Representative Modified Universal Growth Curve for Lemur Subset") 
curve(repmodel_BMRbad, 0, 31, add = TRUE,col='red')



#add columns 
LCATsub$tauBMRrep= dftauBMR$vectauBMR

LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep, colour=Name)) +
  geom_line() + #geom_line connects the data points for each lemur together 
  geom_point() + stat_function(fun=repmodel_BMRbad) + ggtitle("Modified Universal Growth Curve for Lemur Subset Using Poorly Fitting Betas") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight")

#multiple ways to plot 
LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep, colour=Name))  +
  geom_point() + stat_function(fun=repmodel_BMRbad, aes(colour = "The Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve For Lemur Subset Using Poorly Fitting Betas") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight")

LCATsub %>%
  ggplot(aes(x=tauBMRrep, y=Rrep))  +
  geom_point() + stat_function(fun=repmodel_BMRbad, aes(colour = "Modified Universal Growth Curve"), show.legend = TRUE) + ggtitle("Modified Universal Growth Curve for Lemur Subset Using Poorly Fitting Betas") +
  xlab("Dimensionless Time Using BMR") + ylab("Dimensionless Weight") + labs(color='Legend') 






#plotting all RMR individuals for reference 
#we see that each lemur has a unique growth curve
par(mfrow=c(4,2)) 
plot(dfSAPPHO$tauRMR_sapp, dfSAPPHO$r_sapp,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
SAPPHO")
curve(modelsapp_RMR,0, 178, add=TRUE)
plot(dfjones$tauRMR_jones, dfjones$r_jones,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Jones")
curve(modeljones_RMR,0, 175, add=TRUE)

plot(dfjustine$tauRMR_justine, dfjustine$r_justine,xlab="Dimensionless
Time using RMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For Justine")
curve(modeljustine_RMR,0, 340, add=TRUE)

plot(dflilah$tauRMR_lilah, dflilah$r_lilah,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Lilah")
curve(modellilah_RMR,0, 333, add=TRUE)

plot(dfnicaea$tauRMR_nicaea, dfnicaea$r_nicaea,xlab="Dimensionless Time
using RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Nicaea")
curve(modelnicaea_RMR,0, 185, add=TRUE)

plot(dfnikos$tauRMR_nikos, dfnikos$r_nikos,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Nikos")
curve(modelnikos_RMR,0, 140, add=TRUE)

plot(dfonyx$tauRMR_onyx, dfonyx$r_onyx,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Onyx")
curve(modelonyx_RMR,0, 130, add=TRUE)

plot(dfpersephone$tauRMR_persephone, dfpersephone$r_persephone,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Persephone")
curve(modelpersephone_RMR,0, 325, add=TRUE)


plot(dfPHOTIUS$tauRMR_PHOTIUS, dfPHOTIUS$r_PHOTIUS,xlab="Dimensionless
Time using RMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For PHOTIUS")
curve(modelPHOTIUS_RMR,0, 170, add=TRUE)


plot(dfSierraMist$tauRMR_SierraMist, dfSierraMist$r_SierraMist,xlab="Dimensionless Time using
RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Sierra Mist")
curve(modelSierraMist_RMR,0, 255, add=TRUE)

plot(dfSophia$tauRMR_Sophia, dfSophia$r_Sophia,xlab="Dimensionless Time
using RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Sophia")
curve(modelSophia_RMR,0, 345, add=TRUE)

plot(dfStewart$tauRMR_Stewart, dfStewart$r_Stewart,xlab="Dimensionless
Time using RMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For Stewart")
curve(modelStewart_RMR,0, 345, add=TRUE)

par(mfrow=c(6,3)) 
plot(dfTellus$tauRMR_Tellus, dfTellus$r_Tellus,xlab="Dimensionless Time
using RMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Tellus")
curve(modelTellus_RMR,0, 345, add=TRUE)

plot(dfAlena$tauRMR, dfAlena$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alena")
curve(modelAlena_RMR, 0, 327, add=TRUE)

plot(dfAlex$tauRMR, dfAlex$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alexander")
curve(modelAlex_RMR,0, 200, add=TRUE)

plot(dfArt$tauRMR, dfArt$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Artemesia")
curve(modelArt_RMR, 0, 130, add=TRUE)


plot(dfBeri$tauRMR, dfBeri$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Berisades")
curve(modelBeri_RMR, 0, 350, add=TRUE)

plot(dfBrig$tauRMR, dfBrig$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Brigitta")
curve(modelBrig_RMR, 0, 185, add=TRUE)


plot(dfCebes$tauRMR, dfCebes$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Cebes")
curve(modelCebes_RMR, 0, 140, add = TRUE)


plot(dfDor$tauRMR, dfDor$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Dorieus")
curve(modelDor_RMR, 0, 440, add = TRUE)

plot(dfFern$tauRMR, dfFern$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Fern")
curve(modelFern_RMR, 0, 147, add = TRUE)


plot(dfGinger$tauRMR, dfGinger$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ginger")
curve(modelGinger_RMR, 0, 177, add = TRUE)

plot(dfGretl$tauRMR, dfGretl$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Gretl")
curve(modelGretl_RMR, 0, 160, add = TRUE)


plot(dfHero$tauRMR, dfHero$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Herodotus") 
curve(modelHero_RMR, 0, 140, add = TRUE)

plot(dfHib$tauRMR, dfHib$r,
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Hibernia")
curve(modelHib_RMR, 0, 229, add = TRUE)

plot(dfIvy$tauRMR, dfIvy$r, 
     xlab = "Dimensionless Time using RMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ivy") 
curve(modelIvy_RMR, 0, 172, add = TRUE)





#plotting all BMR indivudals for reference 
par(mfrow=c(6,2)) 
plot(dfSAPPHO$tauBMR_sapp, dfSAPPHO$r_sapp,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
SAPPHO")
curve(modelsapp_BMR2,0, 14, add=TRUE)
plot(dfjones$tauBMR_jones, dfjones$r_jones,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Jones")
curve(modeljones_BMR2,0, 15, add=TRUE)

plot(dfjustine$tauBMR_justine, dfjustine$r_justine,xlab="Dimensionless
Time using BMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For Justine")
curve(modeljustine_BMR2,0, 25, add=TRUE)

plot(dflilah$tauBMR_lilah, dflilah$r_lilah,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Lilah")
curve(modellilah_BMR2,0, 25, add=TRUE)

plot(dfnicaea$tauBMR_nicaea, dfnicaea$r_nicaea,xlab="Dimensionless Time
using BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Nicaea")
curve(modelnicaea_BMR2,0, 15, add=TRUE)

plot(dfnikos$tauBMR_nikos, dfnikos$r_nikos,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Nikos")
curve(modelnikos_BMR2,0, 11, add=TRUE)

plot(dfonyx$tauBMR_onyx, dfonyx$r_onyx,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Onyx")
curve(modelonyx_BMR2,0, 11, add=TRUE)

plot(dfpersephone$tauBMR_persephone, dfpersephone$r_persephone,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Persephone")
curve(modelpersephone_BMR2,0, 25, add=TRUE)


plot(dfPHOTIUS$tauBMR_PHOTIUS, dfPHOTIUS$r_PHOTIUS,xlab="Dimensionless
Time using BMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For PHOTIUS")
curve(modelPHOTIUS_BMR2,0, 13, add=TRUE)


plot(dfSierraMist$tauBMR_SierraMist, dfSierraMist$r_SierraMist,xlab="Dimensionless Time using
BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Sierra Mist")
curve(modelSierraMist_BMR2,0, 20, add=TRUE)

plot(dfSophia$tauBMR_Sophia, dfSophia$r_Sophia,xlab="Dimensionless Time
using BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Sophia")
curve(modelSophia_BMR2,0, 25, add=TRUE)

plot(dfStewart$tauBMR_Stewart, dfStewart$r_Stewart,xlab="Dimensionless
Time using BMR",ylab="Dimensionless Weight",main="Modified Universal
Curve For Stewart")
curve(modelStewart_BMR2,0, 14, add=TRUE)

par(mfrow=c(6,3)) 
plot(dfTellus$tauBMR_Tellus, dfTellus$r_Tellus,xlab="Dimensionless Time
using BMR",ylab="Dimensionless Weight",main="Modified Universal Curve For
Tellus")
curve(modelTellus_BMR2,0, 25, add=TRUE)

plot(dfAlena$tauBMR, dfAlena$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alena")
curve(modelAlena_BMR2, 0, 25, add=TRUE)

plot(dfAlex$tauBMR, dfAlex$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Alexander")
curve(modelAlex_BMR2,0, 16, add=TRUE)

plot(dfArt$tauBMR, dfArt$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Artemesia")
curve(modelArt_BMR2, 0, 12, add=TRUE)


plot(dfBeri$tauBMR, dfBeri$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Berisades")
curve(modelBeri_BMR2, 0, 27, add=TRUE)

plot(dfBrig$tauBMR, dfBrig$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Brigitta")
curve(modelBrig_BMR2, 0, 15, add=TRUE)


plot(dfCebes$tauBMR, dfCebes$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Cebes")
curve(modelCebes_BMR2, 0, 12, add = TRUE)


plot(dfDor$tauBMR, dfDor$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Dorieus")
curve(modelDor_BMR2, 0, 33, add = TRUE)

plot(dfFern$tauBMR, dfFern$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Fern")
curve(modelFern_BMR2, 0, 13, add = TRUE)


plot(dfGinger$tauBMR, dfGinger$r,
     xlab = "Dimensionless Time using 2",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ginger")
curve(modelGinger_BMR2, 0, 15, add = TRUE)

plot(dfGretl$tauBMR, dfGretl$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Gretl")
curve(modelGretl_BMR2, 0, 15, add = TRUE)


plot(dfHero$tauBMR, dfHero$r, 
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Herodotus") 
curve(modelHero_BMR2, 0, 15, add = TRUE)

plot(dfHib$tauBMR, dfHib$r,
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Hibernia")
curve(modelHib_BMR2, 0, 17, add = TRUE)

plot(dfIvy$tauBMR, dfIvy$r, 
     xlab = "Dimensionless Time using BMR",
     ylab = "Dimensionless Weight",
     main = "Modified Universal Curve For Ivy") 
curve(modelIvy_BMR2, 0, 15, add = TRUE)




#LOESS fits !

#####plot whole subset 

#fitst do cross validation to find best estimation for alpha
#write function to do cross validation - found online! (documentation: https://rpubs.com/mengxu/loess_cv)
loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.25, 1, by = 0.05), folds = 5){
  # Do model selection using mean absolute error, which is more robust than squared error.
  mean.abs.error <- numeric(length(span.vals))
  
  # Quantify error for each span, using CV
  loess.model <- function(x, y, span){
    loess(y ~ x, span = span, control=loess.control(surface="direct"))
  }
  
  loess.predict <- function(fit, newdata) {
    predict(fit, newdata = newdata)
  }
  
  span.index <- 0
  for (each.span in span.vals) {
    span.index <- span.index + 1
    y.hat.cv <- crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
    non.empty.indices <- !is.na(y.hat.cv)
    mean.abs.error[span.index] <- mean(abs(y[non.empty.indices] - y.hat.cv[non.empty.indices]))
  }
  
  # find the span which minimizes error
  best.span <- span.vals[which.min(mean.abs.error)]
  
  # fit and return the best model
  best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="direct"))
  return(best.model)
}
model.1 <- loess_wrapper_extrapolate(LCATsub$AgeAtWt_y, LCATsub$Weight_g)
summary(model.1) #find the span value caluclated by cross validation

#now plot LOESS fit of whole subset with best estimation of alpha
LOESSfitsub <- LCATsub %>%
  do(LOESSfitsubb = ggplot() +
       geom_point(data = ., aes(x = AgeAtWt_y, y = Weight_g, alpha =  0.35,colour=Name)) +
       geom_smooth (data = ., aes(x = AgeAtWt_y, y = Weight_g), method = loess, se = TRUE)+
       ggtitle("LOESS Fit for Lemur Subset" )+
       xlab("Age in Years")+
       ylab("Weight in Grams"))
#plot
LOESSfitsub$LOESSfitsubb


##LOESS fit subset dimensionsless info RMR
##fitst do cross validation to find best estimation for alpha
model.2 <- loess_wrapper_extrapolate(LCATsub$tauRMRrep, LCATsub$Rrep)
summary(model.2) #find the span value caluclated by cross validation

#now plot lemur subset with RMR dimensionless data 
LOESSfitsubdim<- LCATsub %>%
  do(LOESSfitdim = ggplot() +
       geom_point(data = ., aes(x = tauRMRrep, y = Rrep, alpha = 0.25, colour=Name)) +
       geom_smooth (data = ., aes(x = tauRMRrep, y = Rrep), method = loess, se = TRUE)+
       ggtitle("LOESS Fit for Lemur Subset using Dimensionless Data" )+
       xlab("Dimensionsless Age with RMR")+
       ylab("Dimensionless Weight")) 


#plot
LOESSfitsubdim$LOESSfitdim



##LOESS fit subset dimensionsless info BMR
model.3 <- loess_wrapper_extrapolate(LCATsub$tauBMRrep, LCATsub$Rrep)
summary(model.3) #find the span value caluclated by cross validation
LOESSfitsubdim2<- LCATsub %>%
  do(LOESSfitdim2 = ggplot() +
       geom_point(data = ., aes(x = tauBMRrep, y = Rrep, alpha = 0.25, colour=Name)) +
       geom_smooth (data = ., aes(x = tauBMRrep, y = Rrep), method = loess, se = TRUE)+
       ggtitle("LOESS Fit for Lemur Subset using Dimensionless Data" )+
       xlab("Dimensionsless Age with BMR")+
       ylab("Dimensionless Weight"))
#plot
LOESSfitsubdim2$LOESSfitdim2



# calculating r-squared for the models!
##r^2 determines how well the model fits the data! Did our models do a good job?
# https://www.statology.org/sst-ssr-sse-in-r/ 

# weight vs age for lemur subset LOESS
#find sse: Sum of Squares Error (SSE) 
# sum of squared differences between predicted data points (ŷi) and observed data points (yi)
sse1 <- sum((fitted(model.1) - LCATsub$Weight_g)^2)
sse1
#find ssr: Sum of Squares Regression (SSR) 
# sum of squared differences between predicted data points (ŷi) and the mean of the response variable(y)
ssr1 <- sum((fitted(model.1) - mean(LCATsub$Weight_g))^2)
ssr1
#find sst: Sum of Squares Total (SST) 
# sum of squared differences between individual data points (yi) and the mean of the response variable (y)
# same as ssr + sse
sst1 <- ssr1 + sse1
sst1
# r-squared
ssr1/sst1 # 0.8386289 - decent

#messing around seeing how well the fit was by hand for that first and last value
#this is just a sanity check on our end 
min(LCATsub$Weight_g)
min(fitted(model.1))
max(LCATsub$Weight_g)
max(fitted(model.1))

# RMR dimensionless data LOESS
#model2 <- loess(Rrep ~ tauRMRrep, data = LCATsub)
#find sse
sse2 <- sum((fitted(model.2) - LCATsub$Rrep)^2)
sse2
#find ssr
ssr2 <- sum((fitted(model.2) - mean(LCATsub$Rrep))^2)
ssr2
#find sst
sst2 <- ssr2 + sse2
sst2
# r-squared
ssr2/sst2 # 0.9777611 - good!

# BMR dimensionless data LOESS
#model3 <- loess(Rrep ~ tauBMRrep, data = LCATsub)
#find sse
sse3 <- sum((fitted(model.3) - LCATsub$Rrep)^2)
sse3
#find ssr
ssr3 <- sum((fitted(model.3) - mean(LCATsub$Rrep))^2)
ssr3
#find sst
sst3 <- ssr3 + sse3
sst3
# r-squared
ssr3/sst3 # 0.9809554 - good! 



# BMR dimensionless using modified growth curve (betas esimated by eye)
#getting fitted values from function
fittedvalsBMRrepmodel = repmodel_BMR(tau=LCATsub$Rrep)
sse4 <- sum((fittedvalsBMRrepmodel - LCATsub$Rrep)^2)
sse4
#find ssr
ssr4 <- sum(((fittedvalsBMRrepmodel) - mean(LCATsub$Rrep))^2)
ssr4
#find sst
sst4 <- ssr4 + sse4
sst4
# r-squared
ssr4/sst4 # 0.4977451 - looks like it smooths out data too much! or too much variation amoung indivudal lemurs (i.e. each lemur is growing differently), LOESS provides better fit.
#can interpret this as lemur having random growthspurts. 




# BMR dimensionless using modified growth curve using poorly fitting betas 
#getting fitted values from function
fittedvalsBMRbadrepmodel = repmodel_BMRbad(tau=LCATsub$Rrep)
#print(fittedvalsBMRbadrepmodel)
sse5 <- sum((fittedvalsBMRbadrepmodel - LCATsub$Rrep)^2)
sse5
#find ssr
ssr5 <- sum(((fittedvalsBMRbadrepmodel) - mean(LCATsub$Rrep))^2)
ssr5
#find sst
sst5 <- ssr5 + sse5
sst5
# r-squared
ssr5/sst5 #0.4943752 - looks like it smooths out data too much! or too much variation amoung indivudal lemurs (i.e. each lemur is growing differently), LOESS provides better fit.

#messing around
min(LCATsub$Rrep)
min(fittedvalsBMRbadrepmodel)
min(LCATsub$tauBMRrep)


# RMR dimensionless using modified growth curve
#getting fitted values from function
fittedvalsRMRrepmodel = repmodel_RMR(tau=LCATsub$Rrep)
#print(fittedvalsBMRbadrepmodel)
sse6 <- sum((fittedvalsRMRrepmodel - LCATsub$Rrep)^2)
sse6
#find ssr
ssr6 <- sum(((fittedvalsRMRrepmodel) - mean(LCATsub$Rrep))^2)
ssr6
#find sst
sst6 <- ssr6 + sse6
sst6
# r-squared
ssr6/sst6 # 0.4943514 - looks like it smooths out data too much! or too much variation amoung indivudal lemurs (i.e. each lemur is growing differently), LOESS provides better fit.


