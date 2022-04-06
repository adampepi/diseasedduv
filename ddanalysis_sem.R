##Data analysis of models for paper

rm(list=ls())

#Need to have old version of piecewiseSEM to use GLMMTMB
#library(devtools)

install_version("piecewiseSEM", version = "1.2.1", repos = "http://cran.us.r-project.org")
library(tidyverse)
library(tidyverse)
library(glmmTMB)
library(performance)
library(ordinal)
library(ggeffects)
library(viridis)
library(car)
library(gridExtra)
library(performance)
library(cowplot)
library(piecewiseSEM)
setwd("~/Documents/Research/Woolly Bear Virus Project/25 sites paper")

sitesdata<-read.csv('sitesdata.csv')
str(sitesdata)

sitesdata$insolation360to180.scaled<-scale(sitesdata$insolation360to180,scale=T,center=T) 
sitesdata$insolation360to180log<-log(sitesdata$insolation360to180)
sitesdata$insolation360to180.scaledlog<-scale(sitesdata$insolation360to180log,scale=T,center=T) 
sitesdata$edd360to180.scaled<-scale(sitesdata$edd360to180,scale=T,center=T) 

sitesdata$edd360to180log<-log(sitesdata$edd360to180)
sitesdata$edd360to180.scaledlog<-scale(sitesdata$edd360to180log,scale=T,center=T) 

sitesdata$last_year_scale<-scale(sitesdata$insolation360to180,scale=T,center=T) 

vsites<-read.csv('ads_site_cleaned.csv')

str(vsites)

vsitesm<-vsites
vsitesm1<-left_join(vsitesm,sitesdata,by="site")
vsitesm1

cats<-read.csv('ads_post_mort_cleaned.csv')
catsm<-cats
catsm1<-left_join(catsm,sitesdata,by="site")
str(catsm1)

ads_post_mort_cleaned<-catsm1

ads_post_mort_no_parasite_cleaned<-catsm1 %>% filter(emerged_parasite=="0"|is.na(emerged_parasite))

str(ads_post_mort_no_parasite_cleaned)

ads_post_mort_no_parasite_cleaned$infection_severity.ord

catseverity<-ads_post_mort_no_parasite_cleaned[!is.na(ads_post_mort_no_parasite_cleaned$infection_severity),]

catseverity$infection_severity.ord<-factor(catseverity$infection_severity,levels=c("uninfected","l","m","h","sh"),ordered = T)

##Prepare data for SEM
semdata<-ads_post_mort_cleaned
semdata<-semdata[!is.na(semdata$death_cause),]
semdata$infection_severity.ord
semdata$natural_death<-0
semdata$natural_death[grepl("natural",semdata$death_cause)]<-1
semdata$infection_severity.ord<-factor(semdata$infection_severity.ord,levels=c("uninfected","l","m","h","sh"),ordered = T)
semdata$infection3<-ifelse(as.numeric(semdata$infection_severity.ord)>3,1,0)
semdata$emerged_parasite[is.na(semdata$emerged_parasite)]<-0
semdata<-semdata[semdata$emerged_parasite==0,]
semdata$severity_beta<-((as.numeric(semdata$infection_severity.ord))/5)-0.1
plot(density(semdata$severity_beta))

##SEM component models
denssem1<-glmmTMB(infection~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),family =betabinomial,data=semdata)
summary(denssem1)

denssem2<-glmmTMB(severity_beta~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+infection+(1|site),family =beta_family(link ='logit'),data=semdata)
summary(denssem2)

survsem1<-glmmTMB(natural_death~severity_beta+infection+(1|site),family =betabinomial,data=semdata)
summary(survsem1)

##Run SEM
vlist = list(
 denssem1,
 denssem2,
 survsem1
)

sem.fit(vlist, semdata)
sem.coefs(vlist,semdata)
rsquared(vlist, semdata)

range(semdata$edd360to180.scaledlog)
range(semdata$severity_beta)
range(semdata$infection)
range(log(semdata$last_year.y))

##needs to run with newer piecewise
coefs(denssem1, standardize.type = "latent.linear")


Beta.glm <- fixef(survsem1)$cond[2]

preds <- predict(survsem1, type = "link")

# Compute sd of error variance using theoretical variances
sd.y.LT <- sqrt(var(preds) + pi^2/3)

# Compute sd of x
sd.x <- sd(semdata$severity_beta)

Beta.glm*sd.x/sd.y.LT


#Density
0.645

#infection
0.34

#EDD
-0.089

##severity
-0.388

##Density->Infection->Severity->Survival
0.645*0.34*-0.388
##EDD->Severity->Survival
-0.089*-0.388

-0.0850884/0.034532

