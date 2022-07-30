##Script for analysis of density dependent viral infection
##of Ranchman's tiger moth
##By Adam Pepi and Vincent Pan
rm(list=ls())
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
setwd("~/Documents/Research/Woolly Bear Virus Project/25 sites paper/diseasedduv")

## Data processing
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
vsitesm2<-vsitesm1[!is.na(vsitesm1$infection_n),]
cats<-read.csv('ads_post_mort_cleaned.csv')
str(catsm)
catsm<-cats
catsm1<-left_join(catsm,sitesdata,by="site")
str(catsm1)

ads_post_mort_cleaned<-catsm1
ads_post_mort_no_parasite_cleaned<-catsm1 %>% filter(emerged_parasite=="0"|is.na(emerged_parasite))
str(ads_post_mort_no_parasite_cleaned)
ads_post_mort_no_parasite_cleaned$infection_severity.ord
catseverity<-ads_post_mort_no_parasite_cleaned[!is.na(ads_post_mort_no_parasite_cleaned$infection_severity),]
catseverity$infection_severity.ord<-factor(catseverity$infection_severity,levels=c("uninfected","l","m","h","sh"),ordered = T)

vsitesm2<-vsitesm1[!is.na(vsitesm1$infection_n),]
vsitesm2$num_nat_loose
vsitesm2$this_year.y
vsitesm2$site
str(vsitesm2)


cor(log(catseverity$last_year.y),catseverity$insolation360to180)
cor(log(catseverity$last_year.y),catseverity$edd360to180.scaledlog)
cor(log(vsitesm2$last_year.y),vsitesm2$insolation360to180)
cor(log(vsitesm2$last_year.y),vsitesm2$edd360to180.scaledlog)
cor(log(catseverity$this_year.y),catseverity$edd360to180.scaledlog)
cor(log(vsitesm2$this_year.y),vsitesm2$insolation360to180)
cor(log(catseverity$this_year.y),log(catseverity$last_year.y))
cor(log(vsitesm2$this_year.y),log(vsitesm2$last_year.y))


##Models of infection
pr.inf.m1<-glmmTMB(infection~log(this_year.y)+log(last_year.y)+(1|site),family=binomial(),data=ads_post_mort_no_parasite_cleaned);summary(pr.inf.m1)
check_collinearity(pr.inf.m1)

pr.inf.m2<-glmmTMB(infection~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),family=binomial(),data=ads_post_mort_no_parasite_cleaned);summary(pr.inf.m2)
check_collinearity(pr.inf.m2)

##influence plot
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
pr.inf.m2_influence <- influence_mixed(pr.inf.m2, groups="site")
car::infIndexPlot(pr.inf.m2_influence)


#### Infection severity

clmm1<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+(1|site),data=catseverity)
summary(clmm1)

check_collinearity(clmm1)

clmm1.1<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+(1|site),data=catseverity%>% filter(infection=='1'))
summary(clmm1.1)
check_collinearity(clmm1.1)


clmm2<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),data=catseverity,control = clmm.control( method= 'ucminf',gradTol = 1e-2, maxIter = 1000, maxLineIter = 1000), nAGQ = 5)
summary(clmm2)

check_collinearity(clmm2)

clmm2.1<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),data=catseverity%>% filter(infection=='1'),control = clmm.control( method= 'ucminf',gradTol = 1e-2, maxIter = 1000, maxLineIter = 1000), nAGQ = 5)
summary(clmm2.1)


#Probability of survival to end of life


survival.ed.m1<-glmmTMB(cbind(num_nat_loose,total_reared-num_nat_loose)~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+emerged_parasite_m,family =betabinomial,data=vsitesm2 );summary(survival.ed.m1)

source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
survival.ed.m1_influence <- influence_mixed(survival.ed.m1, groups="site")
car::infIndexPlot(survival.ed.m1_influence)

check_collinearity(survival.ed.m1)


vsitesm1$num_nat_loose_n<-vsitesm1$num_nat_loose/vsitesm1$total_reared

#Infection plots

predframeinfection <- with(ads_post_mort_no_parasite_cleaned,expand.grid(site=levels(site),last_year.y=seq(min(last_year.y),max(last_year.y), .01),this_year.y=mean(this_year.y),insolation360to180=mean(insolation360to180)))

mydf <- ggpredict(pr.inf.m2, terms = c("last_year.y[exp]"))
str(mydf)

mydf$last_year.y<-mydf$x
mydf$infection<-mydf$predicted
mydfl<-mydf
mydfl$infection<-mydf$conf.low
mydfh<-mydf
mydfh$infection<-mydfh$conf.high


infectionplot<-ggplot(ads_post_mort_no_parasite_cleaned,aes(x=last_year.y,y=infection))+scale_x_log10()+
  geom_count()+geom_line(data=mydf, size=1.5)+geom_line(data=mydfl, lty=2)+geom_line(data=mydfh, lty=2) + xlab("Density last year") + ylab("Proportion infected")+ theme_classic()+ggtitle('a')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))
  
infectionplot

mydf2 <- ggpredict(pr.inf.m2, terms = c("edd360to180.scaledlog"))
str(mydf2)

mydf2$edd360to180.scaledlog<-mydf2$x
mydf2$infection<-mydf2$predicted
mydf2l<-mydf2
mydf2l$infection<-mydf2l$conf.low
str(mydf2l)
mydf2h<-mydf2
mydf2h$infection<-mydf2h$conf.high

ggplot(mydf2,aes(x=x,y=predicted))+geom_line()

infectionplot2<-ggplot(ads_post_mort_no_parasite_cleaned,aes(x=edd360to180.scaledlog,y=infection))+
  geom_count()+geom_line(data=mydf2, size=1.5)+ geom_line(data=mydf2l, lty=2)+geom_line(data=mydf2h, lty=2) +xlab("Log Erythemal Daily Dose (scaled)") + ylab("Proportion infected")+ theme_classic()+ggtitle('b')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

infectionplot2

grid.arrange(infectionplot, infectionplot2, ncol=2)


##Infection severity plots

mydf3 <- ggpredict(clmm2, terms = c("last_year.y"))
str(mydf3)

mydf3$last_year.y<-mydf3$x
mydf3$infection_severity.ord<-mydf3$response.level

mydf3$severity<-mydf3$predicted
mydf3l<-mydf3
mydf3l$severity<-mydf3l$conf.low
mydf3h<-mydf3
mydf3h$severity<-mydf3h$conf.high

ads_post_mort_no_parasite_cleaned$severity<-0

severityplot<-ggplot(mydf3,aes(x=last_year.y,y=severity,colour=infection_severity.ord))+scale_x_log10()+geom_line(data=mydf3, size=1.5)+geom_line(data=mydf3l, lty=2)+geom_line(data=mydf3h, lty=2) + xlab("Density last year") + ylab("Proportion")+ theme_classic()+ scale_colour_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ ggtitle('a')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

severityplot

severityplot1.1<-ggplot(mydf3,aes(x=last_year.y,y=severity,colour=infection_severity.ord))+scale_x_log10()+  geom_area(aes(fill = infection_severity.ord)) + ylab("Proportion")+ xlab("Density last year")+theme_classic()+ scale_colour_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ scale_fill_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ ggtitle('b')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

severityplot1.1

mydf5 <- ggpredict(clmm2, terms = c("edd360to180.scaledlog"))
str(mydf5)

mydf5$edd360to180.scaledlog<-mydf5$x
mydf5$infection_severity.ord<-mydf5$response.level
mydf5$severity<-mydf5$predicted
mydf5l<-mydf5
mydf5l$severity<-mydf5l$conf.low
mydf5h<-mydf5
mydf5h$severity<-mydf5h$conf.high

ads_post_mort_no_parasite_cleaned$severity<-0

severityplot3<-ggplot(mydf5,aes(x=edd360to180.scaledlog,y=severity,colour=infection_severity.ord))+geom_line(data=mydf5, size=1.5)+geom_line(data=mydf5l, lty=2)+geom_line(data=mydf5h, lty=2) + xlab("Log Erythemal Daily Dose (scaled)") + ylab("Proportion")+ theme_classic()+ scale_colour_viridis(name = "Infection severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ ggtitle('c')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

severityplot3

severityplot3.1<-ggplot(mydf5,aes(x=edd360to180.scaledlog,y=severity,colour=infection_severity.ord,fill=infection_severity.ord))+geom_area(aes(fill = infection_severity.ord)) + xlab("Log Erythemal Daily Dose (scaled)") + ylab("Proportion")+ theme_classic()+ scale_colour_viridis(name = "Infection severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ scale_fill_viridis(name = "Infection severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+ ggtitle('d')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))
severityplot3.1


grid.arrange(severityplot,severityplot1.1,severityplot3,severityplot3.1)

##survival to end of life



mydf6 <- ggpredict(survival.ed.m1, terms = c("last_year.y"))
str(mydf6)

mydf6$last_year.y<-mydf6$x
mydf6$surv_ed<-mydf6$predicted
mydf6l<-mydf6
mydf6l$surv_ed<-mydf6$conf.low
mydf6h<-mydf6
mydf6h$surv_ed<-mydf6h$conf.high


str(vsitesm1)

vsitesm1$surv_ed<-vsitesm1$num_nat_loose/vsitesm1$total_reared

survplot<-ggplot(vsitesm1,aes(x=last_year.y,y=surv_ed))+scale_x_log10()+
  geom_point(aes(size=total_reared))+labs(size='N')+geom_line(data=mydf6, size=1.5)+geom_line(data=mydf6l, lty=2)+geom_line(data=mydf6h, lty=2) + xlab("Density last year") + ylab("Survival")+ theme_classic()+ggtitle('a')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

survplot

mydf7 <- ggpredict(survival.ed.m1, terms = c("edd360to180.scaledlog"))
str(mydf7)
mydf7$edd360to180.scaledlog<-mydf7$x
mydf7$surv_ed<-mydf7$predicted
mydf7l<-mydf7
mydf7l$surv_ed<-mydf7l$conf.low
str(mydf7l)
mydf7h<-mydf7
mydf7h$surv_ed<-mydf7h$conf.high


survplot2<-ggplot(vsitesm1,aes(x=edd360to180.scaledlog,y=surv_ed))+
geom_point(aes(size=total_reared))+labs(size='N')+geom_line(data=mydf7, size=1.5)+ geom_line(data=mydf7l, lty=2)+geom_line(data=mydf7h, lty=2) +xlab("Log Erythemal Daily Dose (scaled)") + ylab("Survival")+ theme_classic()+ggtitle('b')+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

survplot2


grid.arrange(survplot,survplot2,ncol=2)


##Plot data

toPlot<-catseverity%>%filter(!is.na(infection_severity)) %>% 
  group_by(infection_severity, site)%>%
  summarise(n =n() ,last_year.y=mean(last_year.y),edd=mean(edd360to180.scaledlog))%>%
  group_by(site)%>%
  mutate(prop = n/sum(n))


toPlot%>%
group_by(last_year.y)%>%
summarise(prop=sum(prop))

toPlot%>%
  group_by(site)%>%
  summarise(d=mean(last_year.y))

toPlot$infection_severity<-factor(toPlot$infection_severity,ordered=T,levels=c("uninfected","l","m","h","sh"))

str(toPlot)
toPlot<-rbind(toPlot,toPlot[1,],toPlot[3,],toPlot[3,],toPlot[12,])


toPlot[62,]
#toPlot[62,1]<-'uninfected'
#toPlot[62,3]<-0
#toPlot[62,6]<-0
toPlot[63,]
toPlot[63,1]<-'uninfected'
toPlot[63,3]<-0
toPlot[63,6]<-0

toPlot[64,]
toPlot[64,1]<-'l'
toPlot[64,3]<-0
toPlot[64,6]<-0

toPlot[65,]
toPlot[65,1]<-'uninfected'
toPlot[65,3]<-0
toPlot[65,6]<-0
toPlot<-toPlot[-12,]


##Figure 4

pp1<-toPlot %>% ggplot()+geom_area(aes(y=prop,x=last_year.y,fill = infection_severity),position="stack")+scale_x_log10()+theme_classic()+geom_vline(xintercept=c(4.11, 0.436, 0.188, 3.78, 0.0298, 0.0364, 0.427, 0.164, 0.150, 0.0743, 0.0755, 0.299, 0.0243, 0.0300, 0.0791), linetype="dotted")+ xlab("Site by density last year") + ylab("Proportion")+ theme_classic()+ scale_fill_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))



pp2<-ggplot(ads_post_mort_no_parasite_cleaned,aes(x=last_year.y,y=infection, group=site))+scale_x_log10()+ stat_summary() + xlab("Density last year") + ylab("Proportion")+ theme_classic()+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))

plot_grid(pp2,pp1,labels=c('a','b'), rel_widths=c(1,1.5))


#Figure S2
toPlot2<-catseverity%>%filter(!is.na(infection_severity)) %>% 
  group_by(infection_severity, stage)%>%
  summarise(n =n() )%>%
  group_by(stage)%>%
  mutate(prop = n/sum(n))


toPlot2
toPlot2$infection_severity<-factor(toPlot2$infection_severity,ordered=T,levels=c("uninfected","l","m","h","sh"))
ggplot(toPlot2,aes(x=stage,y=prop,fill=infection_severity))+geom_bar(stat="identity",position='stack')+ xlab("Stage") + ylab("Proportion")+ theme_classic()+ scale_fill_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)


##Figure S3


ppp1<-toPlot %>% ggplot(aes(y=prop,x=reorder(site,last_year.y),fill=infection_severity))+geom_bar(stat="identity",position='stack')+ xlab("Site by density last year") + ylab("Proportion")+ theme_classic()+ scale_fill_viridis(name = "Infection Severity", labels = c("Uninfected", "Low", "Medium","High","Very High"),discrete=T)+guides(n='none')

ppp1

tp3<-ads_post_mort_no_parasite_cleaned[ads_post_mort_no_parasite_cleaned$site!="SHP",]
tp3<-tp3[tp3$site!='RKP',]
ppp2<-ggplot(tp3,aes(x=reorder(site,last_year.y),y=infection ))+ stat_summary() + xlab("Site by density last year") + ylab("Proportion infected")+ theme_classic()
ppp2

plot_grid(ppp2,ppp1,labels=c('a','b'), rel_widths=c(1,1.3))


## Sensitivity analyses

ads_post_mort_no_parasite_cleaned_s<-ads_post_mort_no_parasite_cleaned[ads_post_mort_no_parasite_cleaned$site!='BMR',]

ads_post_mort_no_parasite_cleaned_s<-ads_post_mort_no_parasite_cleaned_s[ads_post_mort_no_parasite_cleaned_s$site!='FVW',]

##Models of infection
pr.inf.m1s<-glmmTMB(infection~log(this_year.y)+log(last_year.y)+(1|site),family=binomial(),data=ads_post_mort_no_parasite_cleaned_s);summary(pr.inf.m1s)

pr.inf.m2s<-glmmTMB(infection~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),family=binomial(),data=ads_post_mort_no_parasite_cleaned_s);summary(pr.inf.m2s)

##Infection severity

catseverity_s<-catseverity[catseverity$site!='BMR',]
catseverity_s<-catseverity_s[catseverity_s$site!='FVW',]

clmm1s<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+(1|site),data=catseverity_s)
summary(clmm1s)

clmm1.1s<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+(1|site),data=catseverity_s%>% filter(infection=='1'))
summary(clmm1.1s)

clmm2s<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),data=catseverity_s,control = clmm.control( method= 'ucminf',gradTol = 1e-2, maxIter = 1000, maxLineIter = 1000), nAGQ = 5)
summary(clmm2s)

clmm2.1s<-clmm(infection_severity.ord~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+(1|site),data=catseverity_s %>% filter(infection=='1'),control = clmm.control( method= 'ucminf',gradTol = 1e-2, maxIter = 1000, maxLineIter = 1000), nAGQ = 5)
summary(clmm2.1s)


#Probability of survival to end of life

vsitesm1s<-vsitesm1[vsitesm1$site!='BMR',]
vsitesm1s<-vsitesm1s[vsitesm1s$site!='FVW',]


survival.ed.m1s<-glmmTMB(cbind(num_nat_loose,total_reared-num_nat_loose)~log(this_year.y)+log(last_year.y)+edd360to180.scaledlog+emerged_parasite_m,family =betabinomial,data=vsitesm1s ,control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")));summary(survival.ed.m1s)

#Revision of Figure 2

pp2<-ggplot(ads_post_mort_no_parasite_cleaned,aes(x=last_year.y,y=infection, group=site))+scale_x_log10()+ stat_summary() + xlab("Density last year") + ylab("Proportion infected")+ theme_classic()+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))+geom_line(data=mydf,aes(x=last_year.y,y=infection,group=group), size=1.5)+geom_line(data=mydfl,aes(x=last_year.y,y=infection,group=group),  lty=2)+geom_line(data=mydfh,aes(x=last_year.y,y=infection,group=group),  lty=2)+ylim(0,1)

pp2


pp3<-ggplot(ads_post_mort_no_parasite_cleaned,aes(x=edd360to180.scaledlog,y=infection, group=site))+scale_x_log10()+ stat_summary() + xlab("Density last year") + ylab("Proportion infected")+ theme_classic()+theme(text=element_text(size=20), plot.title=element_text(hjust=0, size=20))+geom_line(data=mydf2,aes(x=edd360to180.scaledlog,y=infection,group=group), size=1.5)+geom_line(data=mydf2l,aes(x=edd360to180.scaledlog,y=infection,group=group),  lty=2)+geom_line(data=mydf2h,aes(x=edd360to180.scaledlog,y=infection,group=group),  lty=2)+ylim(0,1)+xlab("Log Erythemal Daily Dose (scaled)")

pp3

plot_grid(pp2,pp3,labels=c('a','b'), rel_widths=c(1,1))

