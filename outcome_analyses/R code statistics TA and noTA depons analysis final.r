###################################################
###### bring together time area closure (TA) and no TA outcomes and analysis with a TA covariate
###################################################
library(glmmTMB)
library(ggplot2)

load("~/results/mortality_all.Rdata")

mortality.all.df<-subset(mortality.df,condition=="symmetric" |condition=="nochange")
mortality.all.df$TA<-"noTA"


load("~/resultsTA/mortality_all.Rdata")
mortality.df$TA<-"TA"

mortality.all.df<-rbind(mortality.all.df,mortality.df)

save(mortality.all.df,file="~/resultsTA/mortality_all_wTAcov.Rdata")


mortality.all.df$fyear<-factor(mortality.all.df$year)
mortality.all.df$fprevalence<-factor(mortality.all.df$prevalence)
 
pop.ana0<-glmmTMB(count~condition+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana1<-glmmTMB(count~fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana2<-glmmTMB(count~condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana3<-glmmTMB(count~condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(count~TA+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana5<-glmmTMB(count~TA+condition+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana6<-glmmTMB(count~TA+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana7<-glmmTMB(count~TA+condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana8<-glmmTMB(count~TA+condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(count~TA*condition+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana10<-glmmTMB(count~TA*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(count~TA*condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df)
pop.ana12<-glmmTMB(count~TA*condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(count~TA*fprevalence+condition+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.all.df,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))

AIC(pop.ana13)

AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)

          # df      AIC
# pop.ana0   5 596054.3
# pop.ana1  14 571377.2
# pop.ana2  15 571231.5
# pop.ana3  25 571157.7
# pop.ana4   5 465224.6
# pop.ana5   6 465078.9
# pop.ana6  15 440401.8
# pop.ana7  16 440256.1
# pop.ana8  26 440182.3
# pop.ana9   7 465051.5
# pop.ana10 25 422201.7
# pop.ana11 17 440228.7
# pop.ana12 47 421900.4
# pop.ana13 26 422055.9


library(ggeffects)
count.pred<-ggpredict(pop.ana12,terms=c("fprevalence","condition","TA"))

count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers and area closure"
count.pred$facet<-factor(count.pred$facet)

tiff(file="~/figures for ms/Figure xx abundance.tiff",width=7, height=5,units="in",res=300)
plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="Abundance", x="simulated pinger prevalence",color="condition mediation",title="")+
theme(legend.position="bottom")+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
dev.off()

###############################################################################
#### away from initial conditions, last 20 years


population.year.20<-mortality.all.df[mortality.all.df$year>19,]
population.year.20$fyear<-factor(population.year.20$year-19)
population.year.20$year<-population.year.20$year-19


pop.ana0<-glmmTMB(count~condition+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
pop.ana1<-glmmTMB(count~fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana2<-glmmTMB(count~condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana3<-glmmTMB(count~condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(count~TA+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
pop.ana5<-glmmTMB(count~TA+condition+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana6<-glmmTMB(count~TA+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana7<-glmmTMB(count~TA+condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana8<-glmmTMB(count~TA+condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
pop.ana9<-glmmTMB(count~TA*condition+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana10<-glmmTMB(count~TA*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(count~TA*condition+fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana12<-glmmTMB(count~TA*condition*fprevalence+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(count~TA*fprevalence+condition+(1|run)+ar1(fyear+0|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)

         # df      AIC
# pop.ana0   5       NA
# pop.ana1  14 309059.5
# pop.ana2  15 308944.3
# pop.ana3  25 308843.1
# pop.ana4   5       NA
# pop.ana5   6 234610.4
# pop.ana6  15 217350.5
# pop.ana7  16 217235.3
# pop.ana8  26       NA
# pop.ana9   7 234569.0
# pop.ana10 25 203827.3
# pop.ana11 17 217215.3
# pop.ana12 47 203549.4
# pop.ana13 26 203712.1


library(ggeffects)
count.pred<-ggpredict(pop.ana12,terms=c("fprevalence","condition","TA"))

count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers and area closure"
count.pred$facet<-factor(count.pred$facet)


tiff(file="~/figures for ms/Figure xx abundance 20yrs.tiff",width=7, height=5,units="in",res=300)
plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="Abundance", x="simulated pinger prevalence",color="condition mediation",title="")+
theme(legend.position="bottom")+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
dev.off()

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
### COD
load("~/results/cod_all.Rdata")


cod.df<-subset(cod.try,condition=="symmetric" |condition=="nochange")
cod.df$TA<-"noTA"


load("~/resultsTA/cod_all.Rdata")
cod.try$TA<-"TA"

cod.df<-rbind(cod.df,cod.try)

save(cod.df,file="~/resultsTA/cod_all_wTAcov.Rdata")


cod.df$fyear<-factor(cod.df$year)
cod.df$fprevalence<-factor(cod.df$prevalence)
 
cod.df$pbycatch<-cbind(cod.df$bycatch,(cod.df$n-cod.df$bycatch))
cod.df$count<-mortality.all.df$count
cod.df$rbycatch<-cbind(cod.df$bycatch,(cod.df$count-cod.df$bycatch))



population.year.20<-cod.df[cod.df$year>19,]
population.year.20$fyear<-factor(population.year.20$year-19)
population.year.20$year<-population.year.20$year-19


pop.ana0<-glmmTMB(rbycatch~condition+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana1<-glmmTMB(rbycatch~fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana2<-glmmTMB(rbycatch~condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana3<-glmmTMB(rbycatch~condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(rbycatch~TA+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana5<-glmmTMB(rbycatch~TA+condition+(1|run), family=binomial, data=population.year.20)
pop.ana6<-glmmTMB(rbycatch~TA+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana7<-glmmTMB(rbycatch~TA+condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana8<-glmmTMB(rbycatch~TA+condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(rbycatch~TA*condition+(1|run), family=binomial, data=population.year.20)
pop.ana10<-glmmTMB(rbycatch~TA*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(rbycatch~TA*condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana12<-glmmTMB(rbycatch~TA*condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(rbycatch~TA*fprevalence+condition+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)

          # df      AIC
# pop.ana0   3 85691.77
# pop.ana1  12 71882.78
# pop.ana2  13 71880.04
# pop.ana3  23 71889.84
# pop.ana4   3 82595.16
# pop.ana5   4 82594.23
# pop.ana6  13 69288.84
# pop.ana7  14 69286.22
# pop.ana8  24 69296.13
# pop.ana9   5 82595.88
# pop.ana10 23 69212.35
# pop.ana11 15 69287.39
# pop.ana12 45 69235.35
# pop.ana13 24 69209.79



library(ggeffects)
count.pred<-ggpredict(pop.ana13,type="fixed",terms=c("fprevalence","condition","TA"))
#plot(count.pred)



count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers and area closure"
count.pred$facet<-factor(count.pred$facet)


tiff(file="~/figures for ms/Figure xx bycatch rate 20yrs.tiff",width=7, height=5,units="in",res=300)
plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="bycatch rate", x="simulated pinger prevalence",color="condition mediation",title="")+
theme(legend.position="bottom")+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
dev.off()

  
##################################################################################################################################################


pop.ana0<-glmmTMB(pbycatch~condition+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana1<-glmmTMB(pbycatch~fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana2<-glmmTMB(pbycatch~condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana3<-glmmTMB(pbycatch~condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(pbycatch~TA+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana5<-glmmTMB(pbycatch~TA+condition+(1|run), family=binomial, data=population.year.20)
pop.ana6<-glmmTMB(pbycatch~TA+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana7<-glmmTMB(pbycatch~TA+condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana8<-glmmTMB(pbycatch~TA+condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(pbycatch~TA*condition+(1|run), family=binomial, data=population.year.20)
pop.ana10<-glmmTMB(pbycatch~TA*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(pbycatch~TA*condition+fprevalence+(1|run), family=binomial, data=population.year.20)
pop.ana12<-glmmTMB(pbycatch~TA*condition*fprevalence+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(pbycatch~TA*fprevalence+condition+(1|run), family=binomial, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)

          # df      AIC
# pop.ana0   3 84921.00
# pop.ana1  12 70591.09
# pop.ana2  13 70587.69
# pop.ana3  23 70597.83
# pop.ana4   3 81864.41
# pop.ana5   4 81862.96
# pop.ana6  13 68005.04
# pop.ana7  14 68001.60
# pop.ana8  24 68012.08
# pop.ana9   5 81864.36
# pop.ana10 23 67943.57
# pop.ana11 15 68002.41
# pop.ana12 45 67965.44
# pop.ana13 24 67940.22


library(ggeffects)
count.pred<-ggpredict(pop.ana13,type="fixed",terms=c("fprevalence","condition","TA"))
#plot(count.pred)



count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers and area closure"
count.pred$facet<-factor(count.pred$facet)


tiff(file="~/figures for ms/Figure xx proportion of bycatch deaths 20yrs.tiff",width=7, height=5,units="in",res=300)
plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="proportion of deaths caused by bycatch", x="simulated pinger prevalence",color="condition mediation",title="")+
theme(legend.position="bottom")+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
dev.off()

  
##################################################################################################################################################




population.year.20<-mortality.all.df[mortality.all.df$year>19,]
population.year.20$fyear<-factor(population.year.20$year-19)
population.year.20$year<-population.year.20$year-19


pop.ana0<-glmmTMB(deaths~condition+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana1<-glmmTMB(deaths~fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana2<-glmmTMB(deaths~condition+fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana3<-glmmTMB(deaths~condition*fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(deaths~TA+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana5<-glmmTMB(deaths~TA+condition+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana6<-glmmTMB(deaths~TA+fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana7<-glmmTMB(deaths~TA+condition+fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana8<-glmmTMB(deaths~TA+condition*fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(deaths~TA*condition+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana10<-glmmTMB(deaths~TA*fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(deaths~TA*condition+fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20)
pop.ana12<-glmmTMB(deaths~TA*condition*fprevalence+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(deaths~TA*fprevalence+condition+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))

pop.ana4<-glmmTMB(deaths~TA+offset(log(count))+(1|run), family=poisson, data=population.year.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13w<-glmmTMB(deaths~TA*fprevalence+condition+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)

AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
 
### retained below
# #without ar1
         # df      AIC
# pop.ana0   3 162746.2
# pop.ana1  12 162755.9
# pop.ana2  13 162757.8
# pop.ana3  23 162769.9
# pop.ana4   3 162693.9
# pop.ana5   4 162695.9
# pop.ana6  13 162704.9
# pop.ana7  14 162706.9
# pop.ana8  24 162719.0
# pop.ana9   5 162697.6
# pop.ana10 23 162717.8
# pop.ana11 15 162708.6
# pop.ana12 45 162747.4
# pop.ana13 24 162719.7


count.pred<-ggpredict(pop.ana4,type="fixed",terms=c("TA"))



count.pred$x<-as.character(count.pred$x)
count.pred$x[count.pred$x=="noTA"]<-"pingers only"
count.pred$x[count.pred$x=="TA"]<-"pingers and area closure"
count.pred$x<-factor(count.pred$x)

count.pred20<-count.pred
attr(count.pred20, "x.axis.labels")<-c("pingers only","pingers & area closure")

tiff(file="~/figures for ms/Figure xx death numbers 20yrs.tiff",width=7, height=5,units="in",res=300)
count.20.plot<-plot(count.pred20,dot.size=1.2,line.size=.5)+
				labs(y="number of deaths", x="scenario",title="",caption="B. Last 20 years")+
				ylim(22.7,26)+
				theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

###########################################################################################



mortality.df0<-mortality.all.df[mortality.all.df$year>0,]



pop.ana0<-glmmTMB(deaths~condition+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana1<-glmmTMB(deaths~fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana2<-glmmTMB(deaths~condition+fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana3<-glmmTMB(deaths~condition*fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(deaths~TA+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana5<-glmmTMB(deaths~TA+condition+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana6<-glmmTMB(deaths~TA+fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana7<-glmmTMB(deaths~TA+condition+fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana8<-glmmTMB(deaths~TA+condition*fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(deaths~TA*condition+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana10<-glmmTMB(deaths~TA*fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(deaths~TA*condition+fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)
pop.ana12<-glmmTMB(deaths~TA*condition*fprevalence+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(deaths~TA*fprevalence+condition+offset(log(count))+(1|run)+ar1(fyear+0|run), family=poisson, data=mortality.df0)


AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
 
          # df      AIC
# pop.ana0   5 323007.5
# pop.ana1  14 322871.9
# pop.ana2  15 322871.4
# pop.ana3  25 322886.0
# pop.ana4   5 321760.5
# pop.ana5   6 321759.9
# pop.ana6  15 321582.2
# pop.ana7  16 321581.5
# pop.ana8  26 321595.9
# pop.ana9   7 321760.7
# pop.ana10 25 321473.3
# pop.ana11 17 321582.3
# pop.ana12 47 321503.5
# pop.ana13 26 321472.6

count.pred<-ggpredict(pop.ana13,type="fixed",terms=c("fprevalence","condition","TA"))
plot(count.pred)



count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers & area closure"
count.pred$facet<-factor(count.pred$facet)

count.plot<-plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="number of deaths", x="simulated pinger prevalence",color="",title="",caption="A. whole time series")+
theme(legend.position="bottom")+
ylim(22.5,25.5)+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
 
 library(grid)
 library(gridExtra)

grid.arrange(count.plot,count.20.plot,ncol=2,widths=c(6,3))

tiff(file="~/figures for ms/Figure xx number of deaths.tiff",width=7, height=5,units="in",res=300)
grid.arrange(count.plot,count.20.plot,ncol=2,widths=c(6,3))
dev.off()

#######################################################################################################################################
######################################################################################################################################
#######################################################################################################################################
######################################################################################################################################
#######################################################################################################################################
######################################################################################################################################
#######################################################################################################################################
######################################################################################################################################
####### lifetime reproductive success




library(stringr)
library(ggplot2)

###let's bring it all together
folder.top<-list.dirs("~/results",full.names=T,recursive=F)

condition<-str_split(folder.top,"/",simplify=TRUE)[,7]
prevalence<-c(0,10,100,20,30,40,50,60,70,80,90) ### careful of the text based order

data.types<-c("Statistics","Lifecycle","Mortality","Reproduction")

folder.low<-list.dirs(folder.top[1],recursive=F,full.names=F)

init.files<-list.files(paste(folder.top[1],folder.low[1],sep="/"),full.names=T)

rs<-read.csv(init.files[grep(data.types[4],init.files)],sep=";",header=T)
rs$prevalence<-prevalence[1]
rs$condition<-condition[1]

rs<-rs[-c(1:dim(rs)[1]),]

for (i in 1:length(folder.top)) { 
	for (j in 1:length(folder.low)) {
	init.files<-list.files(paste(folder.top[i],folder.low[j],sep="/"),full.names=T)
	temp<-read.csv(init.files[grep(data.types[4],init.files)],sep=";",header=T)
	if (dim(temp)[2]==1) {
	temp<-read.csv(init.files[grep(data.types[4],init.files)],sep=",",header=T)
	}
	temp$prevalence<-prevalence[j]
	temp$condition<-condition[i]

	rs<-rbind(rs,temp)
	print(c(i,j))
	flush.console()
	}
}

save(rs,file="~/results/reproduction_all.Rdata")


###let's bring it all together
folder.top<-list.dirs("~/resultsTA",full.names=T,recursive=F)

condition<-str_split(folder.top,"/",simplify=TRUE)[,7]
prevalence<-c(0,10,100,20,30,40,50,60,70,80,90) ### careful of the text based order

data.types<-c("Statistics","Lifecycle","Mortality","Reproduction")

folder.low<-list.dirs(folder.top[1],recursive=F,full.names=F)

init.files<-list.files(paste(folder.top[1],folder.low[1],sep="/"),full.names=T)

rsTA<-read.csv(init.files[grep(data.types[4],init.files)],sep=";",header=T)
rsTA$prevalence<-prevalence[1]
rsTA$condition<-condition[1]

rsTA<-rsTA[-c(1:dim(rsTA)[1]),]

for (i in 1:length(folder.top)) { 
	for (j in 1:length(folder.low)) {
	init.files<-list.files(paste(folder.top[i],folder.low[j],sep="/"),full.names=T)
	temp<-read.csv(init.files[grep(data.types[4],init.files)],sep=";",header=T)
	if (dim(temp)[2]==1) {
	temp<-read.csv(init.files[grep(data.types[4],init.files)],sep=",",header=T)
	}
	temp$prevalence<-prevalence[j]
	temp$condition<-condition[i]

	rsTA<-rbind(rsTA,temp)
	print(c(i,j))
	flush.console()
	}
}

save(rsTA,file="~/resultsTA/reproduction_all.Rdata")



rs<-rs[rs$condition=="nochange" |rs$condition=="symmetric",]

rs$TA<-"noTA"
rsTA$TA<-"TA"
 
rs<-rbind(rs,rsTA)

save(rs,file="~/reproduction_all_noTATA.Rdata")


rs$age<-floor(rs$AgeAtDeath) # number of years available to reproduce (repro attempt max)

rs.0<-subset(rs,age>0)

#and then now we restrict to the last 20 years
rs.0$year<-floor(rs.0$tick/(360*48))

rs.0$fyear<-factor(rs.0$year)
rs.0$fprevalence<-factor(rs.0$prevalence)

rs.0.20<-subset(rs.0,year>19)
rs.0.20$fyear<-factor(rs.0.20$year-19)

#a few options for random effect
# first gut feeling is
#(1|run) + (1|year)
# but past 20 years, the (1|year) could be not really there so just (1|run)
#also actually we should perhaps consider it as (1|run/year)
#model selection supports hypothesis 2

pop.ana00<-glmmTMB(CalvesBorn~1+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana0<-glmmTMB(CalvesBorn~condition+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana1<-glmmTMB(CalvesBorn~fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana2<-glmmTMB(CalvesBorn~condition+fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana3<-glmmTMB(CalvesBorn~condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(CalvesBorn~TA+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana5<-glmmTMB(CalvesBorn~TA+condition+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana6<-glmmTMB(CalvesBorn~TA+fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana7<-glmmTMB(CalvesBorn~TA+condition+fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana8<-glmmTMB(CalvesBorn~TA+condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(CalvesBorn~TA*condition+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana10<-glmmTMB(CalvesBorn~TA*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(CalvesBorn~TA*condition+fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana12<-glmmTMB(CalvesBorn~TA*condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(CalvesBorn~TA*fprevalence+condition+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


pop.ana3<-glmmTMB(CalvesBorn~condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana8<-glmmTMB(CalvesBorn~TA+condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana10<-glmmTMB(CalvesBorn~TA*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana12<-glmmTMB(CalvesBorn~TA*condition*fprevalence+offset(log(age))+(1|run), family=poisson, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana00<-glmmTMB(CalvesBorn~1+offset(log(age))+(1|run), family=poisson, data=rs.0.20)

AIC(pop.ana00,pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
 
 
           # df     AIC
# pop.ana00  2 1870448
# pop.ana0   3 1870450
# pop.ana1  12 1870462
# pop.ana2  13 1870464
# pop.ana3  23 1870477
# pop.ana4   3 1870450
# pop.ana5   4 1870452
# pop.ana6  13 1870464
# pop.ana7  14 1870466
# pop.ana8  24 1870479
# pop.ana9   5 1870454
# pop.ana10 23 1870477
# pop.ana11 15 1870468
# pop.ana12 45 1870511
# pop.ana13 24 1870479

library(ggeffects)
count.pred<-ggpredict(pop.ana00,type="fixed",terms=c("age"))
print(count.pred,digits=4)

# age | Predicted |         95% CI
# --------------------------------
  # 0 |      0.00 | [ 0.00,  0.00]
  # 4 |      1.57 | [ 1.56,  1.57]
  # 8 |      3.13 | [ 3.13,  3.13]
 # 12 |      4.70 | [ 4.69,  4.70]
 # 14 |      5.48 | [ 5.47,  5.49]
 # 18 |      7.04 | [ 7.03,  7.05]
 # 22 |      8.61 | [ 8.60,  8.62]
 # 30 |     11.74 | [11.72, 11.76]


##############################################
####LRS

pop.ana00<-glmmTMB(CalvesWeaned~1+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana0<-glmmTMB(CalvesWeaned~condition+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana1<-glmmTMB(CalvesWeaned~fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana2<-glmmTMB(CalvesWeaned~condition+fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana3<-glmmTMB(CalvesWeaned~condition*fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(CalvesWeaned~TA+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana5<-glmmTMB(CalvesWeaned~TA+condition+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana6<-glmmTMB(CalvesWeaned~TA+fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana7<-glmmTMB(CalvesWeaned~TA+condition+fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana8<-glmmTMB(CalvesWeaned~TA+condition*fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(CalvesWeaned~TA*condition+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana10<-glmmTMB(CalvesWeaned~TA*fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(CalvesWeaned~TA*condition+fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)
pop.ana12<-glmmTMB(CalvesWeaned~TA*condition*fprevalence+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(CalvesWeaned~TA*fprevalence+condition+offset(log(age))+(1|run), family=nbinom1, data=rs.0.20)

AIC(pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
 
#poisson
          # df     AIC
# pop.ana0   3 1814958
# pop.ana1  12 1814968
# pop.ana2  13 1814970
# pop.ana3  23 1814981
# pop.ana4   3 1814938
# pop.ana5   4 1814939
# pop.ana6  13 1814950
# pop.ana7  14 1814952
# pop.ana8  24 1814963
# pop.ana9   5 1814941
# pop.ana10 23 1814957
# pop.ana11 15 1814954
# pop.ana12 45 1814985
# pop.ana13 24 1814958

pop.ana4b<-glmmTMB(CalvesWeaned~TA+offset(log(age))+(1|run), family=nbinom2, data=rs.0.20)
#lets make sure mr poisson is appropriate as we deal with quite a few zeros
AIC(pop.ana4b)
# 1810084

pop.ana4c<-glmmTMB(CalvesWeaned~TA+offset(log(age))+(1|run), family=nbinom2, zi=~fprevalence+offset(log(age)),data=rs.0.20)
AIC(pop.ana4c)

#1810106

#neg. binomial
          # df     AIC
# pop.ana0   4 1802545
# pop.ana1  13 1802557
# pop.ana2  14 1802559
# pop.ana3  24 1802572
# pop.ana4   4 1802532
# pop.ana5   5 1802534
# pop.ana6  14 1802547
# pop.ana7  15 1802549
# pop.ana8  25 1802562
# pop.ana9   6 1802536
# pop.ana10 24 1802556
# pop.ana11 16 1802551
# pop.ana12 46 1802589
# pop.ana13 25 1802558


pop.ana4d<-glmmTMB(CalvesWeaned~TA+(1|run), family=nbinom2, data=rs.0.20)
#LRS is the product of calf production and lifespan
AIC(pop.ana4d)

pop.ana00d<-glmmTMB(CalvesWeaned~1+(1|run), family=nbinom2, data=rs.0.20)
pop.ana0d<-glmmTMB(CalvesWeaned~condition+(1|run), family=nbinom2, data=rs.0.20)
pop.ana1d<-glmmTMB(CalvesWeaned~fprevalence+(1|run), family=nbinom2, data=rs.0.20)
pop.ana2d<-glmmTMB(CalvesWeaned~condition+fprevalence+(1|run), family=nbinom2, data=rs.0.20)
pop.ana3d<-glmmTMB(CalvesWeaned~condition*fprevalence+(1|run), family=nbinom2, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
#pop.ana4d<-glmmTMB(CalvesWeaned~TA+(1|run), family=nbinom2, data=rs.0.20)
pop.ana5d<-glmmTMB(CalvesWeaned~TA+condition+(1|run), family=nbinom2, data=rs.0.20)
pop.ana6d<-glmmTMB(CalvesWeaned~TA+fprevalence+(1|run), family=nbinom2, data=rs.0.20)
pop.ana7d<-glmmTMB(CalvesWeaned~TA+condition+fprevalence+(1|run), family=nbinom2, data=rs.0.20)
pop.ana8d<-glmmTMB(CalvesWeaned~TA+condition*fprevalence+(1|run), family=nbinom2, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9d<-glmmTMB(CalvesWeaned~TA*condition+(1|run), family=nbinom2, data=rs.0.20)
pop.ana10d<-glmmTMB(CalvesWeaned~TA*fprevalence+(1|run), family=nbinom2, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11d<-glmmTMB(CalvesWeaned~TA*condition+fprevalence+(1|run), family=nbinom2, data=rs.0.20)
pop.ana12d<-glmmTMB(CalvesWeaned~TA*condition*fprevalence+(1|run), family=nbinom2, data=rs.0.20,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13d<-glmmTMB(CalvesWeaned~TA*fprevalence+condition+(1|run), family=nbinom2, data=rs.0.20)

AIC(pop.ana00d,pop.ana0d,pop.ana1d,pop.ana2d,pop.ana3d,pop.ana4d,pop.ana5d,pop.ana6d,pop.ana7d,pop.ana8d,pop.ana9d,pop.ana10d,pop.ana11d,pop.ana12d,pop.ana13d)
 
           # df     AIC
# pop.ana00d  3 2319438
# pop.ana0d   4 2319440
# pop.ana1d  13 2319455
# pop.ana2d  14 2319457
# pop.ana3d  24 2319472
# pop.ana4d   4 2319440
# pop.ana5d   5 2319442
# pop.ana6d  14 2319457
# pop.ana7d  15 2319459
# pop.ana8d  25 2319473
# pop.ana9d   6 2319444
# pop.ana10d 24 2319470
# pop.ana11d 16 2319461
# pop.ana12d 46 2319504
# pop.ana13d 25 2319472

weaned.pred<-ggpredict(pop.ana4b,"fixed",term=c("TA","age [8]"))

lrs.pred<-ggpredict(pop.ana00d)

print(lrs.pred,digits=4)




count.pred$facet<-as.character(count.pred$facet)
count.pred$facet[count.pred$facet=="noTA"]<-"pingers only"
count.pred$facet[count.pred$facet=="TA"]<-"pingers and area closure"
count.pred$facet<-factor(count.pred$facet)

count.plot<-plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="lifetime calf production", x="simulated pinger prevalence",color="",title="",caption="A. whole time series")+
theme(legend.position="bottom")+
ylim(22.5,25.5)+
scale_color_discrete(labels=c("no condition effects","non-linear condition effects"))
 
 library(grid)
 library(gridExtra)

grid.arrange(count.plot,count.20.plot,ncol=2,widths=c(6,3))

tiff(file="~/figures for ms/Figure xx number of deaths.tiff",width=7, height=5,units="in",res=300)
grid.arrange(count.plot,count.20.plot,ncol=2,widths=c(6,3))
dev.off()



##############################################
####weaning success

rs.0.20.sub<-subset(rs.0.20,CalvesBorn>0)
rs.0.20.sub$pweaned<-cbind(rs.0.20.sub$CalvesWeaned,(rs.0.20.sub$CalvesBorn-rs.0.20.sub$CalvesWeaned))

pop.ana00<-glmmTMB(pweaned~1+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana0<-glmmTMB(pweaned~condition+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana1<-glmmTMB(pweaned~fprevalence+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana2<-glmmTMB(pweaned~condition+fprevalence+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana3<-glmmTMB(pweaned~condition*fprevalence+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(pweaned~TA+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana5<-glmmTMB(pweaned~TA+condition+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana6<-glmmTMB(pweaned~TA+fprevalence+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana7<-glmmTMB(pweaned~TA+condition+fprevalence+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana8<-glmmTMB(pweaned~TA+condition*fprevalence+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(pweaned~TA*condition+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana10<-glmmTMB(pweaned~TA*fprevalence+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(pweaned~TA*condition+fprevalence+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana12<-glmmTMB(pweaned~TA*condition*fprevalence+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(pweaned~TA*fprevalence+condition+(1|run), family=binomial, data=rs.0.20.sub)


AIC(pop.ana00,pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)

          # df     AIC
# pop.ana00  2 1214821
# pop.ana0   3 1214822
# pop.ana1  12 1214789
# pop.ana2  13 1214790
# pop.ana3  23 1214807
# pop.ana4   3 1214749
# pop.ana5   4 1214751
# pop.ana6  13 1214721
# pop.ana7  14 1214723
# pop.ana8  24 1214739
# pop.ana9   5 1214753
# pop.ana10 23 1214729
# pop.ana11 15 1214725
# pop.ana12 45 1214759
# pop.ana13 24 1214730

count.pred<-ggpredict(pop.ana6,type="fixed",terms=c("fprevalence","TA"))
plot(count.pred)



 count.pred$group<-as.character(count.pred$group)
 count.pred$group[count.pred$group=="noTA"]<-"pingers only"
 count.pred$group[count.pred$group=="TA"]<-"pingers and area closure"
 count.pred$group<-factor(count.pred$group)
attr(count.pred, "legend.labels")<-c("pingers only","pingers & area closure")


count.plot<-plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="proportion of calves weaned", x="simulated pinger prevalence",color="",title="")+
theme(legend.position="bottom")+ylim(.68,.7)

tiff(file="~/figures for ms/Figure xx pweaned 20 years.tiff",width=7, height=5,units="in",res=300)
count.plot
dev.off()

##############################################
####weaning success offset?

rs.0.20.sub<-subset(rs.0.20,CalvesBorn>0)
rs.0.20.sub$pweaned<-cbind(rs.0.20.sub$CalvesWeaned,(rs.0.20.sub$CalvesBorn-rs.0.20.sub$CalvesWeaned))

pop.ana00<-glmmTMB(pweaned~1+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana0<-glmmTMB(pweaned~condition+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana1<-glmmTMB(pweaned~fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana2<-glmmTMB(pweaned~condition+fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana3<-glmmTMB(pweaned~condition*fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(pweaned~TA+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana5<-glmmTMB(pweaned~TA+condition+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana6<-glmmTMB(pweaned~TA+fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana7<-glmmTMB(pweaned~TA+condition+fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana8<-glmmTMB(pweaned~TA+condition*fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(pweaned~TA*condition+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana10<-glmmTMB(pweaned~TA*fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(pweaned~TA*condition+fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)
pop.ana12<-glmmTMB(pweaned~TA*condition*fprevalence+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana13<-glmmTMB(pweaned~TA*fprevalence+condition+offset(log(age))+(1|run), family=binomial, data=rs.0.20.sub)


AIC(pop.ana00,pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
#################################################
          # df      AIC
# pop.ana00  2 975230.0
# pop.ana0   3 975231.8
# pop.ana1  12 975223.3
# pop.ana2  13 975225.1
# pop.ana3  23 975243.6
# pop.ana4   3 975177.1
# pop.ana5   4 975179.0
# pop.ana6  13 975172.1
# pop.ana7  14 975173.9
# pop.ana8  24 975192.4
# pop.ana9   5 975180.9
# pop.ana10 23 975183.1
# pop.ana11 15 975175.9
# pop.ana12 45 975222.4
# pop.ana13 24 975184.9


count.pred<-ggpredict(pop.ana6,type="fixed",terms=c("fprevalence","TA"))
plot(count.pred)



 count.pred$group<-as.character(count.pred$group)
 count.pred$group[count.pred$group=="noTA"]<-"pingers only"
 count.pred$group[count.pred$group=="TA"]<-"pingers and area closure"
 count.pred$group<-factor(count.pred$group)
attr(count.pred, "legend.labels")<-c("pingers only","pingers & area closure")


count.plot<-plot(count.pred,dot.size=1.2,line.size=.5)+labs(y="proportion of calves weaned", x="simulated pinger prevalence",color="",title="")+
theme(legend.position="bottom")+ylim(.68,.7)

#################################################
pop.ana00<-glmmTMB(CalvesWeaned~1+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana0<-glmmTMB(CalvesWeaned~condition+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana1<-glmmTMB(CalvesWeaned~fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana2<-glmmTMB(CalvesWeaned~condition+fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana3<-glmmTMB(CalvesWeaned~condition*fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana4<-glmmTMB(CalvesWeaned~TA+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana5<-glmmTMB(CalvesWeaned~TA+condition+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana6<-glmmTMB(CalvesWeaned~TA+fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana7<-glmmTMB(CalvesWeaned~TA+condition+fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana8<-glmmTMB(CalvesWeaned~TA+condition*fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana9<-glmmTMB(CalvesWeaned~TA*condition+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana10<-glmmTMB(CalvesWeaned~TA*fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
pop.ana11<-glmmTMB(CalvesWeaned~TA*condition+fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana12<-glmmTMB(CalvesWeaned~TA*condition*fprevalence+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)
pop.ana13<-glmmTMB(CalvesWeaned~TA*fprevalence+condition+offset(log(CalvesBorn))+(1|run), family=nbinom2, data=rs.0.20.sub)


AIC(pop.ana00,pop.ana0,pop.ana1,pop.ana2,pop.ana3,pop.ana4,pop.ana5,pop.ana6,pop.ana7,pop.ana8,pop.ana9,pop.ana10,pop.ana11,pop.ana12,pop.ana13)
 


