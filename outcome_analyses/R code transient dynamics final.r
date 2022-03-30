####### R code transient dynamics of bycatch system
 plot3d(population.year[population.year$prevalence==0&population.year$condition=="nochange",]$FoodEnergyLevel,population.year[population.year$prevalence==0&population.year$condition=="nochange",]$PorpoiseEnergyLevel,population.year[population.year$prevalence==0&population.year$condition=="nochange",]$PorpoiseCount,xlab="Food",ylab="condition",zlab="abundance")
population.year$col<-as.numeric

 plot3d(population.year[population.year$prevalence==0|population.year$prevalence==100,]$FoodEnergyLevel,population.year[population.year$prevalence==0|population.year$prevalence==100,]$PorpoiseEnergyLevel,population.year[population.year$prevalence==0|population.year$prevalence==100,]$PorpoiseCount,col=as.numeric(factor(population.year[population.year$prevalence==0|population.year$prevalence==100,]$prevalence)),xlab="Food",ylab="condition",zlab="abundance")

## not much df around food availability once we account for porpoise abundance and energy

########### first determine the states
#let's start with TA runs

load("~/resultsTA/count_all_yearly.Rdata")

run=unique(population.year$run)
prevalence=unique(population.year$prevalence)
condition=unique(population.year$condition)
year=unique(population.year$year)

population.year$dcount<-NA
population.year$denergy<-NA
population.year$dfood<-NA


for (i in 1:length(condition)) {
	for (j in 1:length(prevalence)) {
		for (k in 1:length(run)) {
			for (l in 1:(length(year)-1)) {
				population.year$dcount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$denergy[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$dfood[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
			}
		}
	}
}


####################################################################
####################################################################

library(mgcv)

stable.state<-data.frame(condition=rep(c("nochange","symmetric"),each=11),prevalence=rep(seq(0,100,10),2),Estar=NA,Nstar=NA)

for (j in 1:dim(stable.state)[1]) {
sym.bi.gam<-gam(list(dcount~s(PorpoiseCount,PorpoiseEnergyLevel),denergy~s(PorpoiseCount,PorpoiseEnergyLevel)),family=mvn(d=2),data=population.year[population.year$condition==stable.state$condition[j] & population.year$prevalence==stable.state$prevalence[j],])


sym.newdf<-data.frame(dn=NA,de=NA,PorpoiseCount=rep(seq(100,230,1),each=length(seq(1700,4000,1))),PorpoiseEnergyLevel=rep(seq(1700,4000,1),length(seq(100,230,1))))
sym.newdf[,1:2]<-predict(sym.bi.gam,sym.newdf[,3:4])

isocline.n<-data.frame(abundance=100:230,energy=NA)
for (i in 1:dim(isocline.n)[1]) { 
isocline.n$energy[i]<-sym.newdf$PorpoiseEnergyLevel[sym.newdf$PorpoiseCount==isocline.n$abundance[i] & abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.n$abundance[i] ,]$de)==min(abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.n$abundance[i] ,]$de))]
}

isocline.e<-data.frame(abundance=100:230,energy=NA)
for (i in 1:dim(isocline.e)[1]) { 
isocline.e$energy[i]<-sym.newdf$PorpoiseEnergyLevel[sym.newdf$PorpoiseCount==isocline.e$abundance[i] & abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.e$abundance[i] ,]$dn)==min(abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.e$abundance[i] ,]$dn))]
}

lm.e<-(lm(abundance~energy,data=isocline.e))
lm.n<-(lm(abundance~energy,data=isocline.n))
cm<-rbind(coef(lm.e),coef(lm.n))

stable.state[j,3:4]<-c(-solve(cbind(cm[,2],-1))%*%cm[,1])
print(j)
flush.console()

}

stable.state$treatment<-"TA"
save(stable.state,file="~/resultsTA/stablestate.Rdata")

library(ggplot2)

ggplot(stable.state,aes(x=Nstar,y=Estar,colour=factor(prevalence),group=condition))+
geom_point(aes(shape=condition,size=(prevalence+1)/10))

#################################################################################################################################
#################################################################################################################################


load("~/results/count_all_yearly.Rdata")

run=unique(population.year$run)
prevalence=unique(population.year$prevalence)
condition=unique(population.year$condition)
year=unique(population.year$year)

population.year$dcount<-NA
population.year$denergy<-NA
population.year$dfood<-NA


for (i in 1:length(condition)) {
	for (j in 1:length(prevalence)) {
		for (k in 1:length(run)) {
			for (l in 1:(length(year)-1)) {
				population.year$dcount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$denergy[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$dfood[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
			}
		}
	}
}


####################################################################
####################################################################

library(mgvcv)

stable.state.noTA<-data.frame(condition=rep(c("nochange","symmetric"),each=11),prevalence=rep(seq(0,100,10),2),Estar=NA,Nstar=NA)

for (j in 1:dim(stable.state.noTA)[1]) {
sym.bi.gam<-gam(list(dcount~s(PorpoiseCount,PorpoiseEnergyLevel),denergy~s(PorpoiseCount,PorpoiseEnergyLevel)),family=mvn(d=2),data=population.year[population.year$condition==stable.state$condition[j] & population.year$prevalence==stable.state$prevalence[j],])


sym.newdf<-data.frame(dn=NA,de=NA,PorpoiseCount=rep(seq(100,230,1),each=length(seq(1700,4000,1))),PorpoiseEnergyLevel=rep(seq(1700,4000,1),length(seq(100,230,1))))
sym.newdf[,1:2]<-predict(sym.bi.gam,sym.newdf[,3:4])

isocline.n<-data.frame(abundance=100:230,energy=NA)
for (i in 1:dim(isocline.n)[1]) { 
isocline.n$energy[i]<-sym.newdf$PorpoiseEnergyLevel[sym.newdf$PorpoiseCount==isocline.n$abundance[i] & abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.n$abundance[i] ,]$de)==min(abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.n$abundance[i] ,]$de))]
}

isocline.e<-data.frame(abundance=100:230,energy=NA)
for (i in 1:dim(isocline.e)[1]) { 
isocline.e$energy[i]<-sym.newdf$PorpoiseEnergyLevel[sym.newdf$PorpoiseCount==isocline.e$abundance[i] & abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.e$abundance[i] ,]$dn)==min(abs(sym.newdf[sym.newdf$PorpoiseCount==isocline.e$abundance[i] ,]$dn))]
}

lm.e<-(lm(abundance~energy,data=isocline.e))
lm.n<-(lm(abundance~energy,data=isocline.n))
cm<-rbind(coef(lm.e),coef(lm.n))

stable.state.noTA[j,3:4]<-c(-solve(cbind(cm[,2],-1))%*%cm[,1])
print(j)
flush.console()

}

stable.state.noTA$treatment<-"noTA"
stable.state<-rbind(stable.state,stable.state.noTA)
save(stable.state,file="~/all_stablestate.Rdata")


ggplot(stable.state,aes(x=Nstar,y=Estar,colour=treatment,group=condition))+
geom_point(aes(shape=condition))


ggplot(subset(stable.state,treatment=="noTA"),aes(x=log10(Nstar),y=log10(Estar),group=condition,colour=factor(prevalence)))+
geom_point()+
geom_path()

ggplot(subset(stable.state,treatment=="TA"),aes(x=Nstar,y=Estar,group=condition,colour=factor(prevalence)))+
geom_point()+
geom_path()

population.year$treatment<-"noTA"
save(population.year,file="~/populationyearnoTA.Rdata")

population.year.noTA<-population.year
population.year.noTA<-population.year.noTA[population.year.noTA$condition=="nochange"|population.year.noTA$condition=="symmetric",]
rm(population.year)

#################################################################################################################################
#################################################################################################################################


load("~/resultsTA/count_all_yearly.Rdata")

run=unique(population.year$run)
prevalence=unique(population.year$prevalence)
condition=unique(population.year$condition)
year=unique(population.year$year)

population.year$dcount<-NA
population.year$denergy<-NA
population.year$dfood<-NA


for (i in 1:length(condition)) {
	for (j in 1:length(prevalence)) {
		for (k in 1:length(run)) {
			for (l in 1:(length(year)-1)) {
				population.year$dcount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$denergy[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$PorpoiseEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
				population.year$dfood[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]<-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l+1]]-
				population.year$FoodEnergyLevel[population.year$condition==condition[i]&population.year$prevalence==prevalence[j]&population.year$run==run[k]&population.year$year==year[l]]
				
			}
		}
	}
}

population.year$treatment<-"TA"
population.year<-rbind(population.year,population.year.noTA)

save(population.year,file="~/populationyear_all_sym_noch.Rdata")

####################################################################
####################################################################

load("~/all_stablestate.Rdata")
load("~/populationyear_all_sym_noch.Rdata")
####

#now we have the equilibria, we can look at resilience? let's look along the abundance axis

population.year$distN<-NA
condition<-unique(stable.state$condition)
prevalence<-unique(stable.state$prevalence)
treatment<-unique(stable.state$treatment)


for (i in 1:length(condition)) {
	for (j in 1:length(treatment)) {
		for (k in 1:length(prevalence)) {
		
		population.year$distN[population.year$condition==condition[i]&population.year$treatment==treatment[j]&population.year$prevalence==prevalence[k]]<-
		population.year$PorpoiseCount[population.year$condition==condition[i]&population.year$treatment==treatment[j]&population.year$prevalence==prevalence[k]]-
		stable.state$Nstar[stable.state$condition==condition[i]&stable.state$treatment==treatment[j]&stable.state$prevalence==prevalence[k]]
		
		}
		print(j)
		flush.console()
	}

}

population.year$speed<-sqrt(scale(population.year$dcount)^2+scale(population.year$denergy)^2)


nochnoTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="noTA"),aes(x=distN,y=dcount,colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)+
xlim(-50,60)+ylim(-30,25)

symnoTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="noTA"),aes(x=distN,y=dcount,colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)+
xlim(-50,60)+ylim(-30,25)

nochTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="TA"),aes(x=distN,y=dcount,colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)+
xlim(-50,60)+ylim(-30,25)

symTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="TA"),aes(x=distN,y=dcount,colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)+
xlim(-50,60)+ylim(-30,25)


grid.arrange(nochnoTA,symnoTA,nochTA,symTA,ncol=2,nrow=2)



#time to equilibrium



nochnoTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="noTA"),aes(x=year,y=abs(distN),colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)

symnoTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="noTA"),aes(x=year,y=abs(distN),colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)

nochTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="TA"),aes(x=year,y=abs(distN),colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)

symTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="TA"),aes(x=year,y=abs(distN),colour=factor(prevalence)))+
stat_summary(geom = "line", fun = median)

library(grid)
library(gridExtra)

grid.arrange(nochnoTA,symnoTA,nochTA,symTA,ncol=2,nrow=2)


#system speed


nochnoTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="noTA"),aes(x=distN,y=speed,colour=factor(prevalence)))+
stat_summary(geom = "point", fun = median)

symnoTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="noTA"),aes(x=distN,y=speed,colour=factor(prevalence)))+
stat_summary(geom = "point", fun = median)

nochTA<-ggplot(subset(population.year,condition=="nochange"&treatment=="TA"),aes(x=distN,y=speed,colour=factor(prevalence)))+
stat_summary(geom = "point", fun = median)

symTA<-ggplot(subset(population.year,condition=="symmetric"&treatment=="TA"),aes(x=distN,y=speed,colour=factor(prevalence)))+
stat_summary(geom = "point", fun = median)


grid.arrange(nochnoTA,symnoTA,nochTA,symTA,ncol=2,nrow=2)


glm1<-(glm(speed~poly(distN,2)*prevalence*treatment,data=population.year))

library(glmmTMB)
population.year$fprevalence<-factor(population.year$prevalence)
population.year$speed<-as.numeric(population.year$speed)
library(splines)

mem1<-glmmTMB(speed~poly(distN,2)*prevalence*treatment+(1|run),data=population.year)
mem2<-glmmTMB(speed~poly(distN,2)*prevalence*condition+(1|run),data=population.year)
mem3<-glmmTMB(speed~poly(distN,2)*prevalence*treatment*condition+(1|run),data=population.year)
## categorical prevalence continues to be a more useful representation

mem4<-glmmTMB(speed~poly(distN,2)*fprevalence*treatment*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem5<-glmmTMB(speed~poly(distN,2)*fprevalence*treatment+(1|run),data=population.year,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mem6<-glmmTMB(speed~poly(distN,2)*fprevalence*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem7<-glmmTMB(speed~poly(distN,2)*fprevalence+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem8<-glmmTMB(speed~poly(distN,2)*fprevalence*treatment+poly(distN,2)*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem9<-glmmTMB(speed~poly(distN,2)*fprevalence+poly(distN,2)*treatment+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem10<-glmmTMB(speed~poly(distN,2)*fprevalence+poly(distN,2)*treatment+poly(distN,2)*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem11<-glmmTMB(speed~poly(distN,2)*fprevalence*condition+poly(distN,2)*treatment+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem12<-glmmTMB(speed~poly(distN,2)*fprevalence+poly(distN,2)*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))
mem13<-glmmTMB(speed~poly(distN,2)*fprevalence*treatment+poly(distN,2)*fprevalence*condition+(1|run),data=population.year,control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


AIC(mem4,mem5,mem6,mem7,mem8,mem9,mem10,mem11,mem12,mem13)

      # df      AIC
# mem1  14 122159.1
# mem2  14 123168.3
# mem3  26 122146.4
# mem4 134 122012.5
# mem5  68 122017.9
# mem6  68 123054.0
# mem7  35 123073.4
# mem8  71 122011.2
# mem9  38 122288.9


       # df      AIC
# mem4  134 122012.5
# mem5   68 122017.9
# mem6   68 123054.0
# mem7   35 123073.4
# mem8   71 122011.2
# mem9   38 122288.9
# mem10  41 122289.4
# mem11  71 122275.1
# mem12  38 123067.5


library(ggeffects)
speed.pred<-ggpredict(mem13, term=c("distN [all]","fprevalence [0,10,30,40,80,90,100] ","treatment","condition"))

Anova(mem13)
summary(mem13)
spped.plot<-plot(speed.pred)
spped.plot&labs(y="system speed", x="distance from equilibria (# porpoises)",color="pinger prevalence",title="")&
theme(legend.position="right")&scale_colour_discrete(c('#feebe2','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177'))




speed.pred$facet<-as.character(speed.pred$facet)
speed.pred$facet[speed.pred$facet=="noTA"]<-"pingers only"
speed.pred$facet[speed.pred$facet=="TA"]<-"pingers & area closure"
speed.pred$facet<-factor(speed.pred$facet)

speed.pred$panel<-as.character(speed.pred$panel)
speed.pred$panel[speed.pred$panel=="nochange"]<-"no condition effects"
speed.pred$panel[speed.pred$panel=="symmetric"]<-"non-linear condition effects"
speed.pred$panel<-factor(speed.pred$panel)

speed.plot<-plot(speed.pred,line.size=.3,colors=c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4'))&
				labs(y="system speed", x="distance from equilibria (# porpoises)",color="pinger prevalence",title="")&
				theme(legend.position="right")
				

##
tiff(file="~/figures for ms/Figure xx engineering resilience2.tiff",width=10, height=5,units="in",res=300)
speed.plot
dev.off()
### 
