install.packages("DiagrammeR")
install.packages("data.tree")
install.packages("Rtools")
install.packages("heemod")
install.packages("rriskDistributions")
library("data.tree")
library("heemod")
library("shape")
library("diagram")
library("ggplot2")
library("flexsurv")
library("dplyr")
library("DiagrammeR")
library("rriskDistributions")
library("MASS")

#PSA iterations
psait<-1000

#Generate prob vectors
vmutationprev<-rbeta(psait,101,285) #Pan et al
vrrclop<-(rbeta(psait,10,275))/(rbeta(psait,19,367)) #Pan et al
vclopstroke1<-(rnorm(psait,3.971,0.03))/100 #NICE
vmrdstroke1<-(rnorm(psait,3.971,0.03))/100 #NICE
vasastroke1<-(rnorm(psait,4.201,0.03))/100 #NICE
vclopstrokeextra<-(rnorm(psait,7.323,0.69))/100 #NICE
vmrdstrokeextra<-(rnorm(psait,7.323,0.69))/100 #NICE
vasastrokeextra<-(rnorm(psait,7.323,0.69))/100 #NICE
vdisablingclop<-rbeta(psait,392,506) #PROFESS
vdisablingmrd<-rbeta(psait,413,503) #PROFESS
vdisablingasa<-rbeta(psait,413,503) #PROFESS
vcloptolerance<-1-(rbeta(psait,1069,8971)) # PROFESS drug discontinuation
vmrdasatolerance<-1-(rbeta(psait,1650,8405)) #PROFESS drug discontinuation
vclopminorbleed<-runif(psait,0.00837,0.01023)
vmrdminorbleed<-runif(psait,0.00783,0.00957)
vasaminorbleed<-runif(psait,0.00837,0.01023)
vclopmajorbleed<-runif(psait,0.00369,0.00451)
vmrdmajorbleed<-runif(psait,0.00414,0.00506)
vasamajorbleed<-runif(psait,0.00486,0.00594)
vclopchf<-runif(psait,0.00675,0.00825)
vmrdchf<-runif(psait,0.00567,0.00693)
vasachf<-runif(psait,0.00567,0.00693)

vtestcost<-runif(psait,40,80) #Between ?40 and ?80
vtestresourcecost<-(47*(runif(psait,1,20)/60))+(109*(runif(psait,1,10)/60)) # 1 to 20 minutes nurse time, 1 to 10 mins consultant time
vcostinflator<-rgamma(psait,shape=17504.8,rate=17195.57)^10
vcoststroke<-runif(psait,1517.44,1854.64)
vcostdisablingstroke<-runif(psait,4657.90,5692.98)
vcostfatal<-runif(psait,7890.92,9644.46)
vcostfatalnonstroke<-runif(psait,2002.5,2447.5)
vnondisablestrokeevent<-runif(psait,5768.95,7050.93)
vdisablestrokeevent<-runif(psait,12282.64,15012.12)
vcostclop<-runif(psait,28.8,35.2)
vcostmrd<-runif(psait,73.8,90.2)
vcostasa<-runif(psait,19.8,24.2)
vcostaeclop<-runif(psait,18.09,22.11)
vcostaemrd<-runif(psait,23.56,28.80)
costaeasa<-runif(psait,19.87,24.29)

vqolstroke<-rbeta(psait,45.65694,29.20021) #assumed deviation of 0.1 qol either way based on Post et al and fit beta distribution
vqolminorbleed<-runif(psait,0.00297,0.00363)
vqolmajorbleed<-runif(psait,0.12834,0.15686)
vqolchf<-runif(psait,0.01467,0.01793)

param<- define_parameters(
  ACM=acm[markov_cycle],
  vascdeath1=0.00212*(exp(0.0520*markov_cycle))*0.0686*0.791,
  vascdeath2=0.00212*(exp(0.0520*markov_cycle))*0.0686*1.931,
  vascdeath3=0.00212*(exp(0.0520*markov_cycle))*0.0686*4.398
)

#general parameters
N<-60424 #based on crude incidence of 1.07 (2016) and current UK population
cycle_length<-34
dr<-0.035
acm<-c(0.0119775, 0.013126, 0.0142565, 0.0156355, 0.0172055, 
       0.019497, 0.0216465, 0.023874, 0.0271235, 0.0304265, 0.033783,
       0.0375525, 0.041461, 0.0465925, 0.052426, 0.058713, 0.067142, 
       0.0756975, 0.084887, 0.096117, 0.1078025, 0.1217215, 0.13626, 
       0.15048, 0.167086, 0.184636, 0.203201, 0.222409, 0.246654, 
       0.267129, 0.2882125, 0.306311, 0.340515, 0.363598)
costinflator<-vcostinflator[i]
threshold<-c(5000,10000,20000,30000,50000,100000)

draw<-c(1:psait)
psa<-data.frame(Draw=draw,inccosts=numeric(length(draw)),
                incqual=numeric(length(draw)),ICER=numeric(length(draw)),
                CVOIMP=numeric(length(draw)),NBpos=numeric(length(draw)),
                qalyint=numeric(length(draw)),costint=numeric(length(draw)),
                qalycomp=numeric(length(draw)),costcomp=numeric(length(draw)))
for(i in draw){

#Decision Tree Parameters
mutationprev<-vmutationprev[i]
cloptolerance<-vcloptolerance[i]
mrdasatolerance<-vmrdasatolerance[i]
testcost<-vtestcost[i]
testresourcecost<-vtestresourcecost[i]

#pgxclop parameters
rrclop<-vrrclop[i]

#clop parameters
clopstroke1<-vclopstroke1[i]
clopstrokeextra<-vclopstrokeextra[i]
disablingclop<-vdisablingclop[i]

#MRD parameters
mrdstroke1<-vmrdstroke1[i]
mrdstrokeextra<-vmrdstrokeextra[i]
disablingmrd<-vdisablingmrd[i]

#ASA parameters
asastroke1<-vasastroke1[i]
asastrokeextra<-vasastrokeextra[i]
disablingasa<-vdisablingasa[i]

#costs
coststroke<-vcoststroke[i]*costinflator
costdisablingstroke<-vcostdisablingstroke[i]*costinflator
costfatal<-vcostfatal[i]*costinflator
costfatalnonstroke<-vcostfatalnonstroke[i]*costinflator
nondisablestrokeevent<-vnondisablestrokeevent[i]*costinflator
disablestrokeevent<-vdisablestrokeevent[i]*costinflator
costclop<-vcostclop[i]
costmrd<-vcostmrd[i]
costasa<-vcostasa[i]
costaeclop<-vcostaeclop[i]*costinflator
costaemrd<-vcostmrd[i]*costinflator
costaeasa<-vcostasa[i]*costinflator
  
#qol
qolstroke<-vqolstroke[i]
qolstroke1<-qolstroke-0.174
qolstroke2<-qolstroke1

qolminorbleed<-vqolminorbleed[i]
qolmajorbleed<-vqolmajorbleed[i]
qolchf<-vqolchf[i]

clopaeqol<-(vqolminorbleed[i]*vclopminorbleed[i])+(vqolmajorbleed[i]*vclopmajorbleed[i])+(vqolchf[i]*vclopchf[i])
mrdaeqol<-(vqolminorbleed[i]*vmrdminorbleed[i])+(vqolmajorbleed[i]*vmrdmajorbleed[i])+(vqolchf[i]*vmrdchf[i])
asaaeqol<-(vqolminorbleed[i]*vasaminorbleed[i])+(vqolmajorbleed[i]*vasamajorbleed[i])+(vqolchf[i]*vasachf[i])

#Decision tree

genotype<-Node$new("CYP2C19 Genotyping")
 genotyping<-genotype$AddChild("Genotyping",cost=testcost+testresourcecost)
 nogenotyping<-genotype$AddChild("No Genotyping",cost=0)
  pgxclopidogrel<-genotyping$AddChild("Clopidogrel",p=1-mutationprev)
  pgxmrd<-genotyping$AddChild("MRD+ASA",p=mutationprev)
  clopidogrel<-nogenotyping$AddChild("Clopidogrel")
   pgxcloptol<-pgxclopidogrel$AddChild("Clopidogrel tolerated",p=cloptolerance)
   pgxclopnotol<-pgxclopidogrel$AddChild("MRD +ASA",p=1-cloptolerance)
   pgxmrdtol<-pgxmrd$AddChild("MRD + ASA tolerated",p=mrdasatolerance)
   pgxmrdnotol<-pgxmrd$AddChild("ASA",p=1-mrdasatolerance)
   cloptol<-clopidogrel$AddChild("Clopidogrel tolerated",p=cloptolerance)
   clopnotol<-clopidogrel$AddChild("MRD + ASA",p=1-cloptolerance)
    mrdasatol<-pgxclopnotol$AddChild("MRD + ASA tolerated",p=mrdasatolerance)
    mrdasanotol<-pgxclopnotol$AddChild("ASA",p=1-mrdasatolerance)
    mrdasatol2<-clopnotol$AddChild("MRD + ASA tolerated",p=mrdasatolerance)
    mrdasanotol2<-clopnotol$AddChild("ASA",p=1-mrdasatolerance)

#Produce Vector of Conditional Probabilities and Total Costs
rollbacka<-data.frame(prob=pgxcloptol$Get("p",traversal="ancestor"),cost=pgxcloptol$Get("cost",traversal="ancestor"))
rollbackb<-data.frame(prob=mrdasatol$Get("p",traversal="ancestor"),cost=mrdasatol$Get("cost",traversal="ancestor"))
rollbackc<-data.frame(prob=mrdasanotol$Get("p",traversal="ancestor"),cost=mrdasanotol$Get("cost",traversal="ancestor"))
rollbackd<-data.frame(prob=pgxmrdtol$Get("p",traversal="ancestor"),cost=pgxmrdtol$Get("cost",traversal="ancestor"))
rollbacke<-data.frame(prob=pgxmrdnotol$Get("p",traversal="ancestor"),cost=pgxmrdnotol$Get("cost",traversal="ancestor"))
rollbackf<-data.frame(prob=cloptol$Get("p",traversal="ancestor"),cost=cloptol$Get("cost",traversal="ancestor"))
rollbackg<-data.frame(prob=mrdasatol2$Get("p",traversal="ancestor"),cost=mrdasatol2$Get("cost",traversal="ancestor"))
rollbackh<-data.frame(prob=mrdasanotol2$Get("p",traversal="ancestor"),cost=mrdasanotol2$Get("cost",traversal="ancestor"))

condprobs=c(prod(rollbacka[1:2,1]),prod(rollbackb[1:3,1]),prod(rollbackc[1:3,1]),prod(rollbackd[1:2,1]),prod(rollbacke[1:2,1]),prod(rollbackf[1:1,1]),prod(rollbackg[1:2,1]),prod(rollbackh[1:2,1]))
costtest<-c(prod(rollbacka[3,2]),prod(rollbackb[4,2]),prod(rollbackc[4,2]),prod(rollbackd[3,2]),prod(rollbacke[3,2]),prod(rollbackf[3,2]),prod(rollbackg[4,2]),prod(rollbackh[4,2]))
patientflow<-N*condprobs
totalcosts<-patientflow*testcost
totaltestcost<-sum(totalcosts)

#Markov Models
mat_pgxclop<-define_transition(
  state_names=c("No further stroke","1 additional stroke",">1 additional stroke","Vascular Death","Non-vascular death"), 
  C, rrclop*clopstroke1, 0, vascdeath1, ACM,
  0, C, clopstrokeextra, vascdeath2, ACM,
  0, 0, C, vascdeath3, ACM,
  0, 0, 0, 1, 0,
  0, 0, 0, 0, 1)

mat_clop<-define_transition(
  state_names=c("No further stroke","1 additional stroke",">1 additional stroke","Vascular Death","Non-vascular death"), 
  C, clopstroke1, 0, vascdeath1, ACM,
  0, C, clopstrokeextra, vascdeath2, ACM,
  0, 0, C, vascdeath3, ACM,
  0, 0, 0, 1, 0,
  0, 0, 0, 0, 1)

mat_mrd<-define_transition(
  state_names=c("No further stroke","1 additional stroke",">1 additional stroke","Vascular Death","Non-vascular death"), 
  C, mrdstroke1, 0, vascdeath1, ACM,
  0, C, mrdstrokeextra, vascdeath2, ACM,
  0, 0, C, vascdeath3, ACM,
  0, 0, 0, 1, 0,
  0, 0, 0, 0, 1)

mat_asa<-define_transition(
  state_names=c("No further stroke","1 additional stroke",">1 additional stroke","Vascular Death","Non-vascular death"), 
  C, asastroke1, 0, vascdeath1, ACM,
  0, C, asastrokeextra, vascdeath2, ACM,
  0, 0, C, vascdeath3, ACM,
  0, 0, 0, 1, 0,
  0, 0, 0, 0, 1)

state_nfstroke<-define_state(
  cost_health=discount(dispatch_strategy(
pgxclop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
clop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
mrd=(coststroke*(1-disablingmrd))+(costdisablingstroke*disablingmrd),
asa=(coststroke*(1-disablingasa))+(costdisablingstroke*disablingasa)), dr),
  cost_drugs=discount(dispatch_strategy(
    pgxclop=costclop,
    clop=costclop,
    mrd=costmrd,
    asa=costasa),dr),
  cost_ae=discount(dispatch_strategy(
    pgxclop=costaeclop,
    clop=costaeclop,
    mrd=costaemrd,
    asa=costaeasa
  ),dr),
  cost_total=cost_health+cost_drugs+cost_ae,
  qol=discount(dispatch_strategy(
    pgxclop=qolstroke-clopaeqol,
    clop=qolstroke-clopaeqol,
    mrd=qolstroke-mrdaeqol,
    asa=qolstroke-asaaeqol), dr))

state_stroke1<-define_state(
  cost_health=discount(dispatch_strategy(
    pgxclop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
    clop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
    mrd=(coststroke*(1-disablingmrd))+(costdisablingstroke*disablingmrd),
    asa=(coststroke*(1-disablingasa))+(costdisablingstroke*disablingasa)), dr),
  cost_drugs=discount(dispatch_strategy(
    pgxclop=costclop,
    clop=costclop,
    mrd=costmrd,
    asa=costasa),dr),
  cost_ae=discount(dispatch_strategy(
    pgxclop=costaeclop,
    clop=costaeclop,
    mrd=costaemrd,
    asa=costaeasa
  ),dr),
  cost_total=cost_health+cost_drugs+cost_ae,
  qol=discount(dispatch_strategy(
    pgxclop=qolstroke1-clopaeqol,
    clop=qolstroke1-clopaeqol,
    mrd=qolstroke1-mrdaeqol,
    asa=qolstroke1-asaaeqol), dr))

state_stroke2<-define_state(
  cost_health=discount(dispatch_strategy(
    pgxclop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
    clop=(coststroke*(1-disablingclop))+(costdisablingstroke*disablingclop),
    mrd=(coststroke*(1-disablingmrd))+(costdisablingstroke*disablingmrd),
    asa=(coststroke*(1-disablingasa))+(costdisablingstroke*disablingasa)), dr),
  cost_drugs=discount(dispatch_strategy(
    pgxclop=costclop,
    clop=costclop,
    mrd=costmrd,
    asa=costasa),dr),
  cost_ae=discount(dispatch_strategy(
    pgxclop=costaeclop,
    clop=costaeclop,
    mrd=costaemrd,
    asa=costaeasa
  ),dr),
  cost_total=cost_health+cost_drugs+cost_ae,
  qol=discount(dispatch_strategy(
    pgxclop=qolstroke2-clopaeqol,
    clop=qolstroke2-clopaeqol,
    mrd=qolstroke2-mrdaeqol,
    asa=qolstroke2-asaaeqol), dr))

state_vasculardeath<-define_state(
  cost_health=0,
  cost_drugs=0,
  cost_ae=0,
  cost_total=cost_health+cost_drugs+cost_ae,
  qol=0)

state_death<-define_state(
  cost_health=0,
  cost_drugs=0,
  cost_ae=0,
  cost_total=cost_health+cost_drugs+cost_ae,
  qol=0)

strat_pgxclop<-define_strategy(
  transition=mat_pgxclop,
  "No further stroke"=state_nfstroke,
  "1 additional stroke"=state_stroke1,
  ">1 additional stroke"=state_stroke2,
  "Vascular Death"=state_vasculardeath,
  "Non-vascular death"=state_death)

strat_clop<-define_strategy(
  transition=mat_clop,
  "No further stroke"=state_nfstroke,
  "1 additional stroke"=state_stroke1,
  ">1 additional stroke"=state_stroke2,
  "Vascular Death"=state_vasculardeath,
  "Non-vascular death"=state_death)

strat_mrd<-define_strategy(
  transition=mat_mrd,
  "No further stroke"=state_nfstroke,
  "1 additional stroke"=state_stroke1,
  ">1 additional stroke"=state_stroke2,
  "Vascular Death"=state_vasculardeath,
  "Non-vascular death"=state_death)

strat_asa<-define_strategy(
  transition=mat_asa,
  "No further stroke"=state_nfstroke,
  "1 additional stroke"=state_stroke1,
  ">1 additional stroke"=state_stroke2,
  "Vascular Death"=state_vasculardeath,
  "Non-vascular death"=state_death)

res_modpgxclop<-run_model(
  pgxclop=strat_pgxclop,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[1],0,0,0,0)
)

res_modpgxmrd<-run_model(
  mrd=strat_mrd,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[2]+patientflow[4],0,0,0,0)
)

res_modpgxasa<-run_model(
  asa=strat_asa,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[3]+patientflow[5],0,0,0,0)
)

res_modclop<-run_model(
  clop=strat_clop,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[6],0,0,0,0)
)

res_modmrd<-run_model(
  mrd=strat_mrd,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[7],0,0,0,0)
)

res_modasa<-run_model(
  asa=strat_asa,
  cycles=cycle_length,
  parameters = param,
  cost=cost_total,
  effect=qol,
  init=c(patientflow[8],0,0,0,0)
)

t<-c(1:34)
firststrokes<-data.frame(get_counts(res_modpgxclop)[1:34,4]*rrclop*clopstroke1,get_counts(res_modpgxmrd)[1:34,4]*mrdstroke1,get_counts(res_modpgxasa)[1:34,4]*asastroke1,get_counts(res_modclop)[1:34,4]*clopstroke1,get_counts(res_modmrd)[1:34,4]*mrdstroke1,get_counts(res_modasa)[1:34,4]*asastroke1)
compstroke1<-data.frame(t,firststrokes[,1],firststrokes[,2],firststrokes[,3],firststrokes[,4],firststrokes[,5],firststrokes[,6])
strokecost1<-data.frame(compstroke1[,2]*((nondisablestrokeevent*(1-disablingclop)+(disablestrokeevent*disablingclop))/((1+dr)^t)),compstroke1[,3]*((nondisablestrokeevent*(1-disablingmrd)+(disablestrokeevent*disablingmrd))/((1+dr)^t)),compstroke1[,4]*((nondisablestrokeevent*(1-disablingasa)+(disablestrokeevent*disablingasa))/((1+dr)^t)),compstroke1[,5]*((nondisablestrokeevent*(1-disablingclop)+(disablestrokeevent*disablingclop))/((1+dr)^t)),compstroke1[,6]*((nondisablestrokeevent*(1-disablingmrd)+(disablestrokeevent*disablingmrd))/((1+dr)^t)),compstroke1[,7]*((nondisablestrokeevent*(1-disablingasa)+(disablestrokeevent*disablingasa))/((1+dr)^t)))
pgxclopstrokeevent1<-sum(strokecost1[,1])
pgxmrdstrokeevent1<-sum(strokecost1[,2])
pgxasastrokeevent1<-sum(strokecost1[,3])
clopstrokeevent1<-sum(strokecost1[,4])
mrdstrokeevent1<-sum(strokecost1[,5])
asastrokeevent1<-sum(strokecost1[,6])

secondstrokes<-data.frame(get_counts(res_modpgxclop)[35:68,4]*rrclop*clopstroke1,get_counts(res_modpgxmrd)[35:68,4]*mrdstroke1,get_counts(res_modpgxasa)[35:68,4]*asastroke1,get_counts(res_modclop)[35:68,4]*clopstroke1,get_counts(res_modmrd)[35:68,4]*mrdstroke1,get_counts(res_modasa)[35:68,4]*asastroke1)
compstroke2<-data.frame(t,secondstrokes[,1],secondstrokes[,2],secondstrokes[,3],secondstrokes[,4],secondstrokes[,5],secondstrokes[,6])
strokecost2<-data.frame(compstroke2[,2]*((nondisablestrokeevent*(1-disablingclop)+(disablestrokeevent*disablingclop))/((1+dr)^t)),compstroke2[,3]*((nondisablestrokeevent*(1-disablingmrd)+(disablestrokeevent*disablingmrd))/((1+dr)^t)),compstroke2[,4]*((nondisablestrokeevent*(1-disablingasa)+(disablestrokeevent*disablingasa))/((1+dr)^t)),compstroke2[,5]*((nondisablestrokeevent*(1-disablingclop)+(disablestrokeevent*disablingclop))/((1+dr)^t)),compstroke2[,6]*((nondisablestrokeevent*(1-disablingmrd)+(disablestrokeevent*disablingmrd))/((1+dr)^t)),compstroke2[,7]*((nondisablestrokeevent*(1-disablingasa)+(disablestrokeevent*disablingasa))/((1+dr)^t)))
pgxclopstrokeevent2<-sum(strokecost2[,1])
pgxmrdstrokeevent2<-sum(strokecost2[,2])
pgxasastrokeevent2<-sum(strokecost2[,3])
clopstrokeevent2<-sum(strokecost2[,4])
mrdstrokeevent2<-sum(strokecost2[,5])
asastrokeevent2<-sum(strokecost2[,6])

deaths<-data.frame(get_counts(res_modpgxclop)[103:136,4],get_counts(res_modpgxmrd)[103:136,4],get_counts(res_modpgxasa)[103:136,4],get_counts(res_modclop)[103:136,4],get_counts(res_modmrd)[103:136,4],get_counts(res_modasa)[103:136,4])
compdeaths<-data.frame(deaths[,1],deaths[,2],deaths[,3],deaths[,4],deaths[,5],deaths[,6])
t<-c(1:33)
incdeaths<-data.frame(t,compdeaths[2:34,1]-compdeaths[1:33,1],compdeaths[2:34,2]-compdeaths[1:33,2],compdeaths[2:34,3]-compdeaths[1:33,3],compdeaths[2:34,4]-compdeaths[1:33,4],compdeaths[2:34,5]-compdeaths[1:33,5],compdeaths[2:34,6]-compdeaths[1:33,6])
fatalcost<-data.frame(incdeaths[,2]*((costfatal)/((1+dr)^t)),incdeaths[,3]*(costfatal/((1+dr)^t)),incdeaths[,4]*(costfatal/((1+dr)^t)),incdeaths[,5]*(costfatal/((1+dr)^t)),incdeaths[,6]*(costfatal/((1+dr)^t)),incdeaths[,7]*(costfatal/((1+dr)^t)))
pgxclopfatal<-sum(fatalcost[,1])
pgxmrdfatal<-sum(fatalcost[,2])
pgxasafatal<-sum(fatalcost[,3])
clopfatal<-sum(fatalcost[,4])
mrdfatal<-sum(fatalcost[,5])
asafatal<-sum(fatalcost[,6])

nonstrokedeaths<-data.frame(get_counts(res_modpgxclop)[137:170,4],get_counts(res_modpgxmrd)[137:170,4],get_counts(res_modpgxasa)[137:170,4],get_counts(res_modclop)[137:170,4],get_counts(res_modmrd)[137:170,4],get_counts(res_modasa)[137:170,4])
compnonstrokedeaths<-data.frame(nonstrokedeaths[,1],nonstrokedeaths[,2],nonstrokedeaths[,3],nonstrokedeaths[,4],nonstrokedeaths[,5],nonstrokedeaths[,6])
t<-c(1:33)
incnonstrokedeaths<-data.frame(t,compnonstrokedeaths[2:34,1]-compnonstrokedeaths[1:33,1],compnonstrokedeaths[2:34,2]-compnonstrokedeaths[1:33,2],compnonstrokedeaths[2:34,3]-compnonstrokedeaths[1:33,3],compnonstrokedeaths[2:34,4]-compnonstrokedeaths[1:33,4],compnonstrokedeaths[2:34,5]-compnonstrokedeaths[1:33,5],compnonstrokedeaths[2:34,6]-compnonstrokedeaths[1:33,6])
fatalnonstrokecost<-data.frame(incnonstrokedeaths[,2]*((costfatalnonstroke)/((1+dr)^t)),incnonstrokedeaths[,3]*(costfatalnonstroke/((1+dr)^t)),incnonstrokedeaths[,4]*(costfatalnonstroke/((1+dr)^t)),incnonstrokedeaths[,5]*(costfatalnonstroke/((1+dr)^t)),incnonstrokedeaths[,6]*(costfatalnonstroke/((1+dr)^t)),incnonstrokedeaths[,7]*(costfatalnonstroke/((1+dr)^t)))
pgxclopfatal2<-sum(fatalnonstrokecost[,1])
pgxmrdfatal2<-sum(fatalnonstrokecost[,2])
pgxasafatal2<-sum(fatalnonstrokecost[,3])
clopfatal2<-sum(fatalnonstrokecost[,4])
mrdfatal2<-sum(fatalnonstrokecost[,5])
asafatal2<-sum(fatalnonstrokecost[,6])

summarypgxclop<-summary(res_modpgxclop)
tc_tbpgxclop<-c(summarypgxclop[[1]][4],summarypgxclop[[1]][5])
tc_tbpgxclop<-unlist(tc_tbpgxclop)
tc_tbpgxclop<-c(tc_tbpgxclop[1]+pgxclopfatal+pgxclopfatal2+pgxclopstrokeevent1+pgxclopstrokeevent2,tc_tbpgxclop[2])

summarypgxmrd<-summary(res_modpgxmrd)
tc_tbpgxmrd<-c(summarypgxmrd[[1]][4],summarypgxmrd[[1]][5])
tc_tbpgxmrd<-unlist(tc_tbpgxmrd)
tc_tbpgxmrd<-c(tc_tbpgxmrd[1]+pgxmrdfatal+pgxmrdfatal2+pgxmrdstrokeevent1+pgxmrdstrokeevent2,tc_tbpgxmrd[2])

summarypgxasa<-summary(res_modpgxasa)
tc_tbpgxasa<-c(summarypgxasa[[1]][4],summarypgxasa[[1]][5])
tc_tbpgxasa<-unlist(tc_tbpgxasa)
tc_tbpgxasa<-c(tc_tbpgxasa[1]+pgxasafatal+pgxasafatal2+pgxasastrokeevent1+pgxasastrokeevent2,tc_tbpgxasa[2])

summaryclop<-summary(res_modclop)
tc_tbclop<-c(summaryclop[[1]][4],summaryclop[[1]][5])
tc_tbclop<-unlist(tc_tbclop)
tc_tbclop<-c(tc_tbclop[1]+clopfatal+clopfatal2+clopstrokeevent1+clopstrokeevent2,tc_tbclop[2])

summarymrd<-summary(res_modmrd)
tc_tbmrd<-c(summarymrd[[1]][4],summarymrd[[1]][5])
tc_tbmrd<-unlist(tc_tbmrd)
tc_tbmrd<-c(tc_tbmrd[1]+mrdfatal+mrdfatal2+mrdstrokeevent1+mrdstrokeevent2,tc_tbmrd[2])

summaryasa<-summary(res_modasa)
tc_tbasa<-c(summaryasa[[1]][4],summaryasa[[1]][5])
tc_tbasa<-unlist(tc_tbasa)
tc_tbasa<-c(tc_tbasa[1]+asafatal+asafatal2+asastrokeevent1+asastrokeevent2,tc_tbasa[2])

tc_tbintervention<-tc_tbpgxclop+tc_tbpgxmrd+tc_tbpgxasa
tc_tbintervention<-c(tc_tbintervention[1]+totaltestcost,tc_tbintervention[2])
tc_tbcomparator<-tc_tbclop+tc_tbmrd+tc_tbasa

inctc_tb=tc_tbintervention-tc_tbcomparator
ICER<-inctc_tb[[1]]/inctc_tb[[2]]
CVOIMP<-c((threshold*inctc_tb[2])-inctc_tb[1])
names(CVOIMP)<-threshold

psa[i,2]<-inctc_tb[1]/N
psa[i,3]<-inctc_tb[2]/N
psa[i,4]<-ICER
psa[i,5]<-CVOIMP[3]
psa[i,6]<-tc_tbintervention[[2]]/N
psa[i,7]<-tc_tbintervention[[1]]/N
psa[i,8]<-tc_tbcomparator[[2]]/N
psa[i,9]<-tc_tbcomparator[[1]]/N
print(i)  
}

mean(psa[,6])
mean(psa[,7])
mean(psa[,8])
mean(psa[,9])


plot(psa[,3],psa[,2],type="p",xlim=c(-0.3,0.3),ylim=c(-3000,5000),
     xlab=("Incremental Benefit of Testing (QALYs)"),
     ylab=("Incremental Cost of Testing (?)"))
abline(a=c(0,20000))
abline(v=0,h=0)

psa[,10]<-ifelse(psa[,5]>0,1,0)
mean(psa[,10])
plot(density(psa[,5]))

costsaving<-psa[,2]
d<-density(costsaving)
plot(d)
psa[,11]<-ifelse(psa[,2]<0,1,0)
mean(psa[,11])

loss<-c(ifelse(psa[,5]<0,psa[,5],0))
EVPI<--1*mean(loss)
EVPI/N
EVPI

mean(psa[,3])
mean(psa[,8])
mean(psa[,6])
