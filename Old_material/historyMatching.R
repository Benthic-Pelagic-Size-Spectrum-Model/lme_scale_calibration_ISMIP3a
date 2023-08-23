

# 'History Matching' or rejection algorithm for the benthic-pelagic size spectrum model.

# We want to assess which values of input parameters prodcue implausible model outputs and remove these parameter values in further analyses via the emulator

# We use some basic constraints:


# 1. We want the modelled growth rates to fall within the observed range of empirical ( von Bertalanffy) growth rates for fish

# 2. The observed catch rates should fall within the [min, max] of long-term empirical catch rates [0,300] from Reg's dataset

# 3. The size spectra slopes should be < 0


# If any one of these constraints are NOT met the parameter set is IMPLAUSIBLE
 
# # # # # # # # # # # # # # # # # # # # # # # # #


# # Check empirical growth rates - only need to do this once to get the range of values we need.

# # Read in file for empirical growth rates (yr^-1), exclude the first row because it seems wrong

pauly<-read.csv("Pauly1980_VBparams.csv")[-1,]

t<-seq(1:(365*50))/365

# # which t gives weight at 20% of Winf?

# # calculate weight at t

w<-function(Winf,K,t){return(Winf*(1-exp(-K*(t-0)))^3)}

# calculate growth rate using derivative
dwdt<-function(Winf,K,t){return(Winf*(3*(exp(-K*t-0))*K*(1-exp(-K*(t-0)))^2))}

pauly$maxdwdt <- 0
pauly$w.maxdwdt <- 0

# plot(t,w(pauly[2,4],pauly[2,5],t),typ="l") 

for (i in 1: length(pauly[,1])) {

grora<-dwdt(pauly[i,4],pauly[i,5],t)

weight<-w(pauly[i,4],pauly[i,5],t)

pauly$maxdwdt[i] <-max(grora)

pauly$w.maxdwdt[i] <-weight[which(grora==max(grora))]

}

# points(t,w(pauly[i,4],pauly[i,5],t),typ="l") 

#check plots

# plot(pauly$w.maxdwdt,pauly$maxdwdt,xlab="weight (g)",ylab= "Maximum growth rate (grams per year)",log="xy")
plot(pauly$w.maxdwdt,pauly$maxdwdt/pauly$w.maxdwdt,xlab="weight (g)",ylab= "Maximum relative growth rate (per year)",log="x")

# # negative relationship, relative growth rate maximum value of 7, minimum of 0.05 g per year (let's say 0.2) AT 100 GRAMS!

# max(pauly$maxdwdt/pauly$w.maxdwdt)
# min(pauly$maxdwdt/pauly$w.maxdwdt)
 
# max(pauly$w.maxdwdt)
# min(pauly$w.maxdwdt) 

load(file="~/Desktop/model_outputs/thetas.RData")
thetas<-as.data.frame(thetas)

thetas$implausibility = rep(NA,length(thetas[,1]))

save(thetas,file="~/Desktop/model_outputs/thetasHM.RData")
notrun<-1:1000

files<-list.files(path = "~/Desktop/model_outputs", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
           
for (i in 1:1000) notrun[i]<-ifelse(sum(is.element(files,paste(i,".RData",sep="")))>0,notrun[i],NA)
 
rejection<-function(modelfile=paste("~/Desktop/model_outputs/",num,".RData",sep=""), inputparams = "~/Desktop/model_outputs/thetasHM.RData") {
	
load(paste(modelfile))	
load(paste(inputparams))

# if the model was run set implausibility to 1 and multiply as each condition of "implausibilty" is met	(if model not run implausibilty will be NA)
	
thetas$implausibility[num] <- 1	
	
# 1. growth rates

#com<-output[["res"]]
com<-initial_results
params<-initial_results$param

# get modelled relative growth rates (yr^-1)

t = seq(0,3,by=0.1)
ref.max<-length(params$x)
idx<-params$ref:ref.max
wt<-10^params$x[idx]
ts<-params$Neq

gg = com$GG.u[idx,ts] ## Final growth rate

par(mfrow=c(1,1))

plot(pauly$w.maxdwdt,pauly$maxdwdt/pauly$w.maxdwdt,xlab="weight (g)",ylab= "Maximum relative growth rate (per year)",log="xy")
plot(wt[which(wt>0.1&wt<60000)],gg[which(wt>0.1&wt<60000)],log="xy",xlab="weight (g)",ylab="relative growth rate per (yr)",typ="l")

if (max(gg) > 7 ) thetas$implausibility[num]  <- thetas$implausibility[num]*3

if (max(gg) < 0.05 ) thetas$implausibility[num]  <- thetas$implausibility[num]*3

if (gg[which(wt==10)] < gg[which(wt==100)] ) thetas$implausibility[num] <- thetas$implausibility[num]*3


# 2. catch rates

yield = com$Y.u[idx,ts] ## final yield grams per yr per m3

# calculate catch rate t.km^-2.yr^-1 

# convert to tonnes per km^2, use max observed from Reg's catch database

yield = ((yield/1e6)/params$depth)/1e6

yield <- sum(yield*params$dx)
 
if (yield > 300) thetas$implausibility[num] <-  thetas$implausibility[num]*5


# 3. size spectrum 

# check size spectra have negative slopes

pelben = com$U[idx,100] + com$V[idx,100] ## final density 

plot(wt,pelben, ylab="density",xlab="weight (g)",log="xy",typ="l")
slope<-lm(log10(pelben[which(wt> 1 & wt<5000)])~log10(wt[which(wt> 1 & wt<5000)]))

if (pel[which(wt==0.01)] < pel[which(wt==100)] ) thetas$implausibility[num] <- thetas$implausibility[num]*7


save(thetas,file="~/Desktop/model_outputs/thetasHM.RData")

 }



# run rejection algorithm for all saved output files and add column of "reject" variable to input table
for (num in notrun[!is.na(notrun)]) rejection(modelfile=paste("~/Desktop/model_outputs/",num,".RData",sep=""), inputparams = "~/Desktop/model_outputs/thetasHM.RData")

load(file="~/Desktop/model_outputs/thetasHM.RData")
par(mfrow=c(1,1))
plot(1:1000,thetas$implausibility,xlab="theta",ylab="implausibility",typ="l")
 
 
