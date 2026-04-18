# Meta-Analysis-Course
The code for the meta-analysis workshop

dowload the code in R file named 

"Workshop meta with forrest and funnel plots.R"




alternativelly here is the full code so it can be copied directly to R


#CODE

# install packages
install.packages(c("meta", "metafor", "tidyverse", "lattice"), 
                 repos = "https://cran.at.r-project.org/")

# Standard loading
library(meta)
library(metafor)
library(tidyverse)
library (lattice)

#if necessary to list and remove existing lists of data
# ls()
# rm(list=ls())

# generate data for each study / generisem podatke za svaku studiju

# Number of plots per each treatment
n = 50
# Mean yield in general
yield.mu = 50
# Standard deviations 
yield.sd = 10

yield.muTRT = 55

# XXXX radim bez ovoga, nema promjena, jer se radi samo prinos / Mean changes for yield ommited
# mu.d = 10


# Mean age for trees (plots)
age.mu = 10
# sd of age for individual trees
age.sd = 1

# Fix the seed for random number generation - necessary for consistent result, OBLIGATORY run it
set.seed(123)
####generate funcion for simulation of data

data.generator = function(n,age.mu,age.sd,yield.mu,yield.muTRT,yield.sd,center){
# Data from CTRL
age = rnorm(n, age.mu, age.sd)
prinos = rnorm(n,yield.mu,yield.sd)
dat4CTRL = round(cbind(age,prinos))
# Data from Treatment
age = rnorm(n, age.mu, age.sd)
prinos = rnorm(n,yield.muTRT,yield.sd)
dat4drug = round(cbind(age,prinos))
# Put both data matrice\tilde{}s together
dat=data.frame(rbind(dat4CTRL,dat4drug))
# Make "TRT" as a factor for treatment.
dat$TRT = as.factor(rep(c("CTRL", "BIOSTIM"), each=n))
# Make a "Center" to represent the center number
dat$Center = center
# Return the simulated data
dat
} # end of function


# sad je ovo za prvi istrazivacki centar ali generisano funkcijom ovom odozgo

d1 = data.generator(n,age.mu,age.sd,yield.mu,yield.muTRT, yield.sd, 1)

### generisem podatke za ostale istrazivacke centre

# Data from Center 2
yield.muTRT2 = 53
d2 = data.generator(n,age.mu,age.sd,yield.mu,yield.muTRT2,yield.sd,2)
# Data from Center 3
yield.muTRT3 = 58
d3 = data.generator(n,age.mu,age.sd,yield.mu, yield.muTRT3,yield.sd,3)
# Data from Center 4
yield.muTRT4 = 53
d4 = data.generator(n,age.mu,age.sd,yield.mu, yield.muTRT4,yield.sd,4)
# Data from Center 5
yield.muTRT5 = 54.5
d5 = data.generator(n,age.mu,age.sd,yield.mu, yield.muTRT5,yield.sd,5)



dat = data.frame(rbind(d1,d2,d3,d4,d5))
# Change `Center' from numeric to factor
dat$Center = as.factor(dat$Center)
head (dat)
View(dat)#alternatively just dat


# Call bwplot
print(bwplot(prinos~TRT|Center, data=dat,xlab="TRT",
strip=strip.custom(bg="white"),
ylab="Yield",lwd=3,cex=1.3,pch=20,
type=c("p", "r")))


# Call bwplot
print(bwplot(prinos~Center|TRT, data=dat,xlab="orchard",
strip=strip.custom(bg="white"),
ylab="yield",lwd=3,cex=1.3,pch=20,
type=c("p", "r")))


# Model for Center 1
m.c1 = aov(prinos~TRT, data=dat[dat$Center==1,])
# Print the summary
summary(m.c1)

# Model for Center 2
m.c2 = aov(prinos~TRT, data=dat[dat$Center==2,])
# Print the summary
summary(m.c2)



# Model for Center 3
m.c3 = aov(prinos~TRT, data=dat[dat$Center==3,])
# Print the summary
summary(m.c3)


# Model for Center 4
m.c4 = aov(prinos~TRT, data=dat[dat$Center==4,])
# Print the summary
summary(m.c4)



# Model for Center 5
m.c5 = aov(prinos~TRT, data=dat[dat$Center==5,])
# Print the summary
summary(m.c5)




 # Call 'aov' to fit the 3-way model
 lm1 = aov(prinos~TRT*Center*age, data=dat)
summary (lm1)



# Call 'aov' to fit the reduced model
lm2 = aov(prinos~TRT+Center, data=dat)
summary(lm2)


# Get the study sample size
ndat = aggregate(dat$prinos,
list(Center=dat$Center,TRT = dat$TRT), length)
# Print the study specific sample size
ndat



# Calcuate the means by study
mdat = aggregate(dat$prinos,
list(Center=dat$Center,TRT = dat$TRT), mean)
# Print the means
mdat


# Calculate the standard deviations
sddat = aggregate(dat$prinos,
list(Center=dat$Center,TRT = dat$TRT), sd)
# Print the SDs
sddat


# Call the library
## ako treba ####libPaths( "/R/bb paketi" )

library(metafor)
# Calculate the ESs
esdat = escalc(measure="MD",
n1i= ndat$x[ndat$TRT=="BIOSTIM"],
n2i= ndat$x[ndat$TRT=="CTRL"],
m1i= mdat$x[mdat$TRT=="BIOSTIM"],
m2i= mdat$x[mdat$TRT=="CTRL"],
sd1i= sddat$x[sddat$TRT=="BIOSTIM"],
sd2i= sddat$x[sddat$TRT=="CTRL"], append=T)
rownames(esdat) = ndat$Study[ndat$TRT=="TRT"]
# Print the ES dataframe
esdat



# Calculate the z-values for each study
z = esdat$yi/sqrt(esdat$vi)
# Calculate the p-values for each study
pval.studywise = 2*(1-pnorm(abs(z)))
# Print the p-values
pval.studywise



# Random-effects meta-analysis with DL
meta.MD.DL = rma(yi,vi,measure="MD",method="DL", data=esdat)
# Print the result
meta.MD.DL

forest(meta.MD.DL)
funnel (meta.MD.DL)
