library(curl)
library(ape)
library(ggplot2)
library(caper)
library(picante)

f <- curl("https://raw.githubusercontent.com/difiore/ADA2016/master/KamilarAndCooperData.csv")
primatedata <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
head(primatedata)

# Make row names in the data set equal to species names
row.names(primatedata) = primatedata$Scientific_Name

p <- ggplot(data = primatedata, aes(x = Family, y = log(Body_mass_female_mean)))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + ylab("log(Female Body Mass)")
p
par(mfrow = c(1, 2))
plot(x = primatedata$Body_mass_female_mean, y = primatedata$Brain_Size_Female_Mean)
plot(x = log(primatedata$Body_mass_female_mean), y = log(primatedata$Brain_Size_Female_Mean))

# t <- curl("https://raw.githubusercontent.com/difiore/ADA2016/master/Primatetree.nex")
# primatetree <- read.nexus(t)
t <- "primate_phylog_no_polytomies.txt"
primatetree <- read.tree(t)
str(primatetree)

# Visualize the phylogeny
plot.phylo(primatetree)

# Plot a fancy phylogeny
plot.phylo(primatetree, type="fan", edge.color="deeppink", tip.color="green", cex=0.5)

# Check if there are polytomies in the tree
is.binary.tree(primatetree)

# Randomly resolve polytomies with 0-length branches
primatetree = multi2di(primatetree)

# Identify species in the phylogeny but not in the data
DropFromTree = setdiff(primatetree$tip.label, primatedata$Scientific_Name)

# Identify species in the data but not in the phylogeny
DropFromData = setdiff(primatedata$Scientific_Name, primatetree$tip.label)

# Check number of species to drop from the tree and data respectively
length(DropFromTree)
length(DropFromData)

# Prune tree to only species found in the data
primatetree = drop.tip(primatetree, setdiff(primatetree$tip.label, primatedata$Scientific_Name))

# Prune data to only species in the tree
primatedata = primatedata[!(primatedata$Scientific_Name %in% DropFromData),]

# Check that we have the same number of species in the tree and data
length(primatetree$tip.label)
nrow(primatedata)


############################################################################
## Orinary least squares
############################################################################

# Take a look at the data, first untransformed, then log transformed
par(mfrow = c(1,1))	# plot in a 2x2 window
hist(primatedata$Body_mass_female_mean, col=rainbow(8), main="Adult Female Body Mass")
hist(primatedata$Gestation, col=rainbow(8), main="Gestation Length")
hist(log(primatedata$Body_mass_female_mean), col=rainbow(8), main="log(Adult Female Body Mass)")
hist(log(primatedata$Gestation), col=rainbow(8), main="log(Gestation Length)")
par(mfrow = c(1,1))	# reset graphical paramaters to a 1x1 plot window

# Run linear regression of log gestation length against log adult body mass
model.ols = lm(log(Gestation) ~ log(Body_mass_female_mean), data=primatedata)
summary(model.ols)

# Plot the regression line on a scatterplot
plot(log(Gestation) ~ log(Body_mass_female_mean), data=primatedata)
abline(model.ols)

# Visualize phylogenetic pseudoreplocation
points(log(Gestation[Family=="Cercopithecidae"]) ~ log(Body_mass_female_mean[Family=="Cercopithecidae"]), data=primatedata, col="blue", pch=16)

############################################################################
## Phylogenetic generalized least squares
############################################################################

# Create comparative data object
primate = comparative.data(phy=primatetree, data=primatedata, names.col=Scientific_Name, vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE)

# View dropped species to check they make sense
primate$dropped$tips

# Run PGLS model
model.pgls = pgls(log(Gestation) ~ log(Body_mass_female_mean), data=primate, lambda="ML")
summary(model.pgls)

# Plot the results
plot(log(Gestation) ~ log(Body_mass_female_mean), data=primate$data)
abline(model.pgls)

# Set lower bound of lambda to something slightly larger than 0
# In early versions of caper, you must set bounds for all 3 scaling
# parameters (lambda, kappa, delta) if you use the 'bounds' argument
model.pgls2 = pgls(log(Gestation) ~ log(Body_mass_female_mean), data=primate, lambda="ML", bounds=list(lambda=c(1e-05, 1), kappa=c(0,3), delta=c(0,5)))

############################################################################
## Model diagnostics for PGLS
############################################################################

# View diagnostic plots for the PGLS model
par(mfrow=c(1,1))
plot(model.pgls)

############################################################################
## Estimating lambda for a single variable
############################################################################

# Estimate lambda for a single variable
est.lambda = pgls(log(Gestation)~1, data=primate, lambda="ML")
summary(est.lambda)

############################################################################
## Likelihood profiles for lambda
############################################################################

# View likelihood profile for PGLS model
plot(pgls.profile(model.pgls, "lambda"))

# View likelihood profile for lambda of GestationLen_d
plot(pgls.profile(est.lambda, "lambda"))

# Check the confidence interval for lambda 
pgls.confint(model.pgls, "lambda")$ci.val
pgls.confint(est.lambda, "lambda")$ci.val

############################################################################
## Practice exercises
############################################################################

# 1. Plot social group size against home range size, coloring Cebidae in green
plot(MeanGroupSize ~ HomeRange_km2, data=primatedata, xlab="Home Range Size (2 km)", ylab="Mean Social Group Size", main="My Scatterplot")
points(MeanGroupSize[Family=="Cebidae"] ~ HomeRange_km2[Family=="Cebidae"], data=primatedata, col="green", pch=16)

# 2. Plot the primate tree and zoom in on the Saimiri genus
zoom(primatetree, list(grep("Saimiri", primatetree$tip.label)), subtree=TRUE)

# 3. Run a PGLS for social group size as a function of gestation length. How does this 
#    differ from the result obtained with OLS? Plot the results of OLS and PGLS side by side.
par(mfrow = c(1,1))	# plot in a 2x2 window
hist(primatedata$MeanGroupSize, col=rainbow(8), main="Mean Group Size")
hist(primatedata$Gestation, col=rainbow(8), main="Gestation Length")
hist(log(primatedata$MeanGroupSize), col=rainbow(8), main="log(Mean Group Size)")
hist(log(primatedata$Gestation), col=rainbow(8), main="log(Gestation Length)")

model.ols2 = lm(log(MeanGroupSize) ~ log(Gestation), data=primatedata)
model.pgls2 = pgls(log(MeanGroupSize) ~ log(Gestation), data=primate, lambda="ML")

par(mfrow=c(1,2))
plot(log(MeanGroupSize) ~ log(Gestation), data=primatedata)
abline(model.ols2)
plot(log(MeanGroupSize) ~ log(GestationLen), data=primate$data)
abline(model.pgls2)
par(mfrow=c(1,1))

# 4. Find the estimate and confidence interval of lambda for primate social group size, and
#    plot the likelihood profile
est.lambda2 = pgls(log(MeanGroupSize)~1, data=primate, lambda="ML")
plot(pgls.profile(est.lambda2, "lambda"))

# 5. Create a new data set and phylogeny containing only Cercopithecidae species. Find the
#    estimate and confidence interval of lambda for social group size, and plot the likelihood 
#    profile. Compare the results for all primates.
cerconames_index = grep("Cerco", primatedata$Family)
cerconames_not = primatedata$Scientific_Name[grep("Cerco", primatedata$Family, invert=TRUE)]
cercotree = drop.tip(phy=primatetree, tip=cerconames_not)
cercodata = primatedata[cerconames_index,]
cerco = comparative.data(phy=cercotree, data=cercodata, names.col=Scientific_Name, vcv=TRUE, na.omit=FALSE)
est.lambda3 = pgls(log(MeanGroupSize)~1, data=cerco, lambda="ML")
par(mfrow=c(1,2))
plot(pgls.profile(est.lambda2, "lambda"), main="All Primates")
plot(pgls.profile(est.lambda3, "lambda"), main="Cercopithecidae")
par(mfrow=c(1,1))

est.lambda = pgls(log(Body_mass_female_mean)~1, data=primate, lambda="ML")
summary(est.lambda)

library(phytools)
library(phylobase)
library(phylosignal)

# reorder rows of primate data to be in tree tip label order
primatedata <- primatedata[match(primatetree$tip.label,rownames(primatedata)),]

colnames(primatedata)
vars <- c("Scientific_Name", "Family", "Genus", "Species", "Brain_Size_Species_Mean", "Body_mass_male_mean", "Body_mass_female_mean", "Mass_Dimorphism", "MeanGroupSize", "AdultMales", "AdultFemale", "AdultSexRatio", "InterbirthInterval_d", "Gestation", "WeaningAge_d", "MaxLongevity_m", "LitterSz", "GR_MidRangeLat_dd", "Precip_Mean_mm", "Temp_Mean_degC", "AET_Mean_mm", "PET_Mean_mm", "HomeRange_km2", "DayLength_km", "Territoriality", "Fruit", "Leaves", "Fauna", "Canine_Dimorphism", "Feed", "Move", "Rest", "Social")

# let's limit dataset to variables in K&C 2013
primatedata <- dplyr::select(primatedata, vars)

# Calculate K for different variables using phylosignal from picante package
# We first need to set up a vector of the variable and a tree that has only tips for that variable, with the vector being in tree tip label order
# 

DropTaxa <- primatedata[is.na(primatedata$Body_mass_female_mean),]
nrow(DropTaxa)
NewPrimateData <- primatedata[!(primatedata$Scientific_Name %in% DropTaxa$Scientific_Name),]
NewPrimateTree = drop.tip(primatetree, DropTaxa$Scientific_Name)

row.names(NewPrimateData)
length(NewPrimateTree$tip.label)

est.K <- phylosignal(NewPrimateData$Body_mass_female_mean, NewPrimateTree, reps = 999, checkdata=TRUE)
est.K

####
DropTaxa <- primatedata[is.na(primatedata$Brain_Size_Species_Mean),]
nrow(DropTaxa)
NewPrimateData <- primatedata[!(primatedata$Scientific_Name %in% DropTaxa$Scientific_Name),]
NewPrimateTree = drop.tip(primatetree, DropTaxa$Scientific_Name)

row.names(NewPrimateData)
length(NewPrimateTree$tip.label)

est.K <- phylosignal(NewPrimateData$Brain_Size_Species_Mean, NewPrimateTree, reps = 999, checkdata=TRUE)
est.K

####
DropTaxa <- primatedata[is.na(primatedata$Territoriality),]
nrow(DropTaxa)
NewPrimateData <- primatedata[!(primatedata$Scientific_Name %in% DropTaxa$Scientific_Name),]
NewPrimateTree = drop.tip(primatetree, DropTaxa$Scientific_Name)

row.names(NewPrimateData)
length(NewPrimateTree$tip.label)

est.K <- phylosignal(NewPrimateData$Territoriality, NewPrimateTree, reps = 999, checkdata=TRUE)
est.K

####
DropTaxa <- primatedata[is.na(primatedata$Precip_Mean_mm),]
nrow(DropTaxa)
NewPrimateData <- primatedata[!(primatedata$Scientific_Name %in% DropTaxa$Scientific_Name),]
NewPrimateTree = drop.tip(primatetree, DropTaxa$Scientific_Name)

row.names(NewPrimateData)
length(NewPrimateTree$tip.label)

est.K <- phylosignal(NewPrimateData$Precip_Mean_mm, NewPrimateTree, reps = 999, checkdata=TRUE)
est.K
