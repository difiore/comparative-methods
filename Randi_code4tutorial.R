############################################################################
## Preparations
############################################################################

# Comments are ignored by R, but can remind you what the code does.
# A hashtag at the start of a line tells R the line is a comment.

# Install necessary packages
install.packages("ape")
install.packages("caper")
install.packages("picante")

# Load necessary packages from library
library(ape)
library(caper)
library(picante)

# Get information about a function, e.g., sum
?sum

# Get information about a package, e.g., ape
library(help = ape)

############################################################################
## Reading your data and phylogeny into R
############################################################################

# Read the primate data into R
primatedata = read.table("MYPATH/Primatedata.txt", sep="\t", header=T)

# Alternatively, use file.choose() to find the file without typing the path
primatedata = read.table(file.choose(), sep="\t", header=T)

# Look at the structure of your data file
str(primatedata)

# Look at the column names
names(primatedata)

# Look at the first 4 rows of data
head(primatedata, n=4)

# Look at the entire data set
primatedata

# Read phylogeny into R (can also use file.choose as with primatedata)
primatetree = read.nexus(file.choose())

# Look at basic info about the phylogeny
primatetree

# Look at the structure of the phylogeny
str(primatetree)

# Visualize the phylogeny
plot.phylo(primatetree)

# Plot with reduced tip label size
plot.phylo(primatetree, cex=0.5)	

# Zoom in on Cercopithecus species, without and with the rest of the tree
zoom(primatetree, list(grep("Cercopithecus", primatetree$tip.label)), subtree=FALSE)
zoom(primatetree, list(grep("Cercopithecus", primatetree$tip.label)), subtree=TRUE)

# Drop a species from the tree, and view it
primatetree2 = drop.tip(primatetree, "Aotus_azarae_infulatus")
str(primatetree2)

# Find more plotting options
?plot.phylo

# Plot a fancy phylogeny
plot.phylo(primatetree, type="fan", edge.color="deeppink", tip.color="green", cex=0.5)

# Plot a phylogeny with a timescale
plot.phylo(primatetree, cex=0.5)
axisPhylo()

# Check if there are polytomies in the tree
is.binary.tree(primatetree)

# Randomly resolve polytomies with 0-length branches
primatetree = multi2di(primatetree)

############################################################################
## Manipulating your data and phylogeny in R
############################################################################

# Substitute underscores for spaces in the binomial names
primatedata$Binomial = gsub(" ", "_", primatedata$Binomial)

# Make row names in the data set equal to species names
row.names(primatedata) = primatedata$Binomial

# Identify species in the phylogeny but not in the data
DropFromTree = setdiff(primatetree$tip.label, primatedata$Binomial)

# Identify species in the data but not in the phylogeny
DropFromData = setdiff(primatedata$Binomial, primatetree$tip.label)

# Check number of species to drop from the tree and data respectively
length(DropFromTree)
length(DropFromData)

# Prune tree to only species found in the data
primatetree = drop.tip(primatetree, setdiff(primatetree$tip.label, primatedata$Binomial))

# Prune data to only species in the tree
primatedata = primatedata[-which(rownames(primatedata) == DropFromData),]

# Check that we have the same number of species in the tree and data
length(primatetree$tip.label)
nrow(primatedata)

############################################################################
## Orinary least squares
############################################################################

# Take a look at the data, first untransformed, then log transformed
par(mfrow = c(2,2))	# plot in a 2x2 window
hist(primatedata$AdultBodyMass_g, col=rainbow(8), main="Adult Body Mass")
hist(primatedata$GestationLen_d, col=rainbow(8), main="Gestation Length")
hist(log(primatedata$AdultBodyMass_g), col=rainbow(8), main="log(Adult Body Mass)")
hist(log(primatedata$GestationLen_d), col=rainbow(8), main="log(Gestation Length)")
par(mfrow = c(1,1))	# reset graphical paramaters to a 1x1 plot window

# Run linear regression of log gestation length against log adult body mass
model.ols = lm(log(GestationLen_d) ~ log(AdultBodyMass_g), data=primatedata)
summary(model.ols)

# Plot the regression line on a scatterplot
plot(log(GestationLen_d) ~ log(AdultBodyMass_g), data=primatedata)
abline(model.ols)

# Visualize phylogenetic pseudoreplocation
points(log(GestationLen_d[Family=="Cercopithecidae"]) ~ log(AdultBodyMass_g[Family=="Cercopithecidae"]), data=primatedata, col="blue", pch=16)

############################################################################
## Phylogenetic generalized least squares
############################################################################

# Create comparative data object
primate = comparative.data(phy=primatetree, data=primatedata, names.col=Binomial, vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE)

# View dropped species to check they make sense
primate$dropped$tips

# Run PGLS model
model.pgls = pgls(log(GestationLen_d) ~ log(AdultBodyMass_g), data=primate, lambda="ML")
summary(model.pgls)

# Plot the results
plot(log(GestationLen_d) ~ log(AdultBodyMass_g), data=primate$data)
abline(model.pgls)

# Set lower bound of lambda to something slightly larger than 0
# In early versions of caper, you must set bounds for all 3 scaling
# parameters (lambda, kappa, delta) if you use the 'bounds' argument
model.pgls2 = pgls(log(GestationLen_d) ~ log(AdultBodyMass_g), data=primate, lambda="ML", bounds=list(lambda=c(1e-05, 1), kappa=c(0,3), delta=c(0,5)))

############################################################################
## Model diagnostics for PGLS
############################################################################

# View diagnostic plots for the PGLS model
par(mfrow=c(2,2))
plot(model.pgls)

############################################################################
## Estimating lambda for a single variable
############################################################################

# Estimate lambda for a single variable
est.lambda = pgls(log(GestationLen_d)~1, data=primate, lambda="ML")
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
plot(SocialGroupSize ~ HomeRange_km2, data=primatedata, xlab="Home Range Size (2 km)", ylab="Social Group Size", main="My Scatterplot")
points(SocialGroupSize[Family=="Cebidae"] ~ HomeRange_km2[Family=="Cebidae"], data=primatedata, col="green", pch=16)

# 2. Plot the primate tree and zoom in on the Saimiri genus
zoom(primatetree, list(grep("Saimiri", primatetree$tip.label)), subtree=TRUE)

# 3. Run a PGLS for social group size as a function of gestation length. How does this 
#    differ from the result obtained with OLS? Plot the results of OLS and PGLS side by side.
par(mfrow = c(2,2))	# plot in a 2x2 window
hist(primatedata$SocialGroupSize, col=rainbow(8), main="Adult Body Mass")
hist(primatedata$GestationLen_d, col=rainbow(8), main="Gestation Length")
hist(log(primatedata$SocialGroupSize), col=rainbow(8), main="log(Adult Body Mass)")
hist(log(primatedata$GestationLen_d), col=rainbow(8), main="log(Gestation Length)")

model.ols2 = lm(log(SocialGroupSize) ~ log(GestationLen_d), data=primatedata)
model.pgls2 = pgls(log(SocialGroupSize) ~ log(GestationLen_d), data=primate, lambda="ML")

par(mfrow=c(1,2))
plot(log(SocialGroupSize) ~ log(GestationLen_d), data=primatedata)
abline(model.ols2)
plot(log(SocialGroupSize) ~ log(GestationLen_d), data=primate$data)
abline(model.pgls2)
par(mfrow=c(1,1))

# 4. Find the estimate and confidence interval of lambda for primate social group size, and
#    plot the likelihood profile
est.lambda2 = pgls(log(SocialGroupSize)~1, data=primate, lambda="ML")
plot(pgls.profile(est.lambda2, "lambda"))

# 5. Create a new data set and phylogeny containing only Cercopithecidae species. Find the
#    estimate and confidence interval of lambda for social group size, and plot the likelihood 
#    profile. Compare the results for all primates.
cerconames_index = grep("Cerco", primatedata$Family)
cerconames_not = primatedata$Binomial[grep("Cerco", primatedata$Family, invert=TRUE)]
cercotree = drop.tip(phy=primatetree, tip=cerconames_not)
cercodata = primatedata[cerconames_index,]
cerco = comparative.data(phy=cercotree, data=cercodata, names.col=Binomial, vcv=TRUE, na.omit=FALSE)
est.lambda3 = pgls(log(SocialGroupSize)~1, data=cerco, lambda="ML")
par(mfrow=c(1,2))
plot(pgls.profile(est.lambda2, "lambda"), main="All Primates")
plot(pgls.profile(est.lambda3, "lambda"), main="Cercopithecidae")
par(mfrow=c(1,1))





























