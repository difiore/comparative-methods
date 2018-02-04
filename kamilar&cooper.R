library(curl)
library(ape)
library(ggplot2)
library(caper)
library(picante)

f <- curl("https://raw.githubusercontent.com/difiore/ADA2016/master/KamilarAndCooperData.csv")
primatedata <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
head(primatedata)

p <- ggplot(data = primatedata, aes(x = Family, y = log(Body_mass_female_mean)))
p <- p + geom_boxplot()
p <- p + theme(axis.text.x = element_text(angle = 90))
p <- p + ylab("log(Female Body Mass)")
p
par(mfrow = c(1, 2))
plot(x = Body_mass_female_mean, y = Brain_Size_Female_Mean)
plot(x = log(Body_mass_female_mean), y = log(Brain_Size_Female_Mean))

f <- curl("https://raw.githubusercontent.com/difiore/ADA2016/master/Primatetree.nex")
primatetree <- read.nexus(f)

str(primatetree)

# Visualize the phylogeny
plot.phylo(primatetree)

# Plot a fancy phylogeny
plot.phylo(primatetree, type="fan", edge.color="deeppink", tip.color="green", cex=0.5)


# Check if there are polytomies in the tree
is.binary.tree(primatetree)

# Randomly resolve polytomies with 0-length branches
primatetree = multi2di(primatetree)


# Make row names in the data set equal to species names
row.names(primatedata) = primatedata$Scientific_Name

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
