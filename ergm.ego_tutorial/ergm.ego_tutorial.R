## ----setup, include=FALSE, cache=FALSE----------------------------------------
library(knitr)
knitr::opts_chunk$set(
  cache=FALSE, 
  warning=FALSE, 
  comment=NA, 
  fig.align='center'
)


## ----meta-dev, child = 'common/statnet-dev-team.md'---------------------------




## ----meta-project, child = 'common/statnet-project.md'------------------------




## ----eval=F-------------------------------------------------------------------
## install.packages('ergm.ego')


## ----install-universe, eval=FALSE---------------------------------------------
## install.packages(
##   "ergm.ego",
##   repos = c("https://statnet.r-universe.dev", "https://cloud.r-project.org")
## )


## ----load-packages, cache=FALSE-----------------------------------------------
install.packages('ergm.ego')
install.packages("DBI")

library(egor)
library(ergm.ego)
packageVersion('ergm.ego')


## ----eval=F-------------------------------------------------------------------
## example(sample.egor)


## ----eval=FALSE---------------------------------------------------------------
## library(help='ergm.ego')


## ----eval=FALSE---------------------------------------------------------------
## ?as.egor


## ----eval=FALSE---------------------------------------------------------------
## help('ergm.ego-terms')


## ----eval=FALSE---------------------------------------------------------------
## sessionInfo()


## ----cache=FALSE--------------------------------------------------------------
set.seed(1)


## ----mesa-load----------------------------------------------------------------
data(faux.mesa.high)
mesa <- faux.mesa.high


## ----mesa-plot----------------------------------------------------------------
plot(mesa, vertex.col="Grade")
legend('bottomleft',fill=7:12,legend=paste('Grade',7:12),cex=0.75)


## ----mesa-egor----------------------------------------------------------------
mesa.ego <- as.egor(mesa) 


## ----mesa-egor-look-----------------------------------------------------------
names(mesa.ego) # what are the components of this object?
mesa.ego # shows the dimensions of each component
#View(mesa.ego) # opens the component in the Rstudio source window
class(mesa.ego) # what type of "object" is this?



## ----mesa-egor-dfs------------------------------------------------------------
class(mesa.ego$ego) # and what type of objects are the components?
class(mesa.ego$alter)
class(mesa.ego$aatie)


## ----mesa-egor-egos-----------------------------------------------------------
mesa.ego$ego # first few rows of the ego table


## ----mesa-egor-alters---------------------------------------------------------
mesa.ego$alter # first few rows of the alter table

# ties show up twice, but alter info is linked to .altID
# mesa.ego$alter %>% filter((.altID==1 & .egoID==25) | (.egoID==1 & .altID==25))


## ----mesa-egor-aaties---------------------------------------------------------
mesa.ego$aatie # first few rows of the alter table


## ----mesa-egor-writecsv-------------------------------------------------------
# egos
write.csv(mesa.ego$ego, file="ergm.ego_tutorial/mesa.ego.table.csv", row.names = F)

# alters
write.csv(mesa.ego$alter[,-1], file="ergm.ego_tutorial/mesa.alter.table.csv", row.names = F)


## ----mesa-egor-readcsv--------------------------------------------------------
mesa.egos <- read.csv("ergm.ego_tutorial/mesa.ego.table.csv")
head(mesa.egos)
mesa.alts <- read.csv("ergm.ego_tutorial/mesa.alter.table.csv")
head(mesa.alts)


## ----mesa-myegodata (TO CREATE AN EGOR OBJECT)-----------------------------------------------------------
my.egodata <- egor(egos = mesa.egos, 
                   alters = mesa.alts, 
                   ID.vars = list(ego = ".egoID"))
my.egodata

my.egodata$alter
my.egodata$ego


## ----message=F, eval=F--------------------------------------------------------
## example("egor")


## ----mesa-explore-------------------------------------------------------------
# to reduce typing, we'll pull the ego and alter data frames
egos <- mesa.egos$ego
alters <- mesa.egos$alter

table(egos$Sex) # Distribution of `Sex`
table(egos$Race) # Distribution of `Race`
barplot(table(egos$Grade), 
        main = "Ego grade distribution",
        ylab="frequency")


## ----mesa-explore-race, fig.show="hold"---------------------------------------
layout(matrix(1:2, 1, 2))
barplot(table(egos$Race)/nrow(egos),
        main="Ego Race Distn", ylab="percent",
        ylim = c(0,0.5), las = 3)
barplot(table(alters$Race)/nrow(alters),
        main="Alter Race Distn", ylab="percent",
        ylim = c(0,0.5), las = 3)
layout(1)


## ----mesa-mm------------------------------------------------------------------
# to get the crosstabulated counts of ties:
mixingmatrix(mesa.ego,"Grade")

# contrast with the original network crosstab:
mixingmatrix(mesa, "Grade")



## ----mesa-mm2-----------------------------------------------------------------
# to get the row conditional probabilities:

round(mixingmatrix(mesa.ego, "Grade", rowprob=T), 2)
round(mixingmatrix(mesa.ego, "Race", rowprob=T), 2)


## ----mesa-avgdegrees----------------------------------------------------------
# first, using the original network
network.edgecount(faux.mesa.high)

# compare to `egor`
# note that the ties are double counted, so we need to divide by 2.
nrow(mesa.ego$alter)/2

# mean degree -- here we want to count each "stub", so we don't divide by 2
nrow(mesa.ego$alter)/nrow(mesa.ego$ego)


# overall degree distribution
summary(mesa.ego ~ degree(0:20))

# and stratified by sex
summary(mesa.ego ~ degree(0:13, by="Sex"))


## ----mesa-suffstats-----------------------------------------------------------
summary(mesa.ego ~ degree(0:10), scaleto=100000)
summary(mesa.ego ~ degree(0:10), scaleto=nrow(mesa.ego$ego)*100)


## ----mesa-degrees-counts------------------------------------------------------
# degreedist(mesa.ego, plot=TRUE, prob=FALSE) # bug statnet/ergm.ego#82.
degreedist(mesa.ego, by="Sex", plot=TRUE, prob=FALSE)


## ----mesa-degrees-probs-------------------------------------------------------
degreedist(mesa.ego, by="Sex", plot=TRUE, prob=TRUE)


## ----mesa-degrees-brg, cache=FALSE--------------------------------------------
set.seed(1)
degreedist(mesa.ego, brg=TRUE)


## ----mesa-degrees-brg-sex-----------------------------------------------------
degreedist(mesa.ego, by="Sex", prob=TRUE, brg=TRUE)


## ----eval=FALSE---------------------------------------------------------------
## ?control.ergm.ego


## ----fit-edges, message=FALSE-------------------------------------------------
fit.edges <- ergm.ego(mesa.ego ~ edges)
summary(fit.edges)


## ----fit-edges-popsizes-------------------------------------------------------
names(fit.edges)
fit.edges$ppopsize
fit.edges$popsize


## ----fit-edges-ppopsize1000---------------------------------------------------
summary(ergm.ego(mesa.ego ~ edges, 
                 control = control.ergm.ego(ppopsize=1000)))


## ----fit-edges-mcmcdiag-------------------------------------------------------
mcmc.diagnostics(fit.edges, which ="plots")


## ----fit-edges-gof------------------------------------------------------------
plot(gof(fit.edges, GOF="model"))


## ----fit-edges-gof-degree-----------------------------------------------------
plot(gof(fit.edges, GOF="degree"))



## ----fit-edgesdeg0, message=FALSE---------------------------------------------
set.seed(1)
fit.deg0 <- ergm.ego(mesa.ego ~ edges + degree(0), 
                     control = snctrl(ppopsize=1000))
summary(fit.deg0)


## ----fit-edgesdeg0-mcmcdiag---------------------------------------------------
mcmc.diagnostics(fit.deg0, which = "plots")


## ----fit-edgesdeg0-gofplot----------------------------------------------------
plot(gof(fit.deg0, GOF="model"))
plot(gof(fit.deg0, GOF="degree"))


## ----fit-full, message=FALSE--------------------------------------------------
fit.full <- ergm.ego(mesa.ego ~ edges + degree(0:1) 
                     + nodefactor("Sex")
                     + nodefactor("Race", levels = -LARGEST)
                     + nodefactor("Grade")
                     + nodematch("Sex") 
                     + nodematch("Race") 
                     + nodematch("Grade"))
summary(fit.full)


## ----fit-full-mcmcdiag--------------------------------------------------------
mcmc.diagnostics(fit.full, which = "plots")


## ----fit-full-gofplot---------------------------------------------------------
plot(gof(fit.full, GOF="model"))
plot(gof(fit.full, GOF="degree"))


## ----sim-full-----------------------------------------------------------------
sim.full <- simulate(fit.full)
summary(mesa.ego ~ edges + degree(0:1)
                      + nodefactor("Sex")
                      + nodefactor("Race", levels = -LARGEST)
                      + nodefactor("Grade")
                      + nodematch("Sex") + nodematch("Race") + nodematch("Grade"))

summary(sim.full ~ edges + degree(0:1)
                      + nodefactor("Sex")
                      + nodefactor("Race", levels = -LARGEST)
                      + nodefactor("Grade")
                      + nodematch("Sex") + nodematch("Race") + nodematch("Grade"))
plot(sim.full, vertex.col="Grade")
legend('bottomleft',fill=7:12,legend=paste('Grade',7:12),cex=0.75)


## ----sim-full2----------------------------------------------------------------
sim.full2 <- simulate(fit.full, popsize=network.size(mesa)*2)
summary(mesa~edges + degree(0:1)
                      + nodefactor("Sex")
                      + nodefactor("Race", levels = -LARGEST)
                      + nodefactor("Grade")
                      + nodematch("Sex") + nodematch("Race") + nodematch("Grade"))*2

summary(sim.full2~edges + degree(0:1)
                      + nodefactor("Sex")
                      + nodefactor("Race", levels = -LARGEST)
                      + nodefactor("Grade")
                      + nodematch("Sex") + nodematch("Race") + nodematch("Grade"))



## ----parrecov-data------------------------------------------------------------
data(faux.magnolia.high)
faux.magnolia.high -> fmh
N <- network.size(fmh)


## ----parrecov-fit, message=FALSE----------------------------------------------
fit.ergm <- ergm(fmh ~ degree(0:3) 
                 + nodefactor("Race", levels=TRUE) + nodematch("Race")
                 + nodefactor("Sex") + nodematch("Sex") 
                 + absdiff("Grade"))
round(coef(fit.ergm), 3)


## ----parrecov-census, message=FALSE-------------------------------------------
fmh.ego <- as.egor(fmh)
head(fmh.ego)

egofit <- ergm.ego(fmh.ego ~ degree(0:3) 
                   + nodefactor("Race", levels=TRUE) + nodematch("Race")
                   + nodefactor("Sex") + nodematch("Sex") 
                   + absdiff("Grade"), popsize=N,
				  control = snctrl(ppopsize=N))

# A convenience function.
model.se <- function(fit) sqrt(diag(vcov(fit)))

# Parameters recovered:
coef.compare <- data.frame(
  "NW est" = coef(fit.ergm), 
  "Ego Cen est" = coef(egofit)[-1],
  "diff Z" = (coef(fit.ergm)-coef(egofit)[-1])/model.se(egofit)[-1])

round(coef.compare, 3)


## ----parrecov-mcmcdiag, eval=FALSE--------------------------------------------
## # MCMC diagnostics.
## mcmc.diagnostics(egofit, which="plots")


## ----parrecov-gofplot---------------------------------------------------------
plot(gof(egofit, GOF="model"))


## ----parrecov-gofplot2--------------------------------------------------------
plot(gof(egofit, GOF="degree"))


## ----parrecov-samplelarge-----------------------------------------------------
set.seed(1)
fmh.egosampN <- sample(fmh.ego, N, replace=TRUE)

egofitN <- ergm.ego(fmh.egosampN ~ degree(0:3) 
                    + nodefactor("Race", levels=TRUE) + nodematch("Race") 
                    + nodefactor("Sex") + nodematch("Sex")
                    + absdiff("Grade"),
                    popsize=N)

# compare the coef
coef.compare <- data.frame(
  "NW est" = coef(fit.ergm), 
  "Ego SampN est" = coef(egofitN)[-1],
  "diff Z" = (coef(fit.ergm)-coef(egofitN)[-1])/model.se(egofitN)[-1])

round(coef.compare, 3)
                 
# compare the s.e.'s
se.compare <- data.frame(
  "NW SE" = model.se(fit.ergm), 
  "Ego census SE" =model.se(egofit)[-1], 
  "Ego SampN SE" = model.se(egofitN)[-1])

round(se.compare, 3)


## ----parrecov-samplesmall-----------------------------------------------------
set.seed(0) # Some samples have different sets of alter levels from ego levels.

fmh.egosampN4 <- sample(fmh.ego, round(N/4), replace=TRUE)

egofitN4 <- ergm.ego(fmh.egosampN4 ~ degree(0:3) 
                    + nodefactor("Race", levels=TRUE) + nodematch("Race") 
                    + nodefactor("Sex") + nodematch("Sex")
                    + absdiff("Grade"),
                    popsize=N)

# compare the coef
coef.compare <- data.frame(
  "NW est" = coef(fit.ergm), 
  "Ego SampN4 est" = coef(egofitN4)[-1],
  "diff Z" = (coef(fit.ergm)-coef(egofitN4)[-1])/model.se(egofitN4)[-1])

round(coef.compare, 3)

# compare the s.e.'s
se.compare <- data.frame(
  "NW SE" = model.se(fit.ergm), 
  "Ego census SE" =model.se(egofit)[-1], 
  "Ego SampN SE" = model.se(egofitN)[-1],
  "Ego Samp4 SE" = model.se(egofitN4)[-1])

round(se.compare, 3)


## ----session-info, echo = FALSE-----------------------------------------------
sessioninfo::session_info()

