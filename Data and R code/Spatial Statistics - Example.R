# Spatial Statistics- Example
# sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# RStudio Version 1.4.1717

rm(list=ls())

# load the libraries with all dependencies
if(!require('rgdal')) install.packages('rgdal', dependencies=T); library('rgdal')
if(!require('tmap')) install.packages('tmap', dependencies=T); library('tmap')
if(!require('rgeos')) install.packages('rgeos', dependencies=T); library('rgeos')
if(!require('spdep')) install.packages('spdep', dependencies=T); library('spdep')
if(!require('RColorBrewer')) install.packages('RColorBrewer', dependencies=T); library('RColorBrewer')
if(!require('fields')) install.packages('fields', dependencies=T); library('fields')
if(!require('car')) install.packages('car', dependencies=T); library('car')
# Recently, some functions from spdep package have been moved to the 
# spatialreg package. SO we load it in order to avoid warnings
if(!require('spatialreg')) install.packages('spatialreg', dependencies=T); library('spatialreg')
if(!require('spgwr')) install.packages('spgwr', dependencies=T); library('spgwr')


# Shapefile loading---------------
setwd("../Topic2 - Spatial Statistics/Data and R code")
#setwd("C:/Users/franc/Dropbox/EMOS_LO_fra/Prov01012007")

italy <- readOGR("Prov01012007.shp",stringsAsFactors=F)

# Loading data---------------
silcIT3 <- read.csv("silcIT3.csv", sep=";")
str(silcIT3)
# The dataset contain the following variables obtained from a sample of
# Eusilc 2007 wave in Italy: 
# Prov.Cod: province code (in Italy provinces are NUTS-3 areas)
# income: province mean of household equivalised income
# LogIncome: province mean of log household equivalised income (target variable)
# gender: proportion of female household head (hh)
# employment.rate: employment rate
# isced5: rate of graduated people (education level >= isced5)
# owner: rate of households with an owned house

# Merging data
Silc.Spat <- merge(italy, silcIT3, by.x="COD_PROV", by.y="Prov.Cod")

# Mapping data (choropleth map)
# alternative method using sf
# italy <- read_sf("Prov01012007_WGS84.shp")
# italy <- italy[order(as.numeric(italy$COD_PROV)),]
# Silc.Spat  <- st_as_sf(Silc.Spat)
tm_shape(Silc.Spat) + 
  tm_fill("LogIncome",
          palette = "Reds", 
          style = "quantile", 
          title = "") +
  tm_borders(alpha=.7)  +
  tm_layout("",inner.margins=c(0,0,0,0),
            legend.width=0.8,
            legend.text.size = .45,
            legend.position = c("left","bottom"),
            legend.bg.color = "white",
            legend.stack = "vertical")


# Weight matrix---------------
# poly2nb() builds a neighbours list based on regions with contiguous 
# boundaries, that is sharing one or more boundary point.
neigh<-poly2nb(italy,row.names =italy$COD_PROV,queen=TRUE)
# Plot
plot(italy, border = 'lightgrey')
plot(neigh, coordinates(italy), add=TRUE, col='red')
# the function coordinates() returns coordinates in UTM format

# Rook's case neighbors-list
neigh.ROOK<- poly2nb(italy, queen=FALSE)

# nb2mat() generates a weights matrix. Setting style="W" we obtain
# a row standardised proximity matrix. 
# The option zero.policy=TRUE makes it possible to generate
# a list of weights which takes value 'zero' for observations without neighbours
italy.sp_w<-nb2mat(neigh, glist=NULL, style="W", zero.policy=TRUE)
# To transform the list into an actual matrix W, we can use the function nb2listw()
italy.lw <- nb2listw(neigh, style="W", zero.policy=TRUE)


# We can also set style="B" for the basic binary coding
italy.sp_w2<-nb2mat(neigh, glist=NULL, style="B", zero.policy=TRUE)
# we can compute the number of neighbors for each area
n.neigh <- rowSums(italy.sp_w2)
n.neigh

# Matrix based on distance
# Construct spatial weight 
# k-nearest neighbors criteria.
# k-nearest neighbors for k = 2 using knearneigh(). It will give a class of knn, 
# which is similar to class nb and we need to convert using knn2nb()
italy.nb <- knn2nb(knearneigh(coordinates(italy),k=2))

# Inverse distance weight matrix
distance <- rdist.earth(coordinates(italy),coordinates(italy)) 
# computing the Euclidean distance
diag(distance) <- 0
distance.inv <- ifelse(distance!=0, 1/distance, distance) # 1/d_{ij}
# Standardized inverse weight matrix
distance.inv <- mat2listw(distance.inv, style = "W")


# Moran's I plot---------------
moran.plot(Silc.Spat$LogIncome, italy.lw,
           labels=FALSE,
           xlab="Observed distribution of log mean income by NUTS-3", 
           ylab="Spatially lagged log mean income")
# a linear relationship appears between the log(Income) of one province 
# and that of its neighbourhood

# Moran's I index---------------
Moran.I<-moran.test(x=Silc.Spat$LogIncome,listw=italy.lw,
                    zero.policy=TRUE)
Moran.I
# the positive spatial autocorrellation is significant

# Geary's C---------------
geary.test(Silc.Spat$LogIncome, listw=italy.lw, zero.policy=T)

# Local Moran's I---------------
LocalMoran.I <- localmoran(Silc.Spat$LogIncome, italy.lw, 
                           zero.policy=TRUE)
# NB: standard p-values. Using the argument p.adjust.method is possible
# to adjust p-values
# Add results of local Moran I to the shapefile
moran.map <- cbind(Silc.Spat, LocalMoran.I)

# Values of local Moran's I on Italian provinces
tm_shape(moran.map) +
  tm_fill(col = "Ii",
          style = "quantile",
          palette = brewer.pal("Spectral", n = 5),
          title = "") +
  tm_borders(alpha=.7)  +
  tm_layout("",inner.margins=c(0,0,0,0),
            legend.width=0.8,
            legend.text.size = .45,
            legend.position = c("left","bottom"),
            legend.bg.color = "white",
            legend.stack = "vertical")

# Mapping spatial clusters LISA
# centers the variable of interest around its mean 
Silc.Spat$m.LogIncome <- Silc.Spat$LogIncome - mean(Silc.Spat$LogIncome)     
# create a spatially lagged variable
Silc.Spat$lag_m.LogIncome <- lag.listw(italy.lw, Silc.Spat$m.LogIncome)
# we need to identify if each obs belongs to the high-high, low-low,
# high-low, or low-high quadrants
signif <- 0.05 # significance threshold
quadrant <- vector(mode="numeric",length=nrow(LocalMoran.I))
quadrant[Silc.Spat$m.LogIncome >0 & Silc.Spat$lag_m.LogIncome>0] <- 4
quadrant[Silc.Spat$m.LogIncome <0 & Silc.Spat$lag_m.LogIncome<0] <- 1
quadrant[Silc.Spat$m.LogIncome <0 & Silc.Spat$lag_m.LogIncome>0] <- 2
quadrant[Silc.Spat$m.LogIncome >0 & Silc.Spat$lag_m.LogIncome<0] <- 3
quadrant[LocalMoran.I[,5]>signif] <- 0 # non-significant observations

Silc.Spat$quadrant <- car:::recode(quadrant, "0='non-significant';
                                   1='low-low'; 2='low-high';3='high-low';
                                   4='high-high'")
table(Silc.Spat$quadrant)
# We can see that 37 areas have an adjusted p-values less than 0.05 threshold.
# mapping LISA
colors <- c("red",rgb(0,0.5,1,alpha=0.8),rgb(0.8,1,0.8,alpha=0.3))
tm_shape(Silc.Spat) +
  tm_fill(col = "quadrant",
          palette = colors,
          title = "LISA")+
  tm_borders(alpha=.7)  +
  tm_layout("",inner.margins=c(0,0,0,0),
            legend.width=0.8,
            legend.text.size = .45,
            legend.position = c("left","bottom"),
            legend.bg.color = "white",
            legend.stack = "vertical")
# There is statistically significant moderate clustering in log(Income) in Italy

# adjusted p-value
LocalMoran.I.Bonf <- localmoran(silcIT3$LogIncome, italy.lw, 
                                zero.policy=TRUE, 
                                p.adjust.method = "bonferroni")
table(LocalMoran.I.Bonf[,5]<=signif)

# Global autocorrelation test on a categorical variable---------------
# define a new variable that is equal to 1 if the income of a provice is higher
# than the 0.8*median(income), 0 otherwise
inc.threshold <- 0.8*median(silcIT3$income)
pov.ind <- as.factor(ifelse(silcIT3$income>inc.threshold,1,0))
table(pov.ind)
joincount.test(pov.ind, listw=italy.lw, zero.policy=T)

# Spatial econometrics models---------------
# Aim: modelling the log household equivalised income by NUTS-3 in Italy, 
# using as covariates the variables available in the dataset.
# Using the Moran test we have found that the null hypothesis 
# assuming no spatial autocorrelation shoul be rejected

# Defining the model
model <- LogIncome ~gender+ employment.rate+isced5+owner

# OLS model
mod.ols <- lm(model, data=Silc.Spat)
summary(mod.ols)

# extract residuals to evaluate spatial autocorrelation
res.ols <- residuals(mod.ols)
Silc.Spat$std.res.ols <- (res.ols- mean(res.ols)) / sd(res.ols)
summary(Silc.Spat$std.res.ols)
my_breaks <- c(-4,-3,-2,-1,1,2,3,4)

tm_shape(Silc.Spat) +
  tm_fill("std.res.ols", title = "", style = "fixed", 
          breaks = my_breaks, palette = "-RdBu")+
  tm_borders(alpha=.7)  +
  tm_layout("",inner.margins=c(0,0,0,0),
            legend.width=0.8,
            legend.text.size = .45,
            legend.position = c("right","top"),
            legend.bg.color = "white",
            legend.stack = "vertical")
# It seems that there is a spatial patterning of areas of over- and under-
# prediction. So spatial autocorrelation may be present.

# Moran test to residuals
lm.morantest(mod.ols,italy.lw, zero.policy = T, 
             alternative="two.sided")
# statistically significant value for Moran's I
# NB: lm.morantest() takes into account the fact that the variable under 
# consideration is a residual. The usual Moran's I test statistic does not. 

# In order to decide whether to fit a spatial error or a spatially lagged 
# model we run the Lagrange Multiplier tests (bottom-up approach).
# The LM test statistics do allow a distinction between
# spatial error models and spatial lag models. 

lm.LMtests(mod.ols,italy.lw,test="LMerr")
lm.LMtests(mod.ols,italy.lw,test="LMlag")
# the p-value of both tests us below the .05 level and this means that
# we need to run the robust version of these tests

lm.LMtests(mod.ols,italy.lw,test="RLMerr")
lm.LMtests(mod.ols,italy.lw,test="RLMlag")
# for robust version of the tests, only the one for the lag model 
# is significant: the Lagrange Multiplier tests suggest estimating a
# spatial lag model rather than a spatial error model

# SAR Model
mod.sar<-lagsarlm(model, data=Silc.Spat, italy.lw)
summary(mod.sar)
# the spatial autoregressive parameter (rho), which measures the intensity of
# the spatial interdependency, is highly significant and indicates positive
# spatial dependence. The AIC highlights a gain with respect to the OLS.
# However, the LM test for residual autocorrelation is signficant

# Given the presence of an autoregressive term we fi it SDM and then
# apply a LR test on the common factor hypothesis to choose 
# between SEM and SDM

# SDM Model
mod.sdm<-lagsarlm(model, data=Silc.Spat, italy.lw, type="mixed")
# "Durbin" may be used instead of "mixed
summary(mod.sdm)
mod.sem<-errorsarlm(model, data=Silc.Spat, italy.lw)
summary(mod.sem)

# Common factor hypothesis test
# mod.sdm: Constraint-free model
# mod.sem: Constrained model
FC.test<-LR.sarlm(mod.sdm,mod.sem)
print(FC.test)

# the common factor hypothesis in the SDM model is 
# rejected (p-value = 3.126e-05)

# Interpretation in SDM
summary(mod.sdm)
# looking at the Akaike Information Criterion (AIC) we see that 
# the SDM has a lower AIC of than the standard linear model 
# spatial residual autocorrelation test is rejected at 0.05 (p-value=0.08).
# Rho coefficient is positive and statistically significant. 
# In other words, when the log(Income) in surrounding areas 
# increases, so does the log(Income) in each area, even when 
# we adjust for the other explanatory variables in our model. 
# A direct interpretation of regression parameters is 
# not possible because the effects must take into account the 
# effects of endogenous interaction. It is necessary
# to compute the direct and indirect effects:
impactssdm <- impacts(mod.sdm,listw=italy.lw,R=500)
summary(impactssdm, short=TRUE, zstats=TRUE)
# Empirical confidence intervals are found using 500 simulations 
# from empirical distribution

# Other models
# SLX model
mod.slx<-lmSLX(model, data=Silc.Spat, italy.lw)
summary(mod.slx)

# SAC Model
mod.sac<-sacsarlm(model, data=Silc.Spat, italy.lw)
summary(mod.sac)

# Geographically Weighted Regression (GWR)-------------
# It can be used to identify how locally weighted
# regression coefficients may vary across space
# Get the optimal bandwidth.
# adapt=TRUE: find the proportion between 0 and 1 of observations to 
# include in weighting scheme (k-nearest neighbours) 
# adapt=FALSE: find global bandwidth
# gweight: geographical weighting function, gwr.Gauss() default, 
# or gwr.bisquare()
GWRbandwidth <- gwr.sel(model, data=Silc.Spat, adapt=T) 
mod.gwr <- gwr(model,
               data = Silc.Spat,
               adapt=GWRbandwidth,
               hatmatrix=TRUE,
               se.fit=TRUE) 
#print the results of the model
mod.gwr
# The output from the GWR model reveals how the coefficients vary across 
# the italian provinces. You will see how the global coefficients are 
# the same as the coefficients in the standard linear model. 
# For instance, We can see that the coefficients of gender range 
# from a minimum value of -0.31 (1 unit change in proportion of female 
# head household resulting in a drop in average log(Income) of -0.31) 
# to +1.27 (1 unit change in proportion of female hh resulting in an 
# increase in average log(Income) of +1.27). 
# Map results 
results.gwr <-as.data.frame(mod.gwr$SDF)
names(results.gwr)
gwr.map <- cbind(Silc.Spat, as.matrix(results.gwr))
# spatial distribution of gender
map1 <- tm_shape(gwr.map) + 
  tm_fill("gender",
          n = 5,
          style = "quantile")  +
  tm_layout(frame = FALSE,
            legend.text.size = 0.5,
            legend.title.size = 0.6)
# Coefficients of gender
map2 <- tm_shape(gwr.map) + 
  tm_fill("gender.1",
          n = 5,
          style = "quantile")  +
  tm_layout(frame = FALSE,
            legend.text.size = 0.5,
            legend.title.size = 0.6)
# creates a clear grid
grid.newpage()

# assigns the cell size of the grid, in this case 2 by 2
pushViewport(viewport(layout=grid.layout(1,2)))

# prints a map object into a defined cell   
print(map1, vp=viewport(layout.pos.col = 1, layout.pos.row =1))
print(map2, vp=viewport(layout.pos.col = 2, layout.pos.row =1))

# statistical significance coefficients (0.05%)
sigTest = abs(mod.gwr$SDF$gender)-2*mod.gwr$SDF$gender_se 
sigTest
# If this is greater than zero (i.e. the estimate is more than two 
# standard errors away from zero), i.e. it is statistically significantly 
# at (nearly) the 95% confidence level.
which(sigTest>0)

