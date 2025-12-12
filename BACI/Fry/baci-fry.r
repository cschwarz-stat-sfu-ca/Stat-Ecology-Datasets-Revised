# BACI-fry

# This example is based (loosely) on a consulting project from an 
# Independent Power Producer who was interested in monitoring the 
# effects of an in-stream hydroelectric project. 

# The response variable for the project was the minnow density at 
# different locations in the stream.

# The  monitoring design has the river divided into six segments of 
# which three are upstream of the diversion
# and three are downstream of the diversion. In each segment, several sites 
# have been located where minnow fry 
# congregate. In each of the sites, minnow traps are set for various lengths 
# of time. 

# At the end of the soaking period, the traps are removed and the number of 
# minnows are counted and classified by species.
# The counts are standardized to a common soaking time to adjust for the 
# different soak-time.  
# [This could be done directly in the analysis by using the soak-time as a 
# covariate, but the soak-time data was not available.]

# An initial pre-project monitoring study was run in 2000 to 2002. The 
# project became operational in late 2002, and
# post-project monitoring continued in 2003 and 2004.


# A nice slide presentation on lmer() is found at 
#  http://www.stat.wisc.edu/~bates/PotsdamGLMM/LMMD.pdf


library(ggplot2)
library(lmerTest)
library(emmeans)
library(plyr)

source("../../schwarz.functions.r")
#source('http://www.stat.sfu.ca/~cschwarz/Stat-Ecology-Datasets/schwarz.functions.r')
# Get the BACI power program
source("../baci-power.r")


# Read in the data
##***part001b;
# sink('baci-fry-R-001.txt', split=TRUE)
fry <- read.csv('baci-fry.csv', header=TRUE, as.is=TRUE, strip.white=TRUE)

fry$logfry       <- log(fry$Fry)

fry$Period       <- factor(fry$Period)
fry$Site         <- factor(fry$Site)
fry$SiteClass    <- factor(fry$SiteClass)
fry$Period       <- factor(fry$Period)
fry$YearF        <- factor(fry$Year)

fry[1:12,]
##***part001e;
# sink()

str(fry)


# look at the std dev to mean plot for the original and log() data 
# Compute the averages etc
# sink('baci-fry-R-020.txt', split=TRUE)
##***part020b;
mean.fry <- plyr::ddply(fry, c("Year","YearF","Site","SiteClass","Period"), function(x){
  n <- nrow(x)
  nmiss <- sum(is.na(x$Fry))
  mean.fry <- mean(x$Fry, na.rm=TRUE)
  sd.fry   <- sd(  x$Fry, na.rm=TRUE)
  mean.logfry <- mean(x$logfry, na.rm=TRUE)
  sd.logfry   <- sd  (x$logfry, na.rm=TRUE)
  res<- data.frame(n, nmiss, mean.fry, sd.fry, mean.logfry, sd.logfry)
  })
mean.fry[1:5,]
##***part020e;
# sink()


# Plot the std to mean ratios
##***part021b;
sdmeanplot <- ggplot(data=mean.fry, aes(x=mean.fry, y=sd.fry))+
  ggtitle("SD vs Mean - original data")+
  geom_point()
sdmeanplot
##***part021e;
# ggsave(plot=sdmeanplot, file='baci-fry-R-300-diagnostic.png', h=4, w=6, units="in", dpi=300)


# Plot the std to mean ratios

##***part021logb;
sdmeanplotlog <- ggplot(data=mean.fry, aes(x=mean.logfry, y=sd.logfry))+
  ggtitle("SD vs Mean - logged data")+
  geom_point()
sdmeanplotlog
##***part021loge;
# ggsave(plot=sdmeanplot, file='baci-fry-R-021b.png', h=4, w=6, units="in", dpi=300)


# Plot the mean log-fry over time by site
##***part025b;
meanplot <- ggplot(data=mean.fry, aes(x=Year, y=mean.logfry,
                  group=Site, color=SiteClass, shape=Site))+
  ggtitle("Changes in mean over time")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  geom_vline(xintercept=-0.5+min(mean.fry$Year[as.character(mean.fry$Period) == "After"]))
meanplot
##***part025e;
# ggsave(plot=meanplot, file="baci-fry-R-025.png", h=4, w=6, units="in", dpi=300)


#####################################################################################
# The analysis of the mean(log(Fry))

# sink('baci-fry-R-200-type3.txt', split=TRUE)
##***part200b;
# Notice that we use the YearF rather than Year in the random statements
# Because we are analyzing the mean, the Site:Year term is the residual error
result.mean.lmer <- lmerTest::lmer(mean.logfry ~ Period + SiteClass + 
                      Period:SiteClass+
                      (1|YearF)+(1|Site),
                  data=mean.fry) 

##***part200e;
# sink()

##***part200anovab;
anova(result.mean.lmer, ddf="Kenward-Roger")
##***part200anovae;



summary(result.mean.lmer)


# sink('baci-fry-R-200-vc.txt', split=TRUE)
##***part200vcb;
# Extract the variance components
vc <- VarCorr(result.mean.lmer)
vc
##***part200vce;
# sink()

# To extract the variance components for use in R, create a data frame
vc.df <- as.data.frame(vc)
YearSigma    <- vc.df[grep("YearF"   ,vc.df$grp),"sdcor"]
SiteSigma    <- vc.df[grep("Site"    ,vc.df$grp),"sdcor"]
ResidualSigma<- vc.df[grep("Residual",vc.df$grp),"sdcor"] 





# emmeans after a lm() fit
# sink('baci-fry-R-200LSM-SiteClass.txt', split=TRUE)
##***parts200LSM-SiteClassb;
result.mean.lmer.emmo.S <- emmeans::emmeans(result.mean.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass \n\n")
summary(result.mean.lmer.emmo.S)
##***parts200LSM-SiteClasse;
# sink()

# sink('baci-fry-R-200LSM-Period.txt', split=TRUE)
##***part200LSM-Periodb;
result.mean.lmer.emmo.P <- emmeans::emmeans(result.mean.lmer, ~Period)
cat("\nEstimated marginal means for Period \n\n")
summary(result.mean.lmer.emmo.P)
##***part200LSM-Periode;
# sink()

# sink('baci-fry-R-200LSM-int.txt', split=TRUE)
##***part200LSM-intb;
result.mean.lmer.emmo.SP <- emmeans::emmeans(result.mean.lmer, ~SiteClass:Period)
cat("\nEstimated marginal means \n\n")
summary(result.mean.lmer.emmo.SP)
##***part200LSM-inte;
# sink()

# Estimate the BACI contrast
# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.mean.lmer)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]

# sink("baci-fry-R-200baci.txt", split=TRUE)
##***part200bacib; 
# Estimate the BACI contrast along with a se
summary(contrast(result.mean.lmer.emmo.SP, list(baci=c(1,-1,-1,1))), infer=TRUE)
##***part200bacie;
# sink()




# Check the residuals etc
##***part200diagnosticb;
diagplot <- sf.autoplot.lmer(result.mean.lmer)
plot(diagplot)
##***part200diagnostice;
# ggsave(plot=diagplot, file='baci-fry-R-200-diagnostic.png', h=4, w=6, units="in", dpi=300)





#################################################################################
# The analysis of the individual values
# Results will differ slightly from answers on mean because design is unbalanced

# sink('baci-fry-R-300-type3.txt', split=TRUE)
##***part300b;
# Notice that we use the YearF variable in the random effects
result.lmer <- lmerTest::lmer(logfry ~ Period + SiteClass + Period:SiteClass +
                      (1|YearF)+ (1|Site)+(1|YearF:Site), data=fry)
##***part300e;
# sink()


##***part300anovab;
anova(result.lmer, ddf="Kenward-Roger")
##***part300anovae;


summary(result.lmer)


# sink('baci-fry-R-300-vc.txt', split=TRUE)
##***part300vcb;
# Extract the variance components
vc <- VarCorr(result.lmer)
vc
##***part300vce;
# sink()


# emmeans after a lm() fit
# sink('baci-fry-R-300LSM-SiteClass.txt', split=TRUE)
##***parts300LSM-SiteClassb;
result.lmer.emmo.S <- emmeans::emmeans(result.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass \n\n")
summary(result.lmer.emmo.S)
##***parts300LSM-SiteClasse;
# sink()

# sink('baci-fry-R-300LSM-Period.txt', split=TRUE)
##***part300LSM-Periodb;
result.lmer.emmo.P <- emmeans::emmeans(result.lmer, ~Period)
cat("\nEstimated marginal means for Period \n\n")
summary(result.lmer.emmo.P)
##***part300LSM-Periode;
# sink()

# sink('baci-fry-R-300LSM-int.txt', split=TRUE)
##***part300LSM-intb;
result.lmer.emmo.SP <- emmeans::emmeans(result.lmer, ~SiteClass:Period)
cat("\nEstimated marginal means for SiteClass:Period \n\n")
summary(result.lmer.emmo.SP)
##***part300LSM-inte;
# sink()

# Estimate the BACI contrast
# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.lmer)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]


# sink("baci-fry-R-300baci.txt", split=TRUE)
##***part300bacib; 
# Estimate the BACI contrast along with a se
summary(contrast(result.lmer.emmo.SP, list(baci=c(1,-1,-1,1))), infer=TRUE)
##***part300bacie;
# sink()




# Check the residuals etc
##***part300diagnosticb;
diagplot <- sf.autoplot.lmer(result.lmer)
plot(diagplot)
##***part300diagnostice;
# ggsave(plot=diagplot, file='baci-fry-R-300-diagnostic.png', h=4, w=6, units="in", dpi=300)


###################################################################################################
###################################################################################################
###################################################################################################

# Power of the BACI-fry example

# This illustrates the power computations for a BACI design with multiple treatment
# and control sites, multiple years measured before and after treatment is applied, and
# sub-sampling at each year-site combination.

# The proposed experimental design is a BACI experiment with 
# the  monitoring design has the river divided into sites of which 
# some are upstream of the diversion and some are downstream of the diversion. 
# In each site, several minnow traps are set for various lengths of time. 
#
# At the end of the soaking period, the traps are removed and the number of 
# minnows are counted and classified by species.
# The counts are standardized to a common period of time to adjust for the 
# different soak-times. 
#
# This is run for several years before the project starts and after the project starts

# The analysis was originally done on the log-scale from which
# the variance components were extracted (see baci-fry.r code)
# 
# We want to do a power analysis looking at various aspects of the design
#
# The information we need is:
#     Alpha level (usually .05 or .10)
#     Variance components (these were obtained from the baci-fry.r analysis of the previous data)
#        std_site  -> site to site STANDARD DEVIATION
#        std_year  -> year to year STANDARD DEVIATION
#        std_site_year -> site-year interaction STANDARD DEVIATION
#        std_subsample -> within site sub-sampling STANDARD DEVIATION
#     Sub-sampling sample sizes, i.e. the number of minnow traps (sub-samples) set at each site
#       We allow for different sample sizes in the treatment-time combinations
#        but all sites in that combination must have the same number of minnow traps.
#        It is possible to generalize this; contact me for details.
#        n_IA   -> number of subsamples in Treatment-After  combination
#        n_IB   -> number of subsamples in Treatment-Before combination
#        n_CA   -> number of subsamples in Control-After    combination
#        n_CB   -> number of subsamples in Control-Before   combination
#     Number of segments(sites)
#        We allow for a different number of sites for Treatment and Control areaa
#        ns_I     -> number of treatment sites (i.e. downstram of the project)
#        ns_C     -> number of control sites (i.e. upstream of the project)
#     Number of years of monitoring before/after impact
#        ny_B     -> number of years monitoring before project starts
#        ny_A     -> number of years monitoring after  project starts
#     Marginal Means
#        These are used to form the BACI contrast of interest
#        mu_IA, mu_IB, mu_CA, mu_CB (i.e mean of Treatment-After,
#        Treatment-Before, Control-After, Control_Before)
#        These are chosen based on the size of impact that may be biologically important

# This code was originally created by 
#    Tony Booth
#    Department of Ichthyology and Fisheries Science, 
#    Rhodes University, 
#    PO Box 94, 
#    Grahamstown 6140 SOUTH AFRICA
#    t.booth@ru.ac.za

# The computations are based on the 
#    Stroup, W. W. (1999)
#    Mixed model procedures to assess power, precision, and sample size in the design of experiments.
# paper where "dummy" data is generated and "analyzed" and the resulting F-statistics etc
# provide information needed to compute the power


#-----------------------------------------------



# An illustration of how to find the power for one scenario
baci.power(n_IA=3, n_IB=3, n_CA=3, n_CB=3,
           ns_I=3, ns_C=3, 
           ny_B=3, ny_A=2, 
           mu_IA=5.0, mu_IB=5.5, mu_CA=4.5, mu_CB=4.5, 
           sdYear=0.25, sdSite=0.75, sdSiteYear=0.1, sdResid=0.75)


vc <- as.data.frame(VarCorr(result.lmer))
vc


#find the power under various scenarios
##***part500b;

sdSiteYear= vc[ vc$grp=="YearF:Site",    "sdcor"]
sdResid   = vc[ vc$grp=="Residual",      "sdcor"]
sdYear    = vc[ vc$grp=="YearF",         "sdcor"]
sdSite    = vc[ vc$grp=="Site",          "sdcor"]

cat("Estimated variance components are \n",
        "sdSiteYear ", round(sdSiteYear,2), 
    ";\n sdSite ",     round(sdSite,2),
    ";\n sdYear ",     round(sdYear,2),
    "; \nsdResid ",    round(sdResid,2),"\n")
cat("\n")

scenarios <- expand.grid(n_quad=seq(3,9,3),
                         ns_I = 3, ns_C=3,
                         ny_B =3,
                         ny_A =c(2,3,4),
                         baci_effect=seq(0,0.8,.1),
                         sdYear = sdYear, sdSite=sdSite, sdSiteYear=sdSiteYear, sdResid=sdResid)
head(scenarios)

power <- plyr::adply(scenarios,1,function(x){
  #browser()
  power <- baci.power(
    n_IA=x$n_quad, n_IB=x$n_quad, n_CA=x$n_quad, n_CB=x$n_quad, 
    ns_I=x$ns_I, ns_C=x$ns_C, 
    ny_B=x$ny_B, ny_A=x$ny_A, 
    mu_IA=x$baci_effect, mu_IB=0, mu_CA=0, mu_CB=0, 
    sdYear=x$sdYear, sdSite=x$sdSite, sdSiteYear=x$sdSiteYear, sdResid=x$sdResid)
  power
})

head(power)

cat("\n")
head(power[,c("alpha","n_IA","n_IB","n_CA","n_CB","ny_B","ny_A","ns_I","ns_C","baci","power")])

##***part500e;


# sink("baci-fry-power-R-500.txt", split=TRUE)
power[,c("alpha","n_IA","n_IB","n_CA","n_CB","ny_B","ny_A","ns_I","ns_C","baci","power")]
# sink()         



##***part510b;
power.plot <- ggplot(data=power, aes(x=baci_effect, y=power, color=as.factor(n_quad)))+
  ggtitle("Estimated power",
          subtitle=paste("alpha: ", power$alpha[1]))+
  geom_point()+
  geom_line()+
  ylim(0,1)+
  geom_hline(yintercept=0.8, color="blue")+
  xlab("BACI effect size (log scale)")+
  scale_color_discrete(name="# quadrats")+
  facet_wrap(~ny_A, ncol=2, labeller=label_both)
power.plot
##***part510e;
##*


###########################################################################
###########################################################################
###########################################################################

# What to do when number of sub-samples is very large and the design matrix is large
# which makes the power program is very small


sdSiteYear= vc[ vc$grp=="YearF:Site",    "sdcor"]
sdResid   = vc[ vc$grp=="Residual",      "sdcor"]
sdYear    = vc[ vc$grp=="YearF",         "sdcor"]
sdSite    = vc[ vc$grp=="Site",          "sdcor"]

##***part900b;
cat("Estimated variance components are \n",
    "sdSiteYear ", round(sdSiteYear,2), 
    ";\n sdSite ",     round(sdSite,2),
    ";\n sdYear ",     round(sdYear,2),
    "; \nsdResid ",    round(sdResid,2),"\n")
cat("\n")

# We can make an equivalent baci design with 1 observations/site-year (the average)
# and revised sdSiteYear and sdResid

# for example consider a baci design with n_quad=10 in all site.years

cat("\nPower computed using individual observations\n")
power.indiv <- baci.power(n_IA=10, n_IB=10, n_CA=10, n_CB=10,
           ns_I=3, ns_C=3, 
           ny_B=3, ny_A=2, 
           mu_IA=.4, mu_IB=0, mu_CA=0, mu_CB=0, 
           sdYear=sdYear, sdSite=sdSite, sdSiteYear=sdSiteYear, sdResid=sdResid)
power.indiv


# this has equivalent power when you analyze the "averages" and create revised sdSiteYear and sdResid
sdSiteYear.avg <- sqrt(sdSiteYear^2 + sdResid^2/10)
sdResid.avg    <- 0

cat("\nPower computed using averages\n")
power.avg <- baci.power(n_IA=1, n_IB=1, n_CA=1, n_CB=1,
                          ns_I=3, ns_C=3, 
                          ny_B=3, ny_A=2, 
                          mu_IA=.4, mu_IB=0, mu_CA=0, mu_CB=0, 
                          sdYear=sdYear, sdSite=sdSite, sdSiteYear=sdSiteYear.avg, sdResid=sdResid.avg)
power.avg
##***part900e;


##***part901b;
scenarios <- expand.grid(n_quad=seq(3,9,3),
                         ns_I = 3, ns_C=3,
                         ny_B =3,
                         ny_A =c(2,3,4),
                         baci_effect=seq(0,0.8,.1),
                         sdYear = sdYear, sdSite=sdSite, sdSiteYear=sdSiteYear, sdResid=sdResid)
cat("Part of the scenarios dataset\n")
head(scenarios)

power.avg <- plyr::adply(scenarios,1,function(x){
  #browser()
  # do the short cut
  x$sdSiteYear.avg <- sqrt(x$sdSiteYear^2 + x$sdResid^2/x$n_quad) # Note change here
  x$sdResid.avg    <- 0
  
  power <- baci.power(
    n_IA=1, n_IB=1, n_CA=1, n_CB=1, # notice change here
    ns_I=x$ns_I, ns_C=x$ns_C, 
    ny_B=x$ny_B, ny_A=x$ny_A, 
    mu_IA=x$baci_effect, mu_IB=0, mu_CA=0, mu_CB=0, 
    sdYear=x$sdYear, sdSite=x$sdSite, 
    sdSiteYear=x$sdSiteYear.avg, sdResid=x$sdResid.avg) # note change here
  power
})

cat("\nPart of the power dataset computed on averages\n")
head(power.avg)
##***part901e;

# the power plot on the "averages" is the same as before

##***part910b;
power.plot <- ggplot(data=power.avg, aes(x=baci_effect, y=power, color=as.factor(n_quad)))+
  ggtitle("Estimated power computed using short cut",
          subtitle=paste("alpha: ", power$alpha[1]))+
  geom_point()+
  geom_line()+
  ylim(0,1)+
  geom_hline(yintercept=0.8, color="blue")+
  xlab("BACI effect size (log scale)")+
  scale_color_discrete(name="# quadrats")+
  facet_wrap(~ny_A, ncol=2, labeller=label_both)
power.plot
##***part910e;
