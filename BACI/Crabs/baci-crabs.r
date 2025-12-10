# BACI-crabs

# A simple BACI design was used to assess the impact of cooling water 
# discharge on the Density of shore crabs. 
# The beach near the outlet of the cooling water was sampled 
# using several quadrats before and after the plant started operation, 
# as was a control beach on the other side of the body of water.

# As explained in the, notes, analyze on the log() scale.

library(ggfortify)
library(ggplot2)
library(emmeans)
library(plyr)

source("../../schwarz.functions.r")
source("../baci-power.r")


# Read in the actual data and define the factors of interest
# sink("baci-crabs-R-001.txt")
##***part001b;
crabs <- read.csv("baci-crabs.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)

crabs$trt <- interaction(crabs$SiteClass, crabs$Period)
crabs$SiteClass <- factor(crabs$SiteClass)
crabs$Period    <- factor(crabs$Period, levels=c("Before","After"), ordered=TRUE) # sort correctly
crabs$trt       <- factor(crabs$trt)

crabs$logDensity <- log(crabs$Density) # we will wan to analyze on the log scale

cat("Listing of part of the raw data \n")
crabs[1:10,]
##***part001e;
# sink()

str(crabs)
# Preliminary plot

##***part010b;
prelimplot <- ggplot(data=crabs, aes(x=trt, y=logDensity))+
  ggtitle("Preliminary plot to look for outliers etc")+
  geom_point(position=position_jitter(w=0.1))+
  geom_boxplot(alpha=0.1)
prelimplot
##***part010e;

#ggsave(plot=prelimplot, file="baci-crabs-R-010.png", h=4, w=6, units="in", dpi=300)

 

# Get some simple summary statistics
# sink('baci-crabs-R-020.txt', split=TRUE)
##***part020b;
report <- plyr::ddply(crabs, c("Period","SiteClass"), function(x){
   res <- sf.simple.summary(x, "logDensity", crd=TRUE)
   return(res)
})
report
##***part020e;
# sink()

# Draw a profile plot includeing the raw data
##***part030b;
profileplot <- ggplot(data=report, aes(x=Period, y=mean, group=SiteClass, color=SiteClass))+
  ggtitle("Profile plot of crab density")+
  ylab("logDensity with mean and 95% ci")+
  geom_point(position=position_dodge(w=0.1))+
  geom_line(position=position_dodge(w=0.1))+
  geom_errorbar(aes(ymax=ucl, ymin=lcl), width=0.1, position=position_dodge(w=0.1))+
  geom_point(data=crabs, aes(y=logDensity), position=position_dodge(w=0.2))
profileplot
##***part030e;

#ggsave(plot=profileplot, file="baci-crabs-R-030.png", h=4, w=6, units="in", dpi=300) 

# Fit the linear model, get the anova table, and the usual stuff
# CAUTION!!! Because the design is unbalanced, the default model
# fit by aov gives the WRONG sum of squares and F-tests.
# The default tests are "sequential tests" where terms are added
# in the order specified. You want the marginal tests 
# (which are reported in JMP or SAS)
#
# Read the entry at 
#  http://r-eco-evo.blogspot.com/2007/10/infamous-type-iii-ss-before-i-started_20.html
#
##***part100b;
cat("The sum of squares and F-tests from the anova() below are INCORRECT\n")
cat("in unbalanced data because they are sequential and only adjust for effects\n")
cat("that enter the model prior to the term in question.")
result.lm <- lm( log(Density) ~ SiteClass + Period + SiteClass:Period, data=crabs)
cat("\nAnalysis of variance -- this is NOT CORRECT because design is unbalanced \n")
anova(result.lm)

cat("\n\nUse the Type III tests from the Anova() function from the car package")
cat("\nbut you need to set the treatment contrasts to sum rather than treatment")
cat("\nSee http://r.789695.n4.nabble.com/Type-I-v-s-Type-III-Sum-Of-Squares-in-ANOVA-td1573657.html\n")

result.lm2 <- lm( log(Density) ~ SiteClass + Period + SiteClass:Period, data=crabs,
                  contrasts=list(SiteClass="contr.sum", Period='contr.sum'))
##***part100e;

##***part100effectb;
car::Anova(result.lm2,type=3)
##***part100effecte;



# Check the residuals etc

##***part100diagnosticb;
diagplot <- autoplot(result.lm2)
diagplot
##***part100diagnostice;
#ggsave.ggmultiplot(plotdiag, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot   
#       file='baci-crabs-R-100-diagnostic.png', h=4, w=6, units="in", dpi=300)


# Estimate the marginal means and the various effects
# sink('baci-crabs-R-100LSM-SiteClass.txt', split=TRUE)
##***part100LSM-SiteClassb;
result.emmo.S <- emmeans::emmeans(result.lm2, ~SiteClass)
cat("\n\n Estimated marginal means for SiteClass \n\n")
summary(result.emmo.S, infer=TRUE)
##***part100LSM-SiteClasse;
# sink()

# sink('baci-crabs-R-100LSM-Period.txt', split=TRUE)
##***part100LSM-Periodb;
result.emmo.P <- emmeans::emmeans(result.lm2, ~Period)
cat("\n\n Estimated marginal means for Period \n\n")
summary(result.emmo.P)
##***part100LSM-Periode;
# sink()

# sink('baci-crabs-R-100LSM-int.txt', split=TRUE)
##***part100LSM-intb;
result.emmo.SP <- emmeans::emmeans(result.lm2, ~SiteClass:Period)
cat("Estimated marginal means \n\n")
summary(result.emmo.SP, infer=TRUE)
##***part100LSM-inte;
# sink()

# Estimate the BACI contrast along with a se
# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.lm)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]

# sink("baci-crabs-R-100baci.txt", split=TRUE)
##***part100bacib;
summary(contrast(result.emmo.SP, list(baci=c(1,-1,-1,1))), infer=TRUE)
##***part100bacie;
# sink()




###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------

# Power analysis for simple BACI analysis

# For a simple BACI there is 
#   1 site in each of treatment and control -> ns_T=1, ns_C=1
#   1 year before/after                     -> ny_B=1, ny_A=1
# We cannot separate the site-year interaction from the means so we hope
# for the best and suppose that sdSiteYear=0

# Example of crabs.
# We think that a 20 percentage difference in the differential change is important.
# We arbitarily choose 4 mu values that give us this change
#  E.g. mu_tb =0, mu_ta=.20, mu_cb=0, mu_ca=0
#
# We don't have to worry about the sdSite or the sdYear because these "cancel"
# because every site is measured every year.
# The residual standard deviation is about 5
#
# We examine various scenarios for the number of subsamples at each of the 4 site-year combinations


# Get the BACI power program


# illustration of how the power program work for a specific scenario
sdResid <- summary(result.lm2)$sigma
baci.power(n_TA=5, n_TB=5, n_CA=5, n_CB=5, 
           ns_T=1, ns_C=1, 
           ny_B=1, ny_A=1, 
           mu_TA=.2, mu_TB=0, mu_CA=0, mu_CB=0, 
           sdYear=0, sdSite=0, sdSiteYear=0, sdResid=sdResid)



##***part500b;
scenarios <- expand.grid(n_quad=seq(5,20,2),
                         baci_effect=0.2,
                         sdResid=sdResid)
scenarios

power <- plyr::adply(scenarios,1,function(x){
   #browser()
   power <- baci.power(
            n_TA=x$n_quad, n_TB=x$n_quad, n_CA=x$n_quad, n_CB=x$n_quad, 
            ns_T=1, ns_C=1, 
            ny_B=1, ny_A=1, 
            mu_TA=x$baci_effect, mu_TB=0, mu_CA=0, mu_CB=0, 
            sdYear=0, sdSite=0, sdSiteYear=0, sdResid=x$sdResid)
   power
})

power

cat("\n")
power[,c("alpha","n_TA","n_TB","n_CA","n_CB","baci","power")]
##***part500e;



##***part510b;
power.plot <- ggplot(data=power, aes(x=n_quad, y=power))+
   ggtitle("Estimated power",
            subtitle=paste("alpha: ", power$alpha[1]))+
   geom_point()+
   geom_line()+
   ylim(0,1)+
   geom_hline(yintercept=0.8, color="blue")+
   xlab("Number of quadrats at each site:year")
power.plot
##***part510e;


