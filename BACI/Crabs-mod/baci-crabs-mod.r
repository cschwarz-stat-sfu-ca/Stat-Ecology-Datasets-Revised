# BACI design with multiple controls; 2 factor; interaction;

# A BACI design was used to assess the impact 
#   of cooling water discharge on the density of 
#   shore crabs. 

#   The beach near the outlet of the cooling water 
#   was sampled using several quadrats
#   before and after the plant started operation. 
#   Two control beaches at other locations
#   in the inlet were also sampled. 


library(ggplot2)
library(lmerTest)
library(emmeans)
library(plyr)
library(tidyr)

#source("http://www.stat.sfu.ca/~cschwarz/Stat-Ecology-Datasets/schwarz.functions.r")
source("../../schwarz.functions.r")
source("../baci-power.r")


# Read in the actual data
# sink("baci-crabs-mod-R-001.txt", split=TRUE)
##***part001b;
crabs <- read.csv("baci-crabs-mod.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)

crabs$trt <- interaction(crabs$SiteClass, crabs$Site, crabs$Period)
crabs$Site      <- factor(crabs$Site)
crabs$SiteClass <- factor(crabs$SiteClass)
crabs$Period    <- factor(crabs$Period, levels=c("Before","After"), ordered=TRUE)
crabs$trt       <- factor(crabs$trt)

crabs$logDensity<- log(crabs$Density)

cat("Listing of part of the raw data \n")
crabs[1:10,]
##***part001e;
# sink()

str(crabs)

# Preliminary plot

# Get side-by-side dot plots
##***part010b;
prelimplot <- ggplot(data=crabs, aes(x=trt, y=logDensity))+
  ggtitle("Preliminary plot to look for outliers etc")+
  geom_point(position=position_jitter(w=0.1))+
  geom_boxplot(alpha=0.1)
prelimplot
##***part010e;
#ggsave(plot=prelimplot, file="baci-crabs-mod-R-010.png",
#       h=4, w=6, units="in", dpi=300)



# Get some simple summary statistics
# sink('baci-crabs-mod-R-020.txt', split=TRUE)
##***part020b;
report <- plyr::ddply(crabs, c("Period","SiteClass","Site"), function(x){
   res <- sf.simple.summary(x, "logDensity", crd=TRUE)
   return(res)
})
report
##***part020e;
# sink()

# Draw a profile plot
##***part030b;
profileplot <- ggplot(data=report, aes(x=Period, y=mean, 
                    group=Site, color=SiteClass, shape=Site))+
  ggtitle("Profile plot of crab density")+
  ylab("Density with mean and 95% ci")+
  geom_point(position=position_dodge(w=0.1))+
  geom_line(position=position_dodge(w=0.1))+
  geom_errorbar(aes(ymax=ucl, ymin=lcl), width=0.1, position=position_dodge(w=0.1))+
  geom_point(data=crabs, aes(y=logDensity), position=position_dodge(w=0.2))
profileplot
##***part030e;
#ggsave(plot=profileplot, file="baci-crabs-mod-R-030.png",
#       h=4, w=6, units="in", dpi=300)


# There are several ways in which this BACI design can be analyzed.

#################################################################################
#################################################################################
# Do a t-test on the differences of the averages for each site

# sink('baci-crabs-mod-R-100.txt', split=TRUE)
##***part100b;
# The averages are available in the report dataframe
crabs.avg <- report[,c("Site","SiteClass","Period","mean")]
crabs.avg
##***part100e;
# )

# sink('baci-crabs-mod-R-101.txt', split=TRUE)
##***part101b;
# Reshape the file to get the Before and After measurements on the same line
# 
crabs.site.wide <- tidyr::pivot_wider(crabs.avg,
                                      id_cols=c("Site","SiteClass"),
                                      names_from="Period",
                                      values_from="mean"
                                      )
crabs.site.wide$diff <- crabs.site.wide$After - crabs.site.wide$Before
crabs.site.wide
##***part101e;
# sink()


# do the two sample t-test not assuming equal variances
# Unfortunately, you need at least 2 sites in EACH SiteClass, so this does not work here
# sink('baci-crabs-mod-R-104.txt', split=TRUE)
##***part104b;
result <- try(t.test(diff ~ SiteClass, data=crabs.site.wide),silent=TRUE)
if(class(result)=="try-error")
     {cat("Unable to do unequal variance t-test because of small sample size\n")} else
     { 
  result$diff.in.means <- sum(result$estimate*c(1,-1))
  names(result$diff.in.means)<- "diff.in.means"
  print(result)
  cat("Estimated difference in the means on log scale: ",result$diff.in.means,
      "( SE ", result$stderr, ") \n")
     }
##***part104e;
# sink()


# do the two sample t-test  assuming equal variances
# This only requires at least 2 sites in at least one of the SiteClasses
# sink('baci-crabs-mod-R-105.txt', split=TRUE)
##***part105b;
result <- t.test(diff ~ SiteClass, data=crabs.site.wide, var.equal=TRUE)
result$diff.in.means <- sum(result$estimate*c(1,-1))
names(result$diff.in.means)<- "diff.in.means"
result
cat("Estimated difference in the means on log scale: ",result$diff.in.means,
    "( SE ", result$stderr, ") \n")
##***part105e;
# sink()


#################################################################################
#################################################################################
# Do a Mixed effect linear model on the  averages for each site

crabs.avg <- report[,c("Site","SiteClass","Period","mean")]
crabs.avg

# sink('baci-crabs-mod-R-200-type3.txt', split=TRUE)
##***part200b;
result.avg.lmer <- lmerTest::lmer(mean ~ SiteClass+Period+SiteClass:Period +
                    (1|Site), data=crabs.avg)
##***part200e;

##***part201b;
anova(result.avg.lmer, ddf="Kenward-Roger")
##***part201e;

# sink()



# sink('baci-crabs-mod-R-200-vc.txt', split=TRUE)
##***part200vcb;
# Extract the variance components
vc.avg <- VarCorr(result.avg.lmer)
vc.avg
##***part200vce;
# sink()



# Extract the BACI effect 
# You could get this from the summary table looking at the 
# interaction term, but it is more robust to get this from the
# emmeans
summary(result.avg.lmer)

# emmeans after a lm() fit
# Note that there is a emmeans() function in both the emmeans and lmerTest package
# so we must specify which one we want
# sink('baci-crabs-mod-R-s00LSM-SiteClass.txt', split=TRUE)
##***parts200LSM-SiteClassb;
result.avg.lmer.emmo.S <- emmeans::emmeans(result.avg.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass \n\n")
summary(result.avg.lmer.emmo.S)
##***parts200LSM-SiteClasse;
# sink()

# sink('baci-crabs-mod-R-200LSM-Period.txt', split=TRUE)
##***part200LSM-Periodb;
result.avg.lmer.emmo.P <- emmeans::emmeans(result.avg.lmer, ~Period)
cat("\nEstimated marginal means for Period \n\n")
summary(result.avg.lmer.emmo.P)
##***part200LSM-Periode;
# sink()

# sink('baci-crabs-mod-R-200LSM-int.txt', split=TRUE)
##***part200LSM-intb;
result.avg.lmer.emmo.SP <- emmeans::emmeans(result.avg.lmer, ~SiteClass:Period)
cat("\nEstimated marginal means for SiteClass:Period \n\n")
summary(result.avg.lmer.emmo.SP)
##***part200LSM-inte;
# sink()


# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.avg.lmer)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]


# sink("baci-crabs-mod-R-200baci.txt", split=TRUE)
##***part200bacib; 
# Estimate the BACI contrast along with a se
summary(contrast(result.avg.lmer.emmo.SP, list(baci=c(1,-1,-1,1))), infer=TRUE)
##***part200bacie;
# sink()




# Check the residuals etc
##***part200diagnosticb;
diagplot <- sf.autoplot.lmer(result.avg.lmer)
plot(diagplot)
##***part200diagnostice;
#ggsave.ggmultiplot(plotdiag, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot   
#       file='baci-crabs-mod-R-200-diagnostic.png',
#       h=4, w=6, unit="in", dpi=300)






#################################################################################
#################################################################################
# Do a Mixed effect linear model on the individual values
# Results will differ slightly from above analyeses on the averages because
# the design is not balanced.


# sink('baci-crabs-mod-R-300-type3.txt', split=TRUE)
##***part300b;
result.all.lmer <- lmerTest::lmer(log(Density) ~ SiteClass+Period+SiteClass:Period +
                 (1|Site) + (1|Site:Period), data=crabs)
##***part300e;

##***part301b;
anova(result.all.lmer, ddf="Kenward-Roger")
##***part301e;
# sink()

summary(result.all.lmer)


# sink('baci-crabs-mod-R-300-vc.txt', split=TRUE)
##***part300vcb;
# Extract the variance components
vc <- VarCorr(result.all.lmer)
vc
##***part300vce;
# sink()



# emmeans after a lm() fit
# sink('baci-crabs-mod-R-s300LSM-SiteClass.txt', split=TRUE)
##***parts300LSM-SiteClassb;
result.all.lmer.emmo.S <- emmeans::emmeans(result.all.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass\n\n")
summary(result.all.lmer.emmo.S)
##***parts300LSM-SiteClasse;
# sink()

# sink('baci-crabs-mod-R-300LSM-Period.txt', split=TRUE)
##***part300LSM-Periodb;
result.all.lmer.emmo.P <- emmeans::emmeans(result.all.lmer, ~ Period)
cat("\nEstimated marginal means for Period \n\n")
summary(result.all.lmer.emmo.P)
##***part300LSM-Periode;
# sink()

# sink('baci-crabs-mod-R-300LSM-int.txt', split=TRUE)
##***part300LSM-intb;
result.all.lmer.emmo.SP <- emmeans::emmeans(result.all.lmer, ~SiteClass:Period)
cat("\nEstimated marginal means for SiteClass:Period \n\n")
summary(result.all.lmer.emmo.SP)
##***part300LSM-inte;
# sink()


# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.all.lmer)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]


# sink("baci-crabs-mod-R-300baci.txt", split=TRUE)
##***part300bacib; 
# Estimate the BACI contrast along with a se
summary(contrast(result.all.lmer.emmo.SP, list(baci=c(1,-1,-1,1))), infer=TRUE)
##***part300bacie;
# sink()




# Check the residuals etc
##***part300diagnosticb;
diagplot <- sf.autoplot.lmer(result.all.lmer)
plot(diagplot)
##***part300diagnostice;
#ggsave.ggmultiplot(plotdiag, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot   
#       file='baci-crabs-mod-R-300-diagnostic.png',
#              h=4, w=6, units="in", dpi=300)



###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

# Power analysis for BACI design with multiple sites, but one year before/after

# Example of crab - modifieds.
# Density is about 35 and a BACI difference on the log scale of .2 is important
#    We use mu_TA=.2, mu_TB=0, mu_CB=0, mu_cA=0
#
# We don't have to worry about the sdSite or the sdYear because these "cancel"
# because every site is measured every year.
# The residual standard deviation is obtained from the above it.
#
# We examine various scenarios for the number of sites and number o subsamples 



# An illustration of how to run for single case
baci.power(n_TA=5, n_TB=5, n_CA=5, n_CB=5, 
           ns_T=1, ns_C=2, 
           ny_B=1, ny_A=1, 
           mu_TA=30, mu_TB=35, mu_CA=35, mu_CB=35, 
           sdYear=0, sdSite=3.831, sdSiteYear=1.296, sdResid=3.30)


# get the variance components
vc <- as.data.frame(VarCorr(result.all.lmer))

sdResid   <- vc[ vc$grp=="Residual",    "sdcor"]
sdSite    <- vc[ vc$grp=="Site",        "sdcor"]
sdSiteYear<- vc[ vc$grp=="Site:Period", "sdcor"]

##***part500b;

cat("Estimated variance components are \nsdSiteYear ", round(sdSiteYear,2), 
               ";\n sdSite ", round(sdSite,2),
               "; \nsdResid ",round(sdResid,2),"\n")
cat("\n")

scenarios <- expand.grid(n_quad=seq(5,40,5),
                         ns_T = c(1,2),
                         ns_C = seq(2,10,2),
                         baci_effect=0.2,
                         sdSiteYear=sdSiteYear,
                         sdSite =sdSite,
                         sdResid=sdResid)
cat("Some scenarios\n")
head(scenarios)

power <- plyr::adply(scenarios,1,function(x){
  #browser()
  power <- baci.power(
    n_TA=x$n_quad, n_TB=x$n_quad, n_CA=x$n_quad, n_CB=x$n_quad, 
    ns_T=x$ns_T, ns_C=x$ns_C, 
    ny_B=1, ny_A=1, 
    mu_TA=x$baci_effect, mu_TB=0, mu_CA=0, mu_CB=0, 
    sdYear=0, sdSite=x$sdSite, sdSiteYear=x$sdSiteYear, sdResid=x$sdResid)
  power
})
cat("\nSome of the power computations\n")
head(power)

cat("\n")
head(power[,c("alpha","ns_T","ns_C","n_TA","n_TB","n_CA","n_CB","baci","power")])
##***part500e;



# sink("baci-crabs-mod-power-R-500.txt", split=TRUE) # get a smaller report
power[,c("alpha","ns_T","ns_C","n_TA","n_TB","n_CA","n_CB","baci","power")]
#sink()



##***part510b;
power.plot <- ggplot(data=power, aes(x=n_quad, y=power, color=as.factor(ns_C)))+
  ggtitle("Estimated power",
          subtitle=paste("alpha: ", power$alpha[1]))+
  geom_point()+
  geom_line()+
  ylim(0,1)+
  geom_hline(yintercept=0.8, color="blue")+
  xlab("Number of quadrats at each site:year")+
  scale_color_discrete(name="# Site\nControl")+
  facet_wrap(~ns_T, ncol=2, labeller=label_both)
power.plot
##***part510e;


##***part520b;
# you could also have different sample sizes in impact an control areas.
results <-  baci.power(n_TA=20, n_TB=20, n_CA=5, n_CB=5, 
                            ns_T=1, ns_C=4, 
                            ny_B=1, ny_A=1, 
                            mu_TA=0.2, mu_TB=0, mu_CA=0, mu_CB=0, 
                            sdYear=0, sdSite=sdSite, sdSiteYear=sdSiteYear, sdResid=sdResid)
results
##***part520e;


