# BACI-fish

# This is the example that appears in Smith (2002, Table 6). 
# Samples of fish were taken for a period
# of 12 months before and 13 months after a nuclear
# power plant began operations. 

# The power plant is cooled by water that is drawn from a river.

# When the water exits the plant, its temperature is
# elevated. 
# The concern is that the warmed water will adversely affect the abundance 
# and composition of fish below the plant. 

# [It wasn't clear from Smith (2002) what the response measure is, 
# so I'll arbitrarily assume that it is Counts of fish]


library(ggplot2)
library(lmerTest)
library(emmeans)
library(plyr)

source("../../schwarz.functions.r")
# Get the BACI power program
source("../baci-power.r")

# Read in the actual data
# sink("baci-fish-R-001.txt",split=TRUE)
##***part001b;
fish <- read.csv("baci-fish.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)

fish$trt <- interaction(fish$SiteClass, fish$Site, fish$Period)
fish$Period        <- factor(fish$Period)
fish$SiteClass     <- factor(fish$SiteClass)
fish$Site          <- factor(fish$Site)
fish$SamplingTimeF <- factor(fish$SamplingTime)
fish$trt           <- factor(fish$trt)

fish$logCount      <- log(fish$Count+ .5) # add 1/2 of smallest postive count

fish[1:10,]
##***part001e;
# sink()

str(fish)

# Preliminary plot

# Get plot of series over time
##***part010b;
prelimplot <- ggplot(data=fish, aes(x=SamplingTime, y=logCount,
                          group=Site, color=SiteClass, shape=Site))+
  ggtitle("Fish counts over time")+
  geom_point()+
  geom_line()+
  geom_vline(xintercept=-0.5+min(fish$SamplingTime[as.character(fish$Period) == "After"]))
prelimplot
##***part010e;
#ggsave(plot=prelimplot, file="baci-fish-R-010.png",
#       h=4, w=6, units="in", dpi=300)



# There are several ways in which this BACI design can be analyzed.


###############################################################################
# Do a t-test on the differeces of the averages at each sampling time for each site
# There is only one observtion at each sampling time, so we don't have compute
# the averages here


# sink('baci-fish-R-101.txt',split=TRUE)
##***part101b;
# Convert file from long- to wide-format
# to get the Control and Impact measurements on the same line
# and compute the difference
fish.month <- tidyr::pivot_wider(fish,
                                 id_cols=c("SamplingTime","Period"),
                                 names_from="SiteClass",
                                 values_from="logCount")
fish.month$diff <- fish.month$Impact - fish.month$Control
fish.month[1:10,]
##***part101e;
# sink()


# Plot the difference over time
##***part102b;
plotdiff <- ggplot(data=fish.month, aes(x=SamplingTime, y=diff))+
  ggtitle("Plot of differences over time")+
  ylab("Difference (Impact-Control)")+
  geom_point()+
  geom_line()+
  geom_vline(xintercept=-0.5+min(fish$SamplingTime[as.character(fish$Period) == "After"]))
plotdiff
##***part102e;
#ggsave(plot=plotdiff, file="baci-fish-R-102.png", h=4, w=6, units="in", dpi=300)


# do the two sample t-test not assuming equal variances
# sink('baci-fish-R-104.txt',split=TRUE)
##***part104b;
result <- try(t.test(diff ~ Period, data=fish.month),silent=TRUE)
if(class(result)=="try-error")
     {cat("Unable to do unequal variance t-test because of small sample size\n")} else
     { 
   result$diff.in.means <- sum(result$estimate*c(1,-1))
   names(result$diff.in.means)<- "diff.in.means"
   print(result)
   cat("Mean difference in log(Count) :,",result$diff.in.means," (SE ", result$stderr,")\n")
     }
##***part104e;
# sink()


# do the two sample t-test  assuming equal variances
# sink('baci-fish-R-105.txt',split=TRUE)
##***part105b;
result <- t.test(diff ~ Period, data=fish.month, var.equal=TRUE)
result$diff.in.means <- sum(result$estimate*c(1,-1))
names(result$diff.in.means)<- "diff.in.means"
result
cat("Mean difference in log(Count) :,",result$diff.in.means," (SE ", result$stderr,")\n")
##***part105e;
# sink()




# do the two sample Wilcoxon test
# sink('baci-fish-R-107.txt',split=TRUE)
##***part107b;
result <- wilcox.test(diff ~ Period, data=fish.month, conf.int=TRUE)
result
##***part107e;
# sink()





###############################################################################
# Do a Mixed effect linear model on the individual values


# sink('baci-fish-R-300-type3.txt',split=TRUE)
##***part300b;
# Because there is ONLY one measurement per year, the SamplingTime*Site and
# residual variance are total confounded and cannot be separated. This is 
# the residual term.
result.lmer <- lmerTest::lmer(log(Count+.5) ~ SiteClass+Period+SiteClass:Period +
                                     (1|SamplingTimeF),
	                                    data=fish)
##***part300e;

##***part301b;
anova(result.lmer, ddf="Kenward-Roger")
##***part301e;
# sink()

summary(result.lmer)

# sink('baci-fish-R-300-vc.txt',split=TRUE)
##***part300vcb;
# Extract the variance components
vc <- VarCorr(result.lmer)
vc
##***part300vce;
# sink()


# emmeans after a lm() fit
# sink('baci-fish-R-s300LSM-SiteClass.txt',split=TRUE)
##***parts300LSM-SiteClassb;
result.lmer.emmo.S <- emmeans::emmeans(result.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass \n\n")
summary(result.lmer.emmo.S)
##***parts300LSM-SiteClasse;
# sink()

# sink('baci-fish-R-300LSM-Period.txt',split=TRUE)
##***part300LSM-Periodb;
result.lmer.emmo.P <- emmeans::emmeans(result.lmer, ~Period)
cat("\nEstimated marginal means \n\n")
summary(result.lmer.emmo.P)
##***part300LSM-Periode;
# sink()

# sink('baci-fish-R-300LSM-int.txt',split=TRUE)
##***part300LSM-intb;
result.lmer.emmo.SP <- emmeans::emmeans(result.lmer, ~SiteClass:Period)
cat("\nEstimated marginal means \n\n")
summary(result.lmer.emmo.SP)
##***part300LSM-inte;
# sink()

# Estimate the BACI contrast
# You could look at the entry in the summary table from the model fit, but
# this is dangerous as these entries depend on the contrast matrix.
# It is far safer to the contrast function applied to an emmeans object
temp <- summary(result.lmer)$coefficients # get all the coefficients
temp[grepl("SiteClass",rownames(temp)) & grepl("Period", rownames(temp)),]


# sink("baci-fish-R-300baci.txt",split=TRUE)
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
#ggsave.ggmultiplot(plotdiag, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot  
#       file='baci-fish-R-300-diagnostic.png', h=4, w=6, units="in", dpi=300)



#######################################################################################
#######################################################################################
#######################################################################################

# Power analysis for BACI design with a single control and single impact site, but 
# with multiple years before/after impact

# Example of of the fish analysis.
# Fish counts are quite variable and a BACI difference of .5 on the log scale is important
#    We use mu_TA=20, mu_TB=30, mu_CB=30, mu_cA=30
#
# We have the year-to-year variation (sdYear) 
# The site-to-site variation can be set to any arbitrary value (e.g. 0)
# The residual and site-year variance components are confounded to set
# the sdYearSite and the sdResidual to 0
# Also set the n_TA=1 etc.
#
# We examine various scenarios for the number of sites and number subsamples 




# Illustrate how to get computations for a single scenario
baci.power(n_TA=1, n_TB=1, n_CA=1, n_CB=1, 
           ns_T=1, ns_C=1, 
           ny_B=12, ny_A=13, 
           mu_TA=.2, mu_TB=0, mu_CA=0, mu_CB=0, 
           sdYear=31.71, sdSite=99, sdSiteYear=27.07, sdResid=0)



##***part500b;
# Note that because there is only 1 site, there is NO variance component
# associated with site. Use any value here. In general, the site variance
# component is not important because when you meaure the same sites over time
# the site effects cancel out (like blocking) and so site-to-site variability
# is not important. The site-year variation really is what drives the BACI power.

sdResid = vc[ vc$grp=="Residual",      "sdcor"]
sdYear  = vc[ vc$grp=="SamplingTimeF", "sdcor"]

scenarios <- expand.grid(n_quad=1,
                         n.months_B =seq(12,20,2),
                         n.months_A =seq(13,20,1),
                         baci_effect=0.5,
                         sdYear = sdYear, sdResid=sdResid)
head(scenarios)

power <- plyr::adply(scenarios,1,function(x){
  #browser()
  power <- baci.power(
    n_TA=x$n_quad, n_TB=x$n_quad, n_CA=x$n_quad, n_CB=x$n_quad, 
    ns_T=1, ns_C=1, 
    ny_B=x$n.months_B, ny_A=x$n.months_A, 
    mu_TA=x$baci_effect, mu_TB=0, mu_CA=0, mu_CB=0, 
    sdYear=x$sdYear, sdSite=0, sdSiteYear=x$sdResid, sdResid=0)
  power
})

head(power)

##***part500e;



##***part510b;
power.plot <- ggplot(data=power, aes(x=ny_A, y=power, color=as.factor(ny_B)))+
  ggtitle("Estimated power",
          subtitle=paste("alpha: ", power$alpha[1],"; baci_effect:", power$baci_effect[1]))+
  geom_point()+
  geom_line()+
  ylim(0,1)+
  geom_hline(yintercept=0.8, color="blue")+
  xlab("Number of months of monitoring after impact")+
  scale_color_discrete(name="# months\nbefore")
power.plot
##***part510e;






