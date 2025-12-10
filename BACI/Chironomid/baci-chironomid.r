# BACI-Chironomid

# Taken from Krebs, Ecological Methodology, 2nd Edition. Box 10.3.

# Estimates of chironomid abundance in sediments were taken at one station 
# above and below a pulp mill outflow pipe for 3 years before plant #operation 
# and for 6 years after plant operation.


library(ggplot2)
library(emmeans)
library(lmerTest)
library(plyr)
library(tidyr)

source("../../schwarz.functions.r")

# Read in the actual data
# sink("baci-chironomid-R-001.txt", split=TRUE)
##***part001b;
cat(" BACI design measuring chironomid counts with multiple (paired) yearly measurements before/after \n\n")
chironomid <- read.csv("baci-chironomid.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)
cat("Listing of part of the raw data \n")
head(chironomid)
##***part001e;
# sink()

# The data is NOT in the usual format where there is only one column
# for the response and a separate column indicating if it is a control or 
# impact site. We need to restructure the data
# sink('baci-chironomid-R-301.txt', split=TRUE)
##***part002b;
# We  change the data from wide to long format
chironomid.long <- tidyr::pivot_longer(chironomid, 
                                       cols=c("Control.Site","Impact.Site"),
                                       names_to = "SiteClass",
                                       values_to="Count")
chironomid.long$SiteClass     <-  factor(chironomid.long$SiteClass)
chironomid.long$YearF         <-  factor(chironomid.long$Year)
chironomid.long$Period        <-  factor(chironomid.long$Period)

chironomid.long$logCount      <- log(chironomid.long$Count + .5)

head(chironomid.long)
##***part002e;
# sink()
str(chironomid.long)


# Get plot of series over time
##***part010b;
prelimplot <- ggplot(data=chironomid.long,
                     aes(x=Year, y=logCount, group=SiteClass, color=SiteClass))+
  ggtitle("Fish counts over time")+
  geom_point()+
  geom_line()+
  geom_vline(xintercept=-0.5+min(chironomid.long$Year[as.character(chironomid.long$Period) == "After"]))
prelimplot
##***part010e;
#ggsave(plot=prelimplot, file="baci-chironomid-R-010.png", h=4, w=6, units="in", dpi=300)



# There are several ways in which this BACI design can be analyzed.

###########################################################################
# Do a t-test on the differces of the averages for each year
# Because only one measurement was taken at each site in each year, we
# don't have to first average. We can use the wide format data.


# sink('baci-chironomid-R-101.txt', split=TRUE)
##***part101b;
chironomid$diff <- log(chironomid$Impact.Site+.5) - log(chironomid$Control.Site+.5)
head(chironomid)
##***part101e;
# sink()


# Plot the difference over time
##***part102b;
plotdiff <- ggplot(data=chironomid, aes(x=Year, y=diff))+
  ggtitle("Plot of differences over time")+
  ylab("Difference (Impact-Control)")+
  geom_point()+
  geom_line()+
  geom_vline(xintercept=-0.5+min(chironomid$Year[as.character(chironomid$Period) == "After"]))
plotdiff
##***part102e;
#ggsave(plot=plotdiff, file="baci-chironomid-R-102.png", h=4, w=6, units="in", dpi=300)



# do the two sample t-test not assuming equal variances
# sink('baci-chironomid-R-104.txt', split=TRUE)
##***part104b;
result <- try(t.test(diff ~ Period, data=chironomid),silent=TRUE)
if(class(result)=="try-error")
     {cat("Unable to do unequal variance t-test because of small sample size\n")} else
     { 
   result$diff.in.means <- sum(result$estimate*c(1,-1))
   names(result$diff.in.means)<- "diff.in.means"
   print(result)
   cat("Estimated difference in means :",result$diff.in.means, " (SE ",result$stderr, ")\n")
     }
##***part104e;
# sink()


# do the two sample t-test  assuming equal variances
# sink('baci-chironomid-R-105.txt', split=TRUE)
##***part105b;
result.ev <- t.test(diff ~ Period, data=chironomid, var.equal=TRUE)
result.ev$diff.in.means <- sum(result$estimate*c(1,-1))
names(result.ev$diff.in.means)<- "diff.in.means"
result.ev
cat("Estimated difference in means on log scale:",result.ev$diff.in.means, " (SE ",result.ev$stderr, "\n")
##***part105e;
# sink()




# do the two sample Wilcoxon test
# sink('baci-chironomid-R-107.txt', split=TRUE)
##***part107b;
result <- wilcox.test(diff ~ Period, data=chironomid, conf.int=TRUE)
result
##***part107e;
# sink()





##################################################################
# Do a Mixed effect linear model on the individual values



# Because there is ONLY one measurement per year, the SamplingTime*Site and
# residual variance are total confounded and cannot be separated. This is 
# the residual term.
# sink('baci-chironomid-R-300-type3.txt', split=TRUE)
##***part300b;
result.lmer <- lmerTest::lmer(log(Count+.5) ~ SiteClass+Period+SiteClass:Period 
                              + (1|YearF),
                data=chironomid.long)
##***part300e;

##***part301b;
anova(result.lmer, ddf="Kenward-Roger")
##***part301e;
# sink()

summary(result.lmer)


# sink('baci-chironomid-R-300-vc.txt', split=TRUE)
##***part300vcb;
# Extract the variance components
vc <- VarCorr(result.lmer)
vc
##***part300vce;
# sink()


# emmeans after a lm() fit
# sink('baci-chironomid-R-s300LSM-SiteClass.txt', split=TRUE)
##***parts300LSM-SiteClassb;
result.lmer.emmo.S <- emmeans::emmeans(result.lmer, ~SiteClass)
cat("\nEstimated marginal means for SiteClass \n\n")
summary(result.lmer.emmo.S)
##***parts300LSM-SiteClasse;
# sink()

# sink('baci-chironomid-R-300LSM-Period.txt', split=TRUE)
##***part300LSM-Periodb;
result.lmer.emmo.P <- emmeans::emmeans(result.lmer, ~Period)
cat("\nEstimated marginal means \n\n")
summary(result.lmer.emmo.P)
##***part300LSM-Periode;
# sink()

# sink('baci-chironomid-R-300LSM-int.txt', split=TRUE)
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


# sink("baci-chironomid-R-300baci.txt", split=TRUE)
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
#       file='baci-chironomid-R-300-diagnostic.png',
#       h=6, w=6, units="in", dpi=300)







