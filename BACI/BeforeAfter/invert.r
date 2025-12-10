# Invertebrate data from a simple before/after analysis.

#   Two streams were measured at a small hydro electric   project. At each stream
#   multiple samples were taken on the steam, and the number of invertebrates
#   was counted.

#   After 5 years, the hydro electric plant started up, and an additional 5 years
#   of data were collected.
 
#   Is there evidence of a change in abundance after the plant starts?

#  Note that this is NOT a BACI design, and is a VERY poor substitute for such.
#  The key problem is that changes may have occured over time unrelated
#   to the impact, so a "change" may simple be natural. 

#  Also, in this analysis we have ignored autocorrelation which can be
#   a serious problem in long-term studies. 

library(ggfortify)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(plyr)
library(pwr)

source("../../schwarz.functions.r")
source("../before-after-power-stroup.r")

# Read in the data; declare Period as a factor; create factor for Year
# sink('invert-R-001.txt', split=TRUE)
##***part001b;
invert <- read.csv('invert.csv', header=TRUE, as.is=TRUE, strip.white=TRUE)

# Declare Period as an ordered factor so that sorts in proper order
invert$Period  <- factor(invert$Period, levels=c("Before","After"), ordered=TRUE)
# We also want Year to be categorical variable so make a new variable
invert$YearF   <- factor(invert$Year)
# Stream should also be defined as a factor
invert$StreamF <- factor(invert$Stream)

# Compute the log(density) 
invert$logCount <- log(invert$count)

# It turns out that missing values cause subtle problems in some R
# functions that deal with the fitted object. Sigh... Usually, R
# is sensible enough, but not all package authors did a good job.
cat("Dimensions of full dataset :",dim(invert), "\n")
invert <- invert[complete.cases(invert),]
cat("Dimensions of reduced dataset after removing missing values: ", dim(invert), "\n")
invert[1:12,]
##***part001e;
# sink()

str(invert) # check the structure


# look at the std dev to mean plot 
# Compute the averages etc
# sink('invert-R-020.txt', split=TRUE)
##***part020b;
mean.invert <- plyr::ddply(invert, c("Year","YearF","Stream","StreamF","Period"), function(x){
  # compute sample size, missing values, mean, and sd
  n <- nrow(x)
  n.miss <- sum(is.na(x$count))
  count.mean <- mean(x$logCount, na.rm=TRUE)
  count.sd   <- sd  (x$logCount, na.rm=TRUE)
  res <- c(n=n, n.miss=n.miss, count.mean=count.mean, count.sd=count.sd) # put results together
  res
  })
mean.invert
##***part020e;
# sink()


##***part060b;
# Plot the log(std) vs.  log(mean)
# Because there is no relationship between the mean and the sd,
# this plot shows that indicates that no transformation is needed.

sd.plot <- ggplot(data=mean.invert, aes(x=log(count.mean), y=log(count.sd), 
                                        shape=StreamF, color=StreamF))+
  ggtitle("SD vs. Mean")+
  geom_point(size=4)
sd.plot
##***part060e;

# ggsave(plot=sd.plot, file='invert-R-060-varmeanplot.png', height=4, width=6, units="in", dpi=300)



# Plot the mean count over time by stream
##***part070b;
plot.mean <- ggplot(data=invert, aes(x=Year, y=logCount, 
                    shape=StreamF, group=StreamF, color=StreamF))+
  ggtitle("Preliminary plot of data and means")+
  geom_point(size=2, position=position_dodge(w=0.2))+
  geom_vline(xintercept=.5+max(invert$Year[as.character(invert$Period)=="Before"],na.rm=TRUE))+
  stat_summary(fun="mean", geom="line") # , position=position_dodge(w=0.2))
plot.mean
##***part070e;


# ggsave(plot=plot.mean, file='invert-R-070-trendplot.png', height=4, width=6, units='in', dpi=300)




############################################################################################

# Analysis of Stream 1 alone  to illustrate how to analyze one stream's data

##***part99b;
# Subset the data
stream1 <- invert[ invert$Stream ==1,]
head(stream1)

# The analysis of the mean(Count)) for each year
mean.invert1 <- mean.invert[ mean.invert$Stream ==1, ]
mean.invert1 
##***part99e;

# Test for a period effect using the mean count data
# sink('invert-R-100.txt', split=TRUE)
##***part100b;
result100 <- lm( count.mean ~ Period, data=mean.invert1)
anova(result100)
##***part100e;
# sink()


# Estimate the period effect
# Do NOT use the summary output because it gives odd results
# with ordered factors. Compare the two following models summary output.
# They should give identical estimates but are not even though the p-value is the same
summary(result100) #with ordered factor
summary(lm(count.mean ~ factor(as.character(Period)), data=mean.invert1))
# It is better to formally estimate the difference using the
# emmeans package and use the pairs to estimate the difference
# sink('invert-R-105.txt', split=TRUE)
##***part105b;
# Estimate the marginal means and the differences.
result100.emmo <- emmeans::emmeans(result100, ~Period)
cat("\nEstimated marginal means for each Period \n")
summary(result100.emmo)
cat("Estimated difference between means in Before and After periods\n")
summary(pairs(result100.emmo), infer=TRUE)
##***part105e;
# sink()



# Residual and other plots

##***part106b;
plot.diag100 <- ggplot2::autoplot(result100)
plot.diag100
##***part106e;
# ggsave.ggmultiplot(plot.diag100, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot  
#       file='invert-R-106.png', height=4, width=6, units="in", dpi=300)




#################################################################################
# Analysis of the individual quadrat data for Stream 1.
# 

# sink('invert-R-200-type3.txt', split=TRUE)
##***part200b;
# Don't forget that you need a Year as a factor.
# We created YearF as a factor when we read in the data
# We don't need to transform counts because they are largish.
result200 <- lmerTest::lmer(log(count) ~ Period + (1|YearF), data=stream1)
anova(result200, ddf="Kenward-Roger")
##***part200e;
# sink()

ba.invert.result200.fit <- result200

# Again, don't trust the estimates from the summary table.
summary(result200)


# sink('invert-R-200-vc.txt', split=TRUE)
##***part200-vcb;
# Extract the variance components.
vc <- VarCorr(result200)
vc
##***part200-vce;
# sink()


# emmeans after a lmer() fit
# sink('invert-R-200LSM-Period.txt', split=TRUE)
##***part200LSM-Periodb;
# Estimate the marginal means.
result200.emmo <- emmeans::emmeans(result200, ~Period)
cat("Estimated marginal means\n")
summary(result200.emmo)
##***part200LSM-Periode;
# sink()


# sink("invert-R-200baci.txt", split=TRUE)
##***part200bacib; 
# Estimate the period effect along with a se
cat("Estimated difference between means in Before and After periods\n")
summary(pairs(result200.emmo), infer=TRUE)
##***part200bacie;
# sink()


# Diagnostic plots for mixed models are not as standardized as for
# ordinary linear models so there are several than can be produced.
##***part200diagnosticb;
plot.diag200 <- sf.autoplot.lmer(result200)
plot(plot.diag200)
##***part200diagnostice;

# ggsave.ggmultiplot(plot.diag200, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot   
#       file='invert-R-200-diagnostic.png', height=4, width=6, units="in", dpi=300)






#--------------------------------------------------- Analyze all streams together --------------



# The analysis of the means values



# sink('invert-R-300-type3.txt', split=TRUE)
##***part300b;
# Model for all streams together
# Note that Stream, Period, and Year all need to be factors
# These were defined as such when the data were read in.
result300 <- lmerTest::lmer(count.mean~ Period + (1|Stream) + (1|YearF), data=mean.invert)
anova(result300, ddf='Kenward-Roger')
##***part300e;
# sink()

# Again be careful with the summary table because Period is an ordered factor
summary(result300)


# sink('invert-R-300-vc.txt', split=TRUE)
##***part300-vcb;
# Extract the variance components
vc <-VarCorr(result300)
vc
##***part300-vce;
# sink()


# Estimate the marginal means
# sink('invert-R-300LSM-Period.txt', split=TRUE)
##***part300LSM-Periodb;
# Estimate the marginal means
result300.emmo <- emmeans::emmeans(result300, ~Period)
cat("Estimated marginal means for each Period \n")
summary(result300.emmo)
##***part300LSM-Periode;
# sink()


# sink("invert-R-300baci.txt", split=TRUE)
##***part300bacib;
# Estimate the Period effect
cat("Estimated difference between means in Before and After periods\n")
summary(pairs(result300.emmo), infer=TRUE)
##***part300bacie;
# # sink()




# Check the lmer residuals 
##***part300diagnosticb;
plot.diag300 <- sf.autoplot.lmer(result300)
plot(plot.diag300)
##***part300diagnostice;

# ggsave.ggmultiplot(plot.diag300, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot  
#       file='invert-R-300-diagnostic.png', height=4, width=6, units="in", dpi=300)




#######################################################################################
# Now for the analysis of the entire dataset (whew)


# # sink('invert-R-400-type3.txt', split=TRUE)
##***part400b;
# Year, Period, Stream must all be defined as factors.
# This was done when the data were read in.
result400 <- lmerTest::lmer(log(count) ~ Period + (1|StreamF)+(1|YearF)+(1|StreamF:YearF), data=invert)
anova(result400, ddf="Kenward-Roger")
##***part400e;
# sink()

ba.invert.result400.fit <- result400

# Again be very careful about using the results
# from the summary() because Period is defined as an ordered factor
summary(result400)


# sink('invert-R-400-vc.txt', split=TRUE)
##***part400-vcb;
# Extract the variance components
vc <- VarCorr(result400)
vc
##***part400-vce;
# sink()




# sink('invert-R-400LSM-Period.txt', split=TRUE)  # Note these have the wrong standard error
##***part400LSM-Periodb;
# Estimate the marginal means
result400.emmo <- emmeans::emmeans(result400, ~Period)
cat("Estimated marginal means for each Period \n")
summary(result400.emmo)
##***part400LSM-Periode;
# sink()


# sink("invert-R-400baci.txt", split=TRUE)  # note these have the wrong df
##***part400bacib; 
# Estimate the Period effect
cat("Estimated difference between means in Before and After periods\n")
summary(pairs(result400.emmo), infer=TRUE)
##***part400bacie;
# sink()


# diagnostic plots for the lmer fit

##***part400diagnosticb;
plot.diag400 <- sf.autoplot.lmer(result400)
plot(plot.diag400)
##***part400diagnostice;

# ggsave.ggmultiplot(plot.diag400, # see https://github.com/sinhrks/ggfortify/issues/98 for bug in autoplot   
#       file='invert-R-400-diagnostic.png', height=4, width=6, units="in", dpi=300)


#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

# Power analysis for before vs. after study. 
# Invertebrate data from a simple before/after analysis.

#   Note that this is NOT a BACI design, and is a VERY poor substitute for such.
#   The key problem is that changes may have occured over time unrelated
#   to the impact, so a "change" may simple be natural. 

#   Also, in this analysis we have ignored autocorrelation which can be
#   a serious problem in long-term studies. 

#---------------------------------
# 1. Power analysis for a single stream before/after study 
#      based on the analysis of the averages.
#   Here the residual error is a combination of the year-to-year variation
#   and the quadrat-to-qaudrat variation based on a "average" number of quadrats
#   sampled each year. So the only thing that can change is the number of years before
#   or after impact occurred.

#   The power analysis is exactly the same as that for a two-sample t-test with 
#   the number of years before or after as the sample sizes

delta <- .15 # proportional change to detect

single.stream.avg.sdResid <- summary(result100)$sigma


#  An illustration of how to compute power for a single scenario using the pwr package
pwr.t2n.test(n1=5, n2=5, d=delta/single.stream.avg.sdResid, alternative="two.sided")

#  An illustration of how to use the before.after.power.stroup() function
#  Note that the X.after and X.before should have NO elements in common
#  Not unexpectedly you get the same answers as above

before.after.power.stroup(Shift=delta, X.before=1:5, X.after=6:10, 
                          Process.SD=single.stream.avg.sdResid, 
                          Sampling.SD=0 )


# sink("invert-power-R-100.txt", split=TRUE)
##***part500b;
# Note that d=effect size = difference to detect / standard deviation
power <- plyr::ldply(seq(5,20,1), function(n2){
  power1 <- pwr.t2n.test(n1=5, n2=n2, d=delta/single.stream.avg.sdResid, 
                         alternative="two.sided")$power
  power2 <- before.after.power.stroup(Shift=delta, X.before=1:5, X.after=5+(1:n2), 
                                      Process.SD=single.stream.avg.sdResid,
                                      Sampling.SD=0 )
  res <- c(n2=n2, power.pwr=power1, power.stroup=power2[1,"power.2s"])
  res     
})

cat("\nPower for detecting ", delta," proportional change using one stream with 5 years before and \n")
cat(    "different number of years after; average of 4.5 quadrats/year\n")
power
##***part500e;
# sink()


# Egads - over 20 years post impact are required to detect the proportional change in density !!! 
# Perhaps we can improve things by changing the number of quadrats measured each years.
#   So we now need to separte the quadarat-to-quadrat and year-to-year variation.
#   We need to refer to the analysis on the individual values 




#----------------------------------------------
# 2. Power analysis for a single stream before/after study 
#      after changing the number of quadrats measured each year.

#   We need estimates of the quadrat-to-quadrat variation and the year-to-year variation.

#   The power analysis is the same as for a two-sample t-test with sub-sampling 

delta <- .15 # proportional change to detect

vc <- as.data.frame(VarCorr(ba.invert.result200.fit))
single.stream.indiv.sdResid <- vc[ vc$grp=="Residual",  "sdcor"]
single.stream.indiv.sdYear  <- vc[ vc$grp=="YearF",     "sdcor"]

cat("Single stream - individual values - sdYear:  ", single.stream.indiv.sdYear , "\n")
cat("Single stream - individual values - sdResid: ", single.stream.indiv.sdResid , "\n")

cat("Compare single.stream.avg.sdResid^2", single.stream.avg.sdResid^2, "\n")
cat("to sum of vc from individual analysis ", single.stream.indiv.sdYear^2 + single.stream.indiv.sdResid^2/4.5, "\n")


# We could use the previous method with different estimates of residual variance depending on the number of quadrats
#   measured per year. For example, suppose you could do 10 quadrats/year. 

stddev <- sqrt(single.stream.indiv.sdYear^2 + single.stream.indiv.sdResid^2/10)
pwr.t2n.test(n1=5, n2=5, d=delta/stddev, alternative="two.sided")

# It is much easier to use the before.after.power.stroup function as you 
# don't have to compute the new standard deviation, but rather simply
# specify the process and sampling std deviation along with the approriate 
# number of replcates for X.before and X.after. Note that you don't have
# to have a balanced design with this function.
#  Note that the X.after and X.before should have NO elements in common
before.after.power.stroup(Shift=delta, 
                          X.before=rep(1:5, each=10), 
                          X.after=rep(6:10, each=10), 
                          Process.SD =single.stream.indiv.sdYear, 
                          Sampling.SD=single.stream.indiv.sdResid )





##***part600b;
delta <- .15 # proportional change to detect

vc <- as.data.frame(VarCorr(ba.invert.result200.fit))
single.stream.indiv.sdResid <- vc[ vc$grp=="Residual",  "sdcor"]
single.stream.indiv.sdYear  <- vc[ vc$grp=="YearF",     "sdcor"]

cat("Single stream - individual values - sdYear:  ", single.stream.indiv.sdYear , "\n")
cat("Single stream - individual values - sdResid: ", single.stream.indiv.sdResid , "\n")

stddev <- sqrt(single.stream.indiv.sdYear^2 + single.stream.indiv.sdResid^2/10)

power <- plyr::ldply(seq(5,20,1), function(n2){
  # Note that d=effect size = difference to detect / standard deviation 
  # in the pwr.t2n.test() function
  power1 <- pwr.t2n.test(n1=5, n2=n2, d=delta/stddev, 
                         alternative="two.sided")$power
  power2 <-before.after.power.stroup(Shift=delta, 
                                     X.before=rep(1:5, each=10), 
                                     X.after=5+rep(1:n2, each=10), 
                                     Process.SD =single.stream.indiv.sdYear, 
                                     Sampling.SD=single.stream.indiv.sdResid )
  
  res <- c(n2=n2, power.pwr=power1, power.stroup=power2[1,"power.2s"])
  res     
})
cat("Power for detecting ", delta," proportional change using one stream with 5 years before and \n")
cat("different number of years after; 10 quadrats/year\n")
power
##***part600e;


# Power has improved, but notice that even if you took infinite quadrats/year you can't reduce
#   the residual variance below that of the year-to-year variation.

# Note that the previous power analysis is to detect a difference for THAT particular stream only 
#   and inference would be limited to that single stream only


# -----------------------------------------------------------

# 3. Power analysis for multiple streams before/after study 
#      based on changing the number of streams and number of quadrats measured per year.

# Method of Stroup not readily available in R but send me an email and I'll show you
# how to do this.

# -----------------------------------------------------------


# 4. Power analysis for multiple streams before/after study 
#      based on changing the number of streams and number of quadrats measured per year.
# 
#   Now inference is for the entire set of streams.
#   
#   We need estimates of the stream-to-stream, year-to-year, stream-year interaction variation, and the
#   quadrat-to-quadrat variation.
#  

vc <- as.data.frame(VarCorr(result400))
vc

multiple.stream.indiv.sdResid     <- vc[ vc$grp=="Residual",         "sdcor"]
multiple.stream.indiv.sdYear      <- vc[ vc$grp=="YearF",            "sdcor"]
multiple.stream.indiv.sdSite      <- vc[ vc$grp=="StreamF",          "sdcor"]
multiple.stream.indiv.sdSiteYear  <- vc[ vc$grp=="StreamF:Year",     "sdcor"]


# You can use the Stroup method for power analysis if you want to change the number of years for each site etc
#   or you can use the previous methods if you assume that every site is measured every year and the same
#   number of quadrats is measured in each year-site combination
# I have not implemented this method in R yet.

# Alternatively, if the design is balanced, we can use the previous results by computing
#   the residual variance as
#     year-to-year variation 
#     year-stream interaction variation / n_stream
#     quadrat variance = n_streams*n_quadrates
#   Note that because the same site is measured over time, it "disappears" from the power computation
#   in the same way as blocking causes block effects to disappear.

# Once again notice that the limiting term for the power analysis is the year-to-year variation
#   and that even if you measure thousands of streams, you can't ever reduce the impact of this term.

# Note that for a single stream, the "year-to-year" variance really is a combination of 
# the year-to-year variance and the year-stream interaction variance components as these
# cannot be separated if you have a single stream.

# It appears that multiple-streams have no impact vs. a single stream, but this is misleading. As noted
# earlier, with a single stream, you cannot separate the year-to-year variance and the stream-year
# interaction variance components. With multiple streams, you can separate the two components and 
# reduce the impact of the stream-year interaction by measuring multiple streams.  


#   Additionally, the scope of inference is different. In the case of a single stream, you are interested 
#   ONLY in that stream. In the case of
#   multiple streams, you want to generalize to all streams. In this case, the stream-year interaction
#   variance was estimated to be fzero, so this implies that there is no evidence that individual streams
#   behave differently from each other. If this variance component was non-zero, then you have to deal
#   with the additional variation caused by individual streams behaving differently over time from each other. 

# sink("invert-power-R-800.txt", split=TRUE)
##***part800b;

vc

multiple.stream.indiv.sdResid     <- vc[ vc$grp=="Residual",         "sdcor"]
multiple.stream.indiv.sdYear      <- vc[ vc$grp=="YearF",            "sdcor"]
multiple.stream.indiv.sdSite      <- vc[ vc$grp=="StreamF",          "sdcor"]
multiple.stream.indiv.sdSiteYear  <- vc[ vc$grp=="StreamF:YearF",     "sdcor"]

# Note that d=effect size = difference to detect / standard deviation
scenarios <- expand.grid(n_after=5:20, 
                         n_before =5,
                         sdYear = multiple.stream.indiv.sdYear,
                         sdSite = multiple.stream.indiv.sdSite ,
                         sdSiteYear=multiple.stream.indiv.sdSiteYear,
                         sdResid= multiple.stream.indiv.sdResid,
                         n_stream = 4,
                         n_quad  = 5,
                         alpha = .05,
                         delta = .2)

power <- plyr::adply(scenarios,1,function(x){
  # compute power for one row of the table
  stddev <- sqrt(  x$sdYear^2 + x$sdSiteYear^2/x$n_stream +
                   x$sdResid^2 /x$n_stream/x$n_quad)
  power <- pwr.t2n.test(n1=x$n_before, 
                        n2=x$n_after, 
                        d =x$delta/stddev, 
                        alternative="two.sided")$power
  data.frame(power=power)
})

cat("Power for detecting ", power$delta[1],", proportional change in 4 streams with 5 years before and \n")
cat("different number of years after; 5 quadrats/year\n")
cat("\nNote that power is now for a larger scope of inference \n")
power
##***part800e;
# sink()



