################################# MASTER STATS CODE ###########################################
#Much of this code is based on Dr. Teresa Treat's code from Mixed Effects Modeling in Psychology
#Which is an outstanding course that everyone interested in doing good statistics should take.
#
# The general structures of real-world statistical anlyses performed by the Soundscape & Audiology Research Lab
# are based on code contained in this compilation. The code here is well-suited for statistical analyses of real-world data
# which typically involves repeated sampling. For details, see Oleson et al., 2022 ("Statistical Considertions for 
# Analyzing Ecological Momentary Assessment Data").
#
#Compiled and expanded by Erik Jorgensen, AuD, PhD
#UW-Madison 2024
#Original version: University of Iowa 2019

################################# CONTENTS ##################################################
#1. Read Data
#2. Manipulate Data
#3. Explore Data
#4. Multiple Linear Regression
#5. Mixed Effects with Linear Slopes
#6. Mixed Effects with Quadratic Slopes
#7. Logistic Regression
#8. Ordered Logit with Mixed Effects
#9. Non-Linear Mixed Effects
#10. Model Evaluation
#11. Power Analyses
###############################################################################################

#clear existing variables
rm(list=ls(all=TRUE))

##Update packages
update.packages()

##Load packages
library(nlme)    #Non-linear mixed effects package
library(car)   #Regression package
library(lme4)    #Linear mixed effects package
library(ggplot2)   #Advanced plotting package
library(languageR)  ##Useful functions developed originally for language research
library(emmeans) #pairwise contrasts and marginal means
library(lattice) #plot variable relationships
library(lmerTest) #get p-values for mixed models
library(gplots) #more plotting functions
library(FSA) 
library(MASS)
library(ordinal) #ordered logit
library(pwr) #basic power analyses
library(simr) #power analyses for mixed effects using simulation
library(r.squaredGLMM) #r-squared for mixed models
library(multcomp) #post-hoc tests for LMER
library(tidyverse)

#get data
data = read.table("DATA_FILE", header = TRUE, sep = ",")
data = data 

#omit missing
data <- na.omit(data) 

#Attach data file and view list of variable names.
attach(data)
names(data)

##Check beginning of data.
head(data, n = 10) #view first 10 lines of file

############################### MANIPULATE DATA ################################

#set factors
data$COLUMN_FACTOR <- as.factor(data$COLUMN_FACTOR)
#check factor
is.factor(data$COLUMN_FACTOR)
#alternative way to check factor
class(data$COLUMN_FACTOR)	#another way to check factors

#change default coding from dummy to effect
options(contrasts = c("contr.sum","contr.poly"))

#DUMMY: difference from reference group
#EFFECT: difference from grand mean

#get levels
levels(data$COLUMN_FACTOR)

#Change lables and Levels for factors - DUMMY CODING
data$COLUMN_FACTOR <- factor(data$COLUMN_FACTOR, levels = c(0,1), labels = c("Label 0","Label 1"))

#Change lables and level for factors - EFFECT CODING
data$COLUMN_FACTOR <- factor(data$COLUMN_FACTOR, levels = c(-1,1), labels = c("Label -1","Label 1"))

############################### EXPLORE DATA ####################################

column_mean<-mean(data$column)#get mean

##Descriptives for data set
summary(data)

#Provides correlation of all variables in data set
cor(data)

#Provides correlations for subset of variables in data set.  
predictors <- data[,c("column1", "column2")]
cor(predictors)

#Provides single parametric bivariate correlation 
cor.test(column1, column2)

#Centering variables (around mean)
data$column_c <- data$column - mean(data$column)
summary(data)

#get per subject means for repeated measures
aggregate(x = DATA$COLUMN,    # Specify data column
          by = list(DATA$Subject), # Specify group indicator
          FUN = mean, na.rm=TRUE)  # Specify function (i.e. mean)

############################### MULTIPLE LINEAR REGRESSION ####################################

#FULL DATA SET WITH CENTERED PREDICTORS

##Computes Type III SS, which is used by convention in psychology
model1 <- lm(DEPENENT_VARIABLE ~ PREDICTOR1_CENTERED + PREDICTOR2_CENTERED, data = data)
summary(model1)  #Presents coefficients, standard errors, t values, and p-values
summary.lm(model1)  #Provides linear model output -- same as summary(model1) output
Anova(lm(model1), type="III")  #Provides ANOVA-like output, assuming Type III SS
coef(model1)   #Lists model coefficients/parameters

## "summary.aov(model1)" also provides ANOVA-like output, but assumes Type I SS (rather than Type III SS)

confint(model1, level = .95)  #Provides 95% CIs for model parameters
fitted(model1)   #Provides predicted values for each case
residuals(model1) #Provides residuals for each case
vcov(model1)   #Provides variance-covariance matrix for model parameters

#SUBSET SUBJECT 
data_subX <- subset(data, subject == SUBJECTNUMBERX)
summary(data_subX)

##Linear multiple regression of centered predictors on data for subject X
##Note that you must re-run subsetting code for creation of data_subX before
##you can estimate model2, b/c you computed SI_c AFTER you subsetted the data
model2 <- lm(DEPENENT_VARIABLE ~ PREDICTOR1_CENTERED + PREDICTOR2_CENTERED, data = data_subX)
summary.lm(model2)  
Anova(lm(model2), type="III")  #Provides ANOVA-like output, assuming Type III SS 
coef(model2)

##Check for highly influential observations.
##Command below provides wealth of influence indices: DFBetas, DFFits, Cook's distances, 
##leverage values (hat) etc.
##Cases which are influential with respect to any of these measures are marked with an asterisk. 
influence.measures(model2)  
##Command below indicates whether influence indices clear established cutoffs.
influence.measures(model2)$is.inf  

#FULL MODEL WITH UNCENTERED PREDICTORS
model3 <- lm(DEPENDENT_VARIABLE ~ PREDICTOR1 + PREDICTOR2 + PREDICTOR1:PREDICTOR2, data = data)
summary(model3)  #Presents coefficients, standard errors, t values, and p-values
summary.lm(model3)  #Provides linear model output -- same as summary(model1) output
Anova(lm(model3), type="III")  #Provides ANOVA-like output, assuming Type III SS
coef(model3)   #Lists model coefficients/parameters
vif(model3)

#FULL MODEL WITH CENTERED PREDICTORS
model4 <- lm(DEPENDENT_VARIABLE ~ PREDICTOR1_CENTERED + PREDICTOR2_CENTERED + PREDICTOR1_CENTERED:PREDICTOR2_CENTERED, data = data)
summary(model4)  #Presents coefficients, standard errors, t values, and p-values
summary.lm(model4)  #Provides linear model output -- same as summary(model1) output
Anova(lm(model4), type="III")  #Provides ANOVA-like output, assuming Type III SS
coef(model4)   #Lists model coefficients/parameters
vif(model4)

#REGRESS ALL SUBJECTS
model5 <- lmList(DEPENDENT_VARIABLE ~ PREDICTOR1_CENTERED + PREDICTOR2_CENTERED | subject, data = data)
summary(model5)
coef(model5)

#Get Correlations and Effect Sizes
summary(model1)$r.squared
effect.size(model1)

############################### MIXED EFFECTS WITH ONLY LINEAR SLOPES ###############################################

#Add effects and predictors as needed

#Test random effects
model6 <- lmer(DEPENENT_VARIABLE ~ (1|SUBJECT), data = data) #Variance component analysis: random subject
#calculate Interclass Correlation Coefficients (ICC) from this:
#variance of effect of interest divided by total variance

#Can add random effect for item
model6 <- lmer(DEPENENT_VARIABLE ~ (1+ITEM|SUBJECT), data = data) #Variance component analysis: random subject and item

#add fixed effects
model7 <- lmer(DEPENENT_VARIABLE ~ PREDICTOR1 * PREDICTOR2 + (1|SUBJECT), data = data) #full factorial for 2 predictors, random intercept for subject
summary(model7)  #Presents coefficients, standard errors, t values, and p-values

#add random effect for item
model8 <- lmer(DEPENENT_VARIABLE ~ PREDICTOR1 * PREDICTOR2 + (1+ITEM|SUBJECT), data = data) #full factorial, random subject and item
summary(model8)  #Presents coefficients, standard errors, t values, and p-values

#reduce model appropriately
model9 <- lmer(DEPENENT_VARIABLE ~ PREDICTOR1 + PREDICTOR2 + (1|SUBJECT), data = data) #no interaction term for 2 predictors, random intercept for subject
summary(model9)  #Presents coefficients, standard errors, t values, and p-values

aov(model8, model9) #compare models: if p<0.05, more complex model is needed

#calculate both maringal (fixed) and conditional (fixed + random) r-squared
rsqu1 <- r.squaredGLMM(model1)

#LOOK AT THE DATA
#Spaghetti
g1 <- ggplot(data=data, aes(x=TIME, y=OUTCOME, group = SUBJECT)) + geom_line() + theme_bw()
g2 <- g1 + theme_bw()+scale_x_continuous(name = "NAME_OF_X", breaks = c(0,1,2,3,4,5,6))
g3 <- g2 + scale_y_continuous(name = "NAME_OF_Y")
print(g3)

#Individual Growth Curves - Linear
g1 <- ggplot(data=data, aes(x=PREDICTOR, y=OUTCOME, group = SUBJECT))
g2 <- g1 + geom_point() + stat_smooth(method = "lm", se = FALSE) + facet_wrap(~SUBJECT)
g3 <- g2 + theme_bw()+scale_x_continuous(name = "NAME_OF_X", breaks = c(0,1,2,3,4,5,6))
g4 <- g3 + scale_y_continuous(name = "NAME_OF_Y")
print(g4)

############################### MIXED EFFECTS WITH QUADRATIC SLOPES ###############################################

#Full factorial with quadratic fixed and random effects
model1 <- lmer(OUTCOME ~ time*group + I(time^2)*group + (1 + time + I(time^2)|subject), data = data)
summary(model1)

#Eliminate random quadratic slope effect and compare
model2 <- lmer(OUTCOME ~ time*group + I(time^2)*group + (1 + time|subject), data = data)
summary(model2)
anova(model1, model2)

#Repeat etc.

#LOOK AT DATA
#Individual Growth Curves - Quadratic
g1 <- ggplot(data=data, aes(x=TIME, y=OUTCOME, group = SUBJECT))
g2 <- g1 + geom_point() + stat_smooth(method = "lm", formula = y~poly(x,2), se = FALSE) + facet_wrap(~SUBJECT)
g3 <- g2 + theme_bw()+scale_x_continuous(name = "NAME_OF_X", breaks = c(0,1,2,3,4,5,6))
g4 <- g3 + scale_y_continuous(name = "NAME_OF_Y")
print(g4)


############################ LOGISTIC REGRESSION ######################################

#can specify family = binomal for yes/no or family = poisson for count data
#estimates are in log-odds; must convert to probabilities if you want probabilities!

#model factorial fixed effects with random subject intercept
model1 <- glmer(DEPENDENT_VARIABLE ~ PREDICTOR1 * PREDICTOR2 + (1|SUBJECT), data=data, family = binomial, verbose = TRUE)
summary(model1) #get model output

#model factorial fixed effects with random subject intercept and random item slop
model2 <- glmer(DEPENDENT_VARIABLE ~ PREDICTOR1 * PREDICTOR2 + (1+ITEM|SUBJECT), data=data, family = binomial, verbose = TRUE)
summary(model2) #get model output

exp(fixef(model1)) #compute odds ratios for log odds fixed effects 
plogis(fixef(model1)) #compute inverse logit (probability) from log odds fixed effects

confint(model1,method="Wald") #get confidence intervals in log odds
1/(1+exp(-1*INTEREPT1)) #get CI as probability for intercept

#Caterpillar plot
dotplot(ranef(model1, which = "Participant", condVar = TRUE))

#VIF script for use in GLMER
vif.mer(model1) #check multicollinearity
vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}
#Run vif.mer function 
vif.mer(model1)

#plot predicted values
predict(model1, type="response") # predicted values
hist(predict(model1, type="link")) #log odds
hist(predict(model1, type="response")) #probabilities

#plot residuals
residuals(model1, type="deviance") 
hist(residuals(model1, type="deviance")) # residuals in log odds
hist(residuals(model1, type="response")) # residuals in probabilities

coef(model1) #get individual coefficients
ranef(model1) #get random effects

qqmath(ranef(model1, postVar = TRUE), strip = FALSE)$SUBJECT #qqplot

#transform outcome variable to logit
#NEVER USE RAU
Logit_data <- logit(data)

#VERY COOL PLOT of individual sigmoid functions for subjects
g1 <- ggplot(data = data_plot, aes(x = INDEPENDENT_VARIABLE, y = DEPENDENT_VARIABLE))
g2 <- g1 + geom_point(alpha=.5) + geom_smooth(method = "glm", method.args = list(family = "binomial"))
g3 <- g2 + scale_x_continuous(name = "X AXIS TITLE")
g4 <- g3 + scale_y_discrete(name = "Y AXIS TITLE", limits=c(0,1))
g5 <- g4 + coord_cartesian(ylim =c(-.1,1.1), expand = FALSE) + facet_wrap(~subject)
print(g5)

################################# ORDERED LOGISTIC MIXED EFFECT REGRESSION USING LINKED MODELS #########################

#Ordinal logistic regression using cumulative link with random effect for subject
#Must specify outcome as factor
model1 <- clmm2(OUTCOME ~ PREDICTOR, random=SBUJECT, data=data, Hess=TRUE, nAGQ=10)
summary(model1)

#Plot histograms by proportion
p <- ggplot(data) + stat_count(mapping = aes(x=nz, y=..prop.., group=1))
p + facet_wrap(~condition)

################################ NON LINEAR REGRESSION ###########################################

############################ MAKING GROUPED DATA OBJECT #################################

grouped_data <- groupedData(OUTCOME ~ TIME|Subject, outer = ~ Group, data=data)
labels = list(x = "TIME", y = "OUTCOME")
summary(grouped_data)
plot(grouped_data)
plot(grouped_data, outer = ~ Group)

############################ NLS #################################

#must change a, b, and c start values depending on what your data is!

#a = asymptotic outcome
#b = difference between initial and asymptote
#estimate initial by a + b
#c = rate
#group effect is difference in b - estimate differences in change (difference between initial and asymptote) by 
#subtracting or adding group effect to b

#3 paramter model
model1 <- nls(OUTCOME ~ a-(b*exp(-c*TIME)), data=data,
              start = c(a = 1, b = .3, c = 2), trace=TRUE)
summary(model1)

#2 paramter model: 1, 0.3, 2
model2 <- nls(OUTCOME ~ a-(a*exp(-c*TIME)), data=data,
              start = c(a = 1, c = 2), trace=TRUE)
summary(model2)

anova(model2, model1)

#subset data
group1 <- subset(grouped_data, Group == "group1")
group2 <- subset(grouped_data, Group == "group2")

#nls model for Group subset
model3 <- nls(OUTCOME ~  a-(b*exp(-c*TIME)), data = group1,
              start = c(a = 1, b = .3, c = 2), trace=TRUE)
summary(model3)

#nls model for control
model4 <- nls(OUTCOME ~  a-(b*exp(-c*TIME)), data = group2,
              start = c(a = 1, b = .3, c = 2), trace=TRUE)
summary(model4)

####################################### NLSLIST ##########################################

model5 <- nlsList(OUTCOME ~  a-(b*exp(-c*TIME)), 
                  start = c(a = 1, b = .3, c = 2), data = grouped_data)
summary(model5)
coef(model5)
plot(augPred(model5))

############################ NLME #################################

#full NLME model
model6 <- nlme(OUTCOME ~  a-(b*exp(-c*TIME)), fixed = a+b+c~1, random=pdDiag(a+b+c~1),
               start = c(a = 1, b = .3, c = 2), data = grouped_data, control=list(maxIter=4000, pnlsMaxIter = 100), verbose = TRUE)
summary(model6)
coef(model6)
plot(augPred(model6))

#b,c random effects
model7 <- nlme(OUTCOME ~  a-(b*exp(-c*TIME)), fixed = a+b+c~1, random=pdDiag(b+c~1),
               start = c(a = 1, b = .3, c = 2), data = grouped_data, control=list(maxIter=4000, pnlsMaxIter = 100), verbose = TRUE)
summary(model7)
coef(model7)
plot(augPred(model7))

anova(model7, model6)

#group effect on c, all random effects
model8 <- nlme(OUTCOME ~  a-(b*exp(-c*TIME)), fixed = list(a~1,b~1,c~Group), random=pdDiag(a+b+c~1),
               start = c(a = 1, b = .3, c = 2, rep(.5,1)), data = grouped_data, control=list(maxIter=4000, pnlsMaxIter = 100), verbose = TRUE)
summary(model8)
coef(model8)
plot(augPred(model8))

#graphing results
g1 <- ggplot(data=grouped_data, aes(x=TIME, y = OUTCOME, shape = Group))
g2 <- g1 + stat_summary(fun.data = mean_se, geom="pointrange")
g3 <- g2 + stat_summary(aes(y=fitted(model8), linetype = Group), fun.y=mean, geom="line")
g4 <- g3 + theme_bw() + scale_x_continuous(name = "NAME_OF_X")
g5 <- g4 +scale_y_continuous(name = "OUTCOME")
print(g5)

############################### MODEL EVALUATION/CHECKING ####################################

#DIAGNOSTICS FOR MODEL 1 (FULL MODEL WITH CENTERED PREDICTORS)

##Obtain four diagnostic plots by using plot(model name) command. 
##1) fitted values by unstandardized residuals.  Look for curvature that suggests violation of 
##linearity btw Xs and Y.  You ideally want no pattern at all (white noise).
##2) normal Q-Q (quantile-quantile) plot.  A straight line indicates normality of residuals.
##3) fitted values by square root of abs(standardized residuals).  This provides a positive 
##version of the first graph.  It commonly is used to look for heteroscedasticity. 
##4) leverage by standardized residuals (also lists Cook's distance contours).  Here, you are 
##looking for high-leverage points and points close to or beyond Cook's distance contours.
##Syntax below prints sequence of 4 plots.
plot(model1)
##Better to use syntax below that prints 4 graphs on single page.
layout(matrix(c(1,2,3,4),2,2)) 
plot(model1)

##Conduct multicollinearity check by obtaining Variance Inflation Factor (vif) for each predictor in model.
##VIF for a specific predictor X is (1 / 1-r^2) in context of other predictors.  If r^2 close to 1, then X
##well-predicted by other predictors, and you get a huge VIF (b/c you have a big multicollinearity problem!).
##Typically concerned if single VIF > 10, or average VIF >> (much greater than) 1.
vif(model1)

##Check linearity assumption by examining scatterplot of fitted estimates by model residuals.
plot(fitted(model1),residuals(model1))   #unstandardized residuals
plot(fitted(model1),rstandard(model1))   #standardized residuals

##Check homoscedasticity
plot(fitted(model1), residuals(model1))
plot(fitted(model1), rstandard(model1))

##Take a look at list of all residuals, if you like.
residuals(model1)   #unstandardized residuals
rstandard(model1)   #standardized residuals

##Examining normality of residuals.
hist(residuals(model1))    #unstandardized residuals
hist(rstandard(model1))    #standardized residuals
##Note that you could write these to a tab-delimited text file to figure out problematic cases.
write.table(rstandard(model1), "C:\\Users\\Teresa\\Desktop\\tempfile3.txt", sep="\t") 
##Bonferroni test of standardized residuals
outlierTest(model1)
##Q-Q Plot of residuals
qqnorm(residuals(model1))

##Check for highly influential observations.
##Command below provides wealth of influence indices: DFBetas, DFFits, Cook's distances, 
##leverage values (hat) etc.
##Cases which are influential with respect to any of these measures are marked with an asterisk. 
influence.measures(model1)  
##Command below indicates whether influence indices clear established cutoffs.
influence.measures(model1)$is.inf  

###OTHER

transform_column <- data$column + 1 #change values to column

##Transform variables as needed to deal with skew, outliers, nonlinearities, etc.
data$column_ln <- log(transform_sp)  #natural log
##Other useful functions:  sqrt(variable); abs(variable); exp(variable); asin(variable); 1/variable; variable^3
data$column_sqrt <- sqrt(transform_sp)  #square root
data$column_inverse <- 1/(transform_sp)  #inverse
hist(data$column_ln)
hist(data$column_sqrt)
hist(data$column_inverse)
hist(transform_sp)

mean(data$column_ln)
mean(data$column_sqrt)
mean(data$column_inverse)

####################################### POWER ######################################

## Elements for pwr.t.test##
# n: Number of observation 
# d: Effect size (Cohen's d)
# sig.level: Significane level(alpha) 
# power: Power of test (1-beta)
# type: 'one.sample',"two.sample", or "paired"


#Calculate sample size 
pwr.t.test(n=, d=.60, sig.level = .05, power=.8, type="two.sample") 
pwr.t.test(n=, d=.60, sig.level = .05, power=.9, type="two.sample") 

# simulation and plot power curves
# create user-defined function to generate and analyze data
# this assume that treatment 1 and treatment 2 has the same effect when compared to control
# the sample size N is the total sample size
f_func <- function(simNum, N, d) {  # three groups
  x1 <- rnorm(floor(N/3), 0, 1) # control group
  trt <-   rep(0,  length(x1))
  dat_c <- cbind(x1, trt)
  x2 <- rnorm(floor(N/3), d, 1) # treatment one 
  trt1 <-   rep(1, length(x2))
  dat_t1 <-  cbind(x2, trt1)
  x3 <- rnorm(floor(N/3),d, 1) # treatment two
  trt2 <-   rep(2, length(x3))
  dat_t2 <-  cbind(x3, trt2)
  dat <- data.frame(rbind(dat_c, dat_t1, dat_t2))
  names(dat) <- c('outcome', 'treatment')
  # fit the anova model and extract needed stats
  fit <- aov(outcome ~ treatment, data = dat)
  summary(fit)
  stat <-  summary(fit)[[1]][1,4]  # extract F value
  p <- summary(fit)[[1]][1,5]  # extract p value
  #   summary(fit)[[1]][2,3]  # extract MSE
  return(c(F=stat, p=p, sig=(p < .05)))
} # end of the function below

# varying N and Cohen's d
power_ftest_vary3 <- grid_search(f_func, 
                                 params=list(N=c(25, 50, 100, 250),  # Total sample sizes 
                                             d=c(.58, .68, .77)),        # effect size i.e., d 
                                 n.iter=10000,       # number of replications use about 5000
                                 output='data.frame')   
## Running 120,000 tests...
# table N.test is for each group 
power <- results(power_ftest_vary3) %>%       # we create a summary 
  group_by(N.test, d.test) %>%
  summarise(power=mean(sig))
print(power)
load(ggplot)
ggplot(power, aes(x=N.test, y=power, group=factor(d.test), colour=factor(d.test))) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=.8, linetype="dashed", color = "grey") +
  ylim(c(0, 1)) +
  labs(x='Sample Size', y='Power', colour="Cohen's d") +
  theme_minimal()

#Monte Carlo Simulation for Mixed Effects Power Analysis
#create linear mixed effects model from pilot data
model1 = lmer(OUTCOME ~ predictor1 + predictor2  + (1|subject), data=data)
summary(model1)

#Conduct power analysis for specified fixed effect
#can change simulation number using nsim=
#default for fixed effect is first fixed effect in the model (probably the intercept)
powermodel1 <- powerSim(model1, test=fixed("predictor1"), sim=model1)
print(powermodel1)

#get curves for power levels and sample sizes
p_curve_treat <- powerCurve(model1, test=fcompare(insitu~time), within="subject+predictor1+predictor2", breaks=c(5,10,15,20))
plot(p_curve_treat)

#extend sample size (typically highest level; e.g. subject)
model1_ext_n <- extend(model1, along="subject", n=100)
model1_ext_n

#extend measures number (typically embedded level; e.g. surveys per subject)
model1_ext_measures <- extend(model1_ext_n, within="subject+predictor1+predictor2", n=100)
model1_ext_measures

#power analysis for model with extended subject and extended measures
powermodel_ext_measures <- powerSim(model1_ext_measures, test=fixed("predictor1"),nsim=10)
print(powermodel_ext_measures)

