rm(list = ls())

# Load packages ----------------------------------------------------------------

library(optconerrf)
library(ggplot2)
library(cowplot)

# Functions --------------------------------------------------------------------

### Plot of second stage sample size -------------------------------------------

#Function to create a plot of the second stage sample size for a two sample
#Z-test with variance 1
plot_sample_size <- function(x, range, plotNonMonotoneFunction = TRUE){

  #Create a new non-monotone design
  nonmonotoneDesign <- suppressWarnings(new(
    "TrialDesignOptimalConditionalError",
    alpha = x$alpha,
    alpha1 = x$alpha1,
    alpha0 = x$alpha0,
    conditionalPower = x$conditionalPower,
    conditionalPowerFunction = x$conditionalPowerFunction,
    delta1 = x$delta1,
    firstStageInformation = x$firstStageInformation,
    useInterimEstimate = x$useInterimEstimate,
    likelihoodRatioDistribution = x$likelihoodRatioDistribution,
    deltaLR = x$deltaLR,
    weightsDeltaLR = x$weightsDeltaLR,
    tauLR = x$tauLR,
    kappaLR = x$kappaLR,
    deltaMaxLR = x$deltaMaxLR,
    minimumConditionalError = x$minimumConditionalError,
    maximumConditionalError = x$maximumConditionalError,
    minimumSecondStageInformation = x$minimumSecondStageInformation,
    maximumSecondStageInformation = x$maximumSecondStageInformation,
    levelConstantMinimum = x$levelConstantMinimum,
    levelConstantMaximum = x$levelConstantMaximum,
    ncp1Min = x$ncp1Min,
    ncp1Max = x$ncp1Max,
    enforceMonotonicity = FALSE
  ))

  #Create a sequence of first stage p-values for plotting
  firstStagePValues <- seq(from = range[1], to = range[2], length.out = 1e3)

  #Calculate the second stage sample size for the monotone design
  secondStageSampleSize <- 2 * getSecondStageInformation(
    firstStagePValue = firstStagePValues,
    design = x)

  #Calculate the second stage sample size for the non-monotone design
  nonMonoSecondStageSampleSize <- 2 * getSecondStageInformation(
    firstStagePValue = firstStagePValues,
    design = nonmonotoneDesign)

  #Create plot
  sampleSizePlot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = firstStagePValues,
        y = secondStageSampleSize
      ),
      colour = "black",
      linetype = "solid",
      linewidth = 1.1
    ) +
    ggplot2::labs(x = "First-stage p-value",
                  y = "Second-stage sample size per group") +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
      xintercept = x$alpha0,
      linetype = "dotted",
      col = "red",
      linewidth = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = x$alpha1,
      linetype = "dotted",
      col = "blue",
      linewidth = 0.8
    ) +
    ggplot2::xlim(c(range[1], range[2])
    )

  if(plotNonMonotoneFunction == TRUE){
    sampleSizePlot <- sampleSizePlot + ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = firstStagePValues,
        y = nonMonoSecondStageSampleSize
      ),
      colour = "darkgrey",
      linetype = "dashed",
      linewidth = 1.1
    )
  }

  return(sampleSizePlot)
}

### Implementation of the inverse normal method --------------------------------

#Function to get second-stage level for inverse normal test
findAlpha2 <- function(alpha2, alpha, alpha1, alpha0, weights) {
  correlationMatrix <- matrix(data = c(1, weights[1], weights[1], 1), nrow = 2)
  return(-mvtnorm::pmvnorm(upper = c(qnorm(1-alpha2),stats::qnorm(1-alpha1)),
                           corr = correlationMatrix) +
           mvtnorm::pmvnorm(upper = c(qnorm(1-alpha2), stats::qnorm(1-alpha0)),
                            corr = correlationMatrix) + alpha0 - alpha)
}

#Get Second-stage Level for Inverse Normal Method
getAlpha2InverseNormal <- function(alpha, alpha1, alpha0,
                                   weights = c(sqrt(0.5),sqrt(0.5))) {
  # Solve the helper function findAlpha2
  alpha2 <- stats::uniroot(f = findAlpha2,
                           lower = 0,
                           upper = 1,
                           alpha = alpha,
                           alpha1 = alpha1,
                           alpha0 = alpha0,
                           weights = weights)$root
  return(alpha2)
}

# Function to calculate Inverse Normal Conditional Error
getInverseNormalConditionalError <- function(alpha, alpha0, alpha1,
                                             firstStagePValue,
                                             weights = c(sqrt(0.5),sqrt(0.5)),
                                             alpha2 = NULL) {

  if(is.null(alpha2)) {
    alpha2 <- getAlpha2InverseNormal(alpha = alpha,
                                     alpha1 = alpha1,
                                     alpha0 = alpha0,
                                     weights = weights)
  }

  conditionalError <- NA
  if(firstStagePValue <= alpha1) {
    conditionalError <- 1
  }
  else if(firstStagePValue > alpha0) {
    conditionalError <- 0
  }
  else if(firstStagePValue <= alpha0 && firstStagePValue > alpha1){
    conditionalError <- 1-pnorm((qnorm(1-alpha2)-weights[1]*
                                   qnorm(1-firstStagePValue))/weights[2])
  }
  return(conditionalError)
}

getInverseNormalConditionalError <- Vectorize(FUN = getInverseNormalConditionalError,
                                              vectorize.args = c("firstStagePValue"))

#Function to calculate the expected second stage information for the
#inverse normal method, if the interim estimate is used for cp calculation

getSecondStageInformationIN <- function(p1, al, al0, al1, cp, I1, delta_Min,
                                        weights = c(sqrt(0.5),sqrt(0.5))) {

  effect <-
    pmax(qnorm(1 - p1) / sqrt(I1),
        delta_Min)

  # Calculate conditional error
  conditionalError <- getInverseNormalConditionalError(alpha = al,
                                                       alpha0 = al0,
                                                       alpha1 = al1,
                                                       firstStagePValue = p1,
                                                       weights = weights)

  information <- (optconerrf::getNu(
    alpha = conditionalError,
    conditionalPower = cp
  )) /
    (effect^2)
  return(information)
}

getExpectedSecondStageInformationIN <- function(al, al1, al0, cp, I,
                                                delta_Min, delta,
                                                weights = c(sqrt(0.5),sqrt(0.5))){

  a2 <- getAlpha2InverseNormal(alpha = al, alpha1 = al1, alpha0 = al0)

  integrate <- function(p1){
    optimalCondErr <- getInverseNormalConditionalError(alpha = al,
                                                       alpha1 = al1,
                                                       alpha0 = al0,
                                                       alpha2 = a2,
                                                       firstStagePValue = p1,
                                                       weights = weights)

    Nu <- optconerrf::getNu(alpha=optimalCondErr, conditionalPower = cp)
    del_Min <- pmax(delta_Min, qnorm(1-p1)/sqrt(I))
    res <- Nu*exp(qnorm(1-p1)*sqrt(I)*delta-I*delta^2/2)/del_Min^2
    return(res)
  }
  return(stats::integrate(
    f = integrate,
    lower = alpha1,
    upper = alpha0)$value)
}

#Calculate the overall power for the inverse normal method
getOverallPowerIN <- function(al, al1, al0, cp, I,
                              delta_Min, delta,
                              weights = c(sqrt(0.5),sqrt(0.5)),
                              alternative){
  ncp1 <- alternative * sqrt(I)

  firstStageFutility <- numeric(length(alternative))
  firstStageEfficacy <- numeric(length(alternative))
  overallPower <- numeric(length(alternative))

  for (i in 1:length(alternative)) {
    # Early decision probabilities
    firstStageFutility[i] <- stats::pnorm(
      stats::qnorm(1 - al0) - ncp1[i]
    )
    firstStageEfficacy[i] <- 1 -
      stats::pnorm(stats::qnorm(1 - al1) - ncp1[i])

    # Calculate probability to reject at the second stage for given delta
    secondStageRejection <- function(firstStagePValue) {
      (1 -
         stats::pnorm(
           stats::qnorm(
             1 - getInverseNormalConditionalError(firstStagePValue = firstStagePValue,
                                                  alpha = al,
                                                  alpha0 = al0,
                                                  alpha1 = al1,
                                                  weights = weights)
           ) -
             sqrt(getSecondStageInformationIN(p1 = firstStagePValue,
                                              al = al,
                                              al0 = al0,
                                              al1 = al1,
                                              cp = cp,
                                              I1 = I,
                                              delta_Min = delta_Min)) *
             alternative[i]
         )) *
        exp(qnorm(1 - firstStagePValue) * ncp1[i] - ncp1[i]^2 / 2)
    }

    integral <- stats::integrate(
      f = secondStageRejection,
      lower = al1,
      upper = al0
    )$value

    overallPower[i] <- firstStageEfficacy[i] + integral
  }

  powerResults <- new(
    "PowerResultsOptimalConditionalError",
    alternative = alternative,
    firstStageFutility = firstStageFutility,
    firstStageEfficacy = firstStageEfficacy,
    overallPower = overallPower
  )

  return(powerResults)
}

### Find first stage sample size for overall power control ---------------------

#Define functions that can be used to find the first stage sample size that
#leads to an overall power of 0.8

#base design
overallPower_base <- function(n_1){

  I_1 = 0.5*n_1
  design_base <- optconerrf::getDesignOptimalConditionalErrorFunction(
    alpha = alpha,
    alpha1 = alpha1,
    alpha0 = alpha0,
    conditionalPower = conditionalPower,
    firstStageInformation = I_1,
    useInterimEstimate = TRUE,
    likelihoodRatioDistribution = "fixed",
    deltaLR = deltaLR,
    delta1Min = delta1Min)

  Power <- getOverallPower(design_base, alternative = deltaLR)

  return(Power$overallPower-conditionalPower)
}

#maxlr design
overallPower_maxlr <- function(n_1){

  I_1 = 0.5*n_1
  design_maxlr <- optconerrf::getDesignOptimalConditionalErrorFunction(
    alpha = alpha,
    alpha1 = alpha1,
    alpha0 = alpha0,
    conditionalPower = conditionalPower,
    firstStageInformation = I_1,
    useInterimEstimate = TRUE,
    likelihoodRatioDistribution = "maxlr",
    delta1Min = delta1Min)

  Power <- getOverallPower(design_maxlr, alternative = deltaLR)

  return(Power$overallPower-conditionalPower)
}

#inv. normal equally weighted
overallPower_IN1 <- function(n_1){

  I_1 = 0.5*n_1

  Power <- getOverallPowerIN(al = alpha,
                             al1 = alpha1,
                             al0 = alpha0,
                             cp = conditionalPower,
                             I = I_1,
                             delta_Min = delta1Min,
                             weights = c(sqrt(0.5),sqrt(0.5)),
                             alternative = deltaLR)

  return(Power$overallPower-conditionalPower)
}

#inv. normal different weights
overallPower_IN2 <- function(n_1){

  I_1 = 0.5*n_1

  Power <- getOverallPowerIN(al = alpha,
                             al1 = alpha1,
                             al0 = alpha0,
                             cp = conditionalPower,
                             I = I_1,
                             delta_Min = delta1Min,
                             weights = c(sqrt(1/3),sqrt(2/3)),
                             alternative = deltaLR)

  return(Power$overallPower-conditionalPower)
}

# Results ----------------------------------------------------------------------

### Type I error rate ----------------------------------------------------------
alpha <- 0.05
firstStageInformation_maxlr_non_mono <- 10

design_maxlr_non_mono <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = 0,
  alpha0 = 1,
  conditionalPower = 0.8,
  firstStageInformation = firstStageInformation_maxlr_non_mono,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "maxlr",
  delta1Min = 0.125,
  enforceMonotonicity = FALSE)

tab_type_I_error <- data.frame()

for (alternative in c(0*sqrt(firstStageInformation_maxlr_non_mono),-0.25*sqrt(firstStageInformation_maxlr_non_mono), -0.5*sqrt(firstStageInformation_maxlr_non_mono))){
density <- function(p1){
  exp(qnorm(1-p1)*alternative - alternative^2/2)
}

intergrand_ocef <- function(p1){
  density(p1)*getOptimalConditionalError(design = design_maxlr_non_mono, p1)
}

intergrand_invn1 <- function(p1){
  density(p1)*getInverseNormalConditionalError(alpha = alpha,
                                               alpha1 = 0,
                                               alpha0 = 1,
                                               firstStagePValue = p1,
                                               weights = c(sqrt(0.5),sqrt(0.5)))
}

intergrand_invn2 <- function(p1){
  density(p1)*getInverseNormalConditionalError(alpha = alpha,
                                               alpha1 = 0,
                                               alpha0 = 1,
                                               firstStagePValue = p1,
                                               weights = c(sqrt(1/3),sqrt(2/3)))
}

res_ocef <- integrate(f = intergrand_ocef, lower = 0, upper = 1)$value
res_invn1 <- integrate(f = intergrand_invn1, lower = 0, upper = 1)$value
res_invn2 <- integrate(f = intergrand_invn2, lower = 0, upper = 1)$value

tab_type_I_error <- rbind(tab_type_I_error, c(res_invn1, res_invn2, res_ocef))
}

tab_type_I_error <- cbind(c(0, -0.25, -0.5), tab_type_I_error)

colnames(tab_type_I_error) <- c("alternative", "inv. normal 1", "inv. normal 2", "ocef")

round(tab_type_I_error,3)

### Design parameters ----------------------------------------------------------

#Specify globally parameters of the design
alpha <- 0.05 #Overall type I error
overallPower <- 0.8 #Over all power
deltaLR <- 0.2 #Assumed effect/mean difference under which to calculate the likelihood ratio

#Sample size needed per group in a fixed-size-sample test
n <- ceiling(2*(qnorm(overallPower)- qnorm(alpha))^2/deltaLR^2);n

#Take a third of the sample size at the first stage
n_third <- ceiling(1/3 * n); n_third

#Convert to information (Two sample Z-test with a sample size of 104 per group)
firstStageInformation <- n_third * 0.5

alpha1 <- 0.001 #Early rejection boundary
alpha0 <- 0.5 #Futility boundary
conditionalPower <- 0.8 #conditional power
delta1Min <- 0.125 #delta_0

### design_base ----------------------------------------------------------------

#Create design for the first example
design_base <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = deltaLR,
  delta1Min = delta1Min)

### Level constant and monotonisation constants --------------------------------

#Overview over the design parameters including level constant and monotonisation
#constants
design_base

### Figure 1 -------------------------------------------------------------------

#Plot of the optimal conditional error function and the optimal monotone
#conditional error function

A <- plot(type = 4,
     design_base,
     plotNonMonotoneFunction = TRUE,
     range = c(0,alpha0+0.05))

B <- plot(type = 1,
          design_base,
          plotNonMonotoneFunction = TRUE,
          range = c(0,alpha0+0.05))

#Plot of the corresponding second stage sample sizes
C <- plot_sample_size(design_base, range = c(0, alpha0+0.05))

#Combine the plots to one figure
figure1 <- cowplot::plot_grid(A, B, C, labels = "AUTO", ncol = 3)

### design_direct_mono ---------------------------------------------------------

#Create a design that does not need to be monotonized
design_direct_mono <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = deltaLR,
  delta1Min = 0.2)

### design_constant ------------------------------------------------------------

#Create a design with a constant optimal monotone conditional error function
design_constant <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = 0,
  delta1Min = delta1Min)

### Figure 2 -------------------------------------------------------------------

#Plot of the optimal conditional error function (already monotone)
A <- plot(type = 1,
          design_direct_mono,
          range = c(0,alpha0+0.05))

#Plot of the optimal conditional error function and the optimal monotone
#conditional error function
B <- plot(type = 1,
          design_constant,
          plotNonMonotoneFunction = TRUE,
          range = c(0,alpha0+0.05))

#Combine the plots to one figure
figure2 <- cowplot::plot_grid(A, B, labels = "AUTO")

### design_maxlr ---------------------------------------------------------------

#Design that uses the maximum likelihood ratio
design_maxlr <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "maxlr",
  delta1Min = delta1Min)

### Figure 3 -------------------------------------------------------------------

#Create a sequence of first stage p-values for plotting
firstStagePValues <- seq(from = 0, to = alpha0+0.05, length.out = 1e3)

#Calculate the corresponding optimal monotone conditional errors
optimalConditionalErrors_fixed <- getOptimalConditionalError(
  firstStagePValue = firstStagePValues,
  design = design_base)

optimalConditionalErrors_maxlr <- getOptimalConditionalError(
  firstStagePValue = firstStagePValues,
  design = design_maxlr)

#Calculate conditional errors for the inverse normal method
a2 <- getAlpha2InverseNormal(alpha = alpha, alpha1 = alpha1, alpha0 = alpha0)
ConditionalErrors_InverseNormal <-
  getInverseNormalConditionalError(alpha = alpha, alpha1 = alpha1, alpha0 = alpha0,
                                   alpha2 = a2,
                                   firstStagePValue = firstStagePValues)

#Calculate conditional errors for the inverse normal method weighted 1/3 2/3
a2 <- getAlpha2InverseNormal(alpha = alpha, alpha1 = alpha1, alpha0 = alpha0)
ConditionalErrors_InverseNormal_w <-
  getInverseNormalConditionalError(alpha = alpha, alpha1 = alpha1, alpha0 = alpha0,
                                   alpha2 = a2,
                                   firstStagePValue = firstStagePValues,
                                   weights = c(sqrt(1/3),sqrt(2/3)))

#Create data set used for plotting the three different functions
dat_fig3 <- data.frame(firstStagePValues = rep(firstStagePValues, 4),
                       fun = rep(1:4, each = length(firstStagePValues)),
                       cond_err = c(ConditionalErrors_InverseNormal,
                                    ConditionalErrors_InverseNormal_w,
                                    optimalConditionalErrors_fixed,
                                    optimalConditionalErrors_maxlr))

dat_fig3$fun <- factor(dat_fig3$fun, levels = 1:4, labels = c("inv. normal 1",
                                                              "inv. normal 2",
                                                              "Delta = 0.2",
                                                              "max. likel."))

figure3 <- ggplot2::ggplot(data = dat_fig3, ggplot2::aes(x = firstStagePValues,
                                                         y = cond_err))+
  ggplot2::geom_line(ggplot2::aes(linetype = fun),
                     colour = "black",
                     linewidth = 1.1) +
  ggplot2::labs(x = "First-stage p-value",
                y = "Conditional error function") +
  ggplot2::theme_bw() +
  ggplot2::geom_vline(xintercept = alpha1,
                      linetype = "dotted",
                      col = "blue",
                      linewidth = 0.8) +
  ggplot2::geom_vline(xintercept = alpha0,
                      linetype = "dotted",
                      col = "red",
                      linewidth = 0.8) +
  ggplot2::xlim(c(0, alpha0+0.05))+
  ggplot2::scale_fill_discrete(breaks = c("inv. normal 1", "inv. normal 2", "Delta = 0.2", "max. likel."))+
  ggplot2::scale_linetype_manual(name = "Conditional error function",
                                 values = 1:4,
                                 labels = c(expression(paste("inv. normal (", w[1], " = ", w[2], " = ", sqrt(1/2), ")")),
                                            expression(paste("inv. normal (", w[1], " = ", sqrt(1/3), ", ", w[2], " = ", sqrt(2/3), ")")),
                                            expression(paste("optimal monotone (", Delta, " = 0.2)")),
                                            "optimal monotone (max. likel.)"))+
  ggplot2::theme(legend.text = ggplot2::element_text(size=10),
                 legend.position = "inside",
                 legend.position.inside = c(0.6, 0.8),
                 legend.background = ggplot2::element_rect(fill = "white",
                                                           color = "black"),
                 legend.key.width = ggplot2::unit(3, "line"))


### Maximum sample sizes -------------------------------------------------------

#maximum sample size inverse normal method - equally weighted
ceiling(2*getSecondStageInformationIN(p1 = alpha0,
                                      al = alpha,
                                      al0 = alpha0,
                                      al1 = alpha1,
                                      cp = conditionalPower,
                                      I1 = firstStageInformation,
                                      delta_Min = delta1Min))

#maximum sample size inverse normal method - weighted
ceiling(2*getSecondStageInformationIN(p1 = alpha0,
                                      al = alpha,
                                      al0 = alpha0,
                                      al1 = alpha1,
                                      cp = conditionalPower,
                                      I1 = firstStageInformation,
                                      delta_Min = delta1Min,
                                      weights = c(sqrt(1/3),sqrt(2/3))))

#maximum sample size fixed effect design
ceiling(2*getSecondStageInformation(
  design = design_base,
  firstStagePValue = alpha0
))

#maximum sample size maxlr
ceiling(2*getSecondStageInformation(
  design = design_maxlr,
  firstStagePValue = alpha0
))

### design_base_0 ------------------------------------

design_base_0 <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = 0,
  delta1Min = delta1Min)

### design_base_0.125 ------------------------------------

design_base_0.125 <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = 0.125,
  delta1Min = delta1Min)

### Table expected sample sizes ------------------------------------------------

tab_expected_sample_size <- data.frame()

for (Delta in c(0, 0.125, 0.2)){

Exp_Info_inorm <- getExpectedSecondStageInformationIN(al = alpha,
                                                      al1 = alpha1,
                                                      al0 = alpha0,
                                                      cp = conditionalPower,
                                                      I = firstStageInformation,
                                                      delta_Min = delta1Min,
                                                      delta = Delta)

Exp_Info_inorm_weighted <- getExpectedSecondStageInformationIN(al = alpha,
                                                               al1 = alpha1,
                                                               al0 = alpha0,
                                                               cp = conditionalPower,
                                                               I = firstStageInformation,
                                                               delta_Min = delta1Min,
                                                               delta = Delta,
                                                               weights = c(sqrt(1/3),sqrt(2/3)))

Exp_Info_base_0 <- optconerrf::getExpectedSecondStageInformation(
  design = design_base_0,
  likelihoodRatioDistribution = "fixed",
  deltaLR = Delta)

Exp_Info_base_0.125 <- optconerrf::getExpectedSecondStageInformation(
  design = design_base_0.125,
  likelihoodRatioDistribution = "fixed",
  deltaLR = Delta)

Exp_Info_base <- optconerrf::getExpectedSecondStageInformation(
  design = design_base,
  likelihoodRatioDistribution = "fixed",
  deltaLR = Delta)

Exp_Info_maxlr <- optconerrf::getExpectedSecondStageInformation(
  design = design_maxlr,
  likelihoodRatioDistribution = "fixed",
  deltaLR = Delta)

row <- round(c(Delta, 2*Exp_Info_inorm, 2*Exp_Info_inorm_weighted, 2*Exp_Info_base_0,
               2*Exp_Info_base_0.125, 2*Exp_Info_base, 2*Exp_Info_maxlr), 3)

tab_expected_sample_size <- rbind(tab_expected_sample_size, row)
}

colnames(tab_expected_sample_size) = c("True delta", "Inv. normal 1", "Inv. normal 2", "\u0394 = 0", "\u0394 = 0.125", "\u0394 = 0.2", "Max. likel.")

tab_expected_sample_size

### Overall power (cp = 0.8) ---------------------------------------------------

#Overall Power at deltaLR = 0, deltaLR = 0.125, deltaLR = 0.2
alternative <- c(0, 0.125, 0.2)

power_IN1 <- getOverallPowerIN(al = alpha,
                  al1 = alpha1,
                  al0 = alpha0,
                  cp = conditionalPower,
                  I = firstStageInformation,
                  delta_Min = delta1Min,
                  weights = c(sqrt(0.5),sqrt(0.5)),
                  alternative = alternative)

power_IN2 <- getOverallPowerIN(al = alpha,
                  al1 = alpha1,
                  al0 = alpha0,
                  cp = conditionalPower,
                  I = firstStageInformation,
                  delta_Min = delta1Min,
                  weights = c(sqrt(1/3),sqrt(2/3)),
                  alternative = alternative)

power_base_0 <- getOverallPower(design_base_0, alternative = alternative)
power_base_0.125 <- getOverallPower(design_base_0.125, alternative = alternative)
power_base <- getOverallPower(design_base, alternative = alternative)
power_maxlr <- getOverallPower(design_maxlr, alternative = alternative)

tab_power <- data.frame(alternative, power_IN1$overallPower, power_IN2$overallPower,
                        power_base_0$overallPower, power_base_0.125$overallPower,
                        power_base$overallPower, power_maxlr$overallPower)

colnames(tab_power) = c("True delta", "Inv. normal 1", "Inv. normal 2", "\u0394 = 0", "\u0394 = 0.125", "\u0394 = 0.2", "Max. likel.")

round(tab_power, 3)

### First stage sample size for overall power control --------------------------

#Derive the minimal first stage sample size to reach an overall power of 0.8 at
#deltaLR

n_1_IN1 <- ssanv::uniroot.integer(overallPower_IN1, c(100, 400),
                                  pos.side = TRUE)$root
n_1_IN2 <- ssanv::uniroot.integer(overallPower_IN2, c(100, 400),
                                  pos.side = TRUE)$root
n_1_base <- ssanv::uniroot.integer(overallPower_base, c(100, 400),
                                   pos.side = TRUE)$root
n_1_maxlr <- ssanv::uniroot.integer(overallPower_maxlr, c(100, 400),
                                    pos.side = TRUE)$root

n_1_IN1; n_1_IN2; n_1_base; n_1_maxlr

### design_base_op -------------------------------------------------------------

#Base design that reaches overall power
#(different choice of first stage sample size)

design_base_op <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = 0.5*n_1_base,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "fixed",
  deltaLR = deltaLR,
  delta1Min = delta1Min)

### design_maxlr_op ------------------------------------------------------------

#Maximum likelihood design that reaches overall power
#(different choice of first stage sample size)
design_maxlr_op <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = 0.5*n_1_maxlr,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "maxlr",
  delta1Min = delta1Min)

### Maximum overall sample size ------------------------------------------------

#maximum sample size inverse normal method - equally weighted
ceiling(2*getSecondStageInformationIN(p1 = alpha0,
                                      al = alpha,
                                      al0 = alpha0,
                                      al1 = alpha1,
                                      cp = conditionalPower,
                                      I1 = 0.5 * n_1_IN1,
                                      delta_Min = delta1Min)) + n_1_IN1

#maximum sample size inverse normal method - weighted
ceiling(2*getSecondStageInformationIN(p1 = alpha0,
                                      al = alpha,
                                      al0 = alpha0,
                                      al1 = alpha1,
                                      cp = conditionalPower,
                                      I1 = 0.5 * n_1_IN2,
                                      delta_Min = delta1Min,
                                      weights = c(sqrt(1/3),sqrt(2/3)))) + n_1_IN2

#maximum sample size fixed effect design
ceiling(2*getSecondStageInformation(
  design = design_base_op,
  firstStagePValue = alpha0
)) + n_1_base

#maximum sample size maxlr
ceiling(2*getSecondStageInformation(
  design = design_maxlr_op,
  firstStagePValue = alpha0
)) + n_1_maxlr

### Expected overall sample size -----------------------------------------------

tab_expected_sample_size_ov <- data.frame()

for (Delta in c(0, 0.125, 0.2)){

  Exp_Info_inorm <- getExpectedSecondStageInformationIN(al = alpha,
                                                        al1 = alpha1,
                                                        al0 = alpha0,
                                                        cp = conditionalPower,
                                                        I = 0.5 * n_1_IN1,
                                                        delta_Min = delta1Min,
                                                        delta = Delta,
                                                        weights = c(sqrt(1/2),sqrt(1/2)))

  Exp_Info_inorm_weighted <- getExpectedSecondStageInformationIN(al = alpha,
                                                                 al1 = alpha1,
                                                                 al0 = alpha0,
                                                                 cp = conditionalPower,
                                                                 I = 0.5 * n_1_IN2,
                                                                 delta_Min = delta1Min,
                                                                 delta = Delta,
                                                                 weights = c(sqrt(1/3),sqrt(2/3)))

  Exp_Info_base <- optconerrf::getExpectedSecondStageInformation(
    design = design_base_op,
    likelihoodRatioDistribution = "fixed",
    deltaLR = Delta)

  Exp_Info_maxlr <- optconerrf::getExpectedSecondStageInformation(
    design = design_maxlr_op,
    likelihoodRatioDistribution = "fixed",
    deltaLR = Delta)

  row <- round(c(Delta,
                 2*Exp_Info_inorm + n_1_IN1,
                 2*Exp_Info_inorm_weighted + n_1_IN2,
                 2*Exp_Info_base + n_1_base,
                 2*Exp_Info_maxlr+ n_1_maxlr), 3)

  tab_expected_sample_size_ov <- rbind(tab_expected_sample_size_ov, row)
}

colnames(tab_expected_sample_size_ov) = c("True delta", "Inv. normal 1", "Inv. normal 2",
                                       "Delta = 0.2", "Max. likel.")

tab_expected_sample_size_ov

### design_constraints ---------------------------------------------------------

design_constraints <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1,
  alpha0 = alpha0,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "maxlr",
  delta1Min = delta1Min,
  maximumConditionalError = 0.25,
  maximumSecondStageInformation = 310
)

### Figure 4 -------------------------------------------------------------------

A <- plot(type = 1,
          design_constraints,
          plotNonMonotoneFunction = TRUE,
          range = c(0, alpha0+0.05))

B <- plot_sample_size(design_constraints, range = c(0, alpha0+0.05))

figure4 <- cowplot::plot_grid(A, B, labels = "AUTO")

### design_almost_free_lunch ---------------------------------------------------

alpha1_modified <- 0
alpha0_modified <- 0.9
firstStageInformation_modified <- firstStageInformation

design_almost_free_lunch <- optconerrf::getDesignOptimalConditionalErrorFunction(
  alpha = alpha,
  alpha1 = alpha1_modified,
  alpha0 = alpha0_modified,
  conditionalPower = conditionalPower,
  firstStageInformation = firstStageInformation_modified,
  useInterimEstimate = TRUE,
  likelihoodRatioDistribution = "maxlr",
  delta1Min = delta1Min,
  minimumConditionalError = alpha,
  maximumConditionalError = 0.5
)

### Figure 5 -------------------------------------------------------------------

A <- plot(type = 1,
          design_almost_free_lunch,
          range = c(0, alpha0_modified + 0.05),
          plotNonMonotoneFunction = FALSE)

#Function for the second stage sample size of the fixed design

sample_size_fixed <- function(p_1){
  delta_hat <- pmax(delta1Min, qnorm(1-p_1)/sqrt(firstStageInformation_modified))
  sample_size <- 2*(qnorm(1-alpha) + qnorm(0.8))^2/delta_hat^2
  return(sample_size)
}
sample_size_fixed <- Vectorize(sample_size_fixed)

firstStagePValues <- seq(from = 0, to = alpha0_modified + 0.05, length.out = 1e3)
sample_size_fixed_values <- sample_size_fixed(firstStagePValues)

B <- plot_sample_size(design_almost_free_lunch,
                      range = c(0, alpha0_modified + 0.05),
                      plotNonMonotoneFunction = FALSE)+
     ggplot2::geom_line(mapping = ggplot2::aes(x = firstStagePValues, y = sample_size_fixed_values),
                        linewidth = 1.1, linetype = "dotted")


figure5 <- cowplot::plot_grid(A, B, labels = "AUTO")

### Maximum overall sample size ------------------------------------------------

#Adaptive design
ceiling(2*getSecondStageInformation(
  design = design_almost_free_lunch,
  firstStagePValue = alpha0_modified
))

#Separate pilot study
ceiling(2*(qnorm(1-alpha) + qnorm(0.8))^2/delta1Min^2)

### Expected overall sample size -----------------------------------------------

tab_expected_sample_size <- data.frame()

for (Delta in c(0, 0.125, 0.2)){

  #Optimal adative design
  Exp_sample_size_adaptive <-
    2*optconerrf::getExpectedSecondStageInformation(
      design = design_almost_free_lunch,
      likelihoodRatioDistribution = "fixed",
      deltaLR = Delta)

  #Separate pilot study
  int <- function(p_1){
    delta_hat <- pmax(delta1Min, qnorm(1-p_1)/sqrt(firstStageInformation_modified))
    nu_alpha <- (qnorm(1-0.05)+qnorm(0.8))^2
    density <- exp(qnorm(1-p_1)*sqrt(firstStageInformation_modified)*Delta-
                     firstStageInformation_modified*Delta^2/2)
    return((nu_alpha/delta_hat^2)*density)
  }

  Exp_sample_size_pilot <-
  2*integrate(f = int, lower = alpha1_modified, upper = 1)$value

  row <- round(c(Delta, Exp_sample_size_adaptive, Exp_sample_size_pilot), 3)

  tab_expected_sample_size <- rbind(tab_expected_sample_size, row)
}

colnames(tab_expected_sample_size) = c("True delta", "Adaptive", "Pilot")

tab_expected_sample_size

#Overall power and first-stage futility stopping probability (Optimal adaptive design A)
getOverallPower(design_almost_free_lunch, alternative = c(0, 0.125, 0.2))

#Overall power (separate study design S)
Delta <- 0
secondStageRejection <- function(p_1) {

  ncp1 <- Delta * sqrt(firstStageInformation_modified)
  delta_hat <- pmax(delta1Min, qnorm(1-p_1)/sqrt(firstStageInformation_modified))
  nu_alpha <- (qnorm(1-0.05)+qnorm(0.8))^2
  (1 -
      stats::pnorm(
        stats::qnorm(
          1 - alpha) -
          sqrt(nu_alpha/delta_hat^2) *
          Delta
      )) *
    exp(qnorm(1 - p_1) * ncp1 - ncp1^2 / 2)
}

integral <- stats::integrate(
  f = secondStageRejection,
  lower = alpha1_modified,
  upper = 1
)$value; integral

Delta <- 0.125
secondStageRejection <- function(p_1) {

  ncp1 <- Delta * sqrt(firstStageInformation_modified)
  delta_hat <- pmax(delta1Min, qnorm(1-p_1)/sqrt(firstStageInformation_modified))
  nu_alpha <- (qnorm(1-0.05)+qnorm(0.8))^2
  (1 -
      stats::pnorm(
        stats::qnorm(
          1 - alpha) -
          sqrt(nu_alpha/delta_hat^2) *
          Delta
      )) *
    exp(qnorm(1 - p_1) * ncp1 - ncp1^2 / 2)
}

integral <- stats::integrate(
  f = secondStageRejection,
  lower = alpha1_modified,
  upper = 1
)$value; integral


Delta <- 0.2
secondStageRejection <- function(p_1) {

  ncp1 <- Delta * sqrt(firstStageInformation_modified)
  delta_hat <- pmax(delta1Min, qnorm(1-p_1)/sqrt(firstStageInformation_modified))
  nu_alpha <- (qnorm(1-0.05)+qnorm(0.8))^2
    (1 -
       stats::pnorm(
         stats::qnorm(
           1 - alpha) -
           sqrt(nu_alpha/delta_hat^2) *
           Delta
       )) *
      exp(qnorm(1 - p_1) * ncp1 - ncp1^2 / 2)
  }

integral <- stats::integrate(
  f = secondStageRejection,
  lower = alpha1_modified,
  upper = 1
)$value; integral


### Save plots -----------------------------------------------------------------
ggplot2::ggsave(figure1, filename = "figure1.png",
                units = "px", width = 3750, height = 1200, dpi = 300)

ggplot2::ggsave(figure2, filename = "figure2.png",
                units = "px", width = 2500, height = 1200, dpi = 300)

ggplot2::ggsave(figure3, filename = "figure3.png",
                units = "px", width = 2500, height = 1200, dpi = 300)

ggplot2::ggsave(figure4, filename = "figure4.png",
                units = "px", width = 2500, height = 1200, dpi = 300)

ggplot2::ggsave(figure5, filename = "figure5.png",
                units = "px", width = 2500, height = 1200, dpi = 300)


