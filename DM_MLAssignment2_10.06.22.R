
# Name: Drashti Mehta
# UNC email ID: dmehta12@uncc.edu
# Date: 10/06/22
# 
# ---------------------------------------------------------------------------------------------------------------------------------
#   
#   I am hereby using two methods: 
#   1. Parametric Linear Regression and 
# - I am performing this analysis as it is the most basic analysis and is usually first performed
# - It is by far the easiest way to analyse the data
# - This model is used when we want to check if a variable is influencing other or not
# - For us, we need to check here if Bootstrap (BS) is related to S or not and weither LLD (as a fucntion of S) is affected or not
# - This model doesn't seem to be a ideal model to test the correlation between BS and S
#
# 2. Linear Correlation using BoxCox transformation method.
# - I here used this as it seems to be a good option for transforming part of our data to check if that gives a better result
# - It is used to transform the data into normal distribution
# - We here use lambda as a measure of power to which the data is raised
# - Log y (where y is LLD) is considered instead of y because the y will be only 1 and not transform
# - Here it was giving a p-value of zero
# - Results of model having trnasformed values are printed after new model is obtained
# - A plot is generated here too
# - We get a lot of underestimated and overestimated values
# - If we compare the two plots here, the curve seems to be decreased and the variance is evernly spread than before
#
# QUESTION:
# 
# "Based on your analysis, can Boostrap frequencies be an indirect method of measuring support 
# sensu Hacking (1965) in maximum likelihood phylogenetic analyses? Why?"
# 
# ANSWER:
# 
# Based of the methods I used, it seems like the bootstrap frequencies cannot be a method for measuring support
# because from the plot and analysis we have understood that BS is not proportional to the LLD and hence not to 
# S. Results stay a little overestimated or underestimated in the analysis. They have high frequency of outliers. 
# Linear and non linear models do not have a higher difference in the prediction quality. They leave quite broad
# confidence intervals (around 95%).


#!/usr/local/bin/Rscript

# providing data through Table 3 tsv file
readData <- function () {
  rawData <- read.csv("onlineAppendix3_TableS3.tsv", header = TRUE, sep = "\t")
  rawData <- rawData[ -c(4) ]
  colnames(rawData) <- c("Dataset", "BS", "LLD")
  return(rawData)
}
rawData <- readData()
head(rawData)

# Regression_Analysis -----------------------------------------------------

model_reg <- lm(rawData$BS ~ rawData$LLD)
coef(model_reg)
# plotting regression model
plot(
  x = rawData$BS,
  y = rawData$LLD,
  main = "Regression model: BS and LLD",
  xlab = "Bootstrap values (BS)",
  ylab = "Log-likelihood Differences (LLD)"
)

#.............................BOX-COX TRANSFORMATIONS: 

# loading libraries:
cat("> Loading necessary libraries:\n") # to print to the screen or a file
if (!require("ggplot2")) install.packages("ggplot2"); 
library(ggplot2)
theme_set(theme_classic()) #to set the theme for plots created with GGPLOT2 in this R file
if (!require("MASS")) install.packages("MASS"); 
library(MASS)

# Get fresh data:
cat("> Raw data (head):\n")
rawData <- read.csv("onlineAppendix3_TableS3.tsv", header = TRUE, sep = "\t")
rawData <- rawData[ -c(4) ]
colnames(rawData) <- c("Dataset", "BS", "LLD")
head(rawData)

# TRANSFORMATION:



# BoxCox_Transformation ---------------------------------------------------

# Define the powerTransform function (for Box Cox transformations)
# lambda indicates the power to which all the data should be raised
powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    if (lam1 == 0L) {
      log(y + lam2) # log is considered coz if only y is cannot be transformed and will be 1 always in such cases
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  switch(method,
         boxcox = boxcoxTrans(y, lambda1, lambda2),
         tukey = y^lambda1
  )
}

# one-parameter box-cox transformation:
y <- rawData$LLD
x <- rawData$BS

# Running box-cox transformation - plotting
bc <- boxcox(y ~ x)

# Get lambda
lambda <- bc$x[which.max(bc$y)]
fit <- lm(powerTransform(y, lambda) ~ x) # new model fitted
cat("> Summary of the linear model after Box-Cox transformation of LLD:")
summary(fit)
cat(paste0("> Adjusted R-squared = ", summary(fit)$adj.r.squared, "\n"))
cat(paste0("> Slope = ", fit$coef[[2]], "\n"))
cat(paste0("> Intercept = ", fit$coef[[1]], "\n"))

# Visualize the data in a dotplot
ggplotFit <- function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point(col = alpha("light pink", 0.4)) +
    stat_smooth(method = "lm", col = "light blue") +
    labs(
      subtitle = bquote("Adj." ~ R^2 == .(signif(summary(fit)$adj.r.squared, 5)) ~ "/ Intercept" == .(signif(fit$coef[[1]],5)) ~ "/ Slope" == .(signif(fit$coef[[2]], 5)) ~ "/ P" == .(signif(summary(fit)$coef[2,4], 5))),
      title = "OLS linear regression of concatenated datasets after Box Cox transformation"
    )+
    ylab("Transformed LLD (Box Cox)") +
    xlab("BS")
}

pdf("DM_boxcox1.pdf")

ggplotFit(fit)

garbage <- dev.off()

quit()
