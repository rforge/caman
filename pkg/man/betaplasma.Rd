\name{betaplasma}
\docType{data}
\alias{betaplasma}
\title{Determinants of Beta-Carotene Levels}
\description{
Observational studies have suggested that low dietary intake or low plasma 
concentrations of retinol, beta-carotene, or other carotenoids might be associated 
with increased risk of developing certain types of cancer.  However, relatively few 
studies have investigated the determinants of plasma concentrations of these 
micronutrients. Nierenberg et al. designed a cross-sectional study to investigate the relationship between personal characteristics and dietary factors, and plasma concentrations of 
retinol, beta-carotene and other carotenoids. Study subjects (N = 315) were patients 
who had an elective surgical procedure during a three-year period to biopsy or 
remove a lesion of the lung, colon, breast, skin, ovary or uterus that was found to be non-cancerous. 
}
\usage{betaplasma}
\format{datafile contains 315 observations on 14 variables.  This data set can 
be used to demonstrate multiple 
regression, transformations, categorical variables, 
outliers, pooled tests of significance and model building strategies. 
}
\arguments{
Variable Names in order from left to right:
        \item{AGE}{Age (years)} \cr
        \item{SEX}{Sex (1=Male, 2=Female).} \cr
        \item{SMOKSTAT}{Smoking status (1=Never, 2=Former, 3=Current Smoker)} \cr
        \item{QUETELET}{Quetelet (weight/(height^2))} \cr
        \item{VITUSE}{Vitamin Use (1=Yes, fairly often, 2=Yes, not often, 3=No)} \cr
        \item{CALORIES}{Number of calories consumed per day.} \cr
        \item{FAT}{Grams of fat consumed per day.} \cr
        \item{FIBER}{Grams of fiber consumed per day.} \cr
        \item{ALCOHOL}{Number of alcoholic drinks consumed per week.} \cr
        \item{CHOLESTEROL}{Cholesterol consumed (mg per day).} \cr
        \item{BETADIET}{Dietary beta-carotene consumed (mcg per day).} \cr
        \item{RETDIET}{Dietary retinol consumed (mcg per day)} \cr
        \item{BETAPLASMA}{Plasma beta-carotene (ng/ml)} \cr
        \item{RETPLASMA}{Plasma Retinol (ng/ml)} \cr
        \item{LOGBPL}{log(betaplasma)} \cr
}
\examples{
data(betaplasma)
beta4 <- mixcov(dep="betacaro",fixed=c("chol","sex","bmi"),random="betadiet",data=betaplasma,k=4,family="gaussian")

}
\references{
Schlattmann, P.(2008) \emph{Medical Applications of Finite Mixture Models} Berlin: Springer \cr \cr

These data have not been published yet but a related reference is \cr
Nierenberg DW, Stukel TA, Baron JA, Dain BJ, Greenberg ER.  Determinants of plasma levels of 
beta-carotene and retinol.  \emph{American Journal of Epidemiology} 1989;130:511-521. \cr \cr

The use of the data is kindly permitted Dr. Therese Stukel, Dartmouth Hitchcock Medical Center, USA \cr
The data can also be found at StatLib: http://lib.stat.cmu.edu/index.php

}
\keyword{datasets}
