#' Prostate cancer patient data
#'
#' Data on 475 prostate cancer patients
#' 
#' There are twelve pre-trial covariates measured on each patient,
#' seven may be taken to be continuous, four to be discrete, and one 
#' variable (SG) is an index nearly all of whose values lie between 
#' 7 and 15, and which could be considered either discrete or continuous. 
#' We will treat SG as a continuous variable.
#'
#' A preliminary inspection of the data showed that the sizeof the 
#' primary tumour (SZ) and serum prostatic acid phosphatase (AP) were 
#' both skewed variables. These variables have therefore been transformed. 
#' A square root transformation was used for SZ, and a logarithmic
#' transformation was used for AP to achieve approximate normality. 
#' (As for correlation, skewness over the whole data set does not 
#' necessarily mean skewness within clusters. But when clusters were 
#' formed, within-cluster skewness was observed for these variables.) 
#' 
#' Observations that had missing values in any of the twelve pretreatment 
#' covariates were omitted from furtheranalysis, leaving 475 out of the 
#' original 506 observations available. 
#' 
#' The categorical variable \code{Patient activity} had 4 levels: "Normally 
#' Active", "Bed rest below 50% of waking hours", "Bed rest 50% of waking hours
#' or more", and "Confined to bed". The numbers of the 475 in these groups were
#'  428, 32, 12, and 3. The least active two groups are grouped in our data,
#' giving 3 groups of size 428, 32, and 15.
#'
#' @format A data.frame with 475 rows and 12 columns:
#' \describe{
#' \item{age}{Age in years}
#' \item{wt}{Weight in pounds}
#' \item{pf}{Patient activity}
#' \item{hx}{Family history of cancer}
#' \item{sbp}{Systolic blood pressure}
#' \item{dbp}{Diastolic blood pressure}
#' \item{ekg}{Electrocardiogram code}
#' \item{hg}{Serum haemoglobin}
#' \item{sz}{Size of primary tumour}
#' \item{sg}{Index of tumour stage and histolic grade}
#' \item{ap}{Serum prostatic acid phosphatase}
#' \item{bm}{Bone metastatses}
#' }
#' 
#' @usage 
#' data(cancer.df)
#'
#' @source D.P. Byar and S.B. Green "The choice of treatment for cancer patients
#'   based on covariate information - application to prostate cancer", Bulletin
#'   du Cancer 1980: 67:477--490, reproduced in D.A. Andrews and A.M. Herzberg
#'   "Data: a collection of problems from many fields for the student and
#'   research worker" p.261--274 Springer series in statistics, Springer-Verlag.
#'   New York.
"cancer.df"

#' Contraceptive Method Choice data
#' 
#' This dataset is a subset of the 1987 National Indonesia Contraceptive
#' Prevalence Survey. The cases are 1473 married women who were either not 
#' pregnant or do not know if they were at the time of interview. 
#' 
#' The variables "age" (in years) and "nborn" (ranging from 0 to 16) would
#' normally be treated as continuous; "nborn" is skew and might well
#' be transformed. The remaining 8 variables are categorical. 
#' 
#' The variables "edu", "eduh" and "sol" take values "1,2,3,4", #' they are 
#' ordinal with 1 = low and 4 = high. The variable "husocc" takes the
#' same 4 values, but it is not clear if the order has any significance.
#' 
#' The variables "islam", "working", and "medex" are binary (9,1)-valued with
#' 0=Non-Islam, 1=Islam for "islam";  0=Yes, 1=No for "working"; and 0=Good, 
#' 1=Not good for "medex".
#' 
#' The variable "method" is ternary: 1=No-use, 2=Long-term, 3=Short-term.
#' 
#' @format A data.frame with 1473 rows and 10 columns:
#' \describe{
#' \item{age}{Wife's age}
#' \item{edu}{Wife's education}
#' \item{eduh}{Husband's education}
#' \item{nborn}{Number of children ever born}
#' \item{islam}{Wife's religion}
#' \item{working}{Wife's now working?}
#' \item{husocc}{Husband's occupation}
#' \item{sol}{Standard-of-living index}
#' \item{medex}{Media exposure}
#' \item{method}{Contraceptive method used}
#' }
#' 
#' @usage 
#' data(cmc.df)
#'   
#' @source Tjen-Sien Lim "Contraceptive Method Choice" 1997,
#'    UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. 
#'   Irvine, CA: University of California, School of Information and Computer 
#'   Science.
"cmc.df"

