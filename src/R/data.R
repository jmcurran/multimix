#'Prostate cancer patient data
#'
#'This data set was considered in Byar and Green (1980) and reproduced in
#'Andrews and Herzberg (1985, pp 261– 247).
#'
#'There are 12 pre-trial covariates measured on each patient, seven may be taken
#'to be continuous, four to be discrete, and one variable (\code{sg}) is an
#'index nearly all of whose values lie between 7 and 15, and which could be
#'considered either discrete or continuous. A preliminary inspection of the data
#'showed that the size of the primary tumour (\code{sz}) and serum prostatic
#'acid phosphatase (\code{ap}) were both skewed variables. These variables have
#'therefore been transformed. A square root transformation was used for
#'\code{sz}, and a logarithmic transformation was used for \code{ap} to achieve
#'approximate normality. (As for correlation, skewness over the whole data set
#'does not necessarily mean skewness within clusters but when clusters were
#'formed within-cluster skewness was observed for these variables.) Observations
#'that had missing values in any of the 12 pretreatment covariates were omitted
#'from further analysis, leaving 475 out of the original 506 observations
#'available.


#'
#' @format A data.frame with 475 rows and 12 columns:
#' \describe{
#' \item{age}{Age in years.}
#' \item{wt}{Weight in pounds}
#' \item{pf}{Patient activity.}
#' \item{hx}{Family history of cancer.}
#' \item{sbp}{Systolic blood pressure.}
#' \item{dbp}{Diastolic blood pressure.}
#' \item{ekg}{Electrocardiogram code.}
#' \item{hg}{Serum haemoglobin.}
#' \item{sz}{Size of primary tumour.}
#' \item{sg}{Index of tumour stage and histolic grade. }
#' \item{ap}{Serum prostatic acid phosphatase.}
#' \item{bm}{Bone metastatses}
#' }
#'
#' @source D.P. Byar and S.B. Green "The choice of treatment for cancer patients
#'   based on covariate information - application to prostate cancer", Bulletin
#'   du Cancer 1980: 67:477--490, reproduced in D.A. Andrews and A.M. Herzberg
#'   "Data: a collection of problems from many fields for the student and
#'   research worker" p.261--274 Springer series in statistics, Springer-Verlag.
#'   New York.
"cancer.df"