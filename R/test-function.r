#' Title
#'
#' @return
#' @export
#'
#' @examples
testScript = function(){
  # Use setwd(folder) to change to the folder in which the Multimix functions are.
  
  # source("data_organise.r")
  # source("Functions.r")
  # source("first.Z.to.P.r")
  # source("P.to.Z.r")
  # source("Z.to.P.r")
  # source("mmain.r")
  
  
  coltype = c("numeric", "numeric", "factor", "factor", "numeric", "numeric", "factor", 
  "numeric", "numeric", "numeric", "numeric", "factor")
  cancer = read.table("CANCER11a.DAT", header=FALSE, colClasses=coltype)
  
  Stage = scan(system.file("extdata", "stage.txt",  package = "multimix"))
  
  Stage = Stage - 2
  Init_grp = as.factor(Stage)
  names(cancer)  = c("Age", "Wt", "PF" , "HX", "SBP", "DBP",  "EKG", "HG", "SZ", "SG", "AP", "BM" )
  levels(cancer[,3])  <- c("Active", "Bed49", "Bed51")
  levels(cancer[,4])  <- c("No_hist", "Hist")
  levels(cancer[,7])  <- c("Normal", "Benign", "Rythmic", "H_blocks", "H_strain", "Old_MCI", "New_MCI")
  levels(cancer[,12]) <- c("No_BM", "BM")
  # dframe = transform(cancer, Age=Age/100, Wt=Wt/100, SBP=SBP/10, DBP=DBP/10, HG=HG/100, SG=SG/10)
   dframe = cancer
  # scaling does not seem necessary for the cancer df, and it makes it harder to interpret parameter values.
   
  
  
  cdep <- list(c(5, 6))              # list of MVN cells
  lcdep <- list(c(12, 2, 8), c(3,1)) # list of location cells, discrete first in each cell
  library(mvtnorm)
  D <- data_organise(dframe=dframe, Init_grp=Init_grp, cdep=cdep, lcdep=lcdep)
  zpr <- mmain(D)
  zpr$results
}