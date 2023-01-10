#' S3 method for plotting multimix resuls objects
#'
#' @param x an object of class \code{multimixResults} -- see
#' \code{\link{mmain}} for more information.
#' @param \ldots any other arguments to be passed to \code{plot}. Note that
#' because there are two calls to plot, the \ldots arguments will be passed to 
#' each call, and it is unlikely that this will have the desired effect.
#' 
#' @importFrom graphics abline axis box legend lines par
#' @author James Curran
#' @export
plot.multimixResults = function(x, ...){
  oldPar = par(mfrow = c(2, 1), # two plots, one on each row
               tcl = -0.3, # length of tick marks
               mai = c(0.2, 0.7, 0.2, 0.1), # plot margins (b, l, t, r) in inches
               mgp = c(2, 0.3, 0), # distance of title, axis labels, and axis in lines
               las = 1, # keep tick labels horizontal always
               xaxs = "i", # make the x-axis intersect with the limits of the plot
               yaxs = "i") # make the y-axis intersect with the limits of the plot
  
  numIter = nrow(x$results)
  ll = x$results[, "Log-likelihood"]
  ymin = signif(min(ll))
  ymax = ceiling(max(ll))
  
  plot(1:numIter, ll,
       type = "l",
       lwd = 1.5,
       xlab = "",
       ylab = "Log-likelihood",
       xlim = c(1, numIter),
       ylim = c(ymin, ymax),
       axes = FALSE,
       ...)
  abline(h = max(ll), lty = 2, col = "lightgrey")
  axis(2, at = c(ymin, max(ll), ymax), cex.axis = 0.7)
  box()
  

  par(mai = c(0.7, 0.7, 0, 0.1))
  groupProbs = x$results[, -c(1, ncol(x$results))]
  numGroups = ncol(groupProbs)
  
  plot(1:numIter, groupProbs[ ,1], 
       col = 1, 
       lwd = 1.5,
       type = "l",
       xlab = "Iteration",
       ylab = "Probability",
       ylim = c(0, 1),
       xlim = c(1, numIter),
       axes = FALSE,
       ...)

  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), las = 1, cex.axis = 0.7)
  axis(1, at = c(1, pretty(1:numIter), numIter), cex.axis = 0.7)
  box()
  
  for(i in 2:numGroups){
    lines(1:numIter, groupProbs[, i], col = i, lwd = 2)
  }
  
  legend("topright",
         col = 1:numGroups,
         lty = 1,
         lwd = 1.5,
         legend = paste0("Group ", 1:numGroups),
         cex = 0.7,
         bty = "n")
  
  par(oldPar)
}