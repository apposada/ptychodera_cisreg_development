category_by_quantile <- function(
  x, newvalues = c(1,2,3,4), probs = seq(0,1,0.25), ...
) {
  
  q <- quantile(x, probs = probs)
  
  c <- cut(x, breaks=q, labels=FALSE,include.lowest = TRUE, ...)
  
  d <- data.frame(
    quantile = sort(unique(c)),
    newvalue = newvalues
  )
  
  y <- d$newvalue[match(c,d$quantile)]
  
  return(y)
  
}