# function takes dictionary and string, and returns the original if no match was found
translate_ids <- 
  function(x, dict, value.col = 2, return.missing = TRUE){
    
    dict <- dict[,c(1,value.col)]
    colnames(dict) <- c("key","value")
    dict <- dict[dict$value != "",]
    
    b <- dict[match(x,dict$key),value.col]
    
    if(return.missing == TRUE){
      missing <- which(is.na(b), arr.ind = T)
      b[missing] <- x[missing]
    }
    
    return(b)
    
  }
