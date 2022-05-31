# function to drop slices that sum to zero in specified dimension 

# x is an array of 2 or 3 dimensions (may work for higher dim arrays)

dim.drop0 <- function(x, dropdim) {
  sums <- apply(x, MARGIN = dropdim, sum)
  # returns e.g. for dropdim=1:  x[sums > 0, , ]
  return(eval(parse(text=paste0("x[", paste0(rep(",",dropdim-1), collapse=""), "sums>0", paste0(rep(",", length(dim(x))-dropdim),collapse=""), "]", collapse=""))))
}

array2binom <- function(x) {
  x[x>1] <- 1
  return(x)
}