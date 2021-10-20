# function to create the text equation
lmp_func <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  r2 <- summary(modelobject)$r.squared
  attributes(p) <- NULL
  return(list(p, r2))
}