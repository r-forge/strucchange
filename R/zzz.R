.First.lib <- function(lib, pkg) {
  if(as.numeric(R.Version()$minor) < 7) {
    autoload("confint", "MASS")
  }
}
if(!require(weave)) warning("Could not load package weave")
