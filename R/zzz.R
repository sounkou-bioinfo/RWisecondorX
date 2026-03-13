.onLoad <- function(libname, pkgname) {
  reticulate::py_require(
    c(
      "numpy",
      "pandas",
      "pysam",
      "scipy",
      "scikit-learn<=1.4.2"
    )
  )
}
