# Shared mclust namespace helpers used by both the NIPTeR and
# native WisecondorX layers.

.mclust_fit_in_namespace <- function(data,
                                     G = 2L,
                                     modelNames = NULL,
                                     equalPro = NULL) {
  ns <- asNamespace("mclust")
  env <- list2env(
    list(
      data = data,
      G = G,
      modelNames = modelNames,
      equalPro = equalPro
    ),
    parent = ns
  )

  if (is.null(modelNames) && is.null(equalPro)) {
    stop("At least one of 'modelNames' or 'equalPro' must be specified.",
         call. = FALSE)
  }

  if (is.null(equalPro)) {
    return(evalq(Mclust(data, G = G, modelNames = modelNames), envir = env))
  }

  if (is.null(modelNames)) {
    return(evalq(Mclust(data, G = G, control = emControl(equalPro = equalPro)),
                 envir = env))
  }

  evalq(
    Mclust(data, G = G, modelNames = modelNames,
           control = emControl(equalPro = equalPro)),
    envir = env
  )
}
