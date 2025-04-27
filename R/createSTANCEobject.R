#' @title Create the STANCE object.
#' @param counts original gene expression matrix of dimension n x G,
#' with genes in rows and spots in columns.
#' @param pos the n x 2 dimensional location coordinates matrix.
#' The rownames of pos matrix should match the colnames of counts matrix.
#' @param prop the cell type proportion of each spot (n x K matrix)
#' Rows represent spots and columns represent cell types.
#' The rownames of prop matrix should match the colnames of counts matrix.
#' @param covariates (default NULL) the covariate (if any) design matrix of dimension n x p,
#' modeling gene expression.
#' The rownames of covariates matrix should match the colnames of counts matrix.
#' @return returns STANCE object.
#' @import methods
#' @export


creatSTANCEobject <- function(counts, pos, prop, covariates = NULL){
  if (!is.matrix(counts)) {
    tryCatch({
      counts <- as.matrix(counts)
    }, error = function(e) {
      stop('\'counts\' could not be converted to matrix. Please check that \'counts\' is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.')
    })
  }

  if (!is.matrix(prop)) {
    tryCatch({
      prop <- as.matrix(prop)
    }, error = function(e) {
      stop('\'prop\' could not be converted to matrix. Please check that \'prop\' is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.')
    })
  }

  if (!is.matrix(pos)) {
    tryCatch({
      pos <- as.matrix(pos)
    }, error = function(e) {
      stop('\'pos\' could not be converted to matrix. Please check that \'pos\' is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.')
    })
  }

  genes <- rownames(counts)
  if(is.null(genes)) {
    stop('\'rownames(counts)\' is null!')
  }

  barcodes <- colnames(counts)
  if(is.null(barcodes)) {
    stop('\'colnames(counts)\' is null!')
  }

  spots <- intersect(intersect(barcodes, rownames(pos)), rownames(prop))
  pos.use <- pos[match(spots, rownames(pos)),]
  prop.use <- prop[match(spots, rownames(prop)),]
  counts.use <- counts[,match(spots, colnames(counts))]
  if(!is.null(covariates)){
    covariates.use <- covariates[match(spots, rownames(covariates)),]
  }else{
    covariates.use <- NULL
  }

  object <- methods::new(
    Class = "STANCE",
    gene_expression = counts.use,
    location = pos.use,
    cell_type_compositions = prop.use,
    covariates = covariates.use,
    cell_types = colnames(prop.use)
  )
  return(object)
}
