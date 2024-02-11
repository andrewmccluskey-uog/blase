#' gene_selection_matrix
#'
#' @concept gene-selection
#'
#' @param x Object to generate the gene selection matrix from
#' @param pseudotime_slot= The slot of pseudotime to use for ordering. Defaults to slingPseudotime_1".
#' @param target_matrix_size Matrix size to redimension to. Defaults to 1000.
#' @param n_cores Number of cores for the matrix redimensioning
#'
#' @return A plot of the genes and cells ordered by pseudotime
#' @export
#'
#' @examples
setGeneric(name = "gene_selection_matrix",
           signature = c("x"),
           def = function(x, ...)
           standardGeneric("gene_selection_matrix"))

#' @rdname gene_selection_matrix
#'
#' @concept gene-selection
#'
#' @import methods
#' @import rlang
#' @export
setMethod(
  f = "gene_selection_matrix",
  signature = c(x="Seurat"),
  definition = function(x, waves, genes=c(), pseudotime_slot="slingPseudotime_1", target_matrix_size=1000, n_cores=1){
    rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
    sce = Seurat::as.SingleCellExperiment(x)
    return(gene_selection_matrix(
      sce,
      waves,
      genes=genes,
      pseudotime_slot=pseudotime_slot,
      target_matrix_size=target_matrix_size,
      n_cores=n_cores))
  }
)

#' @rdname gene_selection_matrix
#'
#' @concept gene-selection
#'
#' @import methods
#' @import metR
#' @export
setMethod(
  f = "gene_selection_matrix",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, waves, genes=c(), pseudotime_slot="slingPseudotime_1", target_matrix_size=1000, n_cores=1){

    if ( !any(colnames(x@colData) == pseudotime_slot)) {
      stop(paste0("Pseudotime slot '", pseudotime_slot ,"' does not exist"))
    }

    # TODO removing genes and then sorting messes up the ordering of the graph, I think, because waves$phase has all the genes
    # First we need to subset only the requested genes
    if (length(genes) > 0) {
      # R passes parameters by value not reference so this is safe
      x = x[rownames(x) %in% genes,]
      waves = waves[rownames(waves) %in% genes,]
    }

    # Then get the normalised count matrix
    pseudotime = x@colData[[pseudotime_slot]]
    heatmap_counts = SingleCellExperiment::normcounts(x)[,order(pseudotime)]
    heatmap_counts = heatmap_counts[order(waves$phase),]

    small_heatmap_counts = redim_matrix(
      heatmap_counts,
      target_height = target_matrix_size,
      target_width = target_matrix_size,
      n_core=n_cores
    )

    heatmap_counts_ordered.df <- reshape2::melt(small_heatmap_counts, c("gene", "cell"), value.name = "expression")
    plot = ggplot2::ggplot(data=heatmap_counts_ordered.df,ggplot2::aes(x=cell,y=gene,fill=expression)) +
      ggplot2::geom_tile() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()) +
      ggplot2::scale_fill_gradient(low = "white",
                                   high = "red",
                                   guide = "colorbar")

    return(plot)
  }
)

#' select_genes_by_fourier_method
#'
#' @param x
#' @param ...
#'
#' @concept gene-selection
#'
#' @return
#' @export
#'
#' @examples
setGeneric(name = "select_genes_by_fourier_method",
           signature = c("x"),
           def = function(x, ...)
             standardGeneric("select_genes_by_fourier_method"))

#' select_genes_by_fourier_method
#'
#' @param x Seurat.
#'
#' @concept gene-selection
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "select_genes_by_fourier_method",
  signature = c(x="Seurat"),
  definition = function(x, waves, n_genes=100, n_groups=40, top_n_per_group=1, method="power", force_spread_selection=TRUE, pseudotime_slot="slingPseudotime_1"){
    rlang::check_installed("Seurat", reason = "to handle Seurat objects.")
    sce = Seurat::as.SingleCellExperiment(x)
    return(select_genes_by_fourier_method(sce, waves, n_genes, n_groups, top_n_per_group, method, force_spread_selection, pseudotime_slot))
  }
)

#' select_genes_by_fourier_method
#'
#' @param x SingleCellExperiment.
#'
#' @import metR
#'
#' @concept gene-selection
#'
#' @return
#'
#' @export
#'
#' @examples
setMethod(
  f = "select_genes_by_fourier_method",
  signature = c(x="SingleCellExperiment"),
  definition = function(x, waves, n_genes=100, n_groups=40, top_n_per_group=1, method="power", force_spread_selection=TRUE, pseudotime_slot="slingPseudotime_1"){

    if(method != "power" & method != "amplitude" & method != "r2") {
      stop("Requested method is not valud, must be one of ['power','amplitude','r2']")
    }

    if (force_spread_selection) {

      if (n_genes != 100) {
        warning("n_genes is not used when force_spread_selection==TRUE")
      }

      best_waves_in_spread = data.frame()
      stepsize = max(waves$phase)/n_groups
      for(i in seq(from=0, to=max(waves$phase), length.out=n_groups)) {
        waves_in_block = waves[waves$phase>i-stepsize & waves$phase<i+stepsize,]
        # remove genes in `best_waves_in_spread` from `waves_in_block` and then select best
        best_wave_in_block = waves_in_block[order(-waves_in_block[,method]),][1:top_n_per_group,]
        best_wave_in_block = best_wave_in_block[!(rownames(best_wave_in_block) %in% best_waves_in_spread)]
        best_waves_in_spread = rbind(best_waves_in_spread, best_wave_in_block)
      }
      return(best_waves_in_spread)

    } else {
      if (n_groups != 10 | top_n_per_group != 1) {
        warning("n_groups and top_n_per_group are not used when force_spread_selection==FALSE")
      }

      top_waves = waves[order(-waves$r2),][0:n_genes,]
      return(top_waves)
    }

  }
)


#' get_waves
#'
#' @param sce
#' @param pseudotime
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
get_waves <- function(
    sce,
    pseudotime_slot,
    n_cores
) {

  if ( !any(colnames(sce@colData) == pseudotime_slot)) {
    stop(paste0("Pseudotime slot '", pseudotime_slot ,"' does not exist"))
  }
  pseudotime = sce@colData[[pseudotime_slot]]

  heatmap_counts = SingleCellExperiment::normcounts(sce)[,order(pseudotime)]

  # TODO AM rewrite with https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#single-machine
  waves_list = parallel::mclapply(rownames(heatmap_counts), function (gene) {
    wave = as.data.frame(FitWave(as.matrix(heatmap_counts[gene,]), 1))
    rownames(wave) = c(gene)
    return(wave)
  }, mc.cores=n_cores)
  waves <- do.call("rbind", waves_list)

  waves = as.data.frame(waves)
  waves$phase = waves$phase * (((180/3.141593)/360)*max(pseudotime)) # pseudotime
  waves$gene = rownames(waves)

  # Add power to waves
  waves$total_expression = rowSums(heatmap_counts[rownames(waves),])
  # Bozdech et al. use plus or minus 1/48 (i.e. one bulk either side of the peak, 6.25% window)
  # Here we use plus or minus 5% (i.e. a 10% window)
  five_percent_of_pdt = 0.05*max(waves$phase)
  waves$peak_expression = 0
  for (gene in rownames(waves)) {
    waves[gene,"peak_expression"] = sum(heatmap_counts[gene,pseudotime > waves[gene,]$phase-five_percent_of_pdt & pseudotime < waves[gene,]$phase+five_percent_of_pdt])
    waves[gene,"cellcount_in_peak"] = length(heatmap_counts[gene,pseudotime > waves[gene,]$phase-five_percent_of_pdt & pseudotime < waves[gene,]$phase+five_percent_of_pdt])
  }
  waves$power = waves$amplitude / waves$peak_expression

  return(waves)
}

#' redim_matrix
#'
#' @param mat The matrix to reduce
#' @param target_height target number of rows in the matrix. Default 100.
#' @param target_width target number of columns in the matrix. Default 100.
#' @param summary_func How to reduce the cells. Defaults to the mean.
#' @param n_core Number of cores to use in the redimensioning.
#'
#' @concept gene-selection
#'
#' @return The matrix with reduced dimensions
#'
#' @description
#' Snippet from: https://cansnippet.bioinfo-fr.net/details.php?id=3
redim_matrix <- function(
    mat,
    target_height = 100,
    target_width = 100,
    summary_func = function(x) mean(x, na.rm = TRUE),
    n_core = 1
) {

  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }

  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))

  # complicate way to write a double for loop
  # TODO AM rewrite with https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#single-machine
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, 0.0)
  }, mc.cores = n_core))
}
