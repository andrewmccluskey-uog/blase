#' Atgnat Mapping Result
#'
#' Created by [map_best_bin()]
#'
#' @concept mapping
#'
#' @slot bulk_name The name of the bulk sample being mapped.
#' @slot best_bin The bin that best matched the bulk sample.
#' @slot best_correlation The spearman's rho that the test geneset had between the winning bin and the bulk.
#' @slot top_2_distance The absolute difference between the best and second best mapping buckets. Higher indicates a less doubtful mapping.
#' @slot history A dataframe of the correlation score for each bin.
#'
#' @return A [MappingResult] object
#' @export
#'
#' @seealso [map_best_bin()]
#'
#' @examples
#' counts_matrix <- matrix(c(seq_len(120)/10, seq_len(120)/5), ncol=48, nrow=5)
#' sce <- SingleCellExperiment::SingleCellExperiment(assays=list(
#'   normcounts=counts_matrix, logcounts=log(counts_matrix)))
#' colnames(sce) = seq_len(48)
#' rownames(sce) = as.character(seq_len(5))
#' sce$cell_type = c(rep("celltype_1", 24), rep("celltype_2", 24))
#'
#' sce$pseudotime = seq_len(48)
#' atgnat_data = as.AtgnatData(sce, pseudotime_slot="pseudotime", n_bins=4)
#' atgnat_data@genes = as.character(seq_len(5))
#'
#' bulk_counts = matrix(seq_len(15)*10, ncol=3, nrow=5)
#' colnames(bulk_counts) = c("A", "B", "C")
#' rownames(bulk_counts) = as.character(seq_len(5))
#'
#' # Map to bin
#' result = map_best_bin(atgnat_data, "B", bulk_counts)
#' result
#'
#' # Plot Heatmap
#' plot_mapping_result_heatmap(list(result))
#'
#' # Plot bin
#' sce = scater::runUMAP(sce)
#' sce = assign_pseudotime_bins(sce, pseudotime_slot="pseudotime", n_bins=4)
#' plot_mapping_result(sce, result, group_by_slot="cell_type")
#'
MappingResult = setClass(
  Class = "MappingResult",
  slots = list(
    bulk_name = "ANY",
    best_bin = "numeric",
    best_correlation = "numeric",
    top_2_distance = "numeric",
    history = "data.frame"
  )
)

#' @title Show an MappingResult object
#'
#' @concept mapping
#'
#' @param object an [MappingResult] object
#' @export
#' @inherit MappingResult-class examples
setMethod(f = "show",
          signature = "MappingResult",
          definition = function(object){

            output = c(
              paste0("MappingResult for '", object@bulk_name, "':",
                    " best_bin=", object@best_bin,
                     " correlation=", object@best_correlation,
                     " top_2_distance=", object@top_2_distance),
              paste("\t with history for scores against", nrow(object@history), " bins\n")
            )

            cat(paste(output, collapse = '\n'))
          })
