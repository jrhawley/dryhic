#' Plot a desired region of a contact matrix using genomic coordinates
#'
#' @import data.table
#' @param mat HiC contact matrix (it could be the output of \code{\link{get_contacts_matrix}})
#' @param bins A \code{data.frame} containing chromosome, position and bin ID
#' @param coord A string of the form \code{chr:start-end} specifying the genomic region to plot
#' @param chrom Chromosome to plot. Priority over \code{coord},
#' @param start Starting coordinate. Priority over \code{coord},
#' @param end Ending coordinate. Priority over \code{coord},
#' @param ... Parameters passed onto \code{plot_matrix} function
#' @return returns
#' @export
plot_matrix_region = function(mat, bins, coord,
                chrom = strsplit(coord, "[:-]", perl = TRUE)[[1]][1],
                start = as.integer(strsplit(coord, "[:-]", perl = TRUE)[[1]][2]),
                end = as.integer(strsplit(coord, "[:-]", perl = TRUE)[[1]][3]),
                ...) {
    # check input parameters
    if (is.null(chrom) || is.null(start) || is.null(end)) {
        stop("`chrom`, `start`, and `end` cannot be NULL")
    }
    # ensure start < end
    if (start > end) {
        end = end + start
        start = end - start
        end = end - start
    }

    # get from.id and to.id of genomic position
    bins = as.data.table(bins)
    from_idx = bins[, which(chr == chrom && start >= pos && end > pos)[1]]
    to_idx = bins[, which(chr == chrom && start < pos && end <= pos)[1]]

    # send to plot_matrix
    return(plot_matrix(mat, c(from_idx, to_idx), ...))
}
