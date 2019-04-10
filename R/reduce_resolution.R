#' Reduce resolution (incrrease bin size) of a HiC contact matrix
#'
#' This function takes a contact matrix and returns the corresponding contact matrix with the resolution reducction (increase of bin size)
#' @import magrittr
#' @import Matrix
#' @param mat HiC contact map matrix
#' @param newreso Desired bin size in bp
#' @return A HiC contact matrix with the desired resolution
#' @export
#' @examples
#' plot(0)

reduce_resolution <- function(mat, newreso, oldreso = NA){
    # if oldreso is not specified, infer it from the bin start positions
    if (is.na(oldreso)) {
        oldreso <- rownames(mat) %>%
            # split bin names by ":" (separate into chr and start position)
            strsplit(":") %>%
            # combine into matrix where 1st col = chr, 2nd col = start, row = bin
            do.call(rbind, .) %>%
            as.data.frame %>%
            setNames(c("chr", "pos")) %>%
            mutate(pos = as.integer(pos)) %>%
            group_by(chr) %>%
            arrange(pos) %>%
            # find minimum difference between bin start sites on the same chr
            summarize(m = diff(pos) %>% min) %>%
            ungroup  %>%
            # find minimum difference between bin start sites across all chr
            summarize(m = min(m)) %$%
            # return that smallest value
            m
    }

    # don't adjust if resolutions are the same
    if(oldreso == newreso) return(mat)

    # stop if new resolution is not a multiple of old resolution
    if((newreso %% oldreso) != 0) stop("New resolution should be multiple of the original one")
    
    # extract chromosomes and bin start positions
    ids <- rownames(mat)
    chrs <- gsub(":.*$", "", ids)
    pos <- gsub("^.*:", "", ids) %>% as.integer

    # calculate new start positions from new resolution
    newpos <- as.integer(floor(pos / newreso) * newreso)

    # generate new bin names
    #   these are ordered to match the old bin names, but will be repeated
    #   i.e. defines which new bin each row of mat will be mapped to
    newids <- paste(chrs, newpos, sep = ":")

    # generate new square matrix with newids as col and row names
    #   dim(m) = c(# of new bins, # of old bins)
    m <- fac2sparse(newids)

    # sum read counts in old bins into new bins
    #   dim(newmat) = c(# of new bins, # of new bins)
    newmat <- m %*% (mat %*% t(m))

    # find new bins with the same label as old bins
    bins <- rownames(mat)[rownames(mat) %in% rownames(newmat)]
    # match up new bins to ensure that newmat has the same proper ordering
    #   of bins that mat did
    i <- match(bins, rownames(newmat))

    # force to symmetric matrix since all HiC contact matrices are symmetric
    #   in correct bin ordering
    as(newmat[i,i], class(mat))
}
