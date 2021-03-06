#' Greatest common divisor
#'
#' @param n integer
#' @param m integer
#' @return integer
gcd = function(n, m) {
    ifelse (m == 0, n, gcd(m, n %% m))
}


#' Compute the similarity between raw contact matrices of different resolutions for the same sample
#'
#' @importFrom Matrix norm
#' @param mtx1 Raw contact matrix 1 (dense or sparse)
#' @param mtx2 Raw contact matrix 2 (dense or sparse)
#' @param res1 Resolution (i.e. bin size) of mtx1
#' @param res2 Resolution (i.e. bin size) of mtx2
#' @param method Matrix norm method. Default: Frobenius
#' @return 
matrix_res_similarity = function(mtx1, mtx2, res1, res2, method = "F") {
    # calculate lowest common multiple resolution
    lcm_res = as.integer(res1 * res2 / gcd(res1, res2))
    mtx1_scaled = reduce_resolution(mtx1, newreso = lcm_res, oldreso = res1)
    mtx2_scaled = reduce_resolution(mtx2, newreso = lcm_res, oldreso = res2)
    '%ni%' = Negate("%in%")
    if (dim(mtx1_scaled) != dim(mtx2_scaled)) {
        cat("mtx1 bins not in mtx2:\n")
        print(rownames(mtx1_scaled)[rownames(mtx1_scaled) %ni% rownames(mtx2_scaled)])
        cat("mtx2 bins not in mtx1:\n")
        print(rownames(mtx2_scaled)[rownames(mtx2_scaled) %ni% rownames(mtx1_scaled)])
    }
    return(norm(mtx1_scaled - mtx2_scaled, type = method))
}
