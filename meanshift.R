library(doParallel)
library(foreach)
library(cluster)

cos_similarity_matrix <- function(A, B, k){ return((A %*% t(B)) * k) }
l2_norm <- function(X){ sqrt(rowSums(X^2)) }
l2_normalize <- function(X){ return(X / l2_norm(X)) }

cos_similarity_matrix_with_threshold <- function(A, B, k){
    return(pmax(A %*% t(B) - k, 0))
}

estimate_k <- function(X){
    X <- l2_normalize(X)
    newk <- cos(mean(acos(pmin(1, cos_similarity_matrix(X, X, 1)))) / 2.0)
    return(newk)
}

next_centers <- function(X, C, k){
    stopifnot(ncol(X) == ncol(C))
    #stopifnot(all(dim(X) == dim(C)))
    # W <- cos_similarity_matrix(C, X, k)
    # W <- exp(W - apply(W, 1, max))
    W <- cos_similarity_matrix_with_threshold(C, X, k)
    return(l2_normalize(W %*% X))
}

merge_centers <- function(X, C, merge_threshold){
    stopifnot(all(dim(X) == dim(C)))

    S <- cos_similarity_matrix(C, C, 1) >= merge_threshold
    #not_assigned <- order(rowSums(S >= merge_threshold), decreasing=TRUE)
    not_assigned <- order(rowSums(S), decreasing=TRUE)
    #S <- S[not_assigned, not_assigned]
    
    cluster_labels <- rep(0, length(not_assigned))
    centers <- matrix(0, nrow(X), ncol(X))
    for(i in seq_len(nrow(X))){
        cur_idx <- not_assigned[1]
        # cur_idx <- not_assigned[order(rowSums(S[not_assigned,,drop=FALSE] >= merge_threshold, na.rm=TRUE), decreasing=TRUE)][1]
        centers[i, ] <- C[cur_idx, ]
        #cur_assign_idx <- which(S[cur_idx,] >= merge_threshold)
        cur_assign_idx <- which(S[cur_idx,])
        stopifnot(all(cluster_labels[cur_assign_idx] == 0))
        stopifnot(cur_idx %in% cur_assign_idx)
        cluster_labels[cur_assign_idx] <- i
        S[,cur_assign_idx] <- NA
        not_assigned <- setdiff(not_assigned, cur_assign_idx)
        if(length(not_assigned) == 0){ break }
        
        # cur_idx <- order(rowSums(S))[1]
        # centers[i,] <- C[not_assigned[cur_idx],]
        # cur_flg <- S[cur_idx,]
        # cur_assign_idx <- not_assigned[cur_flg]        
        # stopifnot(all(cluster_labels[cur_assign_idx] == 0))
        # cluster_labels[cur_assign_idx] <- i
        # S <- S[!cur_flg, !cur_flg, drop=FALSE]
        # not_assigned <- not_assigned[!cur_flg]
        # if(length(not_assigned) == 0){ break }
    }
    centers <- centers[seq_len(i),, drop=FALSE]
    stopifnot(all(cluster_labels >= 1))
    stopifnot(nrow(centers) == max(cluster_labels))
    return(list(
        centers=centers,
        labels=cluster_labels
    ))
}

combine_for_shifted_means <- function(...){
    tmp <- list(...)
    return(list(
        Y = do.call(rbind, lapply(tmp, function(x){x$Y})),
        converged = all(sapply(tmp, function(x){x$converged}))
    ))
}

shifted_means <- function(X, max_iter, k, convergence_threshold, n_parallel){
    `%cur_do%` <- ifelse(n_parallel > 1, `%dopar%`, `%do%`)
    registerDoParallel(n_parallel)
    on.exit(stopImplicitCluster())
    tmp <- foreach(
        i = seq_len(nrow(X)), 
        .verbose=FALSE,
        .inorder=TRUE, 
        .multicombine=TRUE,
        .combine=combine_for_shifted_means) %cur_do% {
        
        Y <- X[i, , drop=FALSE]
        for(i in seq_len(max_iter)){
            Y2 <- next_centers(X, Y, k)
            if(max(l2_norm(Y - Y2)) < convergence_threshold){
                return(list(Y=Y, converged=TRUE))
            }
            Y <- Y2
        }
        return(list(Y=Y, converged=FALSE))
    }
    stopifnot(all(dim(tmp$Y) == dim(X)))
    return(tmp)
}

ms_sphere <- function(
    X, 
    max_iter,
    k,
    convergence_threshold,
    merge_threshold,
    n_parallel=detectCores()
){
    stopifnot(nrow(X) > 0)
    stopifnot(ncol(X) > 0)
    stopifnot(max_iter > 0)
    stopifnot(k > 0)
    stopifnot(merge_threshold > 0)

    X <- l2_normalize(X)
    tmp <- shifted_means(X, max_iter, k, convergence_threshold, n_parallel)
    Y <- tmp$Y
    converged <- tmp$converged

    # Y <- X
    # convergec <- FALSE
    # for(i in seq_len(max_iter)){
    #     Y2 <- next_centers(X, Y, k)
    #     if(max(l2_norm(Y - Y2)) < convergence_threshold){
    #         converged <- TRUE
    #         break
    #     }
    #     Y <- Y2
    # }

    tmp <- merge_centers(X, Y, merge_threshold)
    return(list(
        cluster_labels = tmp$labels,
        cluster_centers = tmp$centers,
        shifted_means=Y,
        converged=converged
    ))
}

optimize_silhouette <- function(X, k_seq, ...){
    k_seq <- sort(unique(k_seq))
    dmatrix <- 1 - cos_similarity_matrix(X, X, 1)
    results <- list()
    sres <- list()
    mean_sindex <- c()
    for(i in seq_len(length(k_seq))){
        cur_k <- k_seq[i]
        results[[i]] <- ms_sphere(X, k=cur_k, ...)
        sres[[i]] <- silhouette(results[[i]]$cluster_labels, dmatrix=dmatrix)
        mean_sindex <- c(mean_sindex, mean(sres[[i]][,3]))
    }
    names(results) <- names(k_seq)
    df_sidx <- data.frame(k = k_seq, mean_silhouette=mean_sindex)
    best <- results[[which.max(df_sidx$mean_silhouette)]]
    return(list(
        results    = results,
        silhouette = df_sidx,
        best       = best
    ))
}
