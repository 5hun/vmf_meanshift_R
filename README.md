
# Mean-Shift Clustering on Hypersphere

This repository contains a R script implements von Meses-Fisher (vMF) mean shift clustering [1]. <br> 
vMF mean shift clusteirng is a mean shift clustering on hypersphere (, i.e. this algorithm uses cosine distance).


```R
source("meanshift.R")

library(MeanShift)
library(mclust)
```

    Loading required package: foreach
    Loading required package: iterators
    Loading required package: parallel
    Loading required package: wavethresh
    Loading required package: MASS
    WaveThresh: R wavelet software, release 4.6.8, installed
    
    Copyright Guy Nason and others 1993-2016
    
    Note: nlevels has been renamed to nlevelsWT
    
    Package 'mclust' version 5.3
    Type 'citation("mclust")' for citing this R package in publications.


# Example: Unit Circle

## Data Creation


```R
set.seed(1)
N <- 500
M <- 8
true_labels <- sample(seq_len(M), N, replace=TRUE)
radians <- true_labels / M * 2 * pi + rnorm(N)*0.1
X <- cbind(cos(radians), sin(radians))
```

## Mean Shift clustering


```R
X_n <- t(X) / apply(t(X), 1, sd )
dist_X_n <- dist(t(X_n))
h.cand <- quantile(dist_X_n, seq( 0.05, 0.40, by=0.05 ) )
# Choose parameters via optimizing silhouette index.
suppressMessages({bms_result <- optimize_silhouette(X_n, data.frame(h=h.cand), dist=dist_X_n, cluster_func=function(X, h){ bmsClustering(X, h=h)})})
```

## vMF mean shift clustering


```R
# Choose parameters via optimizing silhouette index.
dmatrix <- 1 - cos_similarity_matrix(X, X, 1)
result <- optimize_silhouette(X, data.frame(k=c(0.5,0.6,0.7,0.8,0.9)),  cluster_func=ms_sphere, dmatrix=dmatrix, max_iter=100, convergence_threshold=1e-5, merge_threshold=0.99)
```

## Results

The results of vMF mean shift and normal mean shift almost the same. vMF is slightly closer to a true cluster.


```R
options(repr.plot.width=8, repr.plot.height=4)
par(mfrow=c(1,2))
plot(X, col=result$best$labels, xlab="1", ylab="2", main=sprintf("vMF mean shift\nadjusted Rand Index = %.3f", adjustedRandIndex(true_labels, result$best$labels)))
plot(X, col=bms_result$best$labels, xlab="1", ylab="2", main=sprintf("Mean Shift\nadjusted Rand Index = %.3f", adjustedRandIndex(true_labels, bms_result$best$labels)))
```


![png](output_10_0.png)


# Example: Unit Sphere

## Data Creation


```R
set.seed(3)
N <- 500
M <- 10
X0 <- l2_normalize(cbind(rnorm(M), rnorm(M), rnorm(M)))
true_labels <- sample(seq_len(M), N, replace=TRUE)
X0 <- X0[true_labels,]
radians_1 <- ifelse(X0[,1] > 0, 0, pi) + atan(X0[,2] / X0[,1]) + rnorm(N) * 0.2
radians_2 <- acos(X0[,3]) + rnorm(N) * 0.2
X <- cbind(
    cos(radians_1) * sin(radians_2),
    sin(radians_1) * sin(radians_2),
    cos(radians_2)
)
```

## Mean shift clustering


```R
X_n <- t(X) / apply(t(X), 1, sd )
dist_X_n <- dist(t(X_n))
h.cand <- quantile(dist_X_n, seq( 0.05, 0.40, by=0.05 ) )
# Choose parameters via optimizing silhouette index.
suppressMessages({bms_result <- optimize_silhouette(X_n, data.frame(h=h.cand), dist=dist_X_n, cluster_func=function(X, h){ bmsClustering(X, h=h)})})
```

## vMF mean shift clustering


```R
# Choose parameters via optimizing Silhouette index.
dmatrix <- 1 - cos_similarity_matrix(X, X, 1)
result <- optimize_silhouette(X, data.frame(k=seq(0.5, 0.99, 0.05)),  dmatrix=dmatrix, cluster_func=ms_sphere, max_iter=100, convergence_threshold=1e-5, merge_threshold=0.99)
```

## Results

The result of vMF mean shift is closer to the true clustering.


```R
options(repr.plot.width=9, repr.plot.height=6)
par(mfrow=c(2,3))

plot(X[,c(1,2)], col=result$best$labels, pch=result$best$labels, xlab="1", ylab="2", 
     main=sprintf("vMF mean shift\nadjusted Rand Index = %.3f", adjustedRandIndex(true_labels, result$best$labels)))
plot(X[,c(2,3)], col=result$best$labels, pch=result$best$labels, xlab="2", ylab="3")
plot(X[,c(1,3)], col=result$best$labels, pch=result$best$labels, xlab="1", ylab="3")

plot(X[,c(1,2)], col=bms_result$best$labels, pch=bms_result$best$labels, xlab="1", ylab="2", 
     main=sprintf("Mean Shift\nadjusted Rand Index = %.3f", adjustedRandIndex(true_labels, bms_result$best$labels)))
plot(X[,c(2,3)], col=bms_result$best$labels, pch=bms_result$best$labels, xlab="2", ylab="3")
plot(X[,c(1,3)], col=bms_result$best$labels, pch=bms_result$best$labels, xlab="1", ylab="3")
```


![png](output_19_0.png)


# Benchmark of Very High Dimension and Large Data


```R
set.seed(1)
N <- 10000
M <- 300
X <- l2_normalize(matrix(rnorm(N*M), N, M))
```


```R
system.time({
    result <- ms_sphere(X, max_iter=100, k=0.9, convergence_threshold=1e-5, merge_threshold=0.99, n_parallel = 8)
})
```


       user  system elapsed 
    490.714  74.078 109.935 


# References

1. T. Kobayashi and N. Otsu, Von Mises-Fisher mean shift for clustering on a hypersphere, ICPR 2010.
