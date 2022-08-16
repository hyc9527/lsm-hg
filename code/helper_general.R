##########################################################
#### Part 0: Basic math and stat functions
##########################################################
## math Function 0.1
norm.l2 = function(x){sqrt(crossprod(x))[1]}

## math Function 0.2
# Calculate area of param.alpha pbdim triangle, p = 1, 2, 3,...
size.triangle = function(index.vec, feat.mat){
    if (length(index.vec) > 3 | length(index.vec)  < 3)
    {
        stop("Input is not an triangle")
    } else {
        if (is.null(dim(feat.mat))){ # 1-dim volume is range of feat.mat
            return( max(feat.mat) - min(feat.mat))
        } else {
            mat = feat.mat[index.vec,]
            a =norm.l2(mat[1,]- mat[2,])
            b = norm.l2(mat[1,]- mat[3,])
            c = norm.l2(mat[2,]- mat[3,])
            if (abs(a) < 1e-10 || abs(b) < 1e-10 || abs(c) < 1e-10){
                return(0)
            } else {
                s = 0.5*(a +b + c)
                return( sqrt(s*(s - a)*(s - b)*(s - c))[1])
            }
        }
    }
}

## math Function 0.3
#  Generate multivariate gaussian
# Dependence: mvnorm() from Rpackage 'mnormt'
get.gaussian.featmat = function(nrow, ncol, mean, var1 = 1, cov12 = 0){
    if(ncol <  1+1e-10 ){
        result = rnorm(nrow, mean = mean[1], sd = sqrt(var1))
    } else{
        sigma.mat = diag(rep(var1,ncol)) + matrix(rep(cov12,ncol^2),ncol,ncol) - diag(rep(cov12,ncol))
        if ( length(mean) < 1 + 1e-10){
            return(mvrnorm(nrow, rep(mean,ncol), sigma.mat))
        }   else {
            if (length(mean) != ncol){stop("mean 's dimension is not correct! ")}
            return(mvrnorm(nrow, mean, sigma.mat))
        }
    }
}

## util Function 0.4
#  Print message  with time
msg <- function(s, num=NULL){
    time <- format(Sys.time(), "%b-%d-%H-%M-%S")
    if (is.null(num)){
        cat(sprintf("%s | %s \n", time, s))
    } else {
        cat(sprintf("%s | %s %.2f\n", time, s, num))
    }

}

## other Function 0.5:  turn.list.to.matrix(pylist:list[i]) -> mat[i]
# Convert a list from python function into a mat in R.  Each cell of a list contains a ROW array of same length.
# Usage:  Input pylist is usually a list output from python function, reticulated in R
turn.list.to.matrix = function(pylist){
    return( matrix(sapply(1:length(pylist), function(i){ as.numeric(pylist[[i]])}), ncol = length(pylist[[1]]), byrow = T))
}

## comb Function 0.6:  subsets_using_reticulate(n:int, k:int) -> mat[int]
# a R wrapper for a python function, which find all k-set out of [n], stored in a matrix of size choose(n,k) by k.
# Dependence: python function sunset_finder(), Rpackage reticulate, a R function turn.list.to.matrix()
#subsets_using_reticulate = function(n ,k ){
#turn.list.to.matrix(subset_finder(x = n,k = k))
#}


## comb Function 0.6':  subsets_using_reticulate(n:int, k:int) -> mat[int]
subsets_using_Rcpp = function(n ,k ){
    find_subset(x = n,k = k)
}



## math Function 0.7
# safe NATRUAL log function: smoothly handles zero and infinity
safelog <- function (x) {
    safelog.f <- function (x)
        if (x == Inf)
            Inf
        else if (x == 0)
            -1e300
        else
            log(x)
    if (length(x) == 1)
        safelog.f(x)
    else
        sapply(x, safelog.f)
}

## math Function 0.8 : assign values to vector
# https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line/7523486#7523486
`%tin%` <- function(x, y) {
    mapply(assign, as.character(substitute(x)[-1]), y,
    MoreArgs = list(envir = parent.frame()))
    invisible()
}

## math Function 0.9
# Logic operator: not belong to
# https://www.r-bloggers.com/the-notin-operator/
`%notin%` <- Negate(`%in%`)

## math Function 0.10: tv.dist(f1:array[int], f2:array[int])-> float
# Calculate Total variational distance of two frequency
# Usage: two input vectors, of not necessarily same length. Then zero(s) will be assigned to missing location(s) in shorter array.
tv.dist = function(f1, f2){
    f1 = f1/sum(f1)
    f2 = f2/sum(f2)
    if (length(f1) > length(f2)){
        #msg('Warning: two distribution vectors are of different lengths when calculating t.v.! ')
        f2 = c(f2, rep(0, length(f1) - length(f2)))
        return(0.5*sum(abs(f1 - f2)))
    } else if (length(f1) < length(f2)){
        #msg('Warning: two distribution vectors are of different lengths when calculating t.v.!')
        f1 = c(f1, rep(0, length(f2) - length(f1)))
        return(0.5*sum(abs(f1 - f2)))
    } else {
        return(0.5*sum(abs(f1 - f2)))
    }
}



## comb Function 0.12:  get.perm_using_reticulate(x:array[int]) -> mat[int]
# Find all permutaions of a vector x
# pyfunction perms_finder(), Rpackage reticulate
get.perm_using_reticulate = function(x){
    turn.list.to.matrix(perms_finder(x))
}


## math Function 0.13
len = function(x){
    length(x)
}






## Global Parameter
# machine precision
machine.precision  =  1e-300
{ # pallete
# Dependency:  function 0.9 '%tin%'
                   c(dark.blue  , light.blue  ,  dark.green, light.green, dark.yellow, light.yellow) %tin%
                c('#1f78b4'  , '#a6cee3'   , '#1b9e77'  , '#b2df8a'  , '#d95f02'  , '#fff7bc'   )

    c(dark.purple, light.purple, dark.red   , light.red  , dark.grey  , light.grey  )   %tin%
                c('#7570b3'  , '#beaed4'   , '#fe4365'  , '#fc9d9a'  , '#363636'  , '#A8A7A7'   )

    # https://medialab.github.io/iwanthue/
    warm.col = c('#db38c3', '#df4d22', '#d140a9', '#d1483b', '#db71b3', '#b44b45', '#8e2e6b', '#e13978', '#98304f', '#dd7798')
    cold.col = c("#3570c5", "#4f83f7", "#0f5eb0", "#4691eb", "#135ac2", "#0f5eb0", "#3a6cd7", "#0f5eb0", "#0f5eb0", "#0f5eb0")
}



####################################################################################
### Part I  Latent space volume model of general hypergraphs : sampling
#####################################################################################

## mod Function 1.1:  exact.volume.cvxhull(subset.vec:array[int], feat.mat: mat[float]) -> float
# calculate EXACT volume for k-node p-dim convex hull ONLY for k = 1 , 2 and 3 and any int p.
exact.volume.cvxhull = function(subset.vec,  feat.mat){
    k = length(subset.vec)
    if (k > 3){
        stop("length of first argument shoudnt be more than 3!")
    }
    if (k == 1){
        result = 0
    } else if( k == 2){
        result = norm.l2(feat.mat[subset.vec[1],] - feat.mat[subset.vec[2],])
    } else {
        result = size.triangle(subset.vec, feat.mat)
    }
    result
}

## mod Function 1.2: volume.cvxhull(subset.vec:array[n], feat.mat:mat[f]) -> float
# Calculating volume for a convex body spaned by k points in R^p , p = 1,2,3.....
# Dependency:  function convhulln() from Rpackage 'geomery'
volume.cvxhull = function(subset.vec,feat.mat){
    k = length(subset.vec)
    if (is.null(dim(feat.mat)[1])){ # 1-dim case, vol is range
        result = max(feat.mat) - min(feat.mat)
    } else if(k < 3 + 1e-10){
            result = exact.volume.cvxhull(subset.vec, feat.mat)
    } else {
            result = 0
            mat.select = feat.mat[subset.vec,]
            try({convhull = convhulln(mat.select , "FA"); result = convhull[[3]]})
    }
    result
}


## mod Function 1.2': volume.cvxhull(subset.vec:array[n], feat.mat:mat[f]) -> float
# Calculate volume FOR k = 4,5,...
# Dependency:  function convhulln() from Rpackage 'geomery'
volume.cvxhull.4pt.plus = function(subset.vec,feat.mat){
    result = 0
    mat.select = feat.mat[subset.vec,]
    try({convhull = convhulln(mat.select , "FA"); result = convhull[[3]]})
    result
}




## mod Function 1.3: dist.shape.cvxhull(subset.mat: mat[int], feat.mat: mat[f], which.kernel:str, param.vec: 2-array[f])
# Calculate distribution of convex hull shape model
# Usage: 'which.kernel' \in ['gamma','gaussian','pareto','window']
# Usage: 'param.vec' includes
#                             "gamma": param.vec = c(\alpha, \beta), see definition @ https://en.wikipedia.org/wiki/Gamma_distribution
#                             "gaussian": param.vec = c(\mu, \sigma)
#                             "pareto": param.vec = c(\alpha, x_m),  see denifinition @ https://en.wikipedia.org/wiki/Pareto_distribution
#                             "window" : param.vec = c(lower.bound , upper.bound)
#                              "identity": param.vec = NA
# Dependency: function volume.cvxhull()
dist.shape.cvxhull = function(subset.mat, feat.mat, which.kernel, param.vec){
    switch(which.kernel,
        "gamma" = { shape.kern = function(d, alpha =  param.vec[1], beta = param.vec[2]){
                        if (d > machine.precision){
                                if(d^(alpha - 1) == Inf){
                                        msg('Inf is coerced to a 1e300 in calculating pareto "shape.kern"')
                                        return(1/machine.precision)
                                } else {
                                        return(d^(alpha-1)*exp(-beta*d))
                                }
                        } else {
                                return(0)
                        }
                    }
        },
        "gaussian" ={   shape.kern = function(d, mu =  param.vec[1], sd = param.vec[2]){
                            if (sd < 10*machine.precision){
                                msg('Warning in calculating gaussian shape.kern: parameter sd is too small.')
                            }
                            if (sd > sqrt(machine.precision)){
                                return( exp(-(d - mu)^2/(2*sd^2)))
                            } else {
                                0
                            }
                            }
        },
        "pareto" = { shape.kern = function(d, alpha =  param.vec[1], beta = param.vec[2]){
                    if (d > beta + machine.precision & alpha > machine.precision & beta > machine.precision){
                        return(d^(-alpha-1))
                    } else {
                        0
                    }
            }
        },
        "window" = { shape.kern = function(d, param.lbd = param.vec[1], param.ubd = param.vec[2]){
                if (param.lbd > param.ubd - machine.precision){
                    stop('window parameters error: lower bound is larger than upper bound. ')
                }
                if (d >  param.lbd + machine.precision & d < param.ubd - 1e-100 ){
                    return(1)
                } else {
                    0
                }
            }
        },
        "identity" = { shape.kern = function(d){
                            d
                        }
        },
        { msg("Shape type is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
            stop()
        }
    )
    n.k = dim(subset.mat)[1]
    dist = sapply(1:n.k, function(i) shape.kern(volume.cvxhull(subset.mat[i,], feat.mat)))
    if(sum(is.na(dist)) > 0.5){ # when kernel function return nan, warning message pops up
        msg('Warning by function "dist.shape.cvxhull": value(s) are  missing in unnormalized distribution')
    }
    if (sum(dist) < machine.precision){ # when kernel function return zero for every k-set, warning message pops up
        msg("Warining by function 'dist.shape.cvxhull' when calculating distribution: zero vector is returned  ")
        return(dist)
    } else {
        return(dist/sum(dist))
    }
}

## mod Function 1.4: sample.from.kset(kset.mat: mat[int], prob:array[float], num.sample:int) -> mat[int]
# Sample multiple k-sets using a discrete probability vector 'prob'
sample.from.kset = function(kset.mat, prob, num.sample){
    n.k = dim(kset.mat)[1]
    if ( n.k != length(prob)){stop("distribution doesn't match subset list !")}
    if ( sum(prob) < machine.precision){
        res = NA
        msg('Error: emmpy is returned by sample.from.kset')
    } else {
        res = sapply(1:num.sample, function(i) kset.mat[sample(n.k, 1, replace = F, prob),])
    }
    return(t(res))
}



# comb Function 2.4
# https://stat.ethz.ch/pipermail/r-help/2011-April/276021.html
rowmatch <- function(A,B) {
    f <- function(...) paste(..., sep=":")
    if(!is.matrix(B)) B <- matrix(B, 1, length(B))
    distance.mat <- do.call("f", as.data.frame(A))
    b <- do.call("f", as.data.frame(B))
    match(b, distance.mat)
}



# mod Function 2.3  turn.kset.sample.to.count(edge.mat:array[int] or mat[int], subset:mat[int])->array[int] ?????? naming error
# Turn k-set samples to  discrete count data Y_i
# Dependence: rowmatch
turn.kset.sample.to.count = function(edges, subset){
    if(is.null(dim(edges))){
        if(length(edges) != length(subset[1,])){
            msg('Error: subset doesnt match with edge.mat!!!')
            stop()
        }
    } else {
        if(length(edges[1,]) != length(subset[1,])){
        msg('Error: subset doesnt match with edge.mat!!!')
        stop()
        }
    }
    count = rep(0, dim(subset)[1])
    for (ind in rowmatch(subset, edges)){
        count[ind] = count[ind]  + 1
    }
    count
}



## mod Function 1.5  sample.from.shape.cvxhull( n:int, k:int,
#                                          subset.mat:mat[i], feat.mat:mat[f], num.sample:int,
#                                          which.shape:str, param.vec:2-arrary[f],
#                                          if.save.sample:logical, save.path:str)
#                 -> list(count.vec:array[int],
#                        feature.mat:mat[f], subset.mat:mat[i],
#                        edge.sampled:mat[int])
# Usage:
#       'which.shape' \in c('gamma','gaussian','pareto','window')
#        save.path = './output'
sample.from.shape.cvxhull = function(n, k = 2, subset.mat = NULL, feat.mat = NULL, num.sample, which.shape, param.vec, if.save.sample, save.path){
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = 100, cov12 = 0)
        #msg('feat.mat is not defined and thus a new feature matrix is generated.')
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
        #msg('subset.mat is not found and thus is calculated.')
    }
    dist  = dist.shape.cvxhull( subset.mat = subset.mat,
                                feat.mat =feat.mat,
                                which.kernel = which.shape,
                                param.vec = param.vec)
    if (sum(dist) > machine.precision){ # dist is a valid distribution: (1) nonnegative (2) sum(dist) = 1
            edges.mat = sample.from.kset(kset.mat = subset.mat, prob = dist, num.sample) # in case of mixture of different length of k-set, k \in [n]
            #edges.mat = turn.list.to.matrix(edges.list) # when edge is of same length
            sim.data = list(count.vec = turn.kset.sample.to.count(edges.mat, subset.mat),
                            feature.mat = feat.mat,
                            subset.mat = subset.mat,
                            edge.sampled = edges.mat)
    } else {
        sim.data = list(count.vec = rep(0,dim(subset.mat)[1]),
                        feature.mat = feat.mat,
                        subset.mat = subset.mat,
                        edge.sampled = NA)
    }
    if(if.save.sample){
        #current.time = format(Sys.time(),"%b-%d-%I")
        save(sim.data, file = paste0(save.path, "samples-from-",which.shape,"-n", n, "-a",param.vec[1],'-b',param.vec[2], ".Rdata"))
    }
    return(sim.data)
}





####################################################################################################################
#### Hypergraph
####################################################################################################################



# Function 2.8
collect.kset.from.authorPaperBinaryAdj = function( max.set.k, authorPaperBinaryAdj.mat){
    res = vector("list", length = max.set.k - 1)
    names(res) = paste(2:max.set.k, 'set', sep = '-')
    for (i in 1:dim(authorPaperBinaryAdj.mat)[2]){
        if( sum(authorPaperBinaryAdj.mat[,i]) > 1 + 1e-1 && sum(authorPaperBinaryAdj.mat[,i]) < max.set.k + 1e-1){
            res[[sum(authorPaperBinaryAdj.mat[,i]) - 1]] = rbind(res[[sum(authorPaperBinaryAdj.mat[,i]) - 1]], which(authorPaperBinaryAdj.mat[,i] > 1e-1))
        }
    }
    res
}


# Function 2.9 Create graph obj from hypergraph
#   Dependance: R package 'igraph'
#          k = 2,3,4 ....
#             (i)  k = 2 graph: convert 2-edge into weighted graph
#              (ii) k = 3 hypergraph:  convert 3-edge to a triangle
#               (iii) k \geq 4 hypergraph: convert k-set to a polygon, rather than a complete clique.
# 							note that degree of node for k \geq 4 is not true, only for better visualization.
turn.coauthorship.to.wt.graph = function(mat.lst, total.author){
    for( which.kset in 1:len(mat.lst)){
        k = dim(mat.lst[[which.kset]])[2]
        if (k == 2){
            path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
            path = path.uniq.wt[,-(k+1)]
            wt = path.uniq.wt[,(k+1)]
            wt.graph <- graph(as.vector(t(path[,-(k+1)])), n = total.author, directed=F)
            E(wt.graph)$weight <- wt
            E(wt.graph)$color <- "gray85"
            E(wt.graph)$width <- 2*E(wt.graph)$weight
            E(wt.graph)$type = paste0('k',k)
            E(wt.graph)$lty = 6
        } else if (k == 3) {
            path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
            path = path.uniq.wt[,-(k+1)]
            path = cbind(path, path[,1]) # turn 3-edge to triangle
            wt = path.uniq.wt[,(k+1)]
            if (!exists('wt.graph')){
                wt.graph = graph.empty(n=total.author, directed=F)
            }
            if (is.null(dim(path))){
                wt.graph = wt.graph + path(as.vector(path), color =  cold.col[1], weight = wt[i],  width = 1* wt[i], type = paste0('k',k), lty = 1)
            } else {
                for ( i in 1:dim(path)[1]){
                wt.graph = wt.graph + path(path[i,], color =  cold.col[1], weight = wt[i],  width = 1* wt[i], type = paste0('k',k), lty = 1)
            }
            }
        } else {
            path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
            path = path.uniq.wt[,-(k+1)]
            wt = path.uniq.wt[,(k+1)]
            if (!exists('wt.graph')){
                wt.graph = graph.empty(n=total.author, directed=F)
            }
            if (is.null(dim(path))){
                wt.graph = wt.graph + path(as.vector(path), color =  cold.col[1], weight = wt[i],  width = 1* wt[i], type = paste0('k',k), lty = 3)
            } else {
                for ( i in 1:dim(path)[1]){
                wt.graph = wt.graph + path(path[i,], color =  cold.col[1], weight = wt[i],  width = 1* wt[i], type = paste0('k',k), lty = 3)
            }
            }
        }
    }
    wt.graph
}

# Function 2.10 Turn a list of k-sets into a projected complete graph object
#   Dependance: R package 'igraph'
#          k = 2,3,4 ....
#               k2 graph: convert 2-edge into weighted graph
#               k3 hypergraph:  convert 3-edge to a triangle
#               k \geq 4 hypergraph: convert k-set to a complete polygon
# Note this function is to calculate deg distribution of a hypergraph, while func 2.19 is to create a wt graph for visualization

turn.coauthorship.to.wt.graph.complete = function(mat.lst, total.author){
    for( which.kset in 1:len(mat.lst)){
		if (is.null(dim(mat.lst[[which.kset]]))){
			k = len(mat.lst[[which.kset]])
		} else{
			k = dim(mat.lst[[which.kset]])[2]
		}
		#print(k)
        if (k == 2){
            path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
            path = path.uniq.wt[,-(k+1)]
            wt = path.uniq.wt[,(k+1)]
            wt.graph <- graph(as.vector(t(path[,-(k+1)])), n = total.author, directed=F)
            E(wt.graph)$weight <- wt
            E(wt.graph)$color <- "gray85"
            E(wt.graph)$width <- 2*E(wt.graph)$weight
            E(wt.graph)$type = paste0('k',k)
            E(wt.graph)$lty = 6
        } else if (k == 3) {
            path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
            path = path.uniq.wt[,-(k+1)]
            path = cbind(path, path[,1]) # turn 3-edge to triangle
            wt = path.uniq.wt[,(k+1)]
            if (!exists('wt.graph')){
                wt.graph = graph.empty(n=total.author, directed=F)
            }
            if (is.null(dim(path))){
                wt.graph = wt.graph + path(as.vector(path), color =  cold.col[1], weight = wt,  width = 1* wt[i], type = paste0('k',k), lty = 1)
            } else {
                for ( i in 1:dim(path)[1]){
                wt.graph = wt.graph + path(path[i,], color =  cold.col[1], weight = wt[i],  width = 1* wt[i], type = paste0('k',k), lty = 1)
            }
            }
        } else { # for k =4,5,6,...
			if (is.null(dim(mat.lst[[which.kset]]))){
				path = mat.lst[[which.kset]]
				wt = 1
			} else{
				path.uniq.wt = find.unique.path.weight(mat.lst[[which.kset]])
				path = path.uniq.wt[,-(k+1)]
				wt = path.uniq.wt[,(k+1)]
			}
			#if (!exists('wt.graph')){
            #    wt.graph = graph.empty(n=total.author, directed=F)
            #}
            if (is.null(dim(path))){ # single path
				polygon = path
				ind.subset.mat = find_subset(len(polygon),2)
				edge.to.add = ind.subset.mat
				for(h in 1: dim(ind.subset.mat)[1]){
					for (j in 1:2){
						edge.to.add[h,j] = polygon[ind.subset.mat[h,j]]
					}
				}
				edge.complete = edge.to.add
				#print(edge.complete)
                wt.graph = graph(as.vector(t(edge.complete)), n = total.author, directed=F)
            } else {
				edge.complete = c()
                for ( i in 1:dim(path)[1]){
					polygon = path[i,]
					ind.subset.mat = find_subset(len(polygon),2)
					edge.to.add = ind.subset.mat
					for(h in 1: dim(ind.subset.mat)[1]){
						for (j in 1:2){
							edge.to.add[h,j] = polygon[ind.subset.mat[h,j]]
						}
					}
					edge.complete = rbind(edge.complete, edge.to.add)
				}
				wt.graph <- graph(as.vector(t(edge.complete)), n = total.author, directed=F)
            }
            }
        }
    wt.graph
}








