## mod Function 1.6 sequential.importance.sample.from.shape.cvxhull (n: int, k:int, feat.mat:mat[f], subset.mat:mat[i],
#                                                                    pi.init: n-array[f], which.node.importance:str, which.kernel:str,
#                                                                    param.vec:2-array[f], num.sample:int, if.save.si.sample:logical, save.path:str)
#                    -> list(count.vec:array[int],
#                            feature.mat:mat[f],
#                            subset.mat:mat[int],
#                            edge.sampled:mat[int])
# Sequential importance sampling algorithm: generate multiple samples from convex hull shape model
# Usaage:
#           pi.init is a discrete prob over [n]
#           which.node.importance \in c('small.chull','near.nb')
#           which.kernel \in c('gamma', 'gaussian','pareto','window')
# Dependency: functions get.gaussian.featmat(), subsets_using_reticulate(), volume.cvxhull()

sequential.importance.sample.from.shape.cvxhull = function(n, k,
                                                            feat.mat = NULL, subset.mat = NULL, feat.var = NULL,
                                                            pi.init = NULL, which.node.importance,
                                                            which.kernel,param.vec,
                                                            num.sample, if.save.si.sample, save.path){
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
        msg('feat.mat is not defined and a feat. mat. is generated.')
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
        msg('subset.mat is not found and a subset. mat. is calculated.')
    }
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
        { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
            stop()
        })
    res = matrix(0, nrow = num.sample, ncol = k, byrow = T)
    if(is.null(pi.init)){
        pi.init = rep(1,n)
    }
    switch(which.node.importance,
        'small.chull' = {
            for(i in 1:num.sample){
            kset = rep(0,k) # a k-set to be drawn
            rem = 1:n # all the authors not being picked up for now
            kset[1] = sample(rem, 1 , prob = pi.init) # pick the first vertex accroding to pi.init
            rem = rem[-kset[1]] # update the available set
            for (j in 2:k){
                prob = sapply(1:length(rem), function(i) shape.kern(volume.cvxhull(c(kset[which(kset > 0.5)], rem[i]), feat.mat)))
                if(sum(prob) < machine.precision){
                    msg("Warning from 'sequential.importance.sample.from.shape.cvxhull()' when calculating sequential distribution: zero vector is returned! A random walk is returned instead")
                    kset[j] = sample(rem, 1, replace = F)
                } else {
                    kset[j] = sample(rem, 1, replace = F, prob = prob)
                }
                rem = rem[-which(rem == kset[j])]
            }
            res[i,] = sort(kset)
            }
        },
        'near.nb' = {
            for(i in 1:num.sample){
                kset = rep(0,k) # a k-set to be drawn
                rem = 1:n # all the authors not being picked up for now
                kset[1] = sample(rem, 1 , prob = pi.init) # pick the first vertex accroding to pi.init
                rem = rem[-kset[1]] # update the available set
                for (j in 2:k){
                    distances = sapply(1:length(rem), function(i) norm.l2(feat.mat[kset[j-1],] - feat.mat[rem[i],])) # distance from current node to all the available nodes
                    if(sum(distances) < machine.precision){
                        msg("Warning from'sequential.importance.sample.from.shape.cvxhull()' when calculating sequential distribution: all nodes are too closed to each other! A random walk is returned instead")
                        kset[j] = sample(rem, 1, replace = F)
                    } else {
                        kset[j] = rem[which(distances == min(distances))]
                    }
                    rem = rem[-which(rem == kset[j])]
                }
                res[i,] = sort(kset)
            }
        },
        { msg("node.importance is  not available! Pick one from c('small.chull','near.nb'). ")
                    stop()
        }
    )
    sim.data = list(count.vec = turn.kset.sample.to.count(res, subset.mat),
                    feature.mat = feat.mat,
                    subset.mat = subset.mat,
                    edge.sampled = res)
    if(if.save.si.sample){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(sim.data, file = paste0(save.path, "sis-samples-", which.kernel,"-n", n, "-a",param.vec[1],'-b',param.vec[2],"-",current.time = format(Sys.time(), "%b-%d-%I-%M-%p"), ".Rdata"))
    }
    return(sim.data)
}




## other Function 1.8: count.vert.freq(edges:mat[int]) -> n-array[int]
# Count occurences for n vertices using samples from cvxhull shape model
# Remark: output is an array with names 1:n
count.vert.freq = function(edges){
    if (is.list(edges)){
        return(table(unlist(edges, use.name = F)))
    } else {
        return(table(edges))
    }
}


## mod Funtion 1.10 kset.freq.from.samples(samples:mat[int], subset.mat:mat[int]) -> array[f]
#  Calculate sample(k-set) frequency using k-set sampls, which is stored in mat.
# Remark： an output vector is of length 'n choose k'.
kset.freq.from.samples = function(samples, subset.mat){
    count = turn.kset.sample.to.count(samples, subset.mat)
    count/sum(count)
}

## mod Function 1.11:  find.pi.sis(k: int, feat.mat: mat[f], dist: array[f]) -> array[f]
# Find optimal initial distribution \Pi0, which is used in the first step of a  seq. import. sampler for shape cvxhull model
# Dependence: function linp() from Rpackage 'limSolve', path.set(), edge.prob(), path.prob()
# To be fixed :
#  (i) k = 2 not working
#  (ii) n =100 , k = 3 : function "pi.zero.sis" takes times  roughly from 10:44 to 11:20, about 45min !!!  Too slow for large n and k
find.optimal.pi.sis= function(k, feat.mat, dist, which.shape, shape.vec){
    if(k < 1+1e-10){ stop("k = 1 is not ready!")}
    n = dim(feat.mat)[1]
    all = subsets_using_reticulate(n,k)
    n.k = length(dist)
    if ( choose(n,k) != n.k) {stop("k-set distribution is not corret!")}
    p.mat = matrix(0, n.k, n)
    for(i in 1 : n.k){
        paper = all[i,]
        for(j in 1: k){
            first = paper[j]
            path = path.set(start = first, path.vec = paper)
            if (k < 2 + 1e-10){
                p.mat[i, first] = path.prob(path = path, feat.mat,which.shape, shape.vec)
            } else {
                prob = 0
                num.path  = dim(path)[1]
                for (l in 1: num.path ){
                    prob =  prob + path.prob(path = path[l,], feat.mat,which.shape, shape.vec)
                    }
                p.mat[i, first]  = prob
            }
        }
    }
    #linp(G = rbind(diag(n),-p.mat), H = c(rep(0,n),-dist), Cost = rep(-1,n))$X
    lsei(A=p.mat, B=dist, E = rep(1,n), F = 1, G= diag(n),H=rep(0,n))$X
}




## mod Function 1.12: path.set(start:int, path.vec:array[int]) -> mat[int]
# Helper in "find.optimal.pi.sis",  return all combinatorial sequential path(s) from a given node and then visits a given set
# Usage: output is of size nrow = (length(path.vec)-1)! \times length(path.vec)
# Usage:
#       Input:  start = 1, path.vec = 1:3
#       Output:
#                  [,1] [,2] [,3]
#             [1,]    1    2    3
#             [2,]    1    3    2
# Dependence: get.perm_using_reticulate(), $notin%
path.set= function(start, path.vec){
    if (start %notin% path.vec){stop("starting node is not in path !")}
    if (length(path.vec) < 2+1e-10){ # k = 2
        return(c(start, path.vec[which(path.vec != start)]))
    } else {
        mat = get.perm_using_reticulate(path.vec[which(path.vec != start)]) # permute all the other nodes than the first node
        #return(t(sapply(1:(length(path.vec)-1), function(i) c(start, mat[i,]))))
        return(cbind(rep(start, dim(mat)[1]), mat))
    }
}



## mod Function 1.13: edge.prob(from.vec: 2-array[int] or int, to.node: int, feat.mat:mat[f]) -> float
# Helper in path.prob(), return the prob of a path P : s.vec --> t.node is traveled by a sequential sampler
# Usage:
#        Input:  from.vec = 1, to.node = 2
#        Output:  \frac{ p:1->2 }{ p:1->2 + p:1->3 + p:1->4}
#    or  Input : from.vec = 12, to.node = 3
#        Output:  \frac{ p:12->3 }{ p:12->3 + p:12->4}
#   and  len(from.vec) < dim(feat.mat) - 1
edge.prob = function(from.vec, to.node, feat.mat, which.shape.func, shape.param.vec){
    if ( to.node %in% from.vec ){
        return( 0 )
    } else {
        switch(which.shape.func,
            "gamma" = { shape.kern = function(d, alpha =  shape.param.vec[1], beta = shape.param.vec[2]){
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
            "gaussian" ={   shape.kern = function(d, mu =  shape.param.vec[1], sd = shape.param.vec[2]){
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
            "pareto" = { shape.kern = function(d, alpha =  shape.param.vec[1], beta = shape.param.vec[2]){
                        if (d > beta + machine.precision & alpha > machine.precision & beta > machine.precision){
                            return(d^(-alpha-1))
                        } else {
                            0
                        }
                }
            },
            "window" = { shape.kern = function(d, param.lbd = shape.param.vec[1], param.ubd = shape.param.vec[2]){
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
            { msg("Shape type is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
                stop()
            }
        )
        candt = (1:dim(feat.mat)[1])[-from.vec]
        num = shape.kern(volume.cvxhull(c(from.vec,to.node), feat.mat))
        den = sum(sapply(candt, function(j) shape.kern(volume.cvxhull(c(from.vec,j), feat.mat)) ))
        num/den
        #exp(-volume.cvxhull(c(from.vec, to.node), feat.mat)) / sum(sapply(candt, function(j) exp(-volume.cvxhull(c(from.vec, j), feat.mat)))) # ????????????????????????????????????????????????çççç
    }
}


## mod Function 1.14 path.prob(path:array[int], feat.mat:mat[f]) -> float ???????????????????????????????????????? optimal
# Helper in find.pi.sis(), calculate the prob that a path is travel sequentially
# Usage:
#        Input path = c(1,2,3), from.node = 1
#        Output (p:1->2->3) := (p:1->2) * (p:12->3)
path.prob = function(path, feat.mat, which.shape, param.vec){
    res = 1
    start = path[1]
    for (i in 2 : length(path)){
        res = res* edge.prob(from.vec = start, to.node = path[i], feat.mat, which.shape.func = which.shape, shape.param.vec = param.vec) #???????????????????????????????????????? optimal
        start = c(start, path[i])
    }
    res
}


## mod Function 1.15: set.init.dist.sis(feat.mat: mat[f], which.rule: str) -> n-array[f] ???????
# Calculate a distribution \Pi_0 over [n], which is used in the first step of  sis function 'sequential.importance.sample.from.shape.cvxhull'
# Usage:
#       which.rule \in  c('uniform','weighted'，'edge','mean','optimal-pi') ##################################################################### ??????????? unfinished
#       'uniform' : uniformly pick a vertex
#       'weighted' : weighted by node's degree
#       'optimal-pi':  implied by true model, will produce true sample but too costly to compute.
#       obselete:
#                 'shortest': find a shortest edges among all n choose 2 edges, and then uniformly pick one node from the edge
#                 'nearest': pick a vertex cloest to sample mean of feature vectors

set.init.dist.sis = function(feat.mat, which.rule){
    if (is.null(dim(feat.mat))){
        n = length(feat.mat)
        p = 1
    } else {
        n = dim(feat.mat)[1]
        p = dim(feat.mat)[2]
    }
    switch(which.rule,
                    "uniform" = {
                        dist =  rep(1,n)
                    },
                    #"optimal.init" = {
                    #    dist= optimal.pi
                    #},
                    "edge" = {
                        distance = matrix(0,n,n)
                        for (i in 1:n){
                            if ( p < 1+1e-10){ # 1-dim feature
                                distance[i,] = sapply(1:n, function(j) norm.l2(feat.mat[i] - feat.mat[j]))
                            }  else { # p-dim feature
                                distance[i,] = sapply(1:n, function(j) norm.l2(feat.mat[i,] - feat.mat[j,]))
                            }
                        }
                        pair = which(distance == min(distance[which(distance > 0)]), arr.ind=TRUE)[1,]
                        dist = rep(0,n)
                        dist[pair] = 1/2
                    },
                    "mean" = {
                        if ( p < 1+1e-10){
                            mu = mean(feat.mat)
                            distance.vec = sapply(1:n, function(j) norm.l2(mu - feat.mat[j]) )
                        } else {
                            mu = colMeans(feat.mat)
                            distance.vec = sapply(1:n, function(j) norm.l2(mu - feat.mat[j,]) )
                        }
                        init.vert = which(distance.vec == min(distance.vec))
                        dist = rep(0,n)
                        dist[init.vert] = 1
                    },
                    { msg("The first step is  not available! pick a method from {uniform, edge, mean}.")
                        stop()
                    }
        )
    dist
}



# mod Function 1.18:  experiment.sis.performance.edge.tv(  n:int, k:int,
#                                                               feat.mat:mat[f], subset.mat:mat[i], true.dist:array[f], which.shape:str, param.vec:2-array[f],
#                                                               rule.vec:array[str], num.sample:int, iter.tag.vec:array[int],
#                                                               if.true.sample:logical, if.show.plot:logical, if.save.si.sample:logical, if.save.result:logical, if.save.plot:logical, save.path:str)
# Experiment compare sis performance on induced edge freq.
# Dependence: find.true.edge.freq(), find.sample.edge.freq()
experiment.sis.performance.edge.tv   = function(n, k,
                                                feat.var = 1 , feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                true.edge.freq = NULL,
                                                init.pi.rules, optimal.init.pi = NULL,
                                                num.sample, iter.tag.vec, if.true.sample,
                                                if.show.plot, if.save.si.sample = T, if.save.experiment.result, if.save.plot, save.path){
    msg('Experiment edge frequency is initiating.')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    if (is.null(true.edge.freq)){
        msg('theoritical induced edge freq is calculating')
        true.freq = find.true.edge.freq(feat.mat, subset.mat, which.shape, param.vec)
        msg('theoritical induced edge freq is done')
    }
    df = data.frame()
    if(if.true.sample){
        sample.mat = sample.from.shape.cvxhull(n = n, k = k,
                                        subset.mat = subset.mat, feat.mat = feat.mat,
                                        num.sample = num.sample, which.shape = which.shape , param.vec = param.vec,
                                        if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
        tv  = sapply(iter.tag.vec, function(num){ tv.dist(  f1 = true.freq,
                                                            f2 = find.sample.edge.freq(n = n, k = k, sample.mat = sample.mat[1:num,])
                                                        )})
        df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("exact.sample")))
    }
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            sample.mat = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            tv  = sapply(iter.tag.vec, function(num){ tv.dist(  f1 = true.freq,
                                                                f2 = find.sample.edge.freq(n = n, k = k, sample.mat = sample.mat[1:num,])
                                                            )})
            df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                sample.mat = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
                tv  = sapply(iter.tag.vec, function(num){ tv.dist(  f1 = true.freq,
                                                                    f2 = find.sample.edge.freq(n = n, k = k, sample.mat = sample.mat[1:num,]))})
                df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    experiment.result = list(df = df, true.freq = true.freq)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-edge-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) + ylim(0, 1)
    cols = c('exact.sample' = dark.red, 'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green, 'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple, 'edge'= light.yellow, 'mean' = dark.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= cols)
    plot <- plot + ylab("tv dist from true edges freq.") + xlab("Iteration") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-error-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment sis sampling error on edge degree is finished. ')
    #if(if.show.plot){
    #    plot
    #}
    return(list(df = df, true.freq = true.freq, plot = plot))
}


# comb Function 1.19 :find.true.vx.freq(subset.mat: mat[int], feat.mat: mat[f], which.shape:str, param.vec: 2-array[f]) -> n-array[f]
# Calculate true vertices frequnce from cvxhull shape distribution. Particularly freq(node i) \propto sum of all weighted out-paths-from-node-i
find.true.vx.freq = function(feat.mat,subset.mat, which.shape, param.vec){
    n = dim(feat.mat)[1]
    k = length(subset.mat[1,])
    switch(which.shape,
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
        { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
            stop()
        }
    )
    res = rep(0,n)
    for (i in 1:n){
        rem.set = (1:n)[-i]
        rem.comb = subsets_using_reticulate(rem.set, k-1)
        n.row = dim(rem.comb)[1]
        all.paths.from.i = cbind(rep(i,n.row), rem.comb)
        all.paths.from.i
        res[i] = sum(sapply(1:n.row, function(i) shape.kern(volume.cvxhull(all.paths.from.i[i,], feat.mat))))
    }
    if(sum(res)<machine.precision){
        msg('warining: true vertex freq are almost zero evrywhere and a uniform distributin is returned!!!')
        return(rep(1/n,n))
    }
    res/sum(res)
}

# comb Function 1.20 :find.true.edge.freq(subset.mat: mat[int], feat.mat: mat[f], which.shape:str, param.vec: 2-array[f]) -> n-array[f]
# Calculate true edge frequncy, induced by cvxhull shape distribution.
find.true.edge.freq = function(feat.mat,subset.mat, which.shape, param.vec){
    n = dim(feat.mat)[1]
    k = length(subset.mat[1,])
    if (k < 2+1e-10){
        return(dist.shape.cvxhull(subset.mat = subset.mat, feat.mat = feat.mat, which.kernel = which.shape, param.vec = param.vec))
    } else {
        switch(which.shape,
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
            { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
                stop()
            }
        )
        n2 = choose(n,2)
        res = rep(0,n2)
        edge.mat = subsets_using_reticulate(n, 2)
        for (i in 1:n2){
            rem.set = (1:n)[-edge.mat[i,]]
            rem.comb = subsets_using_reticulate(rem.set, k-2)
            n.row = dim(rem.comb)[1]
            all.paths.from.edge= cbind(matrix(rep(edge.mat[i,], n.row), n.row, byrow = T), rem.comb)
            res[i] = sum(sapply(1:n.row, function(i) shape.kern(volume.cvxhull(all.paths.from.edge[i,], feat.mat))))
        }
        if(sum(res)<machine.precision){
            msg('warining: true vertex freq are almost zero evrywhere and a uniform distributin is returned!!!')
            return(rep(1/n2,n2))
        }
        res/sum(res)
    }

}



# comb Function 1.21: turn.paths.to.edges(sample:mat[int] or array,total.node:int) -> mat[int]
# Convert k-sets into a graph and then return edgelist of the graph
# Dependence: Rpackage igraph, turn.paths.to.graph.obj().
turn.paths.to.edges = function(sample, total){
    if (!is.matrix(sample)){
        sample = matrix(sample,1)
    }
    if (dim(sample)[2] > 3 - 1e-10){
        g = turn.paths.to.graph.obj(paths = sample, total = total)
        return(as_edgelist(g))
    } else if (dim(sample)[2] < 3 - 1e-10 & dim(sample)[2] > 1 + 1e-10){
        return(sample)
    } else{
        msg('Error from turn.paths.to.edges() : each row must be of length 2 or more.')
        stop()
    }
}



## mod  Function 1.22 find.sample.edge.freq(n:int, k:int, sample.mat:mat[int]) -> array[f]
#  Calculate edge frequency from k-set samples
#  Usage: turn k-set samples into a graph object and returns all edges in a matrix, Then count edges' frequency
# Dependence: turn.paths.to.edges(), turn.kset.sample.to.count()
find.sample.edge.freq = function(n, k, sample.mat){
    sample.edges = turn.paths.to.edges(sample = sample.mat, total = n)
    edge.index = subsets_using_reticulate(n,2)
    turn.kset.sample.to.count(edges = sample.edges, subset = edge.index)
}

# mod Function 1.23 find.sample.vx.degree(n:int, sample:mat[int] or array[int]) -> array[int]
# Dependence: Rpackage igraph
find.sample.vx.degree = function(n, sample){
    if(is.null(dim(sample))){
        sample = matrix(sample,1)
    }
    if(len(sample[1,])<2+1e-10){
        g = turn.edges.to.graph.obj(edges = sample, total = n)
    } else {
        g = turn.paths.to.graph.obj(paths = sample, total = n)
    }
    return(degree(g))
}



# mod Function 1.24 find.sample.vx.avg.degree
# Dependence: Rpackage igraph
find.sample.vx.avg.degree = function(n, sample){
    if(is.null(dim(sample))){
        sample = matrix(sample,1)
    }
    if(len(sample[1,])<2+1e-10){
        g = turn.edges.to.graph.obj(edges = sample, total = n)
    } else {
        g = turn.paths.to.graph.obj(paths = sample, total = n)
    }
    deg = degree_distribution(g)
    return(sum((1:length(deg))* deg))
}

# mod Function 1.24' find.sample.vx.degree.dist
# Dependence: Rpackage igraph
find.sample.vx.degree.dist = function(n, sample){
    if(is.null(dim(sample))){
        sample = matrix(sample,1)
    }
    if(len(sample[1,])<2+1e-10){
        g = turn.edges.to.graph.obj(edges = sample, total = n)
    } else {
        g = turn.paths.to.graph.obj(paths = sample, total = n)
    }
    deg = degree_distribution(g)
    return(deg)
}


# Function 1.25:  experiment.sis.performance.vertex.degree = function( n:int, k:int,
#                                                               feat.mat:mat[f], subset.mat:mat[i],
#                                                               true.vx.freq:array[f], which.shape:str, param.vec:2-array[f],
#                                                               rule.vec:array[str], num.sample:int, iter.tag.vec:array[int],
#                                                               if.true.sample:logical, if.show.plot:logical, if.save.si.sample:logical, if.save.result:logical, if.save.plot:logical, save.path:str)
# Experiment compare sis performance on induced vertex  degree distribution.

experiment.sis.performance.vertex.degree = function(n, k,
                                                feat.var, feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                init.pi.rules, optimal.init.pi = NULL,
                                                num.sample, iter.tag.vec, if.true.sample,
                                                if.show.plot, if.save.si.sample = T, if.save.experiment.result, if.save.plot, save.path){
    msg('Experiment sis sampling error is initiating. ')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    df = data.frame()
    samples = sample.from.shape.cvxhull(n = n, k = k,
                                    subset.mat = subset.mat, feat.mat = feat.mat,
                                    num.sample = min(10*choose(n,k),1e6), which.shape = which.shape , param.vec = param.vec,
                                    if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
    target.dist = find.sample.vx.degree(n, sample = samples)
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                                f2= find.sample.vx.degree(n, sample = samples[1:num,]))})
            df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
                tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                                    f2= find.sample.vx.degree(n, sample = samples[1:num,]))})
                df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    experiment.result = list(df = df, sampled.degree.dist = target.dist)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-vertex-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) + ylim(0, 1)
    cols = c('exact.sample' = dark.red, 'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green, 'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple, 'edge'= light.yellow, 'mean' = dark.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= cols)
    plot <- plot + ylab("t.v. dist from sampled vertex degree.") + xlab("iter") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-performance-vx-deg-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment sis sampling error on vertex deg is finished. ')
    #if(if.show.plot){
    #    plot
    #}
    return(list(df = df, sampled.degree.dist = target.dist, plot = plot))
}



# Function 1.26:  experiment.sis.performance.vertex.avg.degree = function( n:int, k:int,
#                                                               feat.mat:mat[f], subset.mat:mat[i],
#                                                               true.vx.freq:array[f], which.shape:str, param.vec:2-array[f],
#                                                               rule.vec:array[str], num.sample:int, iter.tag.vec:array[int],
#                                                               if.true.sample:logical, if.show.plot:logical, if.save.si.sample:logical, if.save.result:logical, if.save.plot:logical, save.path:str)
# Experiment compare sis performance on induced vertex freq.
# Dependence: count.vert.freq()
experiment.sis.performance.vertex.avg.degree = function(n, k,
                                                feat.var, feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                init.pi.rules, optimal.init.pi = NULL,
                                                num.sample, iter.tag.vec, if.true.sample,
                                                if.show.plot, if.save.si.sample = T, if.save.experiment.result, if.save.plot, save.path){
    msg('Experiment  is initiating. ')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    df = data.frame()
    samples = sample.from.shape.cvxhull(n = n, k = k,
                                    subset.mat = subset.mat, feat.mat = feat.mat,
                                    num.sample = num.sample, which.shape = which.shape , param.vec = param.vec,
                                    if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
    avg.deg = sapply(iter.tag.vec, function(num) find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,], feat.mat = feat.mat, which.shape = which.shape , param.vec = param.vec))
    df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0("exact.sample")))
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            avg.deg = sapply(iter.tag.vec, function(num) find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,],feat.mat = feat.mat, which.shape = which.shape , param.vec = param.vec))
            df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init.optimal-import-',which.importance,'is done.'))
                avg.deg = sapply(iter.tag.vec, function(num) find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,], feat.mat = feat.mat, which.shape = which.shape , param.vec = param.vec))
                df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    experiment.result = list(df = df, sampled.avg.deg = avg.deg)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-vertex-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) #+ ylim(0, 1)
    cols = c('exact.sample' = dark.red, 'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green, 'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple, 'edge'= light.yellow, 'mean' = dark.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= cols)
    plot <- plot + ylab("average vx degree") + xlab("num. of samples") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-performance-vx-deg-dist-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment sis sampling error on vertex deg is finished. ')
    experiment.result$plot = plot
    return(experiment.result)
}



#Function 1.27 snowball.importance.sample.from.nfeat(n : int, k :int,
#                                                num.snb.nb: int,
#                                                feat.mat = NULL, subset.mat = NULL,
#                                                pi.init,
#                                                which.kernel, param.vec,
#                                                num.sample, if.save.snb.sample, save.path)
# Snowball  importance sampling algorithm: generate multiple samples from convex hull shape model
# Usaage:
#           num.sample: consistent to other sequential k-set sampler. snowball round = ceiling(log(num.sample) - log(k) + 1),
#                           k^{snb.round} \approx (k-1)*N so that number of edges are kept same for various samplers to compare
#           num.snb.nb: at each snowball round, 'num.snb.nb' is number of neightbours to talk to for each node being sampled
#           which.kernel \in c('gamma', 'gaussian','pareto','window')
# Dependency: functions get.gaussian.featmat(), subsets_using_reticulate(), volume.cvxhull()
snowball.importance.sample.from.nfeat = function(n, k,
                                                num.snb.nb,
                                                feat.var, feat.mat = NULL, subset.mat = NULL,
                                                pi.init,
                                                which.kernel, param.vec,
                                                num.sample, if.save.snb.sample, save.path){
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
        msg('feat.mat is not defined and a feat. mat. is generated.')
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
        msg('subset.mat is not found and a subset. mat. is calculated.')
    }
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
        { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
            stop()
        })
    if(is.null(pi.init)){
        first = sample(n, 1, prob = pi.init)
    } else{
        first = sample(n,1)
    }
    nodeset.vec = first #  node set
    edgemultiset.mat = list()# edge multiset
    edge.ind = 0
    if (num.sample < k + 1e-10){
        msg('Error from snowball.importance.sample.from.nfeat: num.sample must be larger than k')
    }
    while(edge.ind < num.sample + 1e-10){
        for (cur in nodeset.vec){
            total = 1:n
            wt.edges = sapply(1:length(total), function(i) norm.l2(feat.mat[cur,] - feat.mat[total[i],]))
            wt.edges[cur] = Inf # no self connected edges
            new.snb.nb = sapply(1:num.snb.nb, function(x) which(wt.edges == sort(wt.edges)[x])) # return top nearest nodes, controled by num.snb.nb
            for (node in new.snb.nb){
                edgemultiset.mat[[edge.ind+1]] = sort(c(cur,node))
                edge.ind = edge.ind + 1
            }
            nodeset.vec = unique(c(nodeset.vec, new.snb.nb))
            #print(c(t, edge.ind))
        }
    }
    nodeset.vec = sort(nodeset.vec)
    result = list(edge.sampled = turn.list.to.matrix(edgemultiset.mat), feature.mat = feat.mat, subset.mat = subset.mat, node.visited = nodeset.vec)
        if(if.save.snb.sample){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(result, file = paste0(save.path, "snb-samples-", which.kernel,"-n", n, "-a",param.vec[1],'-b',param.vec[2],"-",current.time = format(Sys.time(), "%b-%d-%I-%M-%p"), ".Rdata"))
    }
    return(result)
}


#Function 1.27' forest.fire.sample.from.nfeat(n : int, k :int,
#                                                ffs.prob : float,
#                                                feat.mat = NULL, subset.mat = NULL,
#                                                pi.init,
#                                                which.kernel, param.vec,
#                                                num.sample, if.save.snb.sample, save.path)
#                                   -> list(edge.sampled: mat[int], feature.mat:mat[int], subset.mat:mat[int], node.visited:array)
# Forest fire importance sampling algorithm: generate multiple samples from convex hull shape model
# Usaage:
#           num.sample: consistent to other sequential k-set sampler. snowball round = ceiling(log(num.sample) - log(k) + 1),
#                           k^{snb.round} \approx (k-1)*N so that number of edges are kept same for various samplers to compare
#           num.snb.nb: at each snowball round, 'num.snb.nb' is number of neightbours to talk to for each node being sampled
#           which.kernel \in c('gamma', 'gaussian','pareto','window')
# Dependency: functions get.gaussian.featmat(), subsets_using_reticulate(), volume.cvxhull()

forest.fire.importance.sample.from.nfeat = function(n, k,
                                        ffs.prob,
                                        feat.var, feat.mat = NULL, subset.mat = NULL,
                                        pi.init,
                                        which.kernel, param.vec,
                                        num.sample, if.save.ffs.sample, save.path){
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
        msg('feat.mat is not defined and a feat. mat. is generated.')
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
        msg('subset.mat is not found and a subset. mat. is calculated.')
    }
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
        { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
            stop()
        })
    if(is.null(pi.init)){
        first = sample(n, 1, prob = pi.init)
    } else{
        first = sample(n,1) # defaul: pick first node uniformly over [n]
    }
    nodeset.vec = first #  node set
    edgemultiset.mat = list()# edges are stroed in list first and then convrt list into a mat, of which a row represents an edge.
    edge.ind = 0
    if (num.sample < k + 1e-10){
        msg('Error from snowball.importance.sample.from.nfeat: num.sample must be larger than k')
    }
    while(edge.ind < num.sample + 1e-10){
        for (cur in nodeset.vec){
            total = 1:n
            wt.edges = sapply(1:length(total), function(i) norm.l2(feat.mat[cur,] - feat.mat[total[i],]))
            wt.edges[cur] = Inf # no self connected edges
            num.ff.nb = 0
            while(num.ff.nb < 1e-10){
                num.ff.nb = min(rgeom(1, prob = ffs.prob), n-1)  # Fores fire sampler determines number of  neighbour(s)
            }
            new.nb = sapply(1:num.ff.nb, function(x) which(wt.edges == sort(wt.edges)[x])) # return top nearest nodes, controled by num.ff.nb
            for (node in new.nb){
                edgemultiset.mat[[edge.ind+1]] = sort(c(cur,node))
                edge.ind = edge.ind + 1
            }
            nodeset.vec = unique(c(nodeset.vec, new.nb))
            #print(c(t, edge.ind))
        }
    }
    nodeset.vec = sort(nodeset.vec)
    result = list(edge.sampled = turn.list.to.matrix(edgemultiset.mat), feature.mat = feat.mat, subset.mat = subset.mat, node.visited = nodeset.vec)
        if(if.save.ffs.sample){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(result, file = paste0(save.path, "snb-samples-", which.kernel,"-n", n, "-a",param.vec[1],'-b',param.vec[2],"-",current.time = format(Sys.time(), "%b-%d-%I-%M-%p"), ".Rdata"))
    }
    return(result)
}








# Function 1.28:  experiment.snb.performance.vertex.freq  = function( n:int, k:int,
#                                                               feat.mat:mat[f], subset.mat:mat[i],
#                                                               true.vx.freq:array[f], which.shape:str, param.vec:2-array[f],
#                                                               rule.vec:array[str], num.sample:int, iter.tag.vec:array[int],
#                                                               if.true.sample:logical, if.show.plot:logical, if.save.si.sample:logical, if.save.result:logical, if.save.plot:logical, save.path:str)
# Experiment compare sis performance on induced vertex freq.
# Dependence: count.vert.freq()
experiment.snb.performance.vertex.freq  = function(  n, k,
                                                        num.snowball.neighbour, forest.fire.prob = NULL,
                                                        feat.var, feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                        true.vx.freq = NULL,
                                                        init.pi.rules, optimal.init.pi = NULL,
                                                        num.sample, iter.tag.vec, if.true.sample,
                                                        if.show.plot, if.save.si.sample = T, if.save.experiment.result, if.save.plot, save.path){
    msg('Experiment sis sampling error is initiating. ')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    if (is.null(true.vx.freq)){
        msg('theoritical induced vertex freq is calculating')
        true.freq = find.true.vx.freq(feat.mat, subset.mat, which.shape, param.vec)
        msg('theoritical induced vertex freq is done')
    }
    df = data.frame()
    if(if.true.sample){
        samples = sample.from.shape.cvxhull(n = n, k = k,
                                        subset.mat = subset.mat, feat.mat = feat.mat,
                                        num.sample = num.sample, which.shape = which.shape , param.vec = param.vec,
                                        if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
        tv  = sapply(iter.tag.vec, function(num){ tv.dist( true.freq, count.vert.freq(samples[1:num,]))})
        df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("exact.sample")))
    }
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            tv  = sapply(iter.tag.vec, function(num){ tv.dist( true.freq, count.vert.freq(samples[1:num,]))})
            df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
                tv  = sapply(iter.tag.vec, function(num){ tv.dist( true.freq, count.vert.freq(samples[1:num,]))})
                df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    edges = snowball.importance.sample.from.nfeat(n = n.demo, k = k.demo, # Snowball sampling
                                                num.snb.nb = num.snowball.neighbour,
                                                feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample.demo, if.save.snb.sample = T, save.path = './output/')$edge.sampled
    msg("Snowball sampling on feature space is done.")
    tv  = sapply(iter.tag.vec, function(num){ tv.dist( true.freq, count.vert.freq(edges[1:num,]))})
    df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0('snowball')))
    if (!is.null(forest.fire.prob)){ # Forest fire sampling
        edges = forest.fire.importance.sample.from.nfeat(n = n.demo, k = k.demo,
                                                ffs.prob = forest.fire.prob,
                                                feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample.demo, if.save.ffs.sample = T, save.path = './output/')$edge.sampled
        msg("Forest fire sampling on feature space is done.")
        tv  = sapply(iter.tag.vec, function(num){ tv.dist( true.freq, count.vert.freq(edges[1:num,]))})
        df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0('forestfire')))
    }
    experiment.result = list(df = df, true.freq = true.freq)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-vertex-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) + ylim(0, 1)
    cols = c('exact.sample' = dark.red, 'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green, 'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple, 'snowball'= light.yellow, 'forestfire' = dark.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= cols)
    plot <- plot + ylab("t.v. dist from true vertex freq.") + xlab("Iteration") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-error-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment sis sampling error on vertex deg is finished. ')
    #if(if.show.plot){
    #    plot
    #}
    return(list(df = df, true.freq = true.freq, plot = plot))
}




# Function 1.29 turn.edges.to.graph.obj(edge: mat[i] or 2-array[i], total:i) -> graph object
# Dependence: Rpackage igrpah
# https://igraph.org/r/doc/add_edges.html
turn.edges.to.graph.obj = function(edges, total){
    if(is.null(dim(edges))){
        edges = matrix(edges, 1)
    }
    if (dim(edges)[2] > 2+1e-10){
        msg("Error from turn.edges.to.graph.obj: length of each row must be TWO.")
        stop()
    } else{
        g <- make_empty_graph(directed = FALSE, n = total) %>% add_edges(c(t(edges)))
        return(g)
    }
}


# Function 1.30 turn.paths.to.graph.obj
#https://igraph.org/r/doc/path.html
turn.paths.to.graph.obj = function(paths, total){
    if(is.null(dim(paths))){
        paths = matrix(paths, 1)
    }
    if (dim(paths)[2] < 2+1e-10){
        msg("Error from turn.paths.to.graph.obj: length of each row must be larger than TWO.")
        stop()
    } else{
        g <- make_empty_graph(directed = FALSE, n = total)
            for (i in 1:dim(paths)[1]){
                g = g + path(paths[i,])
            }
    return(g)
    }
}


# Function 1.31 experiment.snb.performance.vertex.degree
experiment.snb.performance.vertex.degree = function(n, k,
                                                num.snowball.neighbour, forest.fire.prob = NULL,
                                                feat.var, feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                init.pi.rules, optimal.init.pi = NULL,
                                                num.sample, iter.tag.vec, if.true.sample,
                                                if.show.plot, if.save.si.sample = T, if.save.experiment.result, if.save.plot, save.path){
    msg('Experiment snowball sampling error is initiating. ')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    df = data.frame()
    samples = sample.from.shape.cvxhull(n = n, k = k,
                                    subset.mat = subset.mat, feat.mat = feat.mat,
                                    num.sample = min(10*choose(n,k),1e5), which.shape = which.shape , param.vec = param.vec, # Caution 10*choose(n,k)
                                    if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
    target.dist = find.sample.vx.degree(n, sample = samples)
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                                f2= find.sample.vx.degree(n, sample = samples[1:num,]))})
            df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
                tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                                    f2= find.sample.vx.degree(n, sample = samples[1:num,]))})
                df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    edges = snowball.importance.sample.from.nfeat(n = n.demo, k = k.demo,
                                                num.snb.nb = num.snowball.neighbour,
                                                feat.var = feat.var, feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample, if.save.snb.sample = T, save.path = './output/')$edge.sampled
    msg("Snowball sampling on feature space is done.")
    tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                        f2= find.sample.vx.degree(n, sample = edges[1:num,]))})
    df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = 'snowball'))
    if(!is.null(forest.fire.prob)){
        edges = forest.fire.importance.sample.from.nfeat(n = n.demo, k = k.demo,
                                                ffs.prob = forest.fire.prob,
                                                feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample.demo, if.save.ffs.sample = T, save.path = './output/')$edge.sampled
        msg("Forest fire sampling on feature space is done.")
        tv  = sapply(iter.tag.vec, function(num){ tv.dist( f1 = target.dist,
                                                            f2= find.sample.vx.degree(n, sample = edges[1:num,]))})
        df = rbind(df, data.frame(x = iter.tag.vec, y = tv, initial.rule = 'forestfire'))

    }

    experiment.result = list(df = df, sampled.degree.dist = target.dist)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-vertex-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) + ylim(0, 1)
    col.assign = c('exact.sample' = dark.red,
                    'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green,
                    'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple,
                    'snowball'= dark.yellow,  'forestfire' = light.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= col.assign)
    plot <- plot + ylab("t.v. dist from sampled vertex degree") + xlab("iter") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-performance-vx-deg-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment snowball sampling error on vertex deg is finished. ')
    #if(if.show.plot){
    #    plot
    #}
    return(list(df = df, sampled.degree.dist = target.dist, plot = plot))
}






# Function 1.32 find.sample.vx.avg.weighted.degree
# Dependence: Rpackage igraph
find.sample.vx.avg.weighted.degree = function(n, sample.mat, graph.obj = NULL, feat.mat,  which.shape, param.vec){
    edgemat = turn.paths.to.edges(sample = sample.mat, total = n)
    switch(which.shape,
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
            { msg("Kernel is  not available! Pick one from c('gamma','gaussian','pareto','window'). ")
                stop()
            }
    )
    if (is.null(dim(edgemat))){
        return(shape.kern(volume.cvxhull(edgemat, feat.mat)))
    } else {
        n.row = dim(edgemat)[1]
        return(mean(sapply(1:n.row, function(i) shape.kern(volume.cvxhull(edgemat[i,], feat.mat)))))
    }
}

# Function 1.33
experiment.snb.performance.vertex.avg.degree = function(n, k,
                                                        num.snowball.neighbour, forest.fire.prob = NULL,
                                                        feat.var, feat.mat = NULL, subset.mat = NULL, which.shape, param.vec,
                                                        init.pi.rules, optimal.init.pi = NULL,
                                                        num.sample, iter.tag.vec, if.true.sample,
                                                        if.show.plot, if.save.si.sample = T, if.save.experiment.result = T, if.save.plot = T, save.path){
    msg('Experiment snowball sampling performance on avergae vertex weighted degree is initiating. ')
    if(is.null(feat.mat)){
        feat.mat = get.gaussian.featmat(nrow = n, ncol = 2, mean = 0, var1 = feat.var, cov12 = 0)
    }
    if(is.null(subset.mat)){
        subset.mat = subsets_using_reticulate(n, k)
    }
    df = data.frame()
    samples = sample.from.shape.cvxhull(n = n, k = k,
                                    subset.mat = subset.mat, feat.mat = feat.mat,
                                    num.sample = num.sample, which.shape = which.shape , param.vec = param.vec,
                                    if.save.sample = if.save.si.sample, save.path = save.path)$edge.sampled
    avg.deg = sapply(iter.tag.vec, function(num)
                                        find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,], feat.mat = feat.mat,
                                                                            which.shape = which.shape , param.vec = param.vec))
    df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0("exact.sample")))
    for (init.pi in init.pi.rules){
        dist.init = set.init.dist.sis(feat.mat = feat.mat, which.rule = init.pi)
        for (which.importance in c('small.chull','near.nb')){
            samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                    feat.mat = feat.mat, subset.mat = subset.mat,
                                                    pi.init = dist.init, which.node.importance = which.importance,
                                                    which.kernel = which.shape, param.vec = param.vec,
                                                    num.sample = num.sample,
                                                    if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
            msg(paste0('SI sampling: init-',init.pi,', import-',which.importance,'is done.'))
            avg.deg = sapply(iter.tag.vec, function(num) find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,],feat.mat = feat.mat, which.shape = which.shape , param.vec = param.vec))
            df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0(init.pi,"-",which.importance)))
            if(!is.null(optimal.init.pi)){
                msg('Using the given optimal initial pi and run seq. imp-sampling for the model')
                samples = sequential.importance.sample.from.shape.cvxhull(n = n, k = k,
                                                        feat.mat = feat.mat, subset.mat = subset.mat,
                                                        pi.init = optimal.init.pi, which.node.importance = which.importance,
                                                        which.kernel = which.shape, param.vec = param.vec,
                                                        num.sample = num.sample,
                                                        if.save.si.sample = if.save.si.sample, save.path)$edge.sampled
                msg(paste0('SI sampling: init.optimal-import-',which.importance,'is done.'))
                avg.deg = sapply(iter.tag.vec, function(num)
                                                    find.sample.vx.avg.weighted.degree(n, sample.mat = samples[1:num,], feat.mat = feat.mat,
                                                                                        which.shape = which.shape , param.vec = param.vec))
                df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0("optimal-",which.importance)))
            }
        }
    }
    edges = snowball.importance.sample.from.nfeat(n = n.demo, k = k.demo,
                                                num.snb.nb = num.snowball.neighbour,
                                                feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample.demo, if.save.snb.sample = T, save.path = './output/')$edge.sampled
    msg("Snowball sampling on feature space is done.")
    avg.deg  = sapply(iter.tag.vec, function(num){
                                        find.sample.vx.avg.weighted.degree( n, sample.mat = edges[1:num,], feat.mat = feat.mat,
                                                                            which.shape = which.shape , param.vec = param.vec)})
    df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0('snowball')))
    if(!is.null(forest.fire.prob)){
        edges = forest.fire.importance.sample.from.nfeat(n = n.demo, k = k.demo,
                                                ffs.prob = forest.fire.prob,
                                                feat.mat = feat.mat, subset.mat = subset.mat,
                                                pi.init = NULL,
                                                which.kernel = which.shape.demo, param.vec = shape.param.demo,
                                                num.sample = num.sample.demo, if.save.ffs.sample = T, save.path = './output/')$edge.sampled
        msg("Forest fire sampling on feature space is done.")
        avg.deg  = sapply(iter.tag.vec, function(num){
                                            find.sample.vx.avg.weighted.degree( n, sample.mat = edges[1:num,], feat.mat = feat.mat,
                                                                                which.shape = which.shape , param.vec = param.vec)})
        df = rbind(df, data.frame(x = iter.tag.vec, y = avg.deg, initial.rule = paste0('forestfire')))
    }
    experiment.result = list(df = df, sampled.avg.deg = avg.deg)
    if (if.save.experiment.result){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        save(experiment.result, file = paste0(save.path, "Result-experiment-sis-performance-vertex-tv-n", n, "k", k,"-",current.time, ".Rdata"))
    }
    plot <- ggplot(df) #+ ylim(0, 1)
    cols = c('exact.sample' = dark.red,
        'uniform-small.chull' = light.green, 'uniform-near.nb' = dark.green,
        'optimal-small.chull' = light.purple, 'optimal-near.nb' = dark.purple,
        'snowball'= dark.yellow, 'forestfire' = light.yellow)
    plot <- plot+geom_line(aes(x, y , colour= initial.rule)) + scale_colour_manual(values= cols)
    plot <- plot + ylab("average vx weighted degree") + xlab("num. of samples") #+  ggtitle(sprintf('n%i-k%i-a%.1fb%.1f', n, k, param.vec[1], param.vec[2] ))
    plot <- plot + theme_bw() + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.15), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            #legend.position = c(.1, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white")
                            )
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("Plot-Experiment-sis-samlping-performance-vx-deg-dist-n", n, "k", k,'-a', param.vec[1],'b',param.vec[2],'-',current.time ,".pdf"), path = save.path)
        msg(sprintf('A plot is saved in %s', save.path))
    }
    msg('Experiment snowball sampling performance on avergae vertex weighted degree is finished. ')
    experiment.result$plot = plot
    return(experiment.result)
}



# Function 1.34
# Metropolis-Hastings edge sampler
flat.mh.edge.sample.from.nfeat = function(n, k, feat.mat, subset.mat, which.kernel, param.vec, feat.var,
                                            num.sample,
                                            target, max.iter, sample.from.proposal = NULL, proposal, init.edge = NULL ){
    #msg('mh sampling is starting***')
    if (is.null(init.edge)){
        init.edge = 1:2
    }
    res = matrix(0, max.iter, 2)
    res[1,] = init.edge
    rej = 0
    n.k = dim(subset.mat)[1]
        for (t in 2: max.iter){
            res[t,] = res[t-1,]
            edge.tilde = subset.mat[sample(n.k,1),]
            valid = runif(1) < target(edge = edge.tilde, feat.mat = feat.mat) * proposal(edge.tilde, res[t-1,])/(target(edge = res[t-1,], feat.mat = feat.mat) * proposal(res[t-1,], edge.tilde))
            if(valid){
                res[t,] = edge.tilde
            } else {
                rej = rej + 1
            }
            #if (t %% 100 == 0){
            #    msg(sprintf('%.1f percert of sampling is done.',  100*(t/max.iter)))
            #}
        }
    #msg('Flat probit model: mh sampling is done***')
    list(edge.sampled = res[(max.iter-num.sample + 1):max.iter,], rejection.rate = rej/(max.iter-1))
}


# Function 1.35
fill.up.missing.node.with.zero = function(name.list, total){
    res = rep(0,total)
    for (i in 1:len(names(name.list))){
        ind = as.numeric(names(name.list)[i])
        res[ind] = name.list[i]
    }
    res
}


# mod Function 1.36 find.sample.vx.degree.dist
# Calculate  the leading (i.e. largest magnitude) eigenvalue and the corresponding eigenvector is calculated.
# Dependence: Rpackage igraph
find.sample.eigen.decomposition = function(n, sample){
    if(is.null(dim(sample))){
        sample = matrix(sample,1)
    }
    if(len(sample[1,])<2+1e-10){
        g = turn.edges.to.graph.obj(edges = sample, total = n)
    } else {
        g = turn.paths.to.graph.obj(paths = sample, total = n)
    }
    tryCatch(
        expr = {
            spect = spectrum(g, options = c(which = 'LM', maxiter = 9e8))
            return(list(value = spect$values, vector = spect$vectors))
        },
        error = function(err){
            msg(sprintf('Error %s from find.sample.eigen.decomposition() when calculating eigen decomposition.', err))
            return(list(value = NA, vector = NA))
        }#,
        #warning = function(wrn){
        #    msg(sprintf('Warning %s from find.sample.eigen.decomposition() when calculating eigen decomposition.', wrn))
        #    spect = spectrum(g, options = c(which = 'LM', maxiter = 9e8))
        #    return(list(value = spect$values, vector = spect$vectors))
        #}
    )
}



# Function 1.37
find.sample.cluster.coef = function(n, sample){
    if(is.null(dim(sample))){
        sample = matrix(sample,1)
    }
    if(len(sample[1,])<2+1e-10){
        g = turn.edges.to.graph.obj(edges = sample, total = n)
    } else {
        g = turn.paths.to.graph.obj(paths = sample, total = n)
    }
    return(transitivity(g, type = "average"))
}
