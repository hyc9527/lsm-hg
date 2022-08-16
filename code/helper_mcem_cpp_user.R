
#######################################################################################
### MCEM algorithms
#####################################################################################


# Function 3.8
# A helper function in finding MLE of latent feature position in MCEM
fill.dist.matrix = function(pair.dist.vec, total){
    res = matrix(0,total,total)
    res[lower.tri(res, diag=FALSE)] <- pair.dist.vec
    res = t(res)
    res[lower.tri(res, diag=FALSE)] <- pair.dist.vec
    res
}

# Function 3.10
# Find MLE for latent sufficient statistics (t_1,....t_n.k) for convex hull exponential model
# Dependence: package CVXR
mle.ss.cvxhull.exp = function(beta, count.data ){
    tryCatch(
    expr = {
        Y = Variable(length(count.data))
        obj = -beta*sum(count.data*Y) - sum(count.data)*log_sum_exp(-beta*Y)
        constr <- list(Y >= 1e-10)
        cvx.min.prob <- Problem(Maximize(obj),constr)
        result <- solve(cvx.min.prob)
        ss = result[[1]]
        return(ss)
    },
    error = function(err){
        msg(sprintf('%s mle.ss.cvxhull.exp().', err))
        return(NA)
    },
    warning = function(wrn){
            msg(sprintf('%s mle.ss.cvxhull.exp().', wrn))
            return(NA)
        }
    )
}



# Function 3.11
mds.latent.pos.mat.exponential.cvxhull = function(ss.vec, total, k, mu, Sigma){
    if (is.na(sum(ss.vec)) || is.null(ss.vec)){
        #msg('ss.vec in  mds.latent.pos.mat.exponential.cvxhull () must be non-empty value')
        return(NA)
    } else {
        #msg("start setting up mds algorithm on suff stat.")
        dist.mat = fill.dist.matrix(pair.dist.vec = ss.vec, total = total)
        if (min(dist.mat < 0)){
            dist.mat = dist.mat - min(dist.mat)
        }
        tryCatch(
            expr = {
                mds.result = cmdscale(d = dist.mat, k = k)
                #msg("mle of latent feature mat is done.")
                return(mds.result)
            },
            error = function(err){
                msg(sprintf('%s mds.latent.pos.mat.exponential.cvxhull ().', err))
                return(NA)
            }
        )
    }
}

# Function 3.12 center and normalize matrix of size n by p
# adopted from https://rpubs.com/mengxu/procrustes_analysis
standerize.matrix = function(mat){
    mat.centered <- scale(mat, center = TRUE, scale = FALSE)
    mat.size <- norm(mat.centered, type = "F") / (ncol(mat) * nrow(mat))
    mat.normalized <- mat.centered / mat.size
    return(mat.normalized)
}

# Function 3.13
# from https://rpubs.com/mengxu/procrustes_analysis
find.procrustean.transform = function(std.Z0, std.Z){
    svd.results <- svd(std.Z%*% t(std.Z0))
    U <- svd.results$u
    V <- svd.results$v
    T <- V %*% t(U)
    return(T %*% std.Z)
}





# Function 1 initialize_mcem(...)
#     initialize params which are needed in mcem algorihtm
# Return
#       a List of lists, including
#        beta --- init value is random drawn from a gaussian(true beta, sd = random.start.shape.var )
#        Feat ---  init Feature matrix is return by warm.start()
#        mu,  Sigma ---  init is specified by users.
#        likelihood.estep --- init is calculated using beta.init and Feat.init
#        likelihood.mstep --- same as estep
#        rej --- init is 0
# Dependency
#      warm.start.for.Feat()
initialize_mcem = function(default.latent.mat.sd = NULL, if.warm.start = T,
                            true.beta,
                            beta.init, random.start.var = random.start.var.demo,
                            mu.init, Sigma.init,
                            tot.author, all.subset.mat, data.s){
    beta.lst =  Feat.lst=  mu.lst =  Sigma.lst =   likelihood.estep.lst =  likelihood.mstep.lst =  rej.lst =  list()
    if (is.null(beta.init)){
        beta.lst[[1]] =  max(rnorm(1, mean = true.beta, sd = sqrt(random.start.var)), 1e-10) #random start for shape parame in MCEM
    } else{
        beta.lst[[1]] = beta.init
    }
    mu.lst[[1]] = mu.init
    Sigma.lst[[1]] = Sigma.init
    if(is.null(default.latent.mat.sd)){
        default.latent.mat.sd = standerize.matrix(rmnorm(n = tot.author,
                                                        mean =  mu.lst[[1]],
                                                        varcov =  Sigma.lst[[1]]))
    }
    # estimate feature matrix by warm.start.for.Feat(...)
    if(if.warm.start){
            Feat.lst = warm.start.for.Feat(beta.init = beta.lst[[1]],
                                        default.mat = default.latent.mat.sd,
                                        tot.author = tot.author, count.data = data.s)
    } else {
        Feat.lst = default.latent.mat.sd
    }

    likelihood.estep.lst[[1]] = likelihood.mstep.lst[[1]] = log_likelihood_shape_cvxhull_gamma(count = data.s,
                                                                                            Feat =  Feat.lst,
                                                                                            Subset = all.subset.mat,
                                                                                            param = c(1,  beta.lst[[1]]))
    rej.lst[[1]] = 0
    return(list(beta.lst =  beta.lst,  Feat.lst = Feat.lst,
                last.feat.mat = Feat.lst, anchor.feat.mat = Feat.lst,
                mu.lst =  mu.lst,  Sigma.lst =  Sigma.lst,
                likelihood.estep.lst = likelihood.estep.lst, likelihood.mstep.lst = likelihood.mstep.lst,
                rej.lst = rej.lst))
}


# Function 2 warm.start(...)
#   For k=2 graph, find latent features X by maximizing joint likelihood of graph
# Input
#       default.mat is  a prescribed solution when warm start fails
# Return
#      Feat mat X estimated from mds() and CVXR
#      if  CVXR method fails, a normalized random gaussian mat X is returned.
warm.start.for.Feat = function( beta.init,
                                default.mat, tot.author, count.data){
        ss = mle.ss.cvxhull.exp(beta = beta.init, count.data)
        if( is.na(sum(ss))){
            msg("Error: NA value from initial run of  mmle.ss.cvxhull.exp() and a random standarized feature mat is returned instead!")
            X.mat =  default.mat
        } else {
            feat.mle = mds.latent.pos.mat.exponential.cvxhull(ss.vec = ss, total = tot.author, k = 2)  #### ONly for k = 2
            X.mat = standerize.matrix(feat.mle)
        }
    return(X.mat)
}


# Obsolete Function 3 simulated E-step
#  mcem_step_cpp() and its variations are always prefered to mcem_estep()
mcem_estep = function(  em.round,min.sweeps,
                        proposal.individual.feat.vec, Sigma.prop,
                        res.from.last.step,
                        tot.author, all.subset.mat, data.s,
                        if.func.message.suppressed = T){
    rej = 0
    beta.last = res.from.last.step$beta.lst[[em.round-1]]
    if (em.round < 2+1e-1){
            X.mat = res.from.last.step$Feat.lst[[em.round - 1]] # consider when em.round = 1
        } else{
            #mats.from.last = res.from.last.step$Feat.lst[[em.round-1]]
            #X.mat = mats.from.last[[len(mats.from.last)]]  # pick up the last feat.mat from res.from.last.step$Feat.lst
            X.mat = res.from.last.step$last.feat.mat
        }
    X.pivot = X.mat
    feat.sweep.list = vector("list", length = max(min.sweeps, em.round^2))
    prev.log.likelihood.shape.convexhull = res.from.last.step$likelihood.mstep.lst[[em.round-1]] #### ?????????????????????????????????????????????? prev.log.likelihood.shape.convexhull using lkhd with latest X.mat
    for (j in 1: max(min.sweeps, em.round^2) ){ # outter loop to repeat sweeps over all features
        #for (i in sample(1:tot.author, size = tot.author, replace = F)){ # ruffle rows to update
        for (i in 1: tot.author){
            x.old = X.mat[i,]
            x.new = proposal.individual.feat.vec(prev.vec = x.old) #  random walk proposal, i.e. x.new = x.old + N_2(0,\Sigma_{prop})
            X.new = X.mat
            X.new[i,] = x.new
            X.new = standerize.matrix(X.new)  #### It seems too strict !!!!!
            x.new = X.new[i,]
            current.log.likelihood.shape.convexhull = log_likelihood_shape_cvxhull_gamma(data.s,
                                                                                        X.new,
                                                                                        all.subset.mat,
                                                                                        c(1, beta.last))
            log.ratio = current.log.likelihood.shape.convexhull + dmnorm(x = x.new, mean = x.old, varcov = Sigma.prop, log = T)  ####
            log.ratio = log.ratio - prev.log.likelihood.shape.convexhull  - dmnorm(x = x.old, mean = x.old, varcov = Sigma.prop, log = T) #### ?????????????????? f(x.old) = dmnorm(x = x.old, mean = rep(0,p.demo), varcov = Sigma.prop, log = T)
            #log.ratio = log.ratio - prev.log.likelihood.shape.convexhull  - dmnorm(x = x.old, mean = rep(0,2), varcov = diag(rep(var.demo,2)), log = T)
            valid = safelog(runif(1)) < log.ratio
            if (valid){
                X.mat = X.new
                prev.log.likelihood.shape.convexhull = current.log.likelihood.shape.convexhull
            } else {
                rej = rej + 1
            }
        }
        feat.sweep.list[[j]] = find.procrustean.transform(std.Z0 = X.pivot , std.Z = X.mat)
        feat.sweep.list[[j]] = X.new
        if(!if.func.message.suppressed){
            if (j %% 10 == 0){
                elkhd = log.likelihood.shape.convexhull(count.vec = data.s, feat.mat = feat.sweep.list[[j]] ,
                                                        subset.mat = all.subset.mat, which.shape =  'gamma',
                                                        param.vec = c(1, beta.last))
                msg(sprintf('MCEM-%ith outer iteration : no.%i sweep likelihood %f',em.round, j,elkhd))
            }
        }
    }
    if(!if.func.message.suppressed){
        msg(sprintf("MCEM-%ith outer iteration:  E-step is done and M-step is starting ",em.round))
    }
    Feat.lst = feat.sweep.list[ceiling(max(min.sweeps, em.round^2)/2): max(min.sweeps, em.round^2)] # sweep burn-in: throw away roughly fisrt [m/2] feat matrices
    #print(res.Feats.lst[[t]])
    likelihood.estep = mean( sapply(Feat.lst, function(X) log_likelihood_shape_cvxhull_gamma(data.s ,
                                                                                            X,
                                                                                            all.subset.mat,
                                                                                            c(1, beta.last))))
    return(list(Feat.lst = Feat.lst, likelihood.estep = likelihood.estep, rej = rej, last.feat.mat = X.mat))
}


# Function 4 mstep in aggregate mcem
aggr_mcem_mstep = function(em.round, res.from.last.step,
                    mstep.lbd, mstep.ubd,
                    tot.author, all.subset.mat, data.s,
                    if.func.message.suppressed = T){
    beta.last = res.from.last.step$beta.lst[[em.round-1]]
    obj.func = function(x){
        mean(sapply(res.from.last.step$Feat.lst, function(X) -log_likelihood_shape_cvxhull_gamma(data.s,
                                                                                                            X,
                                                                                                            all.subset.mat,
                                                                                                            c(1,x))))}
    m.step.optimization  = optim(par = beta.last , fn = obj.func, method='Brent', lower = mstep.lbd, upper = mstep.ubd)
    beta.est = m.step.optimization$par
    likelihood.mstep = -m.step.optimization$value
    msg(sprintf('%ith round: lkhd.estep = %s, beta = %.3f, lkhd.mstep = %s',
                    em.round,
                    formatC(res.from.last.step$likelihood.estep.lst[[em.round]], format = "e", digits = 5),
                    beta.est,
                    formatC(likelihood.mstep , format = "e", digits = 5)))
    if(!if.func.message.suppressed){
        msg(sprintf("aggr MCEM-%ith outer iteration:  M-step is done",t))
    }
    return(list(beta = beta.est, likelihood.mstep = likelihood.mstep))
}


# Function 5 mstep in adpt mcem
adpt_mcem_mstep = function(em.round, res.from.last.step,
                    mstep.lbd, mstep.ubd,
                    tot.author, all.subset.mat, data.s,
                    if.func.message.suppressed = T){
    beta.last = res.from.last.step$beta.lst[[em.round-1]]
    #lkhd = sapply( res.from.last.step$Feat.lst[[em.round]], function(X) log_likelihood_shape_cvxhull_gamma(data.s,
    #                                                                                                        X,
    #                                                                                                        all.subset.mat,
    #
    #
    #                                                                                                          c(1, beta.last)))
    lkhd = sapply(res.from.last.step$Feat.lst, function(X) log_likelihood_shape_cvxhull_gamma(data.s,
                                                                                                X,
                                                                                                all.subset.mat,
                                                                                                c(1, beta.last)))
    #print(lkhd)
    best.ind = which(lkhd == max(lkhd))
    if (len(best.ind) > 1){
            best.ind = best.ind[1]
    }
    X.best = res.from.last.step$Feat.lst[[best.ind]]
    obj.func = function(x){
        -log_likelihood_shape_cvxhull_gamma(data.s,
                                            X.best,
                                            all.subset.mat,
                                            c(1,x))}
    m.step.optimization  = optim(par = beta.last , fn = obj.func, method='Brent', lower = mstep.lbd, upper = mstep.ubd)
    beta.est = m.step.optimization$par
    likelihood.mstep = -m.step.optimization$value
    msg(sprintf('%ith round: lkhd.estep = %s, beta = %.3f, lkhd.mstep = %s',
                    em.round,
                    formatC(res.from.last.step$likelihood.estep.lst[[em.round]], format = "e", digits = 5),
                    beta.est,
                    formatC(likelihood.mstep , format = "e", digits = 5)))
    if(!if.func.message.suppressed){
        msg(sprintf("adpt MCEM-%ith outer iteration:  M-step is done",t))
    }
    return(list(beta = beta.est, likelihood.mstep = likelihood.mstep))
}




### Function 6 MCEM four variations

# Function 6.1
adpt_mcem_uniform_hg_a1b2c1 = function( num.em.iter,
                                    if.warm.start, default.latent.mat = NULL,
                                    beta.init = NULL, random.start.var, mu.init, Sigma.init, tot.author,
                                    proposal.individual.feat.vec, Sigma.prop,
                                    true.beta,
                                    min.sweeps,
                                    mstep.lbd, mstep.ubd,
                                    global.subset.mat, global.count.vec,
                                    if.func.message.suppressed = F){
    #initialize a list of lists to store mcem result
    res = initialize_mcem(default.latent.mat.sd = default.latent.mat, if.warm.start,
                                true.beta,
                                beta.init, random.start.var, mu.init, Sigma.init,
                                tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
    msg(sprintf(' *** robust MCEM on %i-unif hg is starting',dim(global.subset.mat)[2]))
    msg(sprintf('Initial round: lkhd = %s, beta  = %.3f',
                formatC(res$likelihood.estep.lst[[1]], format = "e", digits = 5),
                res$beta.lst[[1]]))
    for(t in 2:num.em.iter ){
        # E-step
        update.from.estep = mcem_estep_cpp_a1b2c1(t, min.sweeps,
                                            Sigma.prop,
                                            res,
                                            tot.author, global.subset.mat, global.count.vec,
                                            T)
        #msg('estep is done')
        # record update from E-step
        #res$Feat.lst[[t]]  = update.from.estep$Feat.lst # keep all Feat.lst
        res$Feat.lst  = update.from.estep$Feat.lst                                  # keep only last Feat.lst
        #print(res$Feat.lst)
        res$likelihood.estep.lst[[t]] =  update.from.estep$likelihood.estep
        #warm.feat.combine = do.call(rbind, res$Feat.lst[[t]])
        warm.feat.combine = do.call(rbind, res$Feat.lst)
        res$mu.lst[[t]] = colMeans(warm.feat.combine)
        res$Sigma.lst[[t]] = 1/dim(warm.feat.combine)[1]*t(warm.feat.combine -  res$mu.lst[[t]])%*%(warm.feat.combine -  res$mu.lst[[t]])
        res$rej.lst[[t]] = update.from.estep$rej/tot.author/max(min.sweeps, t^2)
        res$last.feat.mat = update.from.estep$last.feat.mat
        # M-step
        #print('m is starting')
        update.from.mstep = adpt_mcem_mstep(em.round = t, res.from.last.step = res,
                                            mstep.lbd, mstep.ubd,
                                            tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
        #print('m is done')
        # record update from M-step
        res$beta.lst[[t]] = update.from.mstep$beta
        res$likelihood.mstep.lst[[t]] = update.from.mstep$likelihood.mstep
    }
    msg(("*** robust MCEM-a1b2c1 is finished"))
    return(res)
}








# Function 6.2
adpt_mcem_uniform_hg_a1b2c3 = function( num.em.iter,
                                    if.warm.start, default.latent.mat = NULL,
                                    beta.init = NULL, random.start.var, mu.init, Sigma.init, tot.author,
                                    proposal.individual.feat.vec, Sigma.prop,
                                    true.beta,
                                    min.sweeps,
                                    mstep.lbd, mstep.ubd,
                                    global.subset.mat, global.count.vec,
                                    if.func.message.suppressed = F){
    #initialize a list of lists to store mcem result
    res = initialize_mcem(default.latent.mat.sd = default.latent.mat, if.warm.start,
                                true.beta,
                                beta.init, random.start.var, mu.init, Sigma.init,
                                tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
    msg(sprintf(' *** robust MCEM on %i-unif hg is starting',dim(global.subset.mat)[2]))
    msg(sprintf('Initial round: lkhd = %s, beta  = %.3f',
                formatC(res$likelihood.estep.lst[[1]], format = "e", digits = 5),
                res$beta.lst[[1]]))
    for(t in 2:num.em.iter ){
        # E-step
        update.from.estep = mcem_estep_cpp_a1b2c3(t, min.sweeps,
                                            Sigma.prop,
                                            res,
                                            tot.author, global.subset.mat, global.count.vec,
                                            T)
        #msg('estep is done')
        # record update from E-step
        #res$Feat.lst[[t]]  = update.from.estep$Feat.lst # keep all Feat.lst
        res$Feat.lst  = update.from.estep$Feat.lst                                  # keep only last Feat.lst
        #print(res$Feat.lst)
        res$likelihood.estep.lst[[t]] =  update.from.estep$likelihood.estep
        #warm.feat.combine = do.call(rbind, res$Feat.lst[[t]])
        warm.feat.combine = do.call(rbind, res$Feat.lst)
        res$mu.lst[[t]] = colMeans(warm.feat.combine)
        res$Sigma.lst[[t]] = 1/dim(warm.feat.combine)[1]*t(warm.feat.combine -  res$mu.lst[[t]])%*%(warm.feat.combine -  res$mu.lst[[t]])
        res$rej.lst[[t]] = update.from.estep$rej/tot.author/max(min.sweeps, t^2)
        res$last.feat.mat = update.from.estep$last.feat.mat
        # M-step
        #print('m is starting')
        update.from.mstep = adpt_mcem_mstep(em.round = t, res.from.last.step = res,
                                            mstep.lbd, mstep.ubd,
                                            tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
        #print('m is done')
        # record update from M-step
        res$beta.lst[[t]] = update.from.mstep$beta
        res$likelihood.mstep.lst[[t]] = update.from.mstep$likelihood.mstep
    }
    msg(("*** robust MCEM-a1b2c3 is finished"))
    return(res)
}


# Function 6.3
aggr_mcem_uniform_hg_a1b2c1 =  function(   num.em.iter,
                                    if.warm.start, default.latent.mat = NULL,
                                    beta.init = NULL, random.start.var, mu.init, Sigma.init, tot.author,
                                    proposal.individual.feat.vec, Sigma.prop,
                                    true.beta,
                                    min.sweeps,
                                    mstep.lbd, mstep.ubd,
                                    global.subset.mat, global.count.vec,
                                    if.func.message.suppressed = F){

    #initialize a list of lists to store mcem result
    res = initialize_mcem(default.latent.mat.sd = default.latent.mat, if.warm.start,
                                true.beta,
                                beta.init = beta.init, random.start.var, mu.init, Sigma.init,
                                tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
    msg(sprintf(' *** aggr MCEM on %i-unif hg is starting',dim(global.subset.mat)[2]))
    msg(sprintf('Initial round: lkhd = %s, beta  = %.3f',
                formatC(res$likelihood.estep.lst[[1]], format = "e", digits = 5),
                res$beta.lst[[1]]))
    for(t in 2:num.em.iter ){
        # E-step

        update.from.estep = mcem_estep_cpp_a1b2c1(t, min.sweeps,
                                    Sigma.prop,
                                    res,
                                    tot.author, global.subset.mat, global.count.vec,
                                    T)
        #msg('estep is done')
        # record update from E-step
        res$Feat.lst  = update.from.estep$Feat.lst
        res$likelihood.estep.lst[[t]] =  update.from.estep$likelihood.estep
        warm.feat.combine = do.call(rbind, res$Feat.lst)
        res$mu.lst[[t]] = colMeans(warm.feat.combine)
        res$Sigma.lst[[t]] = 1/dim(warm.feat.combine)[1]*t(warm.feat.combine -  res$mu.lst[[t]])%*%(warm.feat.combine -  res$mu.lst[[t]])
        res$rej.lst[[t]] = update.from.estep$rej/tot.author/max(min.sweeps, t^2)
        res$last.feat.mat = update.from.estep$last.feat.mat
        # M-step
        update.from.mstep = aggr_mcem_mstep(em.round = t, res.from.last.step = res,
                                            mstep.lbd, mstep.ubd,
                                            tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
        # record update from M-step
        res$beta.lst[[t]] = update.from.mstep$beta
        res$likelihood.mstep.lst[[t]] = update.from.mstep$likelihood.mstep
    }
    msg(("*** aggr MCEM-a1b2c1  is finished"))
    return(res)
}

# Function 6.4
aggr_mcem_uniform_hg_a1b2c3 =  function(   num.em.iter,
                                    if.warm.start, default.latent.mat = NULL,
                                    beta.init = NULL, random.start.var, mu.init, Sigma.init, tot.author,
                                    proposal.individual.feat.vec, Sigma.prop,
                                    true.beta,
                                    min.sweeps,
                                    mstep.lbd, mstep.ubd,
                                    global.subset.mat, global.count.vec,
                                    if.func.message.suppressed = F){

    #initialize a list of lists to store mcem result
    res = initialize_mcem(default.latent.mat.sd = default.latent.mat, if.warm.start,
                                true.beta,
                                beta.init = beta.init, random.start.var, mu.init, Sigma.init,
                                tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
    msg(sprintf(' *** aggr MCEM on %i-unif hg is starting',dim(global.subset.mat)[2]))
    msg(sprintf('Initial round: lkhd = %s, beta  = %.3f',
                formatC(res$likelihood.estep.lst[[1]], format = "e", digits = 5),
                res$beta.lst[[1]]))
    for(t in 2:num.em.iter ){
        # E-step

        update.from.estep = mcem_estep_cpp_a1b2c3(t, min.sweeps,
                                    Sigma.prop,
                                    res,
                                    tot.author, global.subset.mat, global.count.vec,
                                    T)
        #msg('estep is done')
        # record update from E-step
        res$Feat.lst  = update.from.estep$Feat.lst
        res$likelihood.estep.lst[[t]] =  update.from.estep$likelihood.estep
        warm.feat.combine = do.call(rbind, res$Feat.lst)
        res$mu.lst[[t]] = colMeans(warm.feat.combine)
        res$Sigma.lst[[t]] = 1/dim(warm.feat.combine)[1]*t(warm.feat.combine -  res$mu.lst[[t]])%*%(warm.feat.combine -  res$mu.lst[[t]])
        res$rej.lst[[t]] = update.from.estep$rej/tot.author/max(min.sweeps, t^2)
        res$last.feat.mat = update.from.estep$last.feat.mat
        # M-step
        update.from.mstep = aggr_mcem_mstep(em.round = t, res.from.last.step = res,
                                            mstep.lbd, mstep.ubd,
                                            tot.author, all.subset.mat = global.subset.mat, data.s = global.count.vec)
        # record update from M-step
        res$beta.lst[[t]] = update.from.mstep$beta
        res$likelihood.mstep.lst[[t]] = update.from.mstep$likelihood.mstep
    }
    msg(("*** aggr MCEM-a1b2c3 is finished"))
    return(res)
}




################################################################################################################
# Obsolete
################################################################################################################




