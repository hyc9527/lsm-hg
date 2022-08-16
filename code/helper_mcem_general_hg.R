######################################################
# Monte Carlo EM algorithm for general hypergraph
## algorithm type:
## 1: 'rbst-prcr';
## 2: 'rbst-norm';
## 3: 'aggr-prcr';
## 4: 'aggr-norm'.
#######################################################



# Function `run.mcem`
#     Run mcem for general hypergraph.
#        - The algorihtm is a wrapper for four mcems for uniform hypergraph, indiced by arg `which.mcem`
#        - note that all *.demo variables inside the function block  are global variables
run.mcem.general.hypergraph = function(which.mcem=1, isVerboseResult = TRUE, isSaveResult = FALSE)
{
    # mcem type
    mix.mu.est = mix.Sigma.est = mix.feat.est = mix.likelihood = mix.rej=list()
    mix.beta.est = beta.hat = c()
    {
        alg.lst = c('rbst-prcr', 'rbst-norm','aggr-prcr','aggr-norm')
        mcem.chosen.demo <<- alg.lst[which.mcem]                                         #---- Pick an algorithm ---#
        print(c(mcem.chosen.demo,'is chose.'))
    }
    {# mcem iter
        mix.num.em.iter.demo <<- c(20,3,2)                                               # --- set iteration nums for each k-sets ---#
        mix.min.num.sweeps.estep.demo <<- c(10,10,1)  # num of sweeps in E-step
    }
    {# E-step specs
        tausq.prop = 1
        Sigma.prop.demo <<- tausq.prop*diag(p.demo)
        mu.init.demo <<- rep(0,2)
        Sigma.init.demo <<- n.demo*p.demo*diag(p.demo)
        mu.prop = rep(mu.demo, p.demo)
        sample.from.proposal.den = function(prev.vec = c(0,0), mu.vec = mu.prop, Sigma.mat = Sigma.prop.demo){
            prev.vec + rmnorm(1, mu.vec, Sigma.mat)
        }
    }
    {# M-step specs
        mstep.lbd.demo <<- 1e-5 # M-step optimization lower bound
        mstep.ubd.demo <<- 5  # M-step optimization upper bound
        random.start.var.demo <<- 1
        which.kset.to.run.on = mix.k.demo                                                       # ---- pick which k to run ----#
        mix.beta.init.demo <<- rep(0, len(which.kset.to.run.on))
        for (i in 1:len(which.kset.to.run.on)){
            mix.beta.init.demo[i] = max(rnorm(1,mean = mix.true.beta.demo[i], sd = sqrt(random.start.var.demo)), 1e-10)
        }
    }
    # Run on general hypergraph

    msg('#### check mcem is to run on k =',which.kset.to.run.on)
    for(i in 1:len(which.kset.to.run.on)){
        msg(paste0('MCEM on  k = ',mix.k.demo[i]))
        if( i < 1+1e-1){ # run mcem on first uniform hg
            true.beta.demo <<- mix.beta.init.demo[i]
            if (mix.k.demo[i] < 2.1)
            { # warm.start is available in  k=2 case ONLY; otherwise wart.start is not defined.
                wart.start.only.for.k2 = TRUE
            } else {
                wart.start.only.for.k2 = FALSE
            }
            res.subproc = mcem.subproc(which.algorithm = mcem.chosen.demo,
                                test.beta.init =mix.beta.init.demo[i],
                                test.num.em.iter = mix.num.em.iter.demo[i],
                                test.if.warm.start =  wart.start.only.for.k2,
                                test.default.latent.mat = NULL,
                                test.min.sweeps = mix.min.num.sweeps.estep.demo[i],
                                test.subset.mat =  mix.subset.mat.demo[[i]],
                                test.global.count.vec = mix.count.demo[[i]])
        } else{ # run mcem on other uniform hg(s), using result from last uniform hg
            res.subproc = mcem.subproc(which.algorithm = mcem.chosen.demo,
                                test.beta.init = mix.beta.init.demo[i],
                                test.num.em.iter = mix.num.em.iter.demo[i],
                                test.if.warm.start = F,
                                test.default.latent.mat = feat.hat,
                                test.min.sweeps = mix.min.num.sweeps.estep.demo[i],
                                test.subset.mat =  mix.subset.mat.demo[[i]],
                                test.global.count.vec = mix.count.demo[[i]])
        }
        # save res for next  k-uniform
        beta.est.vec = unlist(res.subproc$beta.lst)
        beta.hat = c( beta.hat,beta.est.vec[len(beta.est.vec)])
        mix.beta.est = c(mix.beta.est, beta.est.vec)
        mix.rej = c(mix.rej, res.subproc$rej.lst)
        mix.likelihood = c(mix.likelihood , res.subproc$likelihood.mstep.lst)
        mu.hat =  res.subproc$mu.lst[[len(res.subproc$mu.lst)]]
        Sigma.hat = res.subproc$Sigma.lst[[len(res.subproc$Sigma.lst)]]
        feat.hat = res.subproc$last.feat.mat
    }
    if (isVerboseResult)
    {
        print(paste0('true beta is ', mix.true.beta.demo))
	    print('beta.hat ')
	    print(beta.hat)
	    print('feat.hat')
	    print(dim(feat.hat))
	    print('rej rate along path is :')
	    print(unlist(mix.rej))
        lkhd.lst = unlist(mix.likelihood)
        print('log-likelihood along path is :')
        print(unlist(mix.likelihood))
    }
    if (isSaveResult)
    {
        current.time = format(Sys.time(), "%b-%d-%H-%M-%s")
        random.keys = sample(1:1e3,1)
        save(mix.beta.est, file = paste0("./output/result-",mcem.chosen.demo,"-mcem-feat.hat-",mix.true.beta.demo[1] ,"-",current.time,"-",random.keys,".Rdata"))
        save(mix.beta.est, file = paste0("./output/result-",mcem.chosen.demo,"-mcem-shape-param-b",mix.true.beta.demo[1] ,"-",current.time,"-",random.keys,".Rdata"))
        save(lkhd.lst, file = paste0("./output/result-",mcem.chosen.demo,"-mcem-loglikelihd-b",mix.true.beta.demo[1] ,"-",current.time,"-",random.keys,".Rdata"))
        save(mix.feat.est, file = paste0("./output/result-",mcem.type,"-mcem-feat-b",mix.true.beta.demo[1] ,"-",current.time,".Rdata"))
    }
    msg('*** MCEM for general hypergraph is finished. ')
    return(list(beta.hat=beta.hat, feat.hat=feat.hat))
}








# MCEM for uniform hypergraph
mcem.subproc = function(which.algorithm,
                        test.beta.init,
                        test.num.em.iter,
                        test.if.warm.start = FALSE,
                        test.default.latent.mat,
                        test.min.sweeps, test.subset.mat, test.global.count.vec){
    switch(which.algorithm,
        'rbst-prcr' = {print('MCEM: robust-prcr')
                                        res = adpt_mcem_uniform_hg_a1b2c1(  num.em.iter = test.num.em.iter,
                                                                            if.warm.start = test.if.warm.start, default.latent.mat = test.default.latent.mat,
                                                                            beta.init = test.beta.init, random.start.var = random.start.var.demo,
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = n.demo,
                                                                            proposal.individual.feat.vec = sample.from.proposal.den, Sigma.prop = Sigma.prop.demo,
                                                                            true.beta = true.beta.demo,
                                                                            min.sweeps = test.min.sweeps,
                                                                            mstep.lbd = mstep.lbd.demo, mstep.ubd = mstep.ubd.demo,
                                                                            global.subset.mat =  test.subset.mat , global.count.vec = test.global.count.vec,
                                                                            if.func.message.suppressed = T)
                                        #rej = check.result(res)
                                        #print(rej)
                                        #return(res)
                                        },
        'rbst-norm' = {print('MCEM: robust-norm')
                                        res = adpt_mcem_uniform_hg_a1b2c3(  num.em.iter = test.num.em.iter,
                                                                            if.warm.start = test.if.warm.start, default.latent.mat = test.default.latent.mat,
                                                                            beta.init = test.beta.init, random.start.var = random.start.var.demo,
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = n.demo,
                                                                            proposal.individual.feat.vec = sample.from.proposal.den, Sigma.prop = Sigma.prop.demo,
                                                                            true.beta = true.beta.demo,
                                                                            min.sweeps = test.min.sweeps,
                                                                            mstep.lbd = mstep.lbd.demo, mstep.ubd = mstep.ubd.demo,
                                                                            global.subset.mat =  test.subset.mat , global.count.vec = test.global.count.vec,
                                                                            if.func.message.suppressed = T)
                                        #rej = check.result(res)
                                        #print(rej)
                                        #return(res)
                                        },
        'aggr-prcr' = {print('MCEM: aggr-prcr')
                                        res = aggr_mcem_uniform_hg_a1b2c1(num.em.iter = test.num.em.iter,
                                                                            if.warm.start = test.if.warm.start, default.latent.mat = test.default.latent.mat,
                                                                            beta.init = test.beta.init, random.start.var = random.start.var.demo,
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = n.demo,
                                                                            proposal.individual.feat.vec = sample.from.proposal.den, Sigma.prop = Sigma.prop.demo,
                                                                            true.beta = true.beta.demo,
                                                                            min.sweeps = test.min.sweeps,
                                                                            mstep.lbd = mstep.lbd.demo, mstep.ubd = mstep.ubd.demo,
                                                                            global.subset.mat =  test.subset.mat , global.count.vec = test.global.count.vec,
                                                                            if.func.message.suppressed = T)
                                        #rej = check.result(res)
                                        #print(rej)
                                        #return(res)
                                        },
        'aggr-norm' = {print('MCEM: robust-norm')
                                        res = adpt_mcem_uniform_hg_a1b2c3(  num.em.iter = test.num.em.iter,
                                                                            if.warm.start = test.if.warm.start, default.latent.mat = test.default.latent.mat,
                                                                            beta.init = test.beta.init, random.start.var = random.start.var.demo,
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = n.demo,
                                                                            proposal.individual.feat.vec = sample.from.proposal.den, Sigma.prop = Sigma.prop.demo,
                                                                            true.beta = true.beta.demo,
                                                                            min.sweeps = test.min.sweeps,
                                                                            mstep.lbd = mstep.lbd.demo, mstep.ubd = mstep.ubd.demo,
                                                                            global.subset.mat =  test.subset.mat , global.count.vec = test.global.count.vec,
                                                                            if.func.message.suppressed = T)
                                        #rej = check.result(res)
                                        #print(rej)
                                        #return(res)
                                        },
        {stop('the algorithm you picked up is not available for now.')}
    )
}



plot.lkhd = function(mat){
        bvec=seq(0.01,2,  length.out  = 100)
        log_lkhd = sapply(bvec, function(x) log_likelihood_shape_cvxhull_gamma(mix.count.demo[[1]], mat, mix.subset.mat.demo[[1]], c(1, x)))
        #log.lkhd = sapply(bvec, function(x) log.likelihood.shape.convexhull(mix.count.demo[[1]], feat.mat.demo, mix.subset.mat.demo[[1]], "gamma",  c(1, x)))
        #pdf(file="check_likelihood.pdf")
        plot(bvec, log_lkhd)
        ind = which.max(log_lkhd)
        points(bvec[ind], log_lkhd[ind], col = "red")
        abline(v = bvec[ind], lty = 2, col = 'red')
}



check.result = function(res, if.plot.lkhd = F){
    bvec=seq(0.01,2,  length.out  = 100)
    log_lkhd = sapply(bvec, function(x) log_likelihood_shape_cvxhull_gamma(mix.count.demo[[1]], res$last.feat.mat, mix.subset.mat.demo[[1]], c(1, x)))
    #log.lkhd = sapply(bvec, function(x) log.likelihood.shape.convexhull(mix.count.demo[[1]], feat.mat.demo, mix.subset.mat.demo[[1]], "gamma",  c(1, x)))
    #pdf(file="check_likelihood.pdf")
    if (if.plot.lkhd){
        plot(bvec, log_lkhd)
    }
    ind = which.max(log_lkhd)
    points(bvec[ind], log_lkhd[ind], col = "red")
    abline(v = bvec[ind], lty = 2, col = 'red')
    print('rejection rate is :')
    return(unlist(res$rej.lst))
}
