helper_bayes.R
####################################################################################################################
### Part IV shape convex hull model: Bayesian Inference
####################################################################################################################

### Function 4.1 mh sampler for toy gamma param's post distribution  f(a,b | X)
mh.sampler.gamma.cvxhull = function(sample.from.proposal, log.target, proposal, theta.init, max.iter){
    res = matrix(0, 2 , max.iter)
    res[,1] = theta.init
    rej = 0
    #msg('One chain is starting.')
        for (t in 2: max.iter){
            res[,t] = res[,t-1]
            theta.tilde = sample.from.proposal(theta.prev = res[,t])
            #print(theta.tilde)
            if (sum(theta.tilde > 0) < 1.5 ){
                rej = rej + 1
                next
            } else {
                valid = safelog(runif(1)) < log.target(theta.tilde) - log.target(res[,t])
                #print(log.target(theta.tilde))
                #print(valid)
                if(is.na(valid)){
                    rej = rej + 1
                    next
                } else {
                    if(valid){
                        res[,t] = theta.tilde
                    } else {
                        rej = rej + 1
                    }
                }
            }
            if (t %% 10000 == 0){
                msg(sprintf('%.1f percert of sampling is done.',  100*(t/max.iter)))
            }
            }
        #msg('experiment is finisded.')
        list(sample = res, rejection.rate = rej/(max.iter-1))
}

#### Function 4.2 basic sample plot
# Output: a list of ggplots
plot.samples = function(sample.mat, burn.out){
    if (is.null(dim(sample.mat))){
        p = 1
        sample.list  = list()
            df = data.frame(iter = burn.out, y = sample.mat[burn.out], type = paste0('beta'))
            plot <- ggplot(df, aes(iter, y, group = type)) +
            geom_line(aes(linetype = type, colour = type)) +
            scale_color_manual(name = '' , values = c(light.blue, dark.blue)) +
            scale_linetype_manual(name = '', values = c('solid', 'dashed')) +
            theme_linedraw() +
            labs( x = "iter.", y = " samples ", colour = " " ) +
            theme( #element_line(colour = "grey90" ,linetype = 'dotted'),
                    panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                    panel.grid.minor = element_line(colour = "white", linetype = 'dotted'),
                    panel.border = element_rect(colour = "grey", fill=NA,  linetype = 'dotted'),
                    legend.title = element_blank(),
                    legend.background = element_rect(fill=alpha('grey', 0.3)),
                    legend.box.background = element_rect(colour = "white"),
                    legend.position = c(.9, .9))
            sample.list[[1]] = plot
    } else{
        p = dim(sample.mat)[1]
        sample.list  = list()
        for (j in 1: p){
            df = data.frame(iter = (1:p)[burn.out], y = sample.mat[j,burn.out], type = paste0('theta',j))
            plot <- ggplot(df, aes(iter, y, group = type)) +
            geom_line(aes(linetype = type, colour = type)) +
            scale_color_manual(name = '' , values = c(light.blue, dark.blue)) +
            scale_linetype_manual(name = '', values = c('solid', 'dashed')) +
            theme_linedraw() +
            labs( x = "iter.", y = " samples ", colour = " " ) +
            theme(  #element_line(colour = "grey90" ,linetype = 'dotted'),
                    panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                    panel.grid.minor = element_line(colour = "white", linetype = 'dotted'),
                    panel.border = element_rect(colour = "grey", fill=NA,  linetype = 'dotted'),
                    legend.title = element_blank(),
                    legend.background = element_rect(fill=alpha('grey', 0.3)),
                    legend.box.background = element_rect(colour = "white"),
                    legend.position = c(.9, .9))
            sample.list[[j]] = plot
        }
    }
    return(sample.list)
}

#### Function 4.3 histogram plot of samples
plot.histogram = function(sample.mat, burn.out){
    if (is.null(dim(sample.mat))){
        p = 1
        hist.list  = list()
            df = data.frame(y = sample.mat[burn.out])
            plot = ggplot(df, aes(y, fill = paste0('beta')))
            plot = plot + geom_histogram(bins = 30, color = "white") + scale_fill_manual( name=" ", values= light.blue, labels=paste0('beta'))
            plot = plot + theme_linedraw()
            plot = plot + labs( x = "samples", y = "hist.", colour = " " )
            plot = plot + theme(#element_line(colour = "grey10" , linetype = 'dotted'), #gridline
                                panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                                panel.grid.minor = element_line(colour = "white", linetype = 'dotted'),
                                panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                                legend.title = element_blank(),
                                legend.position = c(.9, .9),
                                legend.background = element_rect(fill=alpha('grey', 0.3)),
                                legend.box.background = element_rect(colour = "white"))
            #plot
            hist.list[[1]] = plot
    } else{
        p = dim(sample.mat)[1]
        hist.list  = list()
        for (j in 1: p){
            df = data.frame(y = sample.mat[j,burn.out])
            plot = ggplot(df, aes(y, fill = paste0('theta ',j)))
            plot = plot + geom_histogram(bins = 30, color = "white") + scale_fill_manual( name=" ", values= light.blue, labels=paste0('theta ',j))
            plot = plot + theme_linedraw()
            plot = plot + labs( x = "samples", y = "hist.", colour = " " )
            plot = plot + theme(#element_line(colour = "grey10" , linetype = 'dotted'), #gridline
                                panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                                panel.grid.minor = element_line(colour = "white", linetype = 'dotted'),
                                panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                                legend.title = element_blank(),
                                legend.position = c(.9, .9),
                                legend.background = element_rect(fill=alpha('grey', 0.3)),
                                legend.box.background = element_rect(colour = "white"))
            #plot
            hist.list[[j]] = plot
        }
    }
    return(hist.list)
}
#plot.histogram(sample.mat = sample.mat, burn.out = burn.out)

#### Function 4.3 acf plot of samples
plot.acf = function(sample.mat, burn.out){
    if (is.null(dim(sample.mat))){
        p = 1
        lag.max = 30000
        ci = 0.95
        acf.list = list()
        data = sample.mat
        list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
        N <- as.numeric(list.acf$n.used)
        df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
        df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
        df1$lag.acf[2] <- 0
        df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
        df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
        df1$acfstd[1] <- 0
        df1 <- select(df1, lag, acf, acfstd)
        plot <- ggplot(data = df1, aes( x = lag, y = acf, fill = paste0(''))) +
                scale_fill_manual( name=" ", values= light.blue , labels=paste0('beta'))+
                geom_line(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), color = dark.green, linetype = 'dashed') +
                geom_line(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), color = dark.green, linetype = 'dashed') +
                geom_col(width = 0.7)
        plot = plot + theme_linedraw()
        plot = plot + labs( x = "lag", y = "acf", colour = " " )
        plot = plot + theme(panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                            panel.grid.minor = element_line(colour = "white", linetype = 'dotted'), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            legend.position = c(.9, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white"))
        acf.list[[1]] = plot
    } else{
        p = dim(sample.mat)[1]
        lag.max = 30000
        ci = 0.95
        acf.list = list()
        for (j in 1: p){
        data = sample.mat[j,]
        list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
        N <- as.numeric(list.acf$n.used)
        df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
        df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
        df1$lag.acf[2] <- 0
        df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
        df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
        df1$acfstd[1] <- 0
        df1 <- select(df1, lag, acf, acfstd)
        plot <- ggplot(data = df1, aes( x = lag, y = acf, fill = paste0('theta ',j))) +
                        scale_fill_manual( name=" ", values= light.blue , labels=paste0('theta ',j))+
                        geom_line(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), color = dark.green, linetype = 'dashed') +
                        geom_line(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), color = dark.green, linetype = 'dashed') +
                        geom_col(width = 0.7)
        plot = plot + theme_linedraw()
        plot = plot + labs( x = "lag", y = "acf", colour = " " )
        plot = plot + theme(panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                            panel.grid.minor = element_line(colour = "white", linetype = 'dotted'), #gridline
                            panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                            legend.title = element_blank(),
                            legend.position = c(.9, .9),
                            legend.background = element_rect(fill=alpha('grey', 0.3)),
                            legend.box.background = element_rect(colour = "white"))
        acf.list[[j]] = plot
    }
    }
    return(acf.list)
}

#### Function 4.4  for every parameter, plot samples, histogram and acf
three.diag.plots.samples = function(sample.mat, burn.out, if.show.plot = F, if.save.plot = T){
    if (is.null(dim(sample.mat))){
            num.param = 1
        } else{
            num.param = dim(sample.mat)[1]
        }
    sample.list = plot.samples(sample.mat = sample.mat, burn.out = burn.out)
    msg('sample plot is done.')
    hist.list = plot.histogram(sample.mat = sample.mat, burn.out = burn.out)
    msg('histogram is done.')
    acf.list = plot.acf(sample.mat = sample.mat, burn.out = burn.out)
    msg('acf is done.')
    mylist = list()
    for (i in 1:(3 * num.param)){
        if (i %% 3 == 1){
            mylist[i] = sample.list[1+ i%/%3]
        } else if (i %% 3 == 2){
            mylist[i] = hist.list[1+ i%/%3]
        } else {
            #print(i)
            mylist[i] = acf.list[i%/%3]
        }
    }
    if(if.show.plot){
        grid.arrange(grobs = mylist, nrow = num.param)
    }
    if(if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        g <- arrangeGrob(grobs = mylist, nrow=2)
        ggsave(filename = paste0("output/three-plots-mh-samples-", current.time, ".pdf"), g)
    }
}


### Function 4.5 plot sample path for multiple chains
# run multiple chain
multi.chain.mh.sampler.gamma.cvxhull= function(max.iter, init.position.list, if.save.result = T){
    k = length(init.position.list)
    multi.chain.list = vector(mode = "list", length = k)
    rej.vec = rep(0,k)
    for (j in 1:k){
        msg(sprintf('the %ith chain is starting***', j))
        res =  mh.sampler.gamma.cvxhull(sample.from.proposal = sample.from.proposal.den,
                                                        log.target = log.target,
                                                        proposal = proposal,
                                                        theta.init = init.positions.list[[j]],
                                                        max.iter = max.iter)
        multi.chain.list[[j]]=res$sample # every list is a 2 by max.iter matrix, of which each row stores one parameter's samples
        rej.vec[j] = res$rejection.rate
        msg(sprintf('the %ith chain is done***', j))
    }
    if(if.save.result){
        current.time = format(Sys.time(),"%b-%d-%H-%M")
        save(multi.chain.list, file = paste0('./output/multicahin-mh-res-',current.time,'.Rdata'))
        msg('multiple chains from mh sampler are saved safely')
    }
    return(list(multi.chain = multi.chain.list, rej = rej.vec))
}

#Function 4.6 Extract a bunred and thined samples for plotting
get.df.from.multichain = function(multi.chain.list, thin.vec){
    maxiter = dim(multi.chain.list[[1]])[2]
    dt = data.frame()
    k = length(multi.chain.list)
    for (j in 1:k ){
        dt = rbind(dt, data.frame(x = (multi.chain.list[[j]])[1,thin.vec] , y = (multi.chain.list[[j]])[2,thin.vec], ind = paste0('path',j)))
    }
    dt
}

# Function 4.7 Plor dataframe from function 4.6
multichain.plot = function(df, if.show.plot, if.save.plot){
    k = length(unique(df[,3]))
    n = length(thin.vec)
    plot <- ggplot(df)
    plot <- plot+ geom_line(aes(x= x, y = y, colour= ind), size= .11)
    if (k > 1+ 1e-10){
        starts = seq(0,k-1)*n + 1
    } else {
        starts = 1
    }
    plot <- plot+ geom_point(data = df[starts,],
                            aes(x =x, y =y, colour = ind),
                            shape = 18,
                            size = 3)
    plot <- plot + ylab("beta") + xlab("alpha")  + theme_linedraw()
    plot <- plot + theme(panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.grid.minor = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                        legend.title = element_blank(),
                        legend.position = c(.9, .9),
                        legend.background = element_rect(fill=alpha('grey', 0.3)),
                        legend.box.background = element_rect(colour = "white"))
    if(if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("./output/multichain-path-mh", current.time, ".pdf"), plot)

    }
    if(if.show.plot){
        plot
    }
}

# Function 4.8
# Input:
#chain: a list of mcmc sample paths as matrices
# warm.up : a seq to be discared from sample path
# num.subchain: num of subchains splitting one path
find.Rhat.from.chain.list = function( chain , warm.up, num.subchain){
    warmed.chain.list = remove.burn.in.from.multi.chain(chain.list = chain, burn.in = warm.up)
    warm.subchain.list = split.chain.list(chain.list = warmed.chain.list, num.subchain.per.chain = num.subchain)
    m = length(warm.subchain.list) # number of chains
    if (is.null(dim(warm.subchain.list[[1]])[1])){
        n = length(warm.subchain.list[[1]]) # length of chain
        param.mat = vector()
        for (i in 1:m){
            param.mat = rbind(param.mat, warm.subchain.list[[i]])
        }
        res = find.Rhat(mat = param.mat)
    } else {
        p = dim(warm.subchain.list[[1]])[1]
        n = dim(warm.subchain.list[[1]])[2] # length of chain
        param.list = vector("list", length = p)
        for (j in 1:p ){
            for (i in 1:m){
                param.list[[j]] = rbind(param.list[[j]], warm.subchain.list[[i]][j,])
            }
        }
        res = sapply(param.list, find.Rhat)
    }
    return(res)
}


# Function 4.9
remove.burn.in.from.multi.chain = function(chain.list, burn.in = NULL){
    k = length(chain.list)
    if (is.null(dim(chain.list[[1]]))){
        n = length(chain.list[[1]])
        if (is.null(burn.in)){
            burn.in = 1 : floor(n/2)
        }
        res = list()
        for (i in 1 : k){
        res[[i]] = chain.list[[i]][-burn.in]
        }
    } else{
        n = dim(chain.list[[1]])[2]
        if (is.null(burn.in)){
            burn.in = 1 : floor(n/2)
        }
        res = list()
        for (i in 1 : k){
        res[[i]] = chain.list[[i]][,-burn.in]
        }
    }
    res
}


# Function 4.10 split matrix into submat
matsplitter<-function(M, r, c) {
    rg <- (row(M)-1)%/%r+1
    cg <- (col(M)-1)%/%c+1
    rci <- (rg-1)*max(cg) + cg
    N <- prod(dim(M))/r/c
    cv <- unlist(lapply(1:N, function(x) M[rci==x]))
    dim(cv)<-c(r,c,N)
    cv
}


# Function 4.11 split mh.chains into multiple sub chains
# Usa.
split.chain.list = function(chain.list, num.subchain.per.chain){
    if (is.null(dim(chain.list[[1]]))){
        p = 1
        num.chain = length(chain.list)
        len.chain = length(chain.list[[1]])
        res = list()
        for (i in 1: num.chain){
            mul.submat = matsplitter(t(as.matrix(chain.list[[i]])), p, floor(len.chain/num.subchain.per.chain))
                for (j in 1:num.subchain.per.chain){
                    res[[(i-1)*num.subchain.per.chain + j]] = mul.submat[,,j]
                }
        }
    } else {
        p = dim(chain.list[[1]])[1]
        num.chain = length(chain.list)
        len.chain = dim(chain.list[[1]])[2]
        res = list()
        for (i in 1: num.chain){
            mul.submat = matsplitter(chain.list[[i]], p, floor(len.chain/num.subchain.per.chain))
                for (j in 1:num.subchain.per.chain){
                    res[[(i-1)*num.subchain.per.chain + j]] = mul.submat[,,j]
                }
        }
    }
    res
}

# Function 4.12
# mat: each row is a subchain for single parameter
find.B = function(mat){
    m = dim(mat)[1]
    n = dim(mat)[2]
    n/(m-1)*sum((rowMeans(mat) - mean(mat))^2)
}


# Function 4.13
find.W = function(mat){
    m = dim(mat)[1]
    n = dim(mat)[2]
    sum((mat - rowMeans(mat))^2)/(n-1)/m
}

# Function 4.14
find.Rhat = function(mat){
    #m = dim(mat)[1]
    n = dim(mat)[2]
    w = find.W(mat)
    b = find.B(mat)
    sqrt((w + (b-w)/n)/w)
}

# Function 4.15
find.Vt = function(v.mat, v.t){
    m = dim(v.mat)[1]
    n = dim(v.mat)[2]
    if (v.t > n - 1 + 1e-10){
        msg(' Waring: in function \'find.Vt\', the lagging periord v.t is larger than its upper bound, thus NA is returned ')
        return(NA)
    } else{
        return(sum((v.mat[,(v.t+1):n] - v.mat[,1:(n-v.t)])^2)/m/(n-v.t))
    }
}

# Function 4.16
# output rhot.hat, defined on (11.7) p286, BDA
find.rhot.hat = function(r.mat, rho.t){
    m = dim(r.mat)[1]
    n = dim(r.mat)[2]
    Vt = find.Vt(v.mat = r.mat, v.t = rho.t)
    if (is.na(Vt)){
        msg('warining: in function \'find.rhot.hat\', rho.t is too large, thus NA is returned')
        return(NA)
    } else {
        W = find.W(mat = r.mat)
        B = find.B(mat = r.mat)
        if ( W > 1e-200 && B > 1e-200){
            return(1 - Vt/(W + (B-W)/n)/2)
        } else {
        msg('warining: in function \'find.rhot.hat\' : W or B is too small, sample rejection rate is too high')
            return(1)
        }
    }
}


# Function 4.17
# Find the first period t such that rho_t + rho_{t+1} < 0, for t = 0,1,...., (num.iter-1);
find.T = function(T.mat, T.rho.vec = NULL){
    #t.mat = mat
    m = dim(T.mat)[1]
    n = dim(T.mat)[2]
    if (n < 3){
        msg('Waring: in function \' find.T\', sample size dim(T.mat)[2] is too small')
    }
    if (is.null(T.rho.vec)){
        T.rho.vec = sapply(seq(0,n-1), function(x) find.rhot.hat(r.mat = T.mat, rho.t = x))
    }
    res = NA
    for (i in seq(from = 2, to = (n-1), by=2)){
        if (T.rho.vec[i] + T.rho.vec[i+1] < 0 ){
            res = i - 1
            break
        }
    }
    if (is.na(res)){
        res = n - 1 # may not accurate !!!!!!!!!
    }
    return(res)
}

# Function 4.18
find.neffhat = function(n.mat){
    m = dim(n.mat)[1]
    n = dim(n.mat)[2]
    rho = sapply(seq(0,n-1), function(x) find.rhot.hat(r.mat = n.mat, rho.t = x))
    T = find.T(T.mat = n.mat, T.rho.vec = rho)
    (m*n)/(1 + 2* sum(rho[1:T]))
}

# Function 4.19
# See def. on p287, BDA
find.neffhat.from.chain.list = function( chain , warm.up, num.subchain){
    warm.list = remove.burn.in.from.multi.chain(chain.list = chain, burn.in = warm.up)
    warm.list = split.chain.list(chain.list = warm.list, num.subchain.per.chain = num.subchain)
    m = length(warm.list) # number of chains
    if (is.null(dim(chain[[1]]))){
        n = length(warm.list[[1]])
        param.mat = vector()
        for (i in 1:m){
            param.mat = rbind(param.mat, warm.list[[i]])
        }
        res = find.neffhat(param.mat)
    } else {
        p = dim(warm.list[[1]])[1]
        param.list = vector("list", length = p)
        for (j in 1:p ){
            for (i in 1:m){
                param.list[[j]] = rbind(param.list[[j]], warm.list[[i]][j,])
            }
        }
        res = sapply(param.list, find.neffhat)
    }
    return(res)
}

# Function 4.20  experiment.find.Rhat
experiment.find.Rhat = function(chain.list, warmup.to.discard, num.subchain.per.chain, if.save.Rhat = T){
    if (is.null(dim(chain.list[[1]]))){
        c(p, n) %tin% c(1, length(chain.list[[1]]))
        ind= seq(max(20, ceiling(max(warmup.to.discard))*1.2), n, by = 50) # index until which R.hat is calculated
        if(min(ind) < 20 ){
            msg('Warning: start calculating neff too early. So reset variable \'ind\' in function \' experiment.find.neffhat\'')
        }
        mat  = matrix(0, p, length(ind))
        msg('Start calculating R.hat ***')
        for (i in 1:length(ind)){
            list.report = lapply(chain.list, function(mat) mat[1:ind[i]])
            mat[i] = find.Rhat.from.chain.list(chain = list.report, warm.up = warmup.to.discard , num.subchain = num.subchain.per.chain)
            print(mat[i])
            msg(sprintf('%.2f perc is completed***',100*i/length(ind)))
        }
        res = list(Rhat.mat = mat, ind.report = ind)
        if (if.save.Rhat){
            current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
            save(res, file = paste0("./output/res.Rhat-",current.time,"-.Rdata"))
            msg('Save Rhat as a matrix in res$Rhat.mat')
        }
    } else {
        c(p, n) %tin% dim(chain.list[[1]])
        if (max(warmup.to.discard)*2 > n-1 ){
            msg('Warining: in function \'experiment.find.Rhat\', warmup.to.discard is too larger')
        }
        ind = seq(max(warmup.to.discard)*2, n, by = 30) # index until which R.hat is calculated
        mat = matrix(0, p, length(ind))
        msg('Start calculating R.hat ***')
        for (i in 1:length(ind)){
            list.report = lapply(chain.list, function(mat) mat[,1:ind[i]])
            mat[,i] = find.Rhat.from.chain.list(chain = list.report, warm.up = warmup.to.discard , num.subchain = num.subchain.per.chain)
            msg(sprintf('%.2f perc is completed***',100*i/length(ind)))
        }
        res = list(Rhat.mat = mat, ind.report = ind)
        if (if.save.Rhat){
            current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
            save(res, file = paste0("./output/res.Rhat-",current.time,"-.Rdata"))
            msg('Save Rhat as a matrix in res$Rhat.mat')
        }
    }
    res
}



# Function 4.21  experiment.find.neffhat
experiment.find.neffhat = function(chain.list, warmup.to.discard, num.subchain.per.chain, if.save.neffhat = T){
    if (is.null(dim(chain.list[[1]]))){
        c(p, n) %tin% c(1, length(chain.list[[1]]))
        ind= seq(max(20, max(warmup.to.discard)*1.2), n, by = 30) # index until which R.hat is calculated
        if(min(ind) < 20 ){
            msg('Warning: start calculating neff too early. So reset variable \'ind\' in function \' experiment.find.neffhat\'')
        }
        mat  = matrix(0, p, length(ind))
        msg('Start calculating neff.hat ***')
        for (i in 1:length(ind)){
            list.report = lapply(chain.list, function(mat) mat[1:ind[i]])
            mat[i] = find.neffhat.from.chain.list(chain = list.report, warm.up = warmup.to.discard , num.subchain = num.subchain.per.chain)
            print(mat[i])
            msg(sprintf('%.2f perc is completed***',100*i/length(ind)))
        }
        res = list(neffhat.mat = mat[1,], ind.report = ind)
        if (if.save.neffhat){
            current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
            save(res, file = paste0("./output/res.neffhat-",current.time,"-.Rdata"))
            msg('Save neffhat as a matrix in res$neffhat.mat')
        }
    } else {
        c(p, n) %tin% dim(chain.list[[1]])
        ind= seq(max(warmup.to.discard)*2, n, by = 30) # index until which R.hat is calculated
        mat  = matrix(0, p, length(ind))
        msg('Start calculating neff.hat ***')
        for (i in 1:length(ind)){
            list.report = lapply(chain.list, function(mat) mat[,1:ind[i]])
            mat[,i] = find.neffhat.from.chain.list(chain = list.report, warm.up = warmup.to.discard , num.subchain = num.subchain.per.chain)
            msg(sprintf('%.2f perc is completed***',100*i/length(ind)))
        }
        res = list(neffhat.mat = mat, ind.report = ind)
        if (if.save.neffhat){
            current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
            save(res, file = paste0("./output/res.neffhat-",current.time,"-.Rdata"))
            msg('Save neffhat as a matrix in res$neffhat.mat')
        }
    }
    return(res)
}


# Funtion 4.22 plot R.hat along iteration
plot.Rhat = function(R.res, if.show.plot, if.save.plot ){
    mat = R.res$Rhat.mat
    ind = R.res$ind.report
    if (is.null(dim(mat)) || dim(mat)[1] == 1){
        c(p, n) %tin% c(1, length(mat))
        dt = data.frame()
        for (i in 1: n){
                dt = rbind(dt, data.frame(x = ind[i] , y = mat[1,i], ind = 'beta'))
        }
    } else {
        c(p, n) %tin% dim(mat)
        dt = data.frame()
        for (i in 1: n){
            for(j in 1:p){
                dt = rbind(dt, data.frame(x = ind[i] , y = mat[j,i], ind = paste0('beta',j)))
            }
        }
    }
    plot <- ggplot(dt)
    plot <- plot+ geom_line(aes(x= x, y = y, colour= ind), size= .11)
    plot <- plot + geom_hline(yintercept = 1, colour = 'grey', linetype = 'dashed')
    plot <- plot + ylab("R.hat") + xlab("iteration")  + theme_linedraw()
    plot <- plot + theme(panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.grid.minor = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                        legend.title = element_blank(),
                        legend.position = c(.9, .9),
                        legend.background = element_rect(fill=alpha('grey', 0.3)),
                        legend.box.background = element_rect(colour = "white"))
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("./output/Plot-Rhat-", current.time,".pdf"), plot)
        msg('R.hat plot is saved')
    }
    if(if.show.plot){
        return(plot)
    }
}


# Funtion 4.23 plot neff.hat along iteration
plot.neffhat = function(nhat.res, if.show.plot = F, if.save.plot = T){
    mat = nhat.res$neffhat.mat
    ind = nhat.res$ind.report
    if (is.null(dim(mat)) || dim(mat)[1] == 1 ){
        c(p, n) %tin% c(1, length(mat))
        dt = data.frame()
        for (i in 1: n){
                dt = rbind(dt, data.frame(x = ind[i] , y = mat[i], ind = paste0('beta')))
        }
    } else {
        c(p, n) %tin% dim(mat)
        dt = data.frame()
        for (i in 1: n){
            for(j in 1:p){
                dt = rbind(dt, data.frame(x = ind[i] , y = mat[j,i], ind = paste0('beta',j)))
            }
        }
    }

    plot <- ggplot(dt)
    plot <- plot+ geom_line(aes(x= x, y = y, colour= ind), size= .11)
    plot <- plot + ylab("neff.hat") + xlab("iteration")  + theme_linedraw()
    plot <- plot + theme(panel.grid.major = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.grid.minor = element_line(colour = "white", linetype = 'dotted'), #gridline
                        panel.border = element_rect(colour = "grey", fill=NA,linetype = 'dotted'),
                        legend.title = element_blank(),
                        legend.position = c(.1, .9),
                        legend.background = element_rect(fill=alpha('grey', 0.3)),
                        legend.box.background = element_rect(colour = "white"))
    if (if.save.plot){
        current.time = format(Sys.time(), "%b-%d-%I-%M-%p")
        ggsave(filename = paste0("output/Plot-neffhat-", current.time,".pdf"), plot)
        msg('neff.hat plot is saved')
    }
    if(if.show.plot){
        return(plot)
    }
}



# Function 3.14
# Collect multiple single-chain rdata files and turn them  into a  matrix or dataframe for later use
combine.chains = function(keyword, data.path = './output/'){
    setwd(data.path) # go to data folder
    files <- list.files(pattern = keyword)
    if(identical(files, character(0))){
        msg('Error from combine.chains(): keyword doesnt match any file names ')
    }
    input.lst <- lapply(files, function(x) {
        load(file = x)
        get(ls()[ls()!= "filename"])
    })
    setwd('..') # come back to parent folder
    res.mat = vector()
    for (i in 1:length(input.lst)){
        res.mat = rbind(res.mat, input.lst[[i]][[1]])
    }
    res.df = data.frame()
    for (i in 1: dim(res.mat)[1]){
        res.df = rbind(res.df, data.frame(x = 1:dim(res.mat)[2] , y = res.mat[i,], ind = paste0('path',i) ))
    }
    return(list(mat = res.mat, df = res.df))
}



# Function 3.15
# Helper function for ploting ROC curve
turn.count.to.onevsall = function(count){
    res = matrix(0, sum(count), len(count))
    row.ind = 0
    for (i in 1:len(count)){
        num.rep = count[i]
        if (num.rep > 1e-10){
            for (j in 1:num.rep){
                res[(row.ind+j), i] = 1
            }
            row.ind = row.ind + count[i]
        }
    }
    res
}


# Function 3.16  find.unique.path.weight(mat: mat[int]) -> mat[int]
# Input: every row of  matrix is a k-set
# Ouput: unique rows of input matrix with its occurence number(count number)
# Dependency: function 'unique' from R pakcage 'data.table'
find.unique.path.weight = function(mat){
    k = dim(mat)[2]
    n = dim(mat)[1]
    res = unique(mat)
    colnames(res) = paste('node',1:k,sep='_')
    wt = rep(0,dim(res)[1])
    for (i in 1: len(wt)){
        count = 0
        for (j in 1: n){
            if(sum(res[i,] == mat[j,])>k - 1e-1){
                count = count + 1
            }
        }
        wt[i] = count
    }
    return(cbind(res,wt))
}

