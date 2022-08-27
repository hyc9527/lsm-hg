######################################################
# Analytics on Star War IV: New Hope
# General hypergraph
    # node set: actors
    # edge set: co-occurence of at least two actors in one scence
#######################################################

# Load dataset
    # Input:
        # csv file: new_hope_adjMat.csv, new_hope_actor_lst.csv
    # Output:
        # int num.author, number of nodes
        # int num.paper, number of hyperedges
        # mat authorPaperBiadj, of shape n(num. of nodes) by m (num. of hyoeredges)
load.dataset.starwar = function(){
    # Import starwar dataset
    authorPaperBiadj  <<- read.csv(file='../data/starwar_actor_scene_adjmat.csv', header = F)            #--- author-paper-adjacency-mat ---#
    author.real.names.vec = read.csv(file = '../data/starwar_actor_lst.csv', header = F)   #--- author name
    author.real.names.vec.demo <<- as.vector(author.real.names.vec[,1])
    num.author <<- dim(authorPaperBiadj )[1]
    num.paper <<- dim(authorPaperBiadj )[2]
    colnames(authorPaperBiadj ) = paste0("p", 1:num.paper)
    rownames(authorPaperBiadj ) = paste0("a", 1:num.author)
    # Summary on star war hypergraph
    print('###### Summary on dataset Star War IV ########')
    print('Number of actors (nodes),  number of scenes (hyperedges) are respectively')
    print(dim(authorPaperBiadj))
    print('Main actors are sorted in decreasing order of nodal degree : ')
    deg=rowSums(authorPaperBiadj)
    names(deg) = author.real.names.vec.demo
    print(sort(deg,decreasing=T))
    print(paste0('The busiest actor: ' , author.real.names.vec.demo[which(deg == max(deg))], ', who shows up in ', max(deg), ' scenes!'))
    print('Histogram on k(num of actors per scene):')
    print(table(colSums(authorPaperBiadj))) # num of authors per paper histogram in original dataset

}



# Convert nodeEdgeBiadj to uniform hyperedge set and the corresponding count vec for each k
    # Input:
        # int max.k.for.mcem
    # Output:
        # int max.set.k.demo
        # int max.set.k.for.mcem.demo
        # list mix.hyperedge.demo: a list of matrices of shape n_k by k, for k = 2,3...,max.set.for.mcem.demo
        # list mix.cout.demo: a list of count vectors of length n choose k
convert.adjmat.to.hyperedgeset = function(max.k.for.mcem=3){
    print('###### Proprocessing dateset for mcem ########')
    max.set.k.demo <<- 1+len(table(colSums(authorPaperBiadj)))    #-- upper bound for original starwar hg --#
    max.set.k.for.mcem.demo <<- max.k.for.mcem                    #-- set upper bound for general hg in mcem  --#
    print(sprintf('set up max k = %i in general hg for viz', max.set.k.demo ))
    print(sprintf('set up max k = %i in general hg for mcem', max.set.k.for.mcem.demo ))
    mix.hyperedge.demo <<- collect.kset.from.authorPaperBinaryAdj(max.set.k = max.set.k.for.mcem.demo,
                                                                authorPaperBinaryAdj.mat = authorPaperBiadj)
    print('general hypergraph dataset for mcem: k-sets count:')
    print(sapply(mix.hyperedge.demo, function(x)dim(x)[1]))
    mix.k.demo <<- 2:max.set.k.for.mcem.demo # which hyperedge set to run mcem on
    mix.true.beta.demo <<- c(1e-1, 1e-3, 1e-5) # shape parameter
    mix.count.demo <<- list()
    mix.subset.mat.demo <<- list()
    #load("69-choose-k234.Rdata") # mix.subset.mat is precomputed and save as variable "mix.subset.mat"
    msg('totol numbe of author is ',num.author) # be careful mix.subset.mat is correct file to read in
    for (i in 1:len(mix.k.demo)){
            k = mix.k.demo[i]
            mix.subset.mat.demo[[i]] <<- find_subset(num.author, k)  # Better precomputed. May be slow for large n and k#
            mix.count.demo[[i]] <<- turn.kset.sample.to.count(mix.hyperedge.demo[[i]] , mix.subset.mat.demo[[i]] )
    }
}









# Function `run.mcem`
#     Run mcem for general hypergraph.
#        - The algorihtm is a wrapper for four mcems for uniform hypergraph, indiced by arg `which.mcem`
##          - algorithm type:
##              1: 'rbst-prcr';
##              2: 'rbst-norm';
##              3: 'aggr-prcr';
##              4: 'aggr-norm'.
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
        mix.num.em.iter.demo <<- c(20,20,10,10,10)                                       # --- set iteration nums for AOS dataset k = 2,3,...,6 ---#
        mix.min.num.sweeps.estep.demo <<- c(10,10,10,10,10)                              #--- num of sweeps in E-step  ---#
    }
    {# E-step specs
        p.demo <<- 2
        mu.demo <<- 0
        tausq.prop = 1
        Sigma.prop.demo <<- tausq.prop*diag(p.demo)
        mu.init.demo <<- rep(0,2)
        Sigma.init.demo <<-  num.author*p.demo*diag(p.demo)
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
    lkhd.lst = unlist(mix.likelihood)
    if (isVerboseResult)
    {
        print(paste0('true beta is ', mix.true.beta.demo))
	    print('beta.hat ')
	    print(beta.hat)
	    print('feat.hat')
	    print(dim(feat.hat))
	    print('rej rate along path is :')
	    print(unlist(mix.rej))
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
        save(mix.feat.est, file = paste0("./output/result-",mcem.chosen.demo,"-mcem-feat-b",mix.true.beta.demo[1] ,"-",current.time,".Rdata"))
    }
    msg('*** MCEM for general hypergraph is finished. ')
    return(list(beta.hat=beta.hat, feat.hat=feat.hat))
}





# Viz starwar:  hypergraph on posterior feature
plot.starwar.hg.posterior.feat = function(isSavePlot=FALSE, isShowPlot=TRUE){
    # create a igraph obj
    vertex.color.demo = light.grey
    vertex.label.color.demo = dark.grey
    #vertex.size.demo = 3
    vertex.label.locs.demo <- rep(1, num.author)
    vertex.label.locs.demo[10] = 0.5
    #author.real.names.vec
    #vertex.label.locs.demo[3] = 1 # luke
    #vertex.label.locs.demo[7] = -0.5 # ben
    #vertex.label.locs.demo[9] = -1  # han
    vertex.label.dist.demo = 1.1
    which.kset.demo = 3:4
    # k2 set
    gg = graph(edges = t(mix.hyperedge.demo[[1]]), n = num.author , directed = F)
    #plot(gg)
    if (!is.null(dim(mix.hyperedge.demo[[1]]))){
        E(gg)$type = 1:dim(mix.hyperedge.demo[[1]])[1]
        count.hyperedge = dim(mix.hyperedge.demo[[1]])[1]
    } else{
        E(gg)$type = 1
        count.hyperedge = 1
        E(gg)$lty=1
    }
    for (j in 3: max.set.k.for.mcem.demo){# k = 3,4,...
        #j = 3
        paths.mat = mix.hyperedge.demo[[j-1]]
        if (is.null(dim(paths.mat))){
            head.node = paths.mat[1]
            tail.node = path.mat[len(paths.mat)]
            gg = gg + path(as.vector(paths.mat), weight = 1, type = 1)#, lty = 1)
            gg = add_edges(gg, c(head.node, tail.node),type = 1)#,lty=1)
            count.hyperedge  = count.hyperedge + 1
        } else {
            for (i in 1:dim(paths.mat)[1]){
                #print(i)
                #i = 37
                head.node = paths.mat[i,1]
                tail.node = paths.mat[i,len(paths.mat[i,])]
                gg = gg + path(paths.mat[i,], weight = 1, type = count.hyperedge + i)#, lty = 1)
                gg = add_edges(gg, c(head.node, tail.node), type = count.hyperedge + i)#, lty = 1)
            }
            count.hyperedge = count.hyperedge + dim(paths.mat)[1]
        }
    }
    V(gg)$color = vertex.color.demo
    V(gg)$label = author.real.names.vec.demo
    V(gg)$label.dist = vertex.label.dist.demo
    V(gg)$label.color = vertex.label.color.demo
    V(gg)$size = 0.1*degree(gg) + 3

        if(isSavePlot){
            ggpubr::ggexport(gg, filename='test_starwar_hg_post_feat.png')
            png("test_starwar4_posterior_viz.png", width = 4, height = 4, units = 'in', res = 300)
            par(mar = c(.2, 1, .2, 1)) #bottom, left, top, and right
            plot(gg,
                #layout = layout_with_graphopt(gg), #7 #  force-directed layout algorithm, that scales relatively well to large graphs.
                #layout = layout_with_gem(gg), # 8!!!!  the GEM force-directed layout algorithm.
                #layout = layout_with_dh(gg), # 5 The Davidson-Harel layout algorithm
                #layout = layout_on_sphere, # 6
                #layout = layout_nicely,   # 6
                #layout = layout_as_tree, # 6
                #layout = layout_as_star, 6
                #layout =layout_on_grid, # NO
                layout = res$feat.hat,
                axes = F,
                #vertex.color = vertex.color.demo,
                vertex.frame.color="white",
                #vertex.size = vertex.size.demo,
                vertex.shape = "circle",
                vertex.label.degree=vertex.label.locs.demo,
                vertex.label=author.real.names.vec.demo, vertex.label.dist = vertex.label.dist.demo,
                vertex.label.color = vertex.label.color.demo, vertex.label.cex = .5,
                edge.color = E(gg)$type,  edge.width = .2,
                edge.arrow.size = 0, edge.curved = curve_multiple(gg, start = .2))
                box(lty = 3 , col = light.grey
            )
            dev.off()
        }
    if(isShowPlot){
        plot(gg,
            #layout = layout_with_graphopt(gg), #7 #  force-directed layout algorithm, that scales relatively well to large graphs.
            #layout = layout_with_gem(gg), # 8!!!!  the GEM force-directed layout algorithm.
            #layout = layout_with_dh(gg), # 5 The Davidson-Harel layout algorithm
            #layout = layout_on_sphere, # 6
            #layout = layout_nicely,   # 6
            #layout = layout_as_tree, # 6
            #layout = layout_as_star, 6
            #layout =layout_on_grid, # NO
            layout = res$feat.hat,
            axes = F,
            #vertex.color = vertex.color.demo,
            vertex.frame.color="white",
            #vertex.size = vertex.size.demo,
            vertex.shape = "circle",
            vertex.label.degree=vertex.label.locs.demo,
            vertex.label=author.real.names.vec.demo, vertex.label.dist = vertex.label.dist.demo,
            vertex.label.color = vertex.label.color.demo, vertex.label.cex = .5,
            edge.color = E(gg)$type,  edge.width = .2,
            edge.arrow.size = 0, edge.curved = curve_multiple(gg, start = .2))
        box(lty = 3 , col = light.grey )
    }
}






# Viz starwar: predictive prediction

plot.starwar.predictive.degree = function(which.nodes.to.compare = c(1,3), isVerbose = FALSE, isSavePlot=FALSE, isShowPlot=TRUE){
    mix.deg.dist = list()
    #print(table(colSums(authorPaperBiadj))) # num of authors per paper histogram in original dataset
	max.iter.demo = c(choose(num.author,2), choose(num.author,3))      #--------- num of rep -------------#
	#max.iter.demo = c(10,10,10)
    mix.num.sample.demo = c(11,24)                                                         #--------- size of hg -------------#
	which.kset.to.sim = c(2,3)                                                             #----------- which kset to simulate
    mix.subset.mat.demo = list()
    for (which.k in 1:len(which.kset.to.sim )){
        #which.k = 1
		k = which.kset.to.sim[which.k]
		num.sample.demo = mix.num.sample.demo[which.k]
		p.demo = 2
		mu.demo = 0
		var.demo = 1
		mix.subset.mat.demo[[which.k]] = find_subset(num.author, k)
		shape.vec = c('gamma','pareto','gaussian', 'window')
		which.shape.demo = shape.vec[1]
		save.path = "./output/"
		node.i.deg = rep(0, max.iter.demo[which.k])
		deg.dist.k = c()
		for (iter in 1:max.iter.demo[which.k]){
			sim.data =sample.from.shape.cvxhull(n = num.author, k = k,
												subset.mat = mix.subset.mat.demo[[which.k]],
												feat.mat = res$feat.hat,
												num.sample = num.sample.demo,
												which.shape = which.shape.demo, param.vec = c(1, (res$beta.hat)[which.k]),
												if.save.sample = F, save.path = save.path )
            # start count degree(node) = num of hyperedge(s) containing it
            edge.mat = sim.data$edge.sampled
            deg.hg.k = rep(0, num.author)
            for (which.vert in 1:num.author){
                    count = 0
                    for (which.row in 1:dim(edge.mat)[1]){
                        hyperedge = edge.mat[which.row,]
                        if ( which.vert %in% hyperedge){
                            count = count + 1
                        }
                    }
                    deg.hg.k[which.vert] = count
            }
            deg.dist.k = rbind(deg.dist.k, as.vector(deg.hg.k))
            if (isVerbose){
			    msg(sprintf('%i-set: %i th simulation is done', k, iter))
            }
        }
        mix.deg.dist[[which.k]] = deg.dist.k
        if (isVerbose){
            msg(sprintf('experiemnt: prior predictive degree distribution on simulation is finished.'))
        }
    }
    for (which.node in 1:len(which.nodes.to.compare)){
        node.i = which.nodes.to.compare[which.node]
        node.i.mix.deg = list()
        for (i in 1:len(which.kset.to.sim)){
            deg.hg = mix.deg.dist[[i]]
            node.i.mix.deg[[i]] = deg.hg[, node.i]
        }
        d1 = data.frame(count = node.i.mix.deg[[1]])
        d2 = data.frame(count = node.i.mix.deg[[2]])
        #d3 = data.frame(count = node.i.mix.deg[[3]])
        df <- dplyr::bind_rows(list(k2=d1,k3=d2),.id="group")
        plt = ggplot(df, aes(count, fill = group)) + geom_histogram(aes(y=..density..), alpha=0.9, bins=30, position = 'identity')
        #plt = plt + scale_x_continuous(limits = c(0, 15), breaks=c(0,5,10,15)) + scale_y_continuous(limits = c(0, 0.5), breaks=c(0,0.1,0.2,.3,.4,.5))
        plt = plt + scale_fill_manual(name=" ", values=c(light.green, light.blue),labels=c("k=2","k=3"))
        plt <- plt + xlab("degree") + ylab("density") + theme_linedraw()
        plt <- plt + theme(panel.grid.major = element_line(colour = "grey", linetype = 'dotted',size = 0.2), #gridline
                            panel.grid.minor = element_line(colour = "grey", linetype = 'dotted',size = 0.2), #gridline
                            panel.border = element_rect(colour = "white", fill=NA,linetype = 'dotted'),
                            legend.position = c(0.9, 0.88))
        plot_name = paste0('node_degree_plot_',which.node)
        assign(plot_name, plt)
        print(paste0('node ', author.real.names.vec.demo[node.i], "'s plot is done."))
    }
    # Combine several ggplots into ONE plot.
    figure = ggpubr::ggarrange(node_degree_plot_1  + ggpubr::rremove("ylab") + ggpubr::rremove("xlab"),
                                node_degree_plot_2  + ggpubr::rremove("ylab") + ggpubr::rremove("xlab"),
                        nrow=1,
                        ncol=2,
                        labels=author.real.names.vec.demo[which.nodes.to.compare],
                        font.label = list(size = 12, color = "black", face = "plain", family = NULL, position = "top"),
                        hjust = -2,
                        common.legend = T,
                        legend = "top",
                        align = "hv"
                    )
    figure = figure + theme(plot.margin = margin(1,0.3,1,0.3, "cm"))   #bottom, left, top, and right.
    figure = annotate_figure(figure,
                    left = text_grob("density", rot = 90, hjust = -0.2, vjust = 0.6, size = 14),
                    bottom = text_grob("degree", size = 14, hjust = 0.5, vjust = -2)
    )
    if(isShowPlot){
        print(figure)
    }
    if(isSavePlot){
        ggpubr::ggexport(figure, filename='test_starwar_pred_deg.png')
    }
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
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = num.author,
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
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = num.author,
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
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = num.author,
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
                                                                            mu.init = mu.init.demo, Sigma.init=Sigma.init.demo, tot.author = num.author,
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























































# Obsolete

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
