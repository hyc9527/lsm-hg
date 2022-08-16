
# Helpers :
#       run experiment to compare shape parameter on  cluster coefficient(cc) for simulated hypergraph



# Find edge set incident at vertex u
find.edge.at.vert = function(u, edge.mat){
   res = c()
   for (j in 1:dim(edge.mat)[1]){
       if (u %in% edge.mat[j,]){
           res = c(res, j)
       }
   }
   return(res)
}


# Calculate deg distribution of hg
find.all.vert.deg.hg = function(n,edge.mat){# degree of distribution
   deg.dist = rep(0, n)
   for (which.vert in 1:n){
           count = 0
           for (which.row in 1:dim(edge.mat)[1]){
               hyperedge = edge.mat[which.row,]
               if ( which.vert %in% hyperedge){
                   count = count + 1
               }
           }
           deg.dist[which.vert] = count
  }
   return(deg.dist)
}


# Find all vertices' neighbourgood
find.all.vert.nbhd.hg = function(n, edge.mat){# degree of distribution
   res = list()
   for (which.vert in 1:n){
           #print(which.vert)
           nbhd.vec = c()
           for (which.row in 1:dim(edge.mat)[1]){
               hyperedge = edge.mat[which.row,]
               if ( which.vert %in% hyperedge){
                   nbhd.vec = rbind(nbhd.vec, hyperedge[hyperedge != which.vert])
               }
           }
           #print(nbhd.vec)
           #print(nbhd.vec)
           if (is.null(nbhd.vec)){
               res[[which.vert]]  = NA
               #print(paste0('node ',which.vert, " is of deg ZERO."))
           } else{
           res[[which.vert]] = sort(unique(nbhd.vec))
           }
           #print(sort(unique(nbhd.vec)))
           #print(len(res))
   }
   return(res)
}


# Find all vertices' incident edge(s)
find.all.vert.edgeset.hg = function(n, edge.mat){# degree of distribution
   res = list()
   for (which.vert in 1:n){
           #print(which.vert)
           node.i.edge = c()
           for (which.row in 1:dim(edge.mat)[1]){
               hyperedge = edge.mat[which.row,]
               if ( which.vert %in% hyperedge){
                   node.i.edge = rbind(node.i.edge, hyperedge)
               }
           }
           if (is.null(node.i.edge)){
               res[[which.vert]]  = NA
               #print(paste0('node ',which.vert, " is of deg ZERO."))
           } else{
           res[[which.vert]] = node.i.edge
           }
           #print(node.i.edge)
   }
   return(res)
}


# Find |E(u) Union E(v)|
find.edge.union.size = function(mat1, mat2){
   #if (is.null(dim(mat1))){mat1 = as.matrix(mat1,1)}
   #if (is.null(dim(mat2))){mat1 = as.matrix(mat2,1)}
   return (dim(unique(rbind(mat1, mat2)))[1])
}


# Find |E(u) intersection E(v)|
find.edge.intersec.size = function(mat1, mat2){
   return (sum(duplicated(rbind(mat1, mat2), fromLast = TRUE)))
}


# Calculate cc for vert v and u
clust.coef.pair = function (umat,vmat){
   num = find.edge.intersec.size(umat, vmat)
   dem = find.edge.union.size(umat, vmat)
   return(num/dem)
}


# Calcualte cc for hg
clust.coef.hg = function(n, hyperedge.mat){
   incident.edge.lst = find.all.vert.edgeset.hg(n, hyperedge.mat)
   cc = 0
   for (i in 1:n){
       i.edgemat = incident.edge.lst[[i]]
       i.nbhd = unique(as.vector(i.edgemat))
       i.nbhd = i.nbhd[i.nbhd!=i]
       if (! is.na(i.nbhd[1])){
           cc.i = 0
           for (which.nb in 1:len(i.nbhd)){
               j = i.nbhd[which.nb]
               #print(j)
               j.edgemat = incident.edge.lst[[j]]
               cc.i = cc.i + clust.coef.pair(i.edgemat,j.edgemat)
           }
           cc = cc + cc.i
       }
   }
   return(cc/n)
}



# Plot
plot.clust.coef = function(isShow.plot=TRUE, isSave.plot=TRUE)
 {# plot clust coeff multiplots
        {
            cc.dt = data.frame()
            for (j in 1:length(true.beta.vec.demo )){
                cc.dt = rbind(cc.dt, data.frame(x =  num.sample.vec.demo , y = cc.hg.gamma[j,], beta = factor(true.beta.vec.demo[j])))
                # print(paste0('beta = ', true.beta.vec.demo[j],' is done.'))
            }
            cc.dt
        }
        {
            beta.col.vec.demo = c(dark.green, dark.blue, light.purple)
            plot <- ggplot(cc.dt)  #+ scale_y_continuous(name = "global cluster coefficient")
            plot <- plot+ geom_line(aes(x= x, y = y, color = beta), size= 1.5)
            plot <- plot + scale_color_manual(values=beta.col.vec.demo )
            plot <- plot + geom_point(aes(x= x, y = y, color = beta, shape = beta), size= 4)
            plot <- plot + xlab("num. hyperedges") + ylab("global cluster coefficient") + theme_linedraw()
            plot <- plot + theme(
                                axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                                panel.grid.major = element_line(colour = "white", linetype = 'dotted',size = 0.02), #gridline
                                panel.grid.minor = element_line(colour = "white", linetype = 'dotted',size = 0.02), #gridline
                                panel.border = element_rect(colour = "white", fill=NA,linetype = 'dotted',size = 1),
                                text=element_text(size=10),
                                legend.position = c(.1, .9),
                                plot.background = element_rect(
                                                                fill = "white",
                                                                colour = "white",
                                                                size = .01))
            cc.plot = plot
            assign("gg4", ggplotGrob(cc.plot))
        }
        ############################################
        for (which.beta in 1:len(true.beta.vec.demo)){# plot three graph for each beta = 0.1, 0.5 ,1
            edge.mat = sample.gamma.lst[[which.beta]]
            { # turn igraph to dataframe for ggplot
                {   # viz on graph
                    feat.all.dt = as.data.frame(feat.all.demo)
                    feat.all.dt$id = 1:n.demo
                    edge.dt = as.data.frame(edge.mat)
                    names(edge.dt) = c('from', 'to')
                    coord = matrix(0,dim(edge.mat)[1],4)
                    for (i in 1:(dim(edge.mat)[1])){
                        coord[i,1] = feat.all.demo[edge.mat[i,1],1]
                        coord[i,2] = feat.all.demo[edge.mat[i,1],2]
                        coord[i,3] = feat.all.demo[edge.mat[i,2],1]
                        coord[i,4] = feat.all.demo[edge.mat[i,2],2]
                    }
                    edge.dt$from_v1 = coord[,1]
                    edge.dt$from_v2 = coord[,2]
                    edge.dt$to_v1 = coord[,3]
                    edge.dt$to_v2 = coord[,4]
                    plt = ggplot()
                    xc=yc=0
                    r = radius.demo+2.5
                    plt =  plt + annotate("path",x=xc+r*cos(seq(0,2*pi,length.out=100)), y=yc+r*sin(seq(0,2*pi,length.out=100)), lty = 2, size = 0.3, color = light.grey)
                    plt = plt + geom_segment(data=edge.dt,aes(x=from_v1, xend = to_v1, y=from_v2,yend = to_v2), colour=beta.col.vec.demo[which.beta],  size= 0.3)
                    plt = plt + geom_point(data=feat.all.dt,aes(x=V1,y=V2),size=.3,colour=dark.grey)
                    plt = plt + geom_point(data=feat.all.dt,aes(x=V1,y=V2),size=.2,colour= 'white')
                    plt = plt + theme_bw()  + labs(x=NULL, y=NULL)
                    plt = plt + theme(
                                axis.text.x=element_blank(), #remove x axis labels
                                axis.ticks.x=element_blank(), #remove x axis ticks
                                axis.text.y=element_blank(),  #remove y axis labels
                                axis.ticks.y=element_blank(),  #remove y axis ticks
                                panel.grid.major = element_line(color = "white"), #gridline
                                panel.grid.minor = element_line(color = "white"), #gridline
                                panel.border = element_blank(),
                                #text=element_text(size=5),
                                #legend.position = c(.3, .85),
                                plot.background = element_rect())
                    plt.name = paste0('gg', which.beta)
                    #print(paste0('gg', which.beta, ' is done.'))
                    assign(plt.name, ggplotGrob(plt))
                }
            }
        }

        { # final plot: three hypergraphs on top + one paths plot at bottom
            plot = ggplot(data.frame(a=1)) + xlim(1, 21) + ylim(1, 32) +
                annotation_custom(gg1, xmin = 1, xmax = 7, ymin = 19, ymax = 29) +      # plot: hg with beta = 0.1
                annotation_custom(gg2, xmin = 7.5, xmax = 13.5, ymin = 19, ymax = 29) +     # plot:hg with beta = 0.5
                annotation_custom(gg3, xmin = 14, xmax = 20, ymin = 19, ymax = 29) +      # plot:hg with beta = 1
                annotation_custom(gg4, xmin = 0, xmax = 20, ymin = 1, ymax = 17)       # plot: coef clustering paths

            #add arrow from hg with beta[0] in green color to somwhere in  plot: coef clustering paths
            # arrow's starting, ending and curvature are determined by argumnet p in bezier(t,p)
            #plot = plot + geom_path(data = as.data.frame(bezier(t = 0:100/100,
            #                                                    p = list(x = c(4, 10, 15, 18.8),
            #                                                            y = c(19.5+1.5, 13, 9, 4.8)))),
            #                        aes(x = V1, y = V2),
            #                        size = 0.1,
            #                        arrow = arrow(length = unit(.01, "npc"), type = "closed"),
            #                        color = light.grey)
            ##add arrow from hg with beta[1] to somwhere in  plot: coef clustering paths
            #plot = plot + geom_path(data = as.data.frame(bezier(t = 0:100/100,
            #                                                    p = list(x = c(11, 15,18.8),
            #                                                            y = c(19.5+1.5, 13, 5.2)))),
            #                        aes(x = V1, y = V2),
            #                        size = 0.1,
            #                        arrow = arrow(length = unit(.01, "npc"), type = "closed"),
            #                        color = light.grey)
            ## add arrow from hg with b[2] to somwhere in  plot: coef clustering paths
            #plot = plot + geom_path(data = as.data.frame(bezier(t = 0:100/100,
            #                                                    p = list(x = c(18, 18.9),
            #                                                            y = c(19.5+1.5, 16.4)))),
            #                        aes(x = V1, y = V2),
            #                        size = 0.1,
            #                        arrow = arrow(length = unit(.01, "npc"), type = "closed"),
            #                        color = light.grey)
            plot = plot + theme(rect = element_blank(),
                                line = element_blank(),
                                text = element_blank(),
                                plot.margin = unit(c(0,0,0,0), "mm"))


            if (isSave.plot)
            {
                ggsave(filename = paste0("test_clust_coef.png"), plot, width = 8, height = 6,  units = "in")
                print("fig is saved in test_clust_coef.png")
            }
            if(isShow.plot)
            {
                plot
            }
        }
}











