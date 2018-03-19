graph.phylo <- function(x, y, phy, bks="FD",color="gray", 
                        direction="rightwards", type="phylogram", 
                        legend="sem dados", position1="topleft", 
                        position2="bottomleft", show.tip.label=TRUE, 
                        label.offset=0.3, cex.tip=3, pch.tip=19, 
                        cex.leg=1, pch.nd=20, cex.nd=1.5, cex.leg.nd=1, 
                        bty="n", col.nd="black", show.node.label=FALSE, 
                        edge.color="black", edge.width=1, 
                        edge.lty=1, font=3, no.margin=TRUE)
{
  
  direction <- match.arg(direction, c("rightwards", "leftwards", "upwards", "downwards"))
  color <- match.arg(color, c("gray", "heat", "rainbow"))
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  position1 <- match.arg(position1, c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))
  position2 <- match.arg(position2, c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))
  if(class(phy)=="phylo"){
    if(class(x[,y])=="numeric")
    { 
      x1 <- x[-which(x==0),y]
      x1 <- na.exclude(x1) 
      if(is.character(bks))
      {
        if(bks=="FD") ## FD
        {
          res <- hist(x1, breaks=nclass.FD(x1))$breaks
        }
        if(bks=="Sturges") #Sturges
        {
          res <- hist(x1)$breaks 
        }
        if(bks=="Scott") #Scott
        {
          res <- hist(x1, breaks=nclass.scott(x1))$breaks
        }}
      if(is.numeric(bks)) #vetor com numeros
      { 
        res <- bks
      }
      x <- data.frame(x, "class"=rep(NA, times=dim(x)[1]))
      for(i in 1:(length(res)-1))
      {
        x$class <- replace(x$class, list=(res[i]<=x[,y] & x[,y]<res[i+1]), length(res)-i)
        x$class <- replace(x$class, list=(range(res)[2]<=x[,y]), 1)
      }
      pas <- c(rep(NA, length(res)))
      for(l in 1:length(res))
      {
        pas[l] <- paste(res[l],res[l+1],sep="-")
      }
      pas <- pas[-length(pas)]
      if(color=="gray")
      {
        plot.phylo(phy, show.tip.label=show.tip.label, label.offset=label.offset, direction=direction, type=type,show.node.label=show.node.label, 
                   edge.color=edge.color, edge.width=edge.width, 
                   edge.lty=edge.lty, font=font, no.margin=no.margin)
        for(i in 1:dim(x)[1])
        {
          tiplabels(tip=i, pch=pch.tip, cex=cex.tip, bg=gray(seq(0.1, 0.9, length=length(pas)))[x[i,"class"]], 
                    col=gray(seq(0.1, 0.9, length=length(pas)))[x[i, "class"]])        
        }
        legend(position1,legend=pas[length(pas):1], col=gray(seq(0.1, 0.9, length=length(pas))), pch=pch.tip, cex=cex.leg, bty=bty)
      }
      if(color=="heat")
      {
        plot.phylo(phy, show.tip.label=show.tip.label, label.offset=label.offset, direction=direction, type=type,show.node.label=show.node.label, 
                   edge.color=edge.color, edge.width=edge.width, 
                   edge.lty=edge.lty, font=font, no.margin=no.margin)
        for(i in 1:dim(x)[1])
        {
          tiplabels(tip=i, pch=pch.tip, cex=cex.tip, bg=heat.colors(length(pas))[x[i,"class"]], 
                    col=heat.colors(length(pas))[x[i,"class"]])
        } 
        legend(position1,legend=pas[length(pas):1], col=heat.colors(length(pas)), pch=pch.tip, cex=cex.leg, bty=bty)
      }
      if(color=="rainbow")
      {
        plot.phylo(phy, show.tip.label=show.tip.label, label.offset=label.offset, direction=direction, type=type,show.node.label=show.node.label, 
                   edge.color=edge.color, edge.width=edge.width, 
                   edge.lty=edge.lty, font=font, no.margin=no.margin)
        for(i in 1:dim(x)[1])
        {
          tiplabels(tip=i, pch=pch.tip, cex=cex.tip, bg=rainbow(length(pas))[x[i,"class"]], 
                    col=rainbow(length(pas))[x[i,"class"]])
        }
        legend(position1,legend=pas[length(pas):1], col=rainbow(length(pas)), pch=pch.tip, cex=cex.leg, bty=bty)
      }
      for (m in which(x[,y]==0))
      {
        tiplabels(tip=m, pch=pch.nd, cex=cex.nd, col=col.nd, bg=col.nd)
      } 
      legend(position2,legend=legend, col=col.nd, pch=pch.nd, cex=cex.leg.nd, bty=bty)
    }
    if(class(x[,y])=="factor")
    { x1 <- which(x[,y]!=0)
      if(color=="rainbow"|color=="heat")
      {
        plot.phylo(phy, show.tip.label=show.tip.label, label.offset=label.offset, direction=direction, type=type,show.node.label=show.node.label, 
                   edge.color=edge.color, edge.width=edge.width, 
                   edge.lty=edge.lty, font=font, no.margin=no.margin)
        for(i in x1)
        {
          tiplabels(tip=i, pch=pch.tip, cex=cex.tip, bg=unique(x[i,y]), 
                    col=unique(x[i,y]))
        }
        legend(position1,legend=unique(x[which(x[,y]!=0),y]), col=match(unique(x[which(x[,y]!=0),y]), names(table(x[,y]))), pch=pch.tip, cex=cex.leg, bty=bty)
      }
      if(color=="gray")
      {
        plot.phylo(phy, show.tip.label=show.tip.label, label.offset=label.offset, direction=direction, type=type,show.node.label=show.node.label, 
                   edge.color=edge.color, edge.width=edge.width, 
                   edge.lty=edge.lty, font=font, no.margin=no.margin)
        for(i in x1)
        {
          tiplabels(tip=i, pch=pch.tip, cex=cex.tip, bg=gray(seq(0.1, 0.9, length=length(unique(x[,y]))))[x[i,y]], 
                    col=gray(seq(0.1, 0.9, length=length(unique(x[,y]))))[x[i,y]])
        }
        legend(position1,legend=unique(x[which(x[,y]!=0),y]), col=gray(seq(0.1,0.9,length=length(unique(x[,y]))))[match(unique(x[which(x[,y]!=0),y]), names(table(x[,y])))], pch=pch.tip, cex=cex.leg, bty=bty)
      }
      for (m in which(x[,y]==0))
      {
        tiplabels(tip=m, pch=pch.nd, cex=cex.nd, col=col.nd, bg=col.nd)
        
      }
      
      legend(x=position2,legend=legend, col=col.nd, pch=pch.nd, cex=cex.leg.nd, bty="n")
    }
  }
  else
  {
    cat("objeto phy não é da classe phylo")
  }
}


