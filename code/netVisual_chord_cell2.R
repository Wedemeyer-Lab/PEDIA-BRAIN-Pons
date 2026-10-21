function (object, signaling = NULL, meta= NULL, group.term = NULL, group.col = NULL ,net = NULL, slot.name = "netP", 
          color.use = NULL, group = NULL, cell.order = NULL, sources.use = NULL, 
          targets.use = NULL, lab.cex = 0.8, small.gap = 1, big.gap = 10, 
          annotationTrackHeight = c(0.03), remove.isolate = FALSE, 
          link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
          reduce = -1, transparency = 0.4, link.border = NA, title.name = NULL, 
          show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, 
          nCol = NULL, thresh = 0.05, ...) 
{
  if (!is.null(signaling)) {
    pairLR.list<- lapply(object.list, function(object) {
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                         key = "pathway_name", matching.exact = T, pair.only = F)
    pairLR$sample_id <- object@meta$sample_id[1]
    pairLR$age.group <- object@meta$age.group[1]
    pairLR$age.name <- object@meta$age.name[1]
    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, ]
    })
    
    prob.list <- lapply(object.list, function(object){
      if (object@meta$age.name[1] == group.term){
        net <- object@net
        return(net$prob)
      }
    })
    pval.list <- lapply(object.list, function(object){
      if (object@meta$age.name[1] == group.term){
        net <- object@net
        return(net$pval)
      }
    })
    rem <- lapply(prob.list, is.null) %>% unlist()
    names(rem) <- NULL
    prob.list<-prob.list[!rem]
    pval.list<-pval.list[!rem]
    
    for (n in dimnames(prob.list[[i]])[[3]]) {
      n.list <- lapply(prob.list,function(pb){pb[,,n]})
    }
    n.df <- row_bind(n.list)
    prob[pval > thresh] <- 0
    
    
    if (length(pairLR.name) > 1) {
      pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
                                           3, sum) != 0]
    } else {
      pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
                                       0]
    }
    if (length(pairLR.name.use) == 0) {
      stop(paste0("There is no significant communication of ", 
                  signaling))
    }else {
      pairLR <- pairLR[pairLR.name.use, ]
    }
    pairLR.name.use <- pairLR.name.use[-which(is.na(pairLR.name.use))]
    nRow <- length(pairLR.name.use)
    prob <- prob[, , pairLR.name.use]
    if (length(dim(prob)) == 2) {
      prob <- replicate(1, prob, simplify = "array")
    }
    if (slot.name == "netP") {
      message("Plot the aggregated cell-cell communication network at the signaling pathway level")
      net <- apply(prob, c(1, 2), sum)
      if (is.null(title.name)) {
        title.name <- paste0(signaling, " signaling pathway network")
      }
      gg <- netVisual_chord_cell_internal(net, color.use = color.use, 
                                          group = group, cell.order = cell.order, sources.use = sources.use, 
                                          targets.use = targets.use, lab.cex = lab.cex, 
                                          small.gap = small.gap, big.gap = big.gap, annotationTrackHeight = annotationTrackHeight, 
                                          remove.isolate = remove.isolate, link.visible = link.visible, 
                                          scale = scale, directional = directional, link.target.prop = link.target.prop, 
                                          reduce = reduce, transparency = transparency, 
                                          link.border = link.border, title.name = title.name, 
                                          show.legend = show.legend, legend.pos.x = legend.pos.x, 
                                          legend.pos.y = legend.pos.y, ...)
    }
    else if (slot.name == "net") {
      message("Plot the cell-cell communication network per each ligand-receptor pair associated with a given signaling pathway")
      if (is.null(nCol)) {
        nCol <- min(length(pairLR.name.use), 2)
      }
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[, , i]
        gg[[i]] <- netVisual_chord_cell_internal(net, 
                                                 color.use = color.use, group = group, cell.order = cell.order, 
                                                 sources.use = sources.use, targets.use = targets.use, 
                                                 lab.cex = lab.cex, small.gap = small.gap, 
                                                 big.gap = big.gap, annotationTrackHeight = annotationTrackHeight, 
                                                 remove.isolate = remove.isolate, link.visible = link.visible, 
                                                 scale = scale, directional = directional, 
                                                 link.target.prop = link.target.prop, reduce = reduce, 
                                                 transparency = transparency, link.border = link.border, 
                                                 title.name = title.name, show.legend = show.legend, 
                                                 legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, 
                                                 ...)
      }
    }
  }
  else if (!is.null(net)) {
    gg <- netVisual_chord_cell_internal(net, color.use = color.use, 
                                        group = group, cell.order = cell.order, sources.use = sources.use, 
                                        targets.use = targets.use, lab.cex = lab.cex, small.gap = small.gap, 
                                        big.gap = big.gap, annotationTrackHeight = annotationTrackHeight, 
                                        remove.isolate = remove.isolate, link.visible = link.visible, 
                                        scale = scale, directional = directional, link.target.prop = link.target.prop, 
                                        reduce = reduce, transparency = transparency, link.border = link.border, 
                                        title.name = title.name, show.legend = show.legend, 
                                        legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, 
                                        ...)
  }
  else {
    stop("Please assign values to either `signaling` or `net`")
  }
  return(gg)
}
