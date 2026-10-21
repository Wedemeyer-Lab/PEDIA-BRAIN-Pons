chord_cell_int2 <- function (net, color.use = NULL, group = NULL, cell.order = NULL, 
          sources.use = NULL, targets.use = NULL, lab.cex = 0.8, small.gap = 1, 
          big.gap = 10, annotationTrackHeight = c(0.03), remove.isolate = FALSE, 
          link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
          reduce = -1, transparency = 0.4, link.border = NA, title.name = NULL, 
          show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, metacolumn = "dev.group",
          ...) 
{
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source", "target")
  } else if (is.data.frame(net)) {
    if (all(c("source", "target", "prob") %in% colnames(net)) == 
        FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source, net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)
  net$sample <- NULL
  net[,metacolumn] <-NULL
  net$pathway_name <-NULL
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  net <- subset(net, prob > 0)
  if (dim(net)[1] <= 0) {
    message("No interaction between those cells")
  }
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source, 
                                                             net$target)))
    if (length(cells.removed) > 0) {

      net.fake <- data.frame(cells.removed, cells.removed, 
                             1e-10 * sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if (nrow(net) > nrow(net.fake)) {
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      scale = TRUE
    }
  }
  df <- net
  cells.use <- union(df$source, df$target)
  order.sector <- cell.levels[cell.levels %in% cells.use]
  if (is.null(color.use)) {
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }
  edge.color <- color.use[as.character(df$source)]
  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  }else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df, order = order.sector,
               col = edge.color, 
               grid.col = grid.col, transparency = transparency, link.border = link.border, 
               directional = directional, direction.type = c("diffHeight", 
                                                             "arrows"), link.arr.type = link.arr.type, annotationTrack = "grid", 
               annotationTrackHeight = annotationTrackHeight, preAllocateTracks = list(track.height = max(strwidth(order.sector))), 
               small.gap = small.gap, big.gap = big.gap, link.visible = link.visible, 
               scale = scale, group = group, link.target.prop = link.target.prop, 
               reduce = reduce)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), 
                                  type = "grid", legend_gp = grid::gpar(fill = grid.col), 
                                  title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, 
                                                        "mm"), y = unit(legend.pos.y, "mm"), just = c("right", 
                                                                                                      "bottom"))
  }
  if (!is.null(title.name)) {
    text(-0, 1.02, title.name, cex = 1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}
