#' Produces a mean-and-stdErr plot across groups of TFs
#'
#' @param d A data.frame with TFs as rows
#' @param vars The variables to plot
#' @param group.by The variable by which to group TFs
#' @param order.by The variable by which to order the groups
#' @param nbreaks Target number of breaks in the x axis scale
#'
#' @returns A ggplot object
doSEMplot <- function(d, vars, group.by="customClass", order.by=vars[1], nbreaks=3){
  d2 <- reshape::melt.data.frame(d, id.vars = group.by, measure.vars = unique(c(vars, order.by)))
  d2 <- d2[which(!is.na(d2$value)),]
  w <- which(d2$variable==order.by)
  d2$group.by <- droplevels(as.factor(d2[[group.by]]))
  d2$byorder <- sapply(split(d2$value[w], d2$group.by[w]), mean)[as.character(d2$group.by)]
  if(!is.null(names(vars)))
    d2$variable <- factor(as.character(d2$variable), as.character(vars), names(vars))
  gl <- lengths(split(d2$value[w], d2$group.by[w], drop=TRUE))
  if((nnas <- sum(is.na(d2$group.by[w])))>0) gl <- c(gl, "NA"=nnas)
  d2$count <- gl[as.character(d2$group.by)]
  d2$count[is.na(d2$count)] <- sum(is.na(d2$group.by[w]))
  midpoint <- round(exp((log1p(min(gl))+log1p(max(gl)))/2)-1L)
  b <- c(min(d2$count), midpoint, max(d2$count))
  d2 <- d2[which(d2$variable %in% c(vars, names(vars))),]
  ggplot(d2, aes(value, reorder(group.by, byorder), color=count)) + 
    stat_summary(fun.data = mean_se) + 
    ggh4x::facet_grid2(cols = vars(variable), scales="free_x", independent = "x") + theme_sleek() +
    scale_color_viridis_c(trans="log1p", breaks=b) +
    labs(x="Importance (% Shapley)", y=NULL, color="# TFs") + 
    scale_x_continuous(breaks=scales::pretty_breaks(nbreaks, min.n=2))
}

myTheme <- function(){
  theme_bw() + theme(plot.tag = element_text(face = "bold"))
}

signedSqrt <- scales::new_transform("signedSqrt", \(x) sign(x)*sqrt(abs(x)), \(x) sign(x)*x^2)

breakStrings <- function(x, minSizeForBreak = 20, lb = "\n"){
  if (is.factor(x)) {
    levels(x) <- breakStrings(levels(x), minSizeForBreak, lb)
    return(x)
  }
  stopifnot(is.character(x))
  sapply(x, minSizeForBreak = minSizeForBreak, lb = lb,
         FUN = function(x, minSizeForBreak, lb) {
    if (nchar(x) <= minSizeForBreak) 
      return(x)
    g <- gregexpr(" ", x)[[1]]
    if (length(g) == 0) 
      return(x)
    if (length(g) == 1 & all(g == -1)) 
      return(x)
    mid <- nchar(x)/2
    mid <- g[order(abs(g - mid))[1]]
    substr(x, mid, mid) <- lb
    return(x)
  })
}