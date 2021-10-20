points_func <- function(ss, cols){
  circos.trackPoints(ss$region, 
                     ss$min_ps.t, 
                     ss$ave_abs_score,
                     col = cols,
                     pch = 16, cex = 0.8)
}