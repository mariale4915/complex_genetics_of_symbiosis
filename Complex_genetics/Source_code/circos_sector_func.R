sector_func <- function(s.pos,e.pos,s.index,t.index, text){
  pos = circlize(c(s.pos, e.pos),c(0,(CELL_META$cell.ylim[2]*2)),
                 sector.index = s.index, track.index = t.index)
  draw.sector(pos[1, "theta"], pos[2, "theta"], pos[1, "rou"], pos[2, "rou"], 
              clock.wise = TRUE, col = "#00FFFF40", border = NA)
  circos.text((s.pos+e.pos)/2, (CELL_META$cell.ylim[2]*2), text, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.6,
              sector.index = s.index, track.index = t.index)
}