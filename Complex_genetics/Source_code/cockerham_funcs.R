# for imperfect correlation (page 88, Cockerham 1963)
crossing = function(Vg1,Vg2,Rg){
  
  out = 2*sqrt(Vg1)*sqrt(Vg2)*(1-Rg)
  
  return(out)
  
}

# for heterogeneous variances (page 88, Cockerham 1963)
hetero = function(Vg1,Vg2){
  
  out = ((sqrt(Vg1)-sqrt(Vg2))^2)
  
  return(out)
  
}

# Calculate imperfect correlation and heterogenenous variance values for each row in data frame
cockerham_input$crossing = sapply(1:nrow(cockerham_input), FUN = function(r){
  Vg1 = cockerham_input$Vg1[r]
  Vg2 = cockerham_input$Vg2[r]
  Rg = cockerham_input$Rg[r]
  
  crossing.val = crossing(Vg1,Vg2,Rg)
  
  return(crossing.val)
})

cockerham_input$hetero = sapply(1:nrow(cockerham_input), FUN = function(r){
  Vg1 = cockerham_input$Vg1[r]
  Vg2 = cockerham_input$Vg2[r]
  
  hetero.val = hetero(Vg1,Vg2)
  
  return(hetero.val)
})
