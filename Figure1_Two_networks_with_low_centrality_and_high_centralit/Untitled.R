## 6. Network

if (channel == "Americas") {
  stock_main = "AXP.UN"
  J = 20} else 
    if (channel == "Europe") {
      stock_main = "ZURN.SW"
      J = 20
    }  

FRM_history = readRDS(paste0(output_path, "/Lambda/FRM_", channel, ".rds"))
#FRM_individ_fixed = read.csv(paste0(output_path, "/Lambda/Fixed/lambdas_fixed_",
#                      date_start_fixed, "_", date_end_fixed, ".csv"), header = TRUE) %>% as.matrix()
FRM_index = read.csv(paste0(output_path, "/Lambda/FRM_", channel, "_index.csv"), header = TRUE)

date_start_fixed = 20240201
date_end_fixed =20240301
N0_fixed_net = which( gsub("-", "", names(FRM_history)) == date_start_fixed)
N1_fixed_net = which( gsub("-", "", names(FRM_history)) == date_end_fixed)

#par(family = ‘sans’)


fig = image_graph(width = 1000, height = 1000, res = 96, bg = "transparent")

t = N0_fixed_net
for (t in N0_fixed_net:N1_fixed_net) {
  adj0 = read.csv(file=paste0(output_path, "/Adj_Matrices/Adj_Matix_", 
                              gsub("-", "", names(FRM_history))[t], ".csv"), 
                  header = TRUE, sep = "," , row.names = 1)
  
  adj0 = as.matrix(adj0)[1:J, 1:J] 
  #For deleting equity in the name
  for (i in 1:J){
    colnames(adj0)[i]= substr(colnames(adj0)[i],1,(nchar(colnames(adj0)[i])-7))
  }
  adj0 = apply(adj0, 2, as.numeric)
  netw1 = graph_from_adjacency_matrix(adj0, mode = "directed", weighted = T)
  V(netw1)$color = ifelse(V(netw1)$name == stock_main, "orange", "lightgrey")
  colors = rep("Gray", alpha.f = .8, length(E(netw1)))
  colors = ifelse(head_of(netw1, E(netw1))$name == stock_main, 'blue', colors) #inflow
  colors = ifelse(tail_of(netw1, E(netw1))$name == stock_main, 'orange', colors) #outflow
  
  # plot(netw1, layout = layout_in_circle, vertex.label = colnames(adj0), edge.width = 0.8, 
  #      edge.color = colors, edge.arrow.size = 0.9, edge.arrow.width = 1 
  #      #vertex.size = 700*FRM_history[t-N0_fixed_net+1, -1]
  #      )
  png(paste0(output_path, "/Network/Network_", t, "_", 
             channel, ".png"), width = 900, height = 600, bg = "transparent")
  plot(netw1, layout = layout_in_circle, vertex.label = colnames(adj0), edge.width = 0.8,
       edge.color = colors, edge.arrow.size = 0.9, edge.arrow.width = 1
       #vertex.size = 700*FRM_history[t-N0_fixed_net+1, -1]
  )
  # title(FRM_index$date[t], sub = paste0("FRM: ", round(FRM_index$frm[t], 5)), cex.lab = 1.15, 
  #       font.lab = 2, cex.sub = 1.15, font.sub = 2)
  title(xlab = paste0(FRM_index$date[t], "\n FRM: ", round(FRM_index$frm[t], 5)), 
        cex.lab = 1.15, font.lab = 2) #, line = -0.5
  dev.off()
  
}
