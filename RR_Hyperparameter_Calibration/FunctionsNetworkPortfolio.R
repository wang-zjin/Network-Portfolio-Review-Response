
p_turnover  = function(weights){
  out<-vector()
  for (i in 2:dim(weights)[1]) {
    # print(i)
    # out[i-1]=sum(abs(weights[i,]-weights[i-1,]))/sum((abs(weights[i,])+abs(weights[i-1,]))>0)
    out[i-1]=sum(abs(weights[i,]-weights[i-1,]))/sum((abs(weights[i,])))
  }
  # print(1)
  out = mean(out)
  return(out)
}

pturnoverDN  = function(weights, rets, freq){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  f = abs(bop - eop)
  out = sum(f)*(1/(length(ep) - 1)) #  
  return(out)
}

mad_ew = function(x){
  a = abs(x - 1/ncol(ret))
  return(a)
}
trans_cost = function(weights, rets, freq, c){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  out = c*row_sums(abs(eop - bop))
  return(out)
} 

calculatePerformanceMeasures = function(start,end){
  collNumbers = vector()
  collNumbers_tc = vector()
  collres = xts()
  weightsNumbers = vector()
  for (stratloop in 1:length(strats)){
    Rb = as.xts(rescoll[,  which(strats %in% c("EqualWeight","EW"))], 
                order.by = as.Date(rownames(rescoll)))
    portfolioret_net = na.omit(rescoll[,stratloop])
    strat_weights = weightscoll[[stratloop]] 
    strat_weights[is.nan(strat_weights)] = 0.0
    portfolioret_net_xts = as.xts(as.matrix(na.omit(rescoll[,stratloop])), 
                                  order.by = as.Date(na.omit(rownames(rescoll))))
    portfolioEquity_net = 1 + cumsum(portfolioret_net)
    cumWealth = tail(portfolioEquity_net, 1)
    firstsignal = start
    rettt = portfolioret_net[firstsignal:end]
    rettt_xts = portfolioret_net_xts[firstsignal:end]
    ret_data = rets_log
    stock_rets = ret_data[firstsignal:end]
    tc = trans_cost(strat_weights[firstsignal:end,], stock_rets, freq, c = transi)
    portfolioret_net_tc = portfolioret_net[firstsignal:(end - 1)] - tc
    portfolioEquity_net_tc = 1 + cumsum(portfolioret_net_tc)
    cumWealth_tc = tail(portfolioEquity_net_tc, 1)
    T = (commonDate[end] - commonDate[firstsignal])/365
    Return.ann = (portfolioEquity_net[end]/portfolioEquity_net[firstsignal - 1])^(1/T) - 1
    Return.ann_tc = (tail(portfolioEquity_net_tc, 1)/portfolioEquity_net_tc[firstsignal - 1])^(1/T) - 1
    Vola.ann = sd(rettt)*sqrt(252);
    Vola.ann_tc = sd(portfolioret_net_tc)*sqrt(252);
    Sharpe.ann = Return.ann/Vola.ann
    Sharpe.ann_tc = Return.ann_tc/Vola.ann_tc
    target_turnover = vector();
    for (i in 2:dim(strat_weights)[1]) {
      target_turnover[i] = sum(abs(matrix(strat_weights[i, ]) - matrix(strat_weights[i - 1,])))/dim(strat_weights)[2] 
    }
    Turnover = mean(na.omit(target_turnover))
    value = portfolioEquity_net 
    ma = unlist(lapply(c(2:length(value)),function(x) max(value[1:x])))
    dddisc = value[-1]/ma - 1
    datums = commonDateR[firstsignal:end]
    num = as.numeric(tail(datums,1)-datums[1])
    PR = as.vector(PainRatio(rettt_xts))
    TurnoverDM  = pturnoverDN(strat_weights[firstsignal:end,], stock_rets, freq)
    Return_annual = as.vector(Return.annualized(rettt_xts, geometric = F))
    AverageDrawdown = as.numeric(AverageDrawdown(rettt_xts))
    Sharpe =  as.numeric(SharpeRatio(rettt_xts))[1]
    StdDev.annualized = as.numeric(StdDev.annualized(rettt_xts))
    collNumbers = cbind(collNumbers, as.vector(c(cumWealth, 100*Sharpe.ann,  
                                                 Turnover, TurnoverDM)))#, 
    collNumbers_tc = cbind(collNumbers_tc, as.vector(c(cumWealth_tc,  
                                                       100*Sharpe.ann_tc, Turnover, TurnoverDM)))
    weightcoll = as.data.frame(strat_weights[first_signal:end])
    weightsNumbers = cbind(weightsNumbers, as.vector(c((mean(apply(weightcoll, 1, min))), mean(apply(weightcoll, 1, max)), 
      mean(apply(weightcoll, 1, sd)), mean(apply(weightcoll, 1, mad_ew)), mean(diff(apply(weightcoll, 1, range))))))
    collNumbers = round(collNumbers,4)
    weightsNumbers = round(weightsNumbers,4)
    res = as.xts(portfolioret_net[first_signal:end],order.by=index(ret)[first_signal:end])
    collres = cbind(collres,res)
    first(index(res))
    last(index(res))
    
  }
  return(list(collNumbers,collres,weightsNumbers, collNumbers_tc))
}

#Functions tangency.portfolio, globalMin.portfolio and efficient.frontier are copied 
#from the package "IntroCompFinR: Introduction to Computational Finance in R" @author Eric Zivot
tangency.portfolio =
  function(er,cov.mat,risk.free, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(er)
    if(risk.free < 0)
      stop("Risk-free rate must be positive")
    er = as.vector(er)
    cov.mat = as.matrix(cov.mat)
    N = length(er)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum variance portfolio
    #
    gmin.port = globalMin.portfolio(er, cov.mat, shorts=shorts)
    if(gmin.port$er < risk.free)
      stop("Risk-free rate greater than avg return on global minimum variance portfolio")
    
    # 
    # compute tangency portfolio
    #
    if(shorts==TRUE){
      cov.mat.inv = solve(cov.mat)
      w.t = cov.mat.inv %*% (er - risk.free) # tangency portfolio
      w.t = as.vector(w.t/sum(w.t))          # normalize weights
    } else if(shorts==FALSE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      er.excess = er - risk.free
      Amat = cbind(er.excess, diag(1,N))
      bvec = c(1, rep(0,N))
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w.t = round(result$solution/sum(result$solution), 6)
    } else {
      stop("Shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    names(w.t) = asset.names
    er.t = crossprod(w.t,er)
    sd.t = sqrt(t(w.t) %*% cov.mat %*% w.t)
    tan.port = list("call" = call,
                     "er" = as.vector(er.t),
                     "sd" = as.vector(sd.t),
                     "weights" = w.t)
    class(tan.port) = "portfolio"
    return(tan.port)
  }



efficient.portfolio =
  function(er, cov.mat, target.return, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(er)
    er = as.vector(er) # assign names if none exist
    N = length(er)
    cov.mat = as.matrix(cov.mat)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semidefinite
    
    #
    # compute efficient portfolio
    #
    if(shorts==TRUE){
      ones = rep(1, N)
      top = cbind(2*cov.mat, er, ones)
      bot = cbind(rbind(er, ones), matrix(0,2,2))
      A = rbind(top, bot)
      b.target = as.matrix(c(rep(0, N), target.return, 1))
      x = solve(A, b.target)
      w = x[1:N]
    } else if(shorts==FALSE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), er, diag(1,N))
      bvec = c(1, target.return, rep(0,N))
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
      w = round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    #
    # compute portfolio expected returns and variance
    #
    names(w) = asset.names
    er.port = crossprod(er,w)
    sd.port = sqrt(w %*% cov.mat %*% w)
    ans = list("call" = call,
                "er" = as.vector(er.port),
                "sd" = as.vector(sd.port),
                "weights" = w) 
    class(ans) = "portfolio"
    return(ans)
  }

efficient.frontier =
  function(er, cov.mat, nport=20, alpha.min=-0.5, alpha.max=1.5, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(er)
    er = as.vector(er)
    N = length(er)
    cov.mat = as.matrix(cov.mat)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    
    #
    # create portfolio names
    #
    port.names = rep("port",nport)
    ns = seq(1,nport)
    port.names = paste(port.names,ns)
    
    #
    # compute global minimum variance portfolio
    #
    cov.mat.inv = solve(cov.mat)
    one.vec = rep(1, N)
    port.gmin = globalMin.portfolio(er, cov.mat, shorts)
    w.gmin = port.gmin$weights
    
    if(shorts==TRUE){
      # compute efficient frontier as convex combinations of two efficient portfolios
      # 1st efficient port: global min var portfolio
      # 2nd efficient port: min var port with ER = max of ER for all assets
      er.max = max(er)
      port.max = efficient.portfolio(er,cov.mat,er.max)
      w.max = port.max$weights    
      a = seq(from=alpha.min,to=alpha.max,length=nport) # convex combinations
      we.mat = a %o% w.gmin + (1-a) %o% w.max	         # rows are efficient portfolios
      er.e = we.mat %*% er							                 # expected returns of efficient portfolios
      er.e = as.vector(er.e)
    } else if(shorts==FALSE){
      we.mat = matrix(0, nrow=nport, ncol=N)
      we.mat[1,] = w.gmin
      we.mat[nport, which.max(er)] = 1
      er.e = as.vector(seq(from=port.gmin$er, to=max(er), length=nport))
      for(i in 2:(nport-1)) 
        we.mat[i,] = efficient.portfolio(er, cov.mat, er.e[i], shorts)$weights
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    names(er.e) = port.names
    cov.e = we.mat %*% cov.mat %*% t(we.mat) # cov mat of efficient portfolios
    sd.e = sqrt(diag(cov.e))					        # std devs of efficient portfolios
    sd.e = as.vector(sd.e)
    names(sd.e) = port.names
    dimnames(we.mat) = list(port.names,asset.names)
    
    # 
    # summarize results
    #
    ans = list("call" = call,
                "er" = er.e,
                "sd" = sd.e,
                "weights" = we.mat)
    class(ans) = "Markowitz"
    ans
  }

# function to calculate recurrence vector
RP.P = 
  function(x,epsilon)
  {
    Recurrence_Matrix = matrix(0, length(x), length(x))
    for (i in 1:length(x)) {
      for (j in 1:length(x)) {
        if(epsilon-abs(x[i]-x[j])>=0) Recurrence_Matrix[i,j]=1
        else Recurrence_Matrix[i,j]=0
      }
    }
    P=matrix(0,length(x)-1,1)
    for (tau in 1:(length(x)-1)) {
      for (j in 1:(length(x)-tau)) {
        P[tau,1]=P[tau,1]+Recurrence_Matrix[j,j+tau]
        }
      P[tau,1]=P[tau,1]/(length(x)-tau)
      }
    P=P-mean(P)
    return(P)
  }

# function to calculate recurrence matrix
RP.Mtx = 
  function(x,epsilon)
  {
    Recurrence_Matrix = matrix(0, length(x), length(x))
    for (i in 1:length(x)) {
      for (j in 1:length(x)) {
        if(epsilon-abs(x[i]-x[j])>=0) Recurrence_Matrix[i,j]=1
        else Recurrence_Matrix[i,j]=0
      }
    }
    return(Recurrence_Matrix)
  }

# Construct network of portfolio based on correlation matrix
network.portfolio = 
  function(returnstd)
  {
    # correlation matrix
    Cormat<-cor(returnstd)                        # correlation matrix
    colnames(Cormat)<-colnames(returnstd)
    rownames(Cormat)<-colnames(returnstd)
    # distance matrix
    # Dist_mat<-sqrt(2-2*Covmat)                    # distance matrix
    # Dist_mat<-sqrt(1-Covmat)                    # distance matrix
    Dist_mat<-1-Cormat                    # distance matrix
    # Dist_mat<-sqrt(2*Covmat+2)                    # distance matrix
    Dist_mat<-as.matrix(Dist_mat)
    Dist_mat[is.nan(Dist_mat)]<-0 
    colnames(Dist_mat)<-colnames(returnstd)
    rownames(Dist_mat)<-colnames(returnstd)
    # construct network
    network=graph_from_adjacency_matrix(Dist_mat,weighted=T,
                                        mode="undirected", diag=F)  
    Edgelist_network<-get.edgelist(network)                           # edges of network
    weight_network<-E(network)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    net_port<- graph_from_data_frame(d=links2_network, directed=F)  # net of whole data
    return(net_port)
  }

# Construct network of portfolio based on correlation matrix
# the difference beween network.portfolio is that this directly use correlation matrix
network.correlation = 
  function(returnstd)
  {
    # correlation matrix
    Cormat<-cor(returnstd)                        # correlation matrix
    colnames(Cormat)<-colnames(returnstd)
    rownames(Cormat)<-colnames(returnstd)
    
    # distance matrix
    # Dist_mat<-sqrt(2-2*Cormat)                    # distance matrix
    # Dist_mat<-sqrt(1-Cormat)                    # distance matrix
    # Dist_mat<-1-Cormat                    # distance matrix
    # Dist_mat<--Cormat                     # distance matrix
    Dist_mat<-Cormat-diag(1,dim(Cormat)[1])                    # distance matrix
    # Dist_mat<-1+Cormat                    # distance matrix
    # Dist_mat<-sqrt(2*Covmat+2)                    # distance matrix
    Dist_mat<-as.matrix(Dist_mat)
    Dist_mat[is.nan(Dist_mat)]<-0
    colnames(Dist_mat)<-colnames(returnstd)
    rownames(Dist_mat)<-colnames(returnstd)
    # construct network
    network=graph_from_adjacency_matrix(Dist_mat,weighted=T,
                                        mode="undirected", diag=F)
    Edgelist_network<-get.edgelist(network)                           # edges of network
    weight_network<-E(network)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    net_port<- graph_from_data_frame(d=links2_network, directed=F)  # net of whole data
    return(net_port)
    
    # network=graph_from_adjacency_matrix(Cormat,weighted=T,
    #                                     mode="undirected", diag=F)  
    # Edgelist_network<-get.edgelist(network)                           # edges of network
    # return(network)
  }

# Construct network of portfolio based on correlation matrix
network.abs = 
  function(returnstd)
  {
    # correlation matrix
    CorMat<-cor(returnstd)                        # correlation matrix
    colnames(CorMat)<-colnames(returnstd)
    rownames(CorMat)<-colnames(returnstd)
    # construct network
    network=graph_from_adjacency_matrix(abs(CorMat)-diag(1,ncol=ncol(CorMat),nrow = nrow(CorMat)),
                                        weighted=T,
                                        mode="undirected", diag=F)  
    Edgelist_network<-get.edgelist(network)                           # edges of network
    weight_network<-E(network)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    net_port<- graph_from_data_frame(d=links2_network, directed=F)  # net of whole data
    return(net_port)
  }
# 
# network.decomposition = 
#   function(returnstd)
#   {
#     # correlation matrix
#     CorMat<-cor(returnstd)                        # correlation matrix
#     colnames(CorMat)<-colnames(returnstd)
#     rownames(CorMat)<-colnames(returnstd)
#     # decomposition into positive and negative matrix
#     Lambda<-CorMat-diag(1,ncol=ncol(CorMat),nrow = nrow(CorMat))
#     Lambda1<-matrix(0,ncol=ncol(CorMat),nrow = nrow(CorMat))
#     Lambda2<-matrix(0,ncol=ncol(CorMat),nrow = nrow(CorMat))
#     for (i in seq(1,nrow(CorMat))) {
#       for (j in seq(1,nrow(CorMat))) {
#         if (CorMat[i,j]>0){Lambda1[i,j]<-CorMat[i,j]}
#         else if (CorMat[i,j]<0){Lambda2[i,j]<-CorMat[i,j]}
#       }
#     }
#     # construct network
#     network.positive=graph_from_adjacency_matrix(Lambda1,
#                                         weighted=T,
#                                         mode="undirected", diag=F)  
#     network.negative=graph_from_adjacency_matrix(Lambda2,
#                                                  weighted=T,
#                                                  mode="undirected", diag=F)  
#     # Edgelist_network<-get.edgelist(network)                           # edges of network
#     # weight_network<-E(network)$weight                                 # weight of network
#     # A<-cbind(Edgelist_network,weight_network)
#     # A<-as.matrix(A)
#     # links2_network<-as.data.frame(A)                                # links of network
#     # colnames(links2_network)<-c("from","to","weight")
#     # net_port<- graph_from_data_frame(d=links2_network, directed=F)  # net of whole data
#     return(network.positive)
#     return(network.negative)
#   }

# Construct network of portfolio based on correlation matrix 
# decompose correlation matrix into an indicator matrix and a matrix Lambda
# Lambda decomposed into a positive matrix and negative matrix
network.posi.naga = 
  function(returnstd)
  {
    # correlation matrix
    CorMat<-cor(returnstd)                        # correlation matrix
    colnames(CorMat)<-colnames(returnstd)
    rownames(CorMat)<-colnames(returnstd)
    # decompose coefficient matrix
    Lambda<- CorMat-diag(1,nrow=nrow(CorMat),ncol = ncol(CorMat))
    Lambda1<- Lambda
    Lambda1[Lambda1<0]<- 0
    Lambda2<- -Lambda
    Lambda2[Lambda2<0]<- 0
    # construct positive coefficient network
    network_positive=graph_from_adjacency_matrix(Lambda1,weighted=T,
                                        mode="undirected", diag=F)  
    Edgelist_network<-get.edgelist(network_positive)                           # edges of network
    weight_network<-E(network_positive)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    network_positive<- graph_from_data_frame(d=links2_network, directed=F)  # net of positive coefficient matrix
    # construct negative coefficient network
    network_negative=graph_from_adjacency_matrix(Lambda2,weighted=T,
                                                 mode="undirected", diag=F)  
    Edgelist_network<-get.edgelist(network_positive)                           # edges of network
    weight_network<-E(network_positive)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    network_negative<- graph_from_data_frame(d=links2_network, directed=F)  # net of negative coefficient matrix
    # output
    network_coefficient<-list("positive"=network_positive,"negative"=network_negative)
    return(network_coefficient)
  }

# Construct Minimum spanning tree of portfolio based on covariance
MST.portfolio = 
  function(returnstd)
  {
    # correlation matrix
    Covmat<-cor(returnstd)                        # correlation matrix
    colnames(Covmat)<-colnames(returnstd)
    rownames(Covmat)<-colnames(returnstd)
    # distance matrix
    Dist_mat<-sqrt(2-2*Covmat)                    # distance matrix
    Dist_mat<-as.matrix(Dist_mat)
    Dist_mat[is.nan(Dist_mat)]<-0 
    colnames(Dist_mat)<-colnames(returnstd)
    rownames(Dist_mat)<-colnames(returnstd)
    # construct network
    network=graph_from_adjacency_matrix(Dist_mat,weighted=T,
                                              mode="undirected", diag=F)  
    Edgelist_network<-get.edgelist(network)                           # edges of network
    weight_network<-E(network)$weight                                 # weight of network
    A<-cbind(Edgelist_network,weight_network)
    A<-as.matrix(A)
    links2_network<-as.data.frame(A)                                # links of network
    colnames(links2_network)<-c("from","to","weight")
    net_port<- graph_from_data_frame(d=links2_network, directed=F)  # net of whole data
    mst_port<- minimum.spanning.tree(net_port)                      # minimum spanning tree
    return(mst_port)
  }

network.globalMin.portfolio =
  function(nc, cov.mat, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(nc)
    nc = as.vector(nc) # assign names if none exist
    cov.mat = as.matrix(cov.mat)
    N = length(nc)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum portfolio
    #
    if(shorts==TRUE){
      cov.mat.inv = solve(cov.mat)
      one.vec = rep(1,N)
      w.gmin = rowSums(cov.mat.inv) / sum(cov.mat.inv)
      w.gmin = as.vector(w.gmin)
    } else if(shorts==FALSE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), diag(1,N))
      bvec = c(1, rep(0,N))
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w.gmin = round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    # if(shorts==TRUE){
    #   cov.mat.inv = solve(cov.mat)
    #   one.vec = rep(1,N)
    #   w.gmin = rowSums(cov.mat.inv) / sum(cov.mat.inv)
    #   w.gmin = as.vector(w.gmin)
    # } else if(shorts==FALSE){
    #   Dmat = 2*cov.mat
    #   dvec = rep.int(0, N)
    #   Amat = cbind(rep(1,N), diag(1,N))
    #   bvec = c(1, rep(0,N))
    #   result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
    #   w.gmin = round(result$solution, 6)
    # } else {
    #   stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    # }
    
    names(w.gmin) = asset.names
    nc.gmin = crossprod(w.gmin,nc)
    sd.gmin = sqrt(t(w.gmin) %*% cov.mat %*% w.gmin)
    gmin.port = list("call" = call,
                     "nc" = as.vector(nc.gmin),
                     "sd" = as.vector(sd.gmin),
                     "weights" = w.gmin)
    class(gmin.port) = "portfolio"
    gmin.port
  }

network.efficient.portfolio =
  function(nc, cov.mat, target.nc, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(nc)
    nc = as.vector(nc) # assign names if none exist
    N = length(nc)
    cov.mat = as.matrix(cov.mat)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semidefinite
    
    #
    # compute efficient portfolio
    #
    if(shorts==TRUE){
      # ones = rep(1, N)
      # top = cbind(2*cov.mat, nc, ones)
      # bot = cbind(rbind(nc, ones), matrix(0,2,2))
      # A = rbind(top, bot)
      # b.target = as.matrix(c(rep(0, N), target.nc, 1))
      # x = solve(A, b.target)
      # w = x[1:N]
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), -nc)
      bvec = cbind(1, -target.nc)
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w = round(result$solution, 6)
    } else if(shorts==FALSE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), -nc, diag(1,N))
      bvec = c(1, -target.nc, rep(0,N))
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w = round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    #
    # compute portfolio expected returns and variance
    #
    names(w) = asset.names
    nc.port = crossprod(nc,w)
    sd.port = sqrt(w %*% cov.mat %*% w)
    ans = list("call" = call,
               "nc" = as.vector(nc.port),
               "sd" = as.vector(sd.port),
               "weights" = w) 
    class(ans) = "portfolio"
    return(ans)
  }

network.3constraint.portfolio =
  function(nc,er, cov.mat, target.nc,target.er, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(nc)
    nc = as.vector(nc) # assign names if none exist
    er = as.vector(er) # assign names if none exist
    N = length(nc)
    cov.mat = as.matrix(cov.mat)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semidefinite
    
    #
    # compute efficient portfolio
    #
    if(shorts==TRUE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), -nc, er)
      bvec = cbind(1, -target.nc, target.er)
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w = round(result$solution, 6)
    } else if(shorts==FALSE){
      Dmat = 2*cov.mat
      dvec = rep.int(0, N)
      Amat = cbind(rep(1,N), -nc, diag(1,N), er)
      bvec = c(1, -target.nc, rep(0,N), target.er)
      result = quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
      w = round(result$solution, 6)
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    #
    # compute portfolio expected returns and variance
    #
    names(w) = asset.names
    nc.port = crossprod(nc,w)
    er.port = crossprod(er,w)
    sd.port = sqrt(w %*% cov.mat %*% w)
    ans = list("call" = call,
               "nc" = as.vector(nc.port),
               "er" = as.vector(er.port),
               "sd" = as.vector(sd.port),
               "weights" = w) 
    class(ans) = "portfolio"
    return(ans)
  }

network.efficient.frontier =
  function(nc, cov.mat, nport=20, alpha.min=-0.5, alpha.max=1.5, shorts=TRUE)
  {
    call = match.call()
    
    #
    # check for valid inputs
    #
    asset.names = names(nc)
    er = as.vector(nc)
    N = length(nc)
    cov.mat = as.matrix(cov.mat)
    if(N != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    
    #
    # create portfolio names
    #
    port.names = rep("port",nport)
    ns = seq(1,nport)
    port.names = paste(port.names,ns)
    
    #
    # compute global minimum variance portfolio
    #
    cov.mat.inv = solve(cov.mat)
    one.vec = rep(1, N)
    port.gmin = network.globalMin.portfolio(nc, cov.mat, shorts)
    w.gmin = port.gmin$weights
    
    if(shorts==TRUE){
      # compute efficient frontier as convex combinations of two efficient portfolios
      # 1st efficient port: global min var portfolio
      # 2nd efficient port: min var port with ER = max of ER for all assets
      nc.min = min(nc)
      port.min = network.efficient.portfolio(nc,cov.mat,nc.min, shorts)
      w.min = port.min$weights    
      a = seq(from=alpha.min,to=alpha.max,length=nport) # convex combinations
      we.mat = a %o% w.gmin + (1-a) %o% w.min	         # rows are efficient portfolios
      nc.e = we.mat %*% nc							                 # expected returns of efficient portfolios
      nc.e = as.vector(nc.e)
    } else if(shorts==FALSE){
      we.mat = matrix(0, nrow=nport, ncol=N)
      we.mat[1,] = w.gmin
      we.mat[nport, which.min(nc)] = 1
      nc.e = as.vector(seq(from=port.gmin$nc, to=min(nc), length=nport))
      for(i in 2:(nport-1)) 
        we.mat[i,] = efficient.portfolio(nc, cov.mat, nc.e[i], shorts)$weights
    } else {
      stop("shorts needs to be logical. For no-shorts, shorts=FALSE.")
    }
    
    names(nc.e) = port.names
    cov.e = we.mat %*% cov.mat %*% t(we.mat) # cov mat of efficient portfolios
    sd.e = sqrt(diag(cov.e))					        # std devs of efficient portfolios
    sd.e = as.vector(sd.e)
    names(sd.e) = port.names
    dimnames(we.mat) = list(port.names,asset.names)
    
    # 
    # summarize results
    #
    ans = list("call" = call,
               "nc" = nc.e,
               "sd" = sd.e,
               "weights" = we.mat)
    class(ans) = "Markowitz"
    ans
  }

##### portfolio using Dantzig type selector with 1 constraints #####
linfun1=function(Sn,b,lambda)
{
  #equivalent to solving   min 1'u
  #such that		u-x>=0
  #						u+x>=0
  #						-hatS x>=-(lambda 1+b)
  #						hatS x>=-lambda 1+b
  #						x^T1=1
  p = length(b)
  a=rep(0,2*p)
  a[c(1,1+p)]=1
  A0=toeplitz(a)
  A0[upper.tri(A0)]=-A0[upper.tri(A0)]
  A1=cbind(matrix(0,p,p),-Sn)
  A2=-A1
  A=rbind(A0,A1,A2)
  rhs=c(rep(0,2*p),c(-b,b)-lambda*rep(1,2*p))
  C=rep(c(1,0),c(p,p))
  #	EE=rep(c(0,1),c(p,p))
  #	FF=1
  solvetheta=linp(#E=EE,F=FF,
    G=A,H=rhs,Cost=C,ispos=FALSE)$X[(p+1):(2*p)]
  # return(solvetheta)
  return(solvetheta)
}
# What's new in "linfun1_1" is that we added a time controler
linfun1_1=function(Sn,b,lambda)
{
  #equivalent to solving   min 1'u
  #such that		u-x>=0
  #						u+x>=0
  #						-hatS x>=-(lambda 1+b)
  #						hatS x>=-lambda 1+b
  #						x^T1=1
  p = length(b)
  a=rep(0,2*p)
  a[c(1,1+p)]=1
  A0=toeplitz(a)
  A0[upper.tri(A0)]=-A0[upper.tri(A0)]
  A1=cbind(matrix(0,p,p),-Sn)
  A2=-A1
  A=rbind(A0,A1,A2)
  rhs=c(rep(0,2*p),c(-b,b)-lambda*rep(1,2*p))
  C=rep(c(1,0),c(p,p))
  
  # Setting time limit for linp
  setTimeLimit(elapsed = 100)  # e.g., 10 seconds max for linp
  
  solvetheta <- tryCatch({
    linp(G = A, H = rhs, Cost = C, ispos = FALSE)$X[(p+1):(2*p)]
  }, error = function(e) {
    message("linp exceeded time limit or encountered an error: ", e$message)
    return(rep(0, length(b)))  # Fallback to a default value if linp fails
  })
  
  # Reset time limit
  setTimeLimit(cpu = Inf, elapsed = Inf)
  
  return(solvetheta)
}
##### portfolio using Dantzig type selector with 2 constraints #####
linfun2=function(Sn,b1,b2,lambda)
{
  #equivalent to solving   min 1'u1 + 1'u1
  #such that		u1-x1>=0
  #						  u1+x1>=0
  #						  u2-x2>=0
  #						  u2+x2>=0
  #				 	 -hatS x1>=-(lambda 1+b1)
  #				    hatS x1>=-lambda 1+b1
  #				 	 -hatS x2>=-(lambda 1+b2)
  #				    hatS x2>=-lambda 1+b2
  #						x^T1=1 (I do not think so)
  a=rep(0,4*p)
  a[c(1,1+2*p)]=1
  A0=toeplitz(a)
  A0[upper.tri(A0)]=-A0[upper.tri(A0)]
  A1=cbind(matrix(0,p,2*p),-Sn,matrix(0,p,p))
  A2=-A1
  A3=cbind(matrix(0,p,3*p),-Sn)
  A4=-A1
  A=rbind(A0,A1,A2,A3,A4)
  rhs=c(rep(0,4*p),c(-b1,b1,-b2,b2)-lambda*rep(1,4*p))
  C=rep(c(1,0),c(2*p,2*p))
  #	EE=rep(c(0,1),c(p,p))
  #	FF=1
  theta=linp(G=A,H=rhs,Cost=C,ispos=FALSE)
  solvetheta=list("theta1"=theta$X[(2*p+1):(3*p)],
                  "theta2"=theta$X[(3*p+1):(4*p)])
  # return(solvetheta)
  return(solvetheta)
}


linfun3=function(Sn,b,lambda,a1=1)
{
  #equivalent to solving   min 1'x
  #such that		x>=0
  #						-hatS x>=-(lambda 1+b)
  #						hatS x>=-lambda 1+b
  #						x(1)=a1
  p <- dim(Sn)[1]
  a=rep(0,p)
  a[1]=1
  A0=diag(1,p,p)
  A1=cbind(-Sn)
  A2=-A1
  A=rbind(A0,A1,A2)
  rhs=c(rep(0,p),c(-b,b)-lambda*rep(1,p))
  C=rep(1,p)
  EE=diag(a)
  a[1]=a1
  FF=a
  #	EE=rep(c(0,1),c(p,p))
  #	FF=1
  solvetheta=linp(E=EE,F=FF,
    G=A,H=rhs,Cost=C,ispos=FALSE)$X
  # return(solvetheta)
  return(solvetheta)
}
linfun3_1=function(Sn,b,lambda,a1=1)
{
  #equivalent to solving   min 1'x
  #such that		x>=0
  #						-hatS x>=-(lambda 1+b)
  #						hatS x>=-lambda 1+b
  #						x(1)=a1
  p <- dim(Sn)[1]
  a=rep(0,p)
  a[1]=1
  A0=diag(1,p,p)
  A1=cbind(-Sn)
  A2=-A1
  A=rbind(A0,A1,A2)
  rhs=c(rep(0,p),c(-b,b)-lambda*rep(1,p))
  C=rep(1,p)
  EE=diag(a)
  a[1]=a1
  FF=a
  #	EE=rep(c(0,1),c(p,p))
  #	FF=1
  solvetheta=linp(E=EE,F=FF,
                  G=A,H=rhs,Cost=C,ispos=FALSE)$X
  
  # Make sure max(solvetheta) is not all 0
  while (max(solvetheta)==0) {
    lambda = lambda + 0.1
    a=rep(0,p)
    a[1]=1
    A0=diag(1,p,p)
    A1=cbind(-Sn)
    A2=-A1
    A=rbind(A0,A1,A2)
    rhs=c(rep(0,p),c(-b,b)-lambda*rep(1,p))
    C=rep(1,p)
    EE=diag(a)
    a[1]=a1
    FF=a
    #	EE=rep(c(0,1),c(p,p))
    #	FF=1
    solvetheta=linp(E=EE,F=FF,
                    G=A,H=rhs,Cost=C,ispos=FALSE)$X
  }
  return(solvetheta)
}
linfun4=function(Sn,b,lambda,a1=1,indx=1)
{
  #equivalent to solving   min 1'x
  #such that		x>=0
  #						-hatS x>=-(lambda 1+b)
  #						hatS x>=-lambda 1+b
  #						x(1)=a1
  
  A0=diag(1,p,p)
  A1=cbind(-Sn)
  A2=-A1
  A=rbind(A0,A1,A2)
  rhs=c(rep(0,p),c(-b,b)-lambda*rep(1,p))
  C=rep(1,p)
  a=rep(0,p)
  a[indx]=1
  EE=diag(a)
  a[indx]=a1
  FF=a
  #	EE=rep(c(0,1),c(p,p))
  #	FF=1
  solvetheta=linp(E=EE,F=FF,
                  G=A,H=rhs,Cost=C,ispos=FALSE)$X
  # return(solvetheta)
  return(solvetheta)
}


Dantzig_eigencent_rollwind_hypertune = function(returnstd, window_size = 500, lambda_max = 1){
  
  tic("Tuning parameter of Dantzig estimation in Eigenvector centrality")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):(window_size+22+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  
  
  
  lmd.EC.Dantzig.list = list()
  for (t in c(1:length(W_in))) {
    n<-dim(W_in[[t]])[1]
    B=100
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    lmd.i=c()
    for (valid.block in 1:n.block) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[t]][train.ind,]
      returnstd.valid=W_in[[t]][valid.ind,]
      mu.train=rep(0,p)
      mu.valid=rep(0,p)
      
      corr.train=cor(returnstd.train)
      corr.train[is.na(corr.train)]=0
      A.train=corr.train-diag(1,p,p)
      corr.valid=cor(returnstd.valid)
      corr.valid[is.na(corr.valid)]=0
      A.valid=corr.valid-diag(1,p,p)
      
      cov.train=A.train-diag(max(eigen(A.train)$value),p,p)
      cov.valid=A.valid-diag(max(eigen(A.valid)$value),p,p)
      lambda.grid=seq(0,lambda_max,length=101)[2:101]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l.lmd=NULL
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        print(i)
        # lmd=0.2
        lin.train=linfun3(cov.train,mu.train,lmd,abs(eigen(cov.valid)$vector[1,1]))
        # max(lin.train)
        # max(cov.train%*%lin.train)
        if(!(all(lin.train==0))){
          error=sum((cov.valid%*%lin.train-mu.valid)^2)
          cv.l.error=c(cv.l.error,error)
          cv.l.lmd=c(cv.l.lmd, lmd)
        }
      }
      lmd.i[valid.block]=min(cv.l.lmd[which(cv.l.error==min(cv.l.error))])
    }
    lmd=mean(lmd.i[lmd.i<Inf])
    # lmd.EC.Dantzig <- 0.682
    lmd.EC.Dantzig <- lmd
    lmd.EC.Dantzig.list[[t]] = lmd.EC.Dantzig
  }
  toc()
  
  save(lmd.EC.Dantzig.list,file=paste0("Dantzig_lambda_rolling_windowsize",window_size,".RData"))
  # load("Dantzig_lambda_rolling_window_125.RData")
  # mean(unlist(lmd.EC.Dantzig.list)) # 0.5085714
  # std(unlist(lmd.EC.Dantzig.list))
  load(paste0("Dantzig_lambda_rolling_windowsize",window_size,".RData"))
  lmd.EC.Dantzig = lmd.EC.Dantzig.list[[1]] 
  
  #   EC_DS <- eigenvector centrality estimated by Dantzig selector
  EC_DS<-list()
  a=c()
  for (t in 1: length(W_in)) {
    print(t)
    EC_DS[[t]] =linfun3(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p),rep(0,p),lambda=lmd.EC.Dantzig,abs(eigen(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p))$vector[1,1]))
    EC_DS[[t]]=EC_DS[[t]]/max(EC_DS[[t]])
    # boxplot(EC_DS[[t]])
    a[t]=sum(EC_DS[[t]]==0)
  }
  a
  eigenvector_centrality = list("eigenvector_absolute_value"=EC_in,"eigenvector_centrality_Dantzig"=EC_DS,
                                "zeors_in_EC_DS"=a,"covariance_matrix"=COV_in,"correlation_matrix"=C_in,
                                "expected_return"=ER_in)
  save(eigenvector_centrality,file = paste0("eigenvector_centrality_WS",window_size,"_20240723.RData"))
}

#### Dantzig selector estimation for eigenvector centrality, rolling window calibrate hyperparameter ####
estimate_rollingwindow_dantzig_lambda_eigenvectorcentrality = function(file_name = "SP500 securities_up_20230306.csv", 
                                                                       window_size=500){
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter of Dantzig Selector for Eigenvector Centrality, using rolling window
  
  tic("Tuning parameter of Dantzig estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  lmd.EC.Dantzig.list = list()
  for (t in c(1:T.windows)) {
    n<-dim(W_in[[t]])[1]
    B=100
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    lmd.i=c()
    for (valid.block in 1:n.block) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[t]][train.ind,]
      returnstd.valid=W_in[[t]][valid.ind,]
      mu.train=rep(0,p)
      mu.valid=rep(0,p)
      
      corr.train=cor(returnstd.train)
      corr.train[is.na(corr.train)]=0
      A.train=corr.train-diag(1,p,p)
      corr.valid=cor(returnstd.valid)
      corr.valid[is.na(corr.valid)]=0
      A.valid=corr.valid-diag(1,p,p)
      
      cov.train=A.train-diag(max(eigen(A.train)$value),p,p)
      cov.valid=A.valid-diag(max(eigen(A.valid)$value),p,p)
      lambda.grid=seq(0,1,length=101)[2:101]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l.lmd=NULL
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        print(i)
        # lmd=0.2
        lin.train=linfun3(cov.train,mu.train,lmd,abs(eigen(cov.valid)$vector[1,1]))
        # max(lin.train)
        # max(cov.train%*%lin.train)
        if(!(all(lin.train==0))){
          error=sum((cov.valid%*%lin.train-mu.valid)^2)
          cv.l.error=c(cv.l.error,error)
          cv.l.lmd=c(cv.l.lmd, lmd)
        }
      }
      lmd.i[valid.block]=min(cv.l.lmd[which(cv.l.error==min(cv.l.error))])
    }
    lmd=mean(lmd.i[lmd.i<Inf])
    # lmd.EC.Dantzig <- 0.682
    lmd.EC.Dantzig <- lmd
    lmd.EC.Dantzig.list[[t]] = lmd.EC.Dantzig
  }
  toc()
  
  save(lmd.EC.Dantzig.list, file=paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
  
}

estimate_rollingwindow_dantzig_eigenvectorcentrality = function(file_name = "SP500 securities_up_20230306.csv", 
                                                                window_size=500){
  
  tic("Estmate Rolling Window Eigenvector Centrality by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter of Dantzig Selector for Eigenvector Centrality, using rolling window
  
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  
  load(paste0("Dantzig_lambda_rolling_window_windowsize",window_size,".RData"))
  
  #   EC_DS <- eigenvector centrality estimated by Dantzig selector
  EC_DS<-list()
  a=c()
  for (t in 1: T.windows) {
    print(t)
    
    lmd.EC.Dantzig = lmd.EC.Dantzig.list[[t]] 
    
    EC_DS[[t]] =linfun3_1(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p),
                          rep(0,p),
                          lambda=lmd.EC.Dantzig,
                          abs(eigen(C_in[[t]]-diag(1,p,p)-diag(max(eigen(C_in[[(t)]])$value),p,p))$vector[1,1])
    )
    EC_DS[[t]]=EC_DS[[t]]/max(EC_DS[[t]])
    a[t]=sum(EC_DS[[t]]==0)
  }
  a
  eigenvector_centrality = list("eigenvector_absolute_value"=EC_in,"eigenvector_centrality_Dantzig"=EC_DS,
                                "zeors_in_EC_DS"=a,"covariance_matrix"=COV_in,"correlation_matrix"=C_in,
                                "expected_return"=ER_in)
  save(eigenvector_centrality,file = paste0("eigenvector_centrality_WS",window_size,"_20250117.RData"))
  
  toc()
}


#### Dantzig selector estimation for portfolio parameters ####

##### CV for tuning lambdas in portfolio solutions

##### CV for tuning lambda in estimation Sigma^-1 1 using the first "windowsize"(500/750/1000) data ######

estimate_dantzig_lambda1_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Use for recurrence
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/5)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:5) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    mu.train=rep(1,p)
    mu.valid=rep(1,p)
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    time_threshold = 300
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      cat("Iteration", i, "of", l.lambda, "\n")
      
      # Set a time limit for the computation
      tryCatch({
        setTimeLimit(elapsed = time_threshold, transient = TRUE)
        lin.train <- linfun1(cov.train, mu.train, lmd)
        if (!all(lin.train == 0)) {
          error <- sum((cov.valid %*% lin.train - mu.valid)^2)
          cv.l.error <- c(cv.l.error, error)
          cv.l <- c(cv.l, lmd)
        }
      }, error = function(e) {
        cat("Iteration", i, "exceeded time limit. Skipping to next iteration.\n")
      }, finally = {
        # Reset the time limit after each iteration
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      })
    }
      
      # # Measure the time taken for each iteration
      # time_taken <- system.time({
      #   lin.train=linfun1(cov.train,mu.train,lmd)
      #   if(!(all(lin.train==0))){
      #     error=sum((cov.valid%*%lin.train-mu.valid)^2)
      #     cv.l.error=c(cv.l.error,error)
      #     cv.l=c(cv.l,lmd)
      #   }
      # })
      # 
      # # If the time taken exceeds the threshold, skip to the next iteration
      # if(time_taken["elapsed"] > time_threshold){
      #   cat("Iteration", i, "took too long. Skipping to next iteration.\n")
      #   next
      # }
     
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd1=mean(lmd.i[lmd.i<Inf])
  
  
  dantzig_lambda1_portfolio_parameter = list("lmd1"=lmd1)
  save(dantzig_lambda1_portfolio_parameter,file = paste0("dantzig_lambda1_portfolio_parameter_WS",window_size,"_20250117.RData"))
  
  toc()
}

estimate_dantzig_lambda1_portfolio_parameter_simplified = function(returnstd, window_size){
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/5)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:5) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    mu.train=rep(1,p)
    mu.valid=rep(1,p)
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    time_threshold = 300
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      cat("Iteration", i, "of", l.lambda, "\n")
      
      lin.train <- linfun1(cov.train, mu.train, lmd)
      if (!all(lin.train == 0)) {
        error <- sum((cov.valid %*% lin.train - mu.valid)^2)
        cv.l.error <- c(cv.l.error, error)
        cv.l <- c(cv.l, lmd)
      }
    }
    
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd1=mean(lmd.i[lmd.i<Inf])
  return(lmd1)
}



##### CV for tuning lambda in estimation Sigma^-1 mu, using the first "windowsize"(500/750/1000) data ######

estimate_dantzig_lambda2_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/5)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:5) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    mu.train=colMeans(returnstd.train)
    mu.valid=colMeans(returnstd.valid)
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      
      
      lin.train=linfun1(cov.train,mu.train,lmd)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l=c(cv.l,lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd2=mean(lmd.i[lmd.i<Inf])
  
  dantzig_lambda2_portfolio_parameter = list("lmd2"=lmd2)
  save(dantzig_lambda2_portfolio_parameter,file = paste0("dantzig_lambda2_portfolio_parameter_WS",window_size,"_20250117.RData"))
  
  toc()
}

estimate_dantzig_lambda2_portfolio_parameter_simplified = function(returnstd, window_size){
  
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/3)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:3) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    mu.train=colMeans(returnstd.train)
    mu.valid=colMeans(returnstd.valid)
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      
      lin.train=linfun1(cov.train,mu.train,lmd)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l=c(cv.l,lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd2=mean(lmd.i[lmd.i<Inf])
  return(lmd2)
}

##### CV for tuning lambda in estimation Sigma^-1 phi, using the first "windowsize"(500/750/1000) data ######
estimate_dantzig_lambda3_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Load eigenvector centrality
  load(paste0("eigenvector_centrality_WS",window_size,"_20250117.RData"))
  EC_DS = eigenvector_centrality$eigenvector_centrality_Dantzig
  
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/5)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:5) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    mu.train=EC_DS[[1]]
    mu.valid=EC_DS[[1]]
    lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      # lmd=0.1
      lin.train=linfun1(cov.train,mu.train,lmd)
      # sum(lin.train==0)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l=c(cv.l,lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd3=mean(lmd.i[lmd.i<Inf])
  
  dantzig_lambda3_portfolio_parameter = list("lmd3"=lmd3)
  save(dantzig_lambda3_portfolio_parameter,file = paste0("dantzig_lambda3_portfolio_parameter_WS",window_size,"_20250117.RData"))
  
  toc()
  
}

estimate_dantzig_lambda3_portfolio_parameter = function(returnstd, window_size){
  
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/3)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  lmd.i=c()
  for (valid.block in 1:3) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    mu.train=EC_DS[[1]]
    mu.valid=EC_DS[[1]]
    lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
    l.lambda=length(lambda.grid)
    cv.l.error=NULL
    cv.l=NULL
    for(i in 1:l.lambda){
      lmd=lambda.grid[i]
      print(i)
      # lmd=0.1
      lin.train=linfun1(cov.train,mu.train,lmd)
      # sum(lin.train==0)
      if(!(all(lin.train==0))){
        error=sum((cov.valid%*%lin.train-mu.valid)^2)
        cv.l.error=c(cv.l.error,error)
        cv.l=c(cv.l,lmd)
      }
    }
    lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
  }
  lmd3=mean(lmd.i[lmd.i<Inf])
  return(lmd3)
}

##### CV for tuning parameter in glasso, using the first "windowsize"(500/750/1000) data #####

estimate_glasso_rho_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV for tuning parameter in glasso")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  n<-dim(returnstd[1:window_size,])[1]
  B=floor(window_size/5)
  n.block=floor(n/B)
  block.start=1+(0:(n.block-1))*B
  # valid.block=sort(sample(1:n.block,floor(n.block/4)))
  rho.i=c()
  for (valid.block in 1:5) {
    valid.ind=NULL
    for(k in valid.block){
      valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
    }
    n.valid=length(valid.ind)
    train.ind=setdiff(1:n,valid.ind)
    n.train=length(train.ind)
    returnstd.train=returnstd[1:window_size,][train.ind,]
    returnstd.valid=returnstd[1:window_size,][valid.ind,]
    mu.train=colMeans(returnstd.train)
    mu.valid=colMeans(returnstd.valid)
    cov.train=cov(returnstd.train)
    cov.valid=cov(returnstd.valid)
    rho.grid=seq(0,0.8,length=101)[2:101]
    l.rho=length(rho.grid)
    cv.rho.error=NULL
    cv.rho=NULL
    for (i in 1:l.rho){
      rho=rho.grid[i]
      prec.glasso=glasso(cov.train,rho=rho)$wi
      error=sum((prec.glasso%*%cov.valid-diag(rep(1,p)))^2)
      cv.rho.error=c(cv.rho.error,error)
      cv.rho=c(cv.rho,rho)
      print(i)
    }
    rho.i[valid.block]=cv.rho[which(cv.rho.error==min(cv.rho.error))]
  }
  rho=mean(rho.i[rho.i<Inf])
  
  glasso_rho_portfolio_parameter = list("rho"=rho)
  save(glasso_rho_portfolio_parameter,file = paste0("glasso_rho_portfolio_parameter_WS",window_size,"_20250117.RData"))
  
  toc()
}

##### Estimate portfolio parameters #####

estimate_portfolio_parameter = function(EC_file_name = "eigenvector_centrality_WS750_20250117.RData", 
                                        DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS750_20250117.RData", 
                                        DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS750_20250117.RData", 
                                        DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS750_20250117.RData", 
                                        window_size=750){
  
  tic("Estimate Portfolio Parameters rho1, rho2, rho3 by Dantzig-type selector")
  
  load(EC_file_name)
  EC_in=eigenvector_centrality$eigenvector_absolute_value
  EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
  ER_in=eigenvector_centrality$expected_return
  C_in=eigenvector_centrality$correlation_matrix
  COV_in=eigenvector_centrality$covariance_matrix
  
  load(DS_lambda1_file_name)
  load(DS_lambda2_file_name)
  load(DS_lambda3_file_name)
  lmd1 = dantzig_lambda1_portfolio_parameter$lmd1
  lmd2 = dantzig_lambda2_portfolio_parameter$lmd2
  lmd3 = dantzig_lambda3_portfolio_parameter$lmd3
  
  rho1<-list()
  rho2<-list()
  rho3<-list()
  p <- dim(COV_in[[1]])[1]
  for(t in 1: length(COV_in)){
    print(t)
    ptm<-proc.time()
    ## compute global minimum variance portfolio ##
    rho1[[t]] =linfun1(COV_in[[t]],rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
    # print('rho1')
    rho2[[t]] =linfun1(COV_in[[t]],ER_in[[t]],lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
    # print('rho2')
    rho3[[t]] =linfun1(COV_in[[t]],EC_DS[[t]],lambda=lmd3) # lambda <= 0.1 will lead to be infeasible
    ptm<-proc.time()-ptm
    print(ptm)
  }
  rho<-list("rho1"=rho1,
              "rho2"=rho2,
              "rho3"=rho3)
  save(rho,file=paste0("rho_Dantzig_portfolio_WS",window_size,"_20250125.RData"))
  
  toc()
}

##### CV for tuning lambda in estimation Sigma^-1 1 using the rolling "windowsize"(500/750/1000) data ######

estimate_rollingwindow_dantzig_lambda1_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", 
                                                                      window_size=375){
  
  tic("CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter lambda1 of Dantzig Selector for portfolio parameter phi_1 estimation, using rolling window
  
  tic("Tuning parameter of Dantzig estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  lmd1.Dantzig.list = list()
  for (t in c(1:T.windows)) {
    
    # Use for recurrence
    n<-dim(W_in[[(t)]])[1]
    B=floor(window_size/3) # 3-fold cross validation
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    # valid.block=sort(sample(1:n.block,floor(n.block/4)))
    lmd.i=c()
    for (valid.block in 1:3) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[(t)]][train.ind,]
      returnstd.valid=W_in[[(t)]][valid.ind,]
      mu.train=rep(1,p)
      mu.valid=rep(1,p)
      cov.train=cov(returnstd.train)
      cov.valid=cov(returnstd.valid)
      lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l=NULL
      time_threshold = 300
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        cat("Iteration", i, "of", l.lambda, "\n")
        
        # Set a time limit for the computation
        tryCatch({
          setTimeLimit(elapsed = time_threshold, transient = TRUE)
          lin.train <- linfun1(cov.train, mu.train, lmd)
          if (!all(lin.train == 0)) {
            error <- sum((cov.valid %*% lin.train - mu.valid)^2)
            cv.l.error <- c(cv.l.error, error)
            cv.l <- c(cv.l, lmd)
          }
        }, error = function(e) {
          cat("Iteration", i, "exceeded time limit. Skipping to next iteration.\n")
        }, finally = {
          # Reset the time limit after each iteration
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        })
      }
      
      # # Measure the time taken for each iteration
      # time_taken <- system.time({
      #   lin.train=linfun1(cov.train,mu.train,lmd)
      #   if(!(all(lin.train==0))){
      #     error=sum((cov.valid%*%lin.train-mu.valid)^2)
      #     cv.l.error=c(cv.l.error,error)
      #     cv.l=c(cv.l,lmd)
      #   }
      # })
      # 
      # # If the time taken exceeds the threshold, skip to the next iteration
      # if(time_taken["elapsed"] > time_threshold){
      #   cat("Iteration", i, "took too long. Skipping to next iteration.\n")
      #   next
      # }
      
      lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
    }
    lmd1=mean(lmd.i[lmd.i<Inf])
    lmd1.Dantzig.list[[t]] <- lmd1
  }
  
  
  dantzig_lambda1_portfolio_parameter = list("lmd1.list"=lmd1.Dantzig.list)
  save(dantzig_lambda1_portfolio_parameter,file = paste0("dantzig_lambda1_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))
  
  toc()
}

# This version controls time for each  optimization lin1_1
estimate_rollingwindow_dantzig_lambda1_portfolio_parameter_timecontrol = function(file_name = "SP500 securities_up_20230306.csv", 
                                                                      window_size=375){
  
  tic("CV hyperparameter lambda1 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter lambda1 of Dantzig Selector for portfolio parameter phi_1 estimation, using rolling window
  
  tic("Tuning parameter of Dantzig estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  lmd1.Dantzig.list = list()
  for (t in c(1:T.windows)) {
    
    # Use for recurrence
    n<-dim(W_in[[(t)]])[1]
    B=floor(window_size/3) # 3-fold cross validation
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    # valid.block=sort(sample(1:n.block,floor(n.block/4)))
    lmd.i=c()
    for (valid.block in 1:3) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[(t)]][train.ind,]
      returnstd.valid=W_in[[(t)]][valid.ind,]
      mu.train=rep(1,p)
      mu.valid=rep(1,p)
      cov.train=cov(returnstd.train)
      cov.valid=cov(returnstd.valid)
      lambda.grid=seq(0.1, max(abs(mu.train)),length=101)[2:100]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l=NULL
      time_threshold = 1000
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        cat("Iteration", i, "of", l.lambda, "\n")
        
        lin.train <- linfun1_1(cov.train, mu.train, lmd)
        if (!all(lin.train == 0)) {
          error <- sum((cov.valid %*% lin.train - mu.valid)^2)
          cv.l.error <- c(cv.l.error, error)
          cv.l <- c(cv.l, lmd)
        }
      }
      
      lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
    }
    lmd1=mean(lmd.i[lmd.i<Inf])
    lmd1.Dantzig.list[[t]] <- lmd1
  }
  
  
  dantzig_lambda1_portfolio_parameter = list("lmd1.list"=lmd1.Dantzig.list)
  save(dantzig_lambda1_portfolio_parameter,file = paste0("dantzig_lambda1_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))
  
  toc()
}

##### CV for tuning lambda in estimation Sigma^-1 mu, using the rolling "windowsize"(500/750/1000) data ######

estimate_rollingwindow_dantzig_lambda2_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda2 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter lambda2 of Dantzig Selector for portfolio parameter phi_2 estimation, using rolling window
  
  tic("Tuning parameter of Dantzig estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  lmd2.Dantzig.list = list()
  for (t in c(1:T.windows)) {
    
    n<-dim(W_in[[(t)]])[1]
    B=floor(window_size/3) # 3-fold cross validation
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    # valid.block=sort(sample(1:n.block,floor(n.block/4)))
    lmd.i=c()
    for (valid.block in 1:3) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[(t)]][train.ind,]
      returnstd.valid=W_in[[(t)]][valid.ind,]
      mu.train=colMeans(returnstd.train)
      mu.valid=colMeans(returnstd.valid)
      cov.train=cov(returnstd.train)
      cov.valid=cov(returnstd.valid)
      lambda.grid=seq(min(max(abs(mu.train))/100,min(abs(mu.train))), 0.01,length=101)[2:100]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l=NULL
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        print(i)
        
        lin.train=linfun1(cov.train,mu.train,lmd)
        if(!(all(lin.train==0))){
          error=sum((cov.valid%*%lin.train-mu.valid)^2)
          cv.l.error=c(cv.l.error,error)
          cv.l=c(cv.l,lmd)
        }
      }
      lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
    }
    lmd2=mean(lmd.i[lmd.i<Inf])
    lmd2.Dantzig.list[[t]] <- lmd2
  }
  
  dantzig_lambda2_portfolio_parameter = list("lmd2.list"=lmd2.Dantzig.list)
  save(dantzig_lambda2_portfolio_parameter,file = paste0("dantzig_lambda2_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))
  
  toc()
}

##### CV for tuning lambda in estimation Sigma^-1 phi, using the rolling "windowsize"(500/750/1000) data ######
estimate_rollingwindow_dantzig_lambda3_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter lambda3 of Dantzig Selector for portfolio parameter phi_1 estimation, using rolling window
  
  tic("Tuning parameter of Dantzig estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
  }
  
  
  # Load eigenvector centrality
  load(paste0("eigenvector_centrality_WS",window_size,"_20250117.RData"))
  EC_DS = eigenvector_centrality$eigenvector_centrality_Dantzig
  
  
  lmd3.Dantzig.list = list()
  for (t in c(1:T.windows)) {
    
    n<-dim(W_in[[(t)]])[1]
    B=floor(window_size/3) # 3-fold cross validation
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    # valid.block=sort(sample(1:n.block,floor(n.block/4)))
    lmd.i=c()
    for (valid.block in 1:3) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[(t)]][train.ind,]
      returnstd.valid=W_in[[(t)]][valid.ind,]
      cov.train=cov(returnstd.train)
      cov.valid=cov(returnstd.valid)
      mu.train=EC_DS[[1]]
      mu.valid=EC_DS[[1]]
      lambda.grid=seq(0.1, max(mu.train),length=101)[1:100]
      l.lambda=length(lambda.grid)
      cv.l.error=NULL
      cv.l=NULL
      for(i in 1:l.lambda){
        lmd=lambda.grid[i]
        print(i)
        # lmd=0.1
        lin.train=linfun1(cov.train,mu.train,lmd)
        # sum(lin.train==0)
        if(!(all(lin.train==0))){
          error=sum((cov.valid%*%lin.train-mu.valid)^2)
          cv.l.error=c(cv.l.error,error)
          cv.l=c(cv.l,lmd)
        }
      }
      lmd.i[valid.block]=min(cv.l[which(cv.l.error==min(cv.l.error))])
    }
    lmd3=mean(lmd.i[lmd.i<Inf])
    lmd3.Dantzig.list[[t]] <- lmd3
  }
  
  dantzig_lambda3_portfolio_parameter = list("lmd3.list"=lmd3.Dantzig.list)
  save(dantzig_lambda3_portfolio_parameter,file = paste0("dantzig_lambda3_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))
  
  toc()
  
}

##### CV for tuning parameter in glasso, using the rolling "windowsize"(500/750/1000) data #####

estimate_rollingwindow_glasso_rho_portfolio_parameter = function(file_name = "SP500 securities_up_20230306.csv", window_size=375){
  
  tic("CV hyperparameter lambda3 for Portfolio Parameter by Dantzig-type selector")
  
  # load data
  prices<-read.csv(file_name)
  ZOO <- zoo(prices[,-1], order.by=as.Date(as.character(prices$Dates), format='%Y-%m-%d'))
  
  #return
  return<- Return.calculate(ZOO, method="log")
  return<- return[-1, ]
  returnstd<-xts(return)
  p=dim(return)[2]
  
  # set label
  node.label=colnames(returnstd)
  names(returnstd) = node.label
  
  # Calibrate the hyperparameter rho of glasso for portfolio construction, using rolling window
  
  tic("Tuning parameter of glasso estimation")
  # rolling window
  W<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W[[(t+1)]]=returnstd[(1+t*22):((window_size+22)+t*22),]
  }
  W_in<-list()
  W_out<-list()
  for(t in 0: (floor((1857-window_size)/22)-1)){
    W_in[[(t+1)]]=W[[t+1]][c(1:window_size),]
    W_out[[(t+1)]]=W[[t+1]][c((window_size+1):(window_size+22)),]
  }
  T.windows<-length(W)
  # correlation matrix, Expected return, covariance matrix
  C_in <- list()
  ER_in <- list()
  COV_in <- list()
  EC_in <- list()
  for(t in 1: length(W_in)){
    C_in[[(t)]] =cor(W_in[[(t)]])
    ER_in[[(t)]] = colMeans(W_in[[(t)]])
    COV_in[[(t)]] = cov(W_in[[(t)]])
    network_port = network.correlation(W_in[[(t)]])
    EC_in[[(t)]] <- eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector
    max(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    min(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
    boxplot(eigen_centrality(network_port,directed = FALSE, scale = TRUE)$vector)
  }
  
  rho.glasso.list = list()
  for (t in c(1:T.windows)) {
    
    n<-dim(W_in[[(t)]])[1]
    B=floor(window_size/3) # 3-fold cross validation
    n.block=floor(n/B)
    block.start=1+(0:(n.block-1))*B
    # valid.block=sort(sample(1:n.block,floor(n.block/4)))
    rho.i=c()
    for (valid.block in 1:3) {
      valid.ind=NULL
      for(k in valid.block){
        valid.ind=c(valid.ind,block.start[k]:(min(block.start[k]+B-1,n)))
      }
      n.valid=length(valid.ind)
      train.ind=setdiff(1:n,valid.ind)
      n.train=length(train.ind)
      returnstd.train=W_in[[(t)]][train.ind,]
      returnstd.valid=W_in[[(t)]][valid.ind,]
      mu.train=colMeans(returnstd.train)
      mu.valid=colMeans(returnstd.valid)
      cov.train=cov(returnstd.train)
      cov.valid=cov(returnstd.valid)
      rho.grid=seq(0,0.8,length=101)[2:101]
      l.rho=length(rho.grid)
      cv.rho.error=NULL
      cv.rho=NULL
      for (i in 1:l.rho){
        rho=rho.grid[i]
        prec.glasso=glasso(cov.train,rho=rho)$wi
        error=sum((prec.glasso%*%cov.valid-diag(rep(1,p)))^2)
        cv.rho.error=c(cv.rho.error,error)
        cv.rho=c(cv.rho,rho)
        print(i)
      }
      rho.i[valid.block]=cv.rho[which(cv.rho.error==min(cv.rho.error))]
    }
    rho=mean(rho.i[rho.i<Inf])
    rho.glasso.list[[t]] <- rho
  }
  
  glasso_rho_portfolio_parameter = list("rho.list"=rho.glasso.list)
  save(glasso_rho_portfolio_parameter,file = paste0("glasso_rho_portfolio_parameter_WS",window_size,"_rollingwindow_20250126.RData"))
  
  toc()
}

##### Estimate rolling window portfolio parameters #####

estimate_rollingwindow_portfolio_parameter = function(EC_file_name = "eigenvector_centrality_WS750_20250117.RData", 
                                        DS_lambda1_file_name = "dantzig_lambda1_portfolio_parameter_WS750_20250117.RData", 
                                        DS_lambda2_file_name = "dantzig_lambda2_portfolio_parameter_WS750_20250117.RData", 
                                        DS_lambda3_file_name = "dantzig_lambda3_portfolio_parameter_WS750_20250117.RData", 
                                        window_size=750){
  
  tic("Estimate Portfolio Parameters rho1, rho2, rho3 by Dantzig-type selector")
  
  load(EC_file_name)
  EC_in=eigenvector_centrality$eigenvector_absolute_value
  EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
  ER_in=eigenvector_centrality$expected_return
  C_in=eigenvector_centrality$correlation_matrix
  COV_in=eigenvector_centrality$covariance_matrix
  
  load(DS_lambda1_file_name)
  load(DS_lambda2_file_name)
  load(DS_lambda3_file_name)
  lmd1.list = dantzig_lambda1_portfolio_parameter$lmd1.list
  lmd2.list = dantzig_lambda2_portfolio_parameter$lmd2.list
  lmd3.list = dantzig_lambda3_portfolio_parameter$lmd3.list
  
  theta1<-list()
  theta2<-list()
  theta3<-list()
  p <- dim(COV_in[[1]])[1]
  for(t in 1: length(COV_in)){
    print(t)
    ptm<-proc.time()
    lmd1 <- lmd1.list[[t]]
    lmd2 <- lmd2.list[[t]]
    lmd3 <- lmd3.list[[t]]
    ## compute global minimum variance portfolio ##
    theta1[[t]] =linfun1(COV_in[[t]],rep(1,p),lambda=lmd1) # lambda <= 0.1 will lead to be infeasible
    # print('theta1')
    theta2[[t]] =linfun1(COV_in[[t]],ER_in[[t]],lambda=lmd2) # lambda <= 0.1 will lead to be infeasible
    # print('theta2')
    theta3[[t]] =linfun1(COV_in[[t]],EC_DS[[t]],lambda=lmd3) # lambda <= 0.1 will lead to be infeasible
    ptm<-proc.time()-ptm
    print(ptm)
  }
  theta<-list("theta1"=theta1,
            "theta2"=theta2,
            "theta3"=theta3)
  save(theta,file=paste0("theta_Dantzig_portfolio_WS",window_size,"_20250126.RData"))
  
  toc()
}
