library(phytools)
library(rBayesianOptimization)
library(castor)

##load tree files

trees<-read.newick("msprime trees constant/Tree_1.newick")
for (i in 2:10){
  trees<-c(trees,read.newick(paste0("msprime trees constant/Tree_",i,".newick")))
}

##get initial parameters

k<-length(G)
n<-G[[1]]$Nnode+1
bt<-matrix(nrow = k,ncol = n)
for (i in 1:k){
  b<-branching.times(G[[i]])
  b<-sort(b,decreasing = F)
  bt[i,]<-c(0,b)
}
T2<-as.vector(bt[,2:n])
Tsort<-c(0,sort(T2,decreasing = F))
Glist<-list(k,n,bt,Tsort)

## define functions

loglike<-function(theta,A,Glist) {
  theta<-10^(theta)
  k<-Glist[[1]]
  n<-Glist[[2]]
  bt<-Glist[[3]]
  Tsort<-Glist[[4]]
  Apts<-A
  Apts[1]<-1
  for (i in 2:length(A)){
    Apts[i]<-sum(A[1:i-1])+1
  }
  Tpts<-Tsort[Apts]
  LogL<-matrix(nrow = k,ncol = n-1)
  for (i in 1:k){
    for (j in 1:(n-1)){
      tstart<-bt[i,j]
      tend<-bt[i,j+1]
      w<-max(which((Tpts<=tstart)))
      u<-c(tstart,Tpts[(tend>Tpts)&(tstart<Tpts)],tend)
      u2<-diff(u)
      #u2<-u[2:(length(u))]-u[1:(length(u)-1)]
      v<-c(theta[w],theta[(tend>Tpts)&(tstart<Tpts)])
      v2<-choose(n+1-j,2)/v
      LogL[i,j]<-log(v2[length(v2)])-sum(v2*u2)
    }
  }
  LogLtot<-sum(LogL)
  return(LogLtot)
}

grouper<-function(eps,Glist){
  Tsort<-Glist[[4]]
  Tsort<-log10(Tsort+1)
  npts<-Glist[[1]]*(Glist[[2]]-1)
  A<-numeric(npts)+1
  u<-diff(Tsort)
  if (max(Tsort)<eps){
    A2<-npts
  } else {
    while (min(u,na.rm = T)<eps) {
      a<-min(which(u<eps))
      if (sum(u[a:npts])<eps){
        b<-sum(u[a:npts])
        u[(a):npts]<-NA
        A[(a):npts]<-0
        a<-max(which(!is.na(u)))
        u[a]<-u[a]+b
        break
      } else {
        i<-1
        while (sum(u[a:(a+i)])<eps) {
          i<-i+1
        }
        u[a]<-sum(u[a:(a+i)])
        u[(a+1):(a+i)]<-NA
        A[(a+1):(a+i)]<-0
      }
    }
    A2<-which(A==1)
    A2<-diff(c(A2,npts+1))
  }
  return(A2)
}

optimfn<-function(eps){
  A<-grouper(eps,Glist)
  m<-length(A)
  k<-Glist[[1]]
  n<-Glist[[2]]
  theta<-numeric(m)+log10(10000)
  optimfn2<-function(theta0,index,theta,A,Glist){
    theta[index]<-theta0
    value<-loglike(theta,A,Glist)
    return(-value)
  }
  for (i in 1:m){
    invisible(capture.output(optim2_out<-optimise(optimfn2,c(0,10),index = i, theta = theta,A = A, Glist  = Glist)))
    theta[i]<-optim2_out$minimum
  }
  LLL<--optim2_out$objective
  optimval<-LLL-m-(m*(m+1))/(k*(n-1)-m-1)
  out<-list(Score = optimval,Pred = 0)
  #out<--optimval
  return(out)
}

optimfn_login<-function(logeps){
  eps<-10^(logeps)
  out<-optimfn(eps)
  return(out)
}

##optimise

optim_out<-BayesianOptimization(optimfn_login,bounds = list(logeps = c(-2,1)),init_points = 10,n_iter = 10)

##evaluate MLE effective population size

epsopt<-10^optim_out$Best_Par
A<-grouper(epsopt,Glist)
m<-length(A)
k<-Glist[[1]]
n<-Glist[[2]]
theta_estim<-numeric(m)+4
optimfn2<-function(theta0,index,theta,A,Glist){
  theta[index]<-theta0
  value<-loglike(theta,A,Glist)
  return(-value)
}
st<-Sys.time()
for (i in 1:m){
  invisible(capture.output(optim2_out<-optimise(optimfn2,c(0,10),index = i, theta = theta_estim,A = A, Glist  = Glist)))
  theta_estim[i]<-optim2_out$minimum
}
ed<-Sys.time()
time_taken<-ed-st
Apts<-A
Apts[1]<-1
for (i in 2:length(A)){
  Apts[i]<-sum(A[1:i-1])+1
}
Tpts<-Tsort[Apts]

x<-numeric(2*length(A))
for (i in 1:(length(A)-1)){
  x[2*i]<-Tpts[i+1]
  x[2*i+1]<-Tpts[i+1]
}
x[1]<-Tpts[1]
x[2*length(A)]<-max(Tsort)

y<-numeric(2*length(A))
for (i in 1:(length(A))){
  y[2*i]<-theta_estim[i]
  y[2*i-1]<-theta_estim[i]
}
plot(x[1:length(y)],y,main = paste("Optimised Skyline for eps=",epsopt),xlab = "time",ylab = "log10(pop)")
lines(x[1:length(y)],y)

x_real<-seq(0,max(x),by=100)
y_real<-log(41000)+numeric(length(x_real))
#y_real<-log(100000)-0.3*x_real
lines(x_real,y_real/log(10),lty=2)

gsp1<-skylineplot.deluxe(trees[[1]])
gsp_step<-stepfun(c(0,gsp1$time),c(0,log10(gsp1$population.size),0))

real_df<-data.frame(x = x_real,y=y_real/log(10))
gsp_df<-data.frame(x = x_real,y = gsp_step(x_real))
mtgsp_df<-data.frame(x = x[1:length(y)],y=y)
real_df$x[1]<-0.1
gsp_df$x[1]<-0.1
mtgsp_df$x[1]<-0.1

for (i in 1:10){
  gsp2<-skylineplot.deluxe(trees[[i]])
  gsp_step2<-stepfun(c(0,gsp2$time),c(0,log10(gsp2$population.size),log10(gsp2$population.size)[length(gsp2$population.size)]))
  gsp_df<-data.frame(x = x_real,y = gsp_step2(x_real))
  gsp_df$x[1]<-0.1
  assign(paste0("gsp_df",i),gsp_df)
}

##plotting

library(ggplot2)
library(scales)
library(gridExtra)

p = ggplot() + 
  geom_line(data = real_df, aes(x = x, y = 10^y,color = "population simulated",linetype = "population simulated"),size = 1.2) +
  geom_line(data = gsp_df1, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1.2) +
  #geom_line(data = gsp_df2, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df3, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df4, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df5, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df6, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df7, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df8, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df9, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  #geom_line(data = gsp_df10, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = mtgsp_df, aes(x = x, y = 10^y,color = "multi-tree generalised skyline",linetype = "multi-tree generalised skyline"),size = 1.2) +
  scale_color_manual("Skyline Method", 
                        values = c("population simulated"="black", "generalised skyline for 1 tree"="red","multi-tree generalised skyline"="blue" )) +
  scale_linetype_manual("", guide = 'none',
                      values = c("population simulated"="dashed", "generalised skyline for 1 tree"="solid","multi-tree generalised skyline"="solid" )) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #coord_cartesian(xlim = c(1, max(x)), ylim = c(1, max(10^(y+1))), expand = FALSE) +
  coord_cartesian(xlim = c(1, max(x)), ylim = c(10^3, max(10^(6))), expand = FALSE) +
  annotation_logticks()  +
  xlab('Time (years ago)') +
  ylab('Effective Population Size') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1), 
        legend.background = element_rect(fill="WHITE", colour = "black",
                                         size=0.5, linetype="solid"),
        legend.text = element_text(size=8,face="bold"),
        legend.title = element_text(size=8,face="bold"))

p2 = ggplot() + 
  geom_line(data = real_df, aes(x = x, y = 10^y,color = "population simulated",linetype = "population simulated"),size = 1.2) +
  geom_line(data = gsp_df1, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df2, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df3, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df4, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df5, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df6, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df7, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df8, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df9, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = gsp_df10, aes(x = x, y = 10^y,color = "generalised skyline for 1 tree",linetype = "generalised skyline for 1 tree"),size = 1) +
  geom_line(data = mtgsp_df, aes(x = x, y = 10^y,color = "multi-tree generalised skyline",linetype = "multi-tree generalised skyline"),size = 1.2) +
  scale_color_manual("Skyline Method", 
                     values = c("population simulated"="black", "generalised skyline for 1 tree"="red","multi-tree generalised skyline"="blue" )) +
  scale_linetype_manual("", guide = 'none',
                        values = c("population simulated"="dashed", "generalised skyline for 1 tree"="solid","multi-tree generalised skyline"="solid" )) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #coord_cartesian(xlim = c(1, max(x)), ylim = c(1, max(10^(y+1))), expand = FALSE) +
  coord_cartesian(xlim = c(1, max(x)), ylim = c(10^3, max(10^(6))), expand = FALSE) +
  annotation_logticks()  +
  xlab('Time (years ago)') +
  ylab('Effective Population Size') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1), 
        legend.background = element_rect(fill="WHITE", colour = "black",
                                         size=0.5, linetype="solid"),
        legend.text = element_text(size=8,face="bold"),
        legend.title = element_text(size=8,face="bold"))



grid.arrange(p,p2,ncol=2)

