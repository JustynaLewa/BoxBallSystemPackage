evolution<-function(k){
  stopifnot(is.numeric(k), is.vector(k), is.numeric(k), length(k)<=30,max(k)<=1, min(k)>=0,length(k)>0,all(k == floor(k)))
  basket<-0
  for (i in 1:length(k)){
    if (k[i]==1){
      k[i]<-0
      basket=basket+1
    } else {
      if (basket>0) {
        k[i]<-1
        basket=basket-1}
    }
  }
  return(k)
}

#\eta->T^{n}\eta
evolutions<-function(k,n){
  if (n>0) {
    for (i in 1:n){
      k<-evolution(k)
    }
  }
  return(k)
}

#plot single configuartion T^{n}\eta
eta<-function(k,n=0){
  stopifnot(is.numeric(k), is.vector(k), is.numeric(k), length(k)<=30,max(k)==1, min(k)==0,length(k)>0,all(k == floor(k)))
  k<-evolutions(k,n)
  x<-seq(from=1, to=length(k),by=1)
  y<-rep(1,length(k))
  plot(x,y,
       xlab="",
       xaxt='n',
       ylab="",
       yaxt='n',
       pch=(k*15)+1,
       axes=F,
       cex=2)
}

#plot dynamics T^{0}\eta, T^{1}\eta,...,T^{n}\eta max n=9
T_eta<-function(k,n=0){
  stopifnot(is.numeric(k), is.vector(k), is.numeric(k), length(k)<=30,max(k)==1, min(k)==0,length(k)>0,all(k == floor(k)))
  x<-seq(from=1, to=length(k),by=1)
  y<-rep((n+1),length(k))
  plot(x,y,
       xlab="",
       xaxt='n',
       ylab="",
       yaxt='n',
       ylim = c(0.5,n+1.5),
       pch=(k*15)+1,
       axes=F,
       cex=2)
  axis(2, at = 1:(n+1), labels = n:0,tick = FALSE,las=1)
  for (j in 1:n){
    k<-evolution(k)
    y<-rep((n-j+1),length(k))
    points(x,y,pch=(k*15)+1,cex=2)
  }
}

#plots walk corresponding to a given configuration
configuration.to.walk.plot<-function(k){
  stopifnot(is.numeric(k), is.vector(k), is.numeric(k), length(k)<=30,max(k)<=1, min(k)>=0,length(k)>0,all(k == floor(k)))
  y<-rep(0,length(k))
  start=0
  for (i in 1:length(k)){
    if (start==0) {
      if(k[i]==1) {start=1
      j=2}
    }
    if (start==1){
      if (k[i]==1){y[j]=y[j-1]+1}
      else {y[j]=y[j-1]-1}
    }
    j=j+1
  }
  x<-seq(from=0, to=length(y)-1,by=1)
  df <- data.frame(x, y)
  ggplot(df, aes(x = x, y = y)) +
    geom_path() +
    labs(x = "",y="")+ scale_x_continuous(breaks=seq(from=0, to=length(y)-1,by=1))
}

idyntify.runs<-function(k){
  stopifnot(is.numeric(k), is.vector(k), is.numeric(k), length(k)<=30,max(k)<=1, min(k)>=0,length(k)>0,all(k == floor(k)),sum(k[which(k==1)])==sum(k[which(k==0)]+1))
  runs<-rep(0,length(k))
  value<-1
  for (i in 2:length(k)){
    if(k[i]==k[i-1]){
      value=value+1
      if (i==length(k)){
        for(l in 0:(value-1)){
          runs[i-l]<-value
        }
      }
    }
    else {
      if(i==length(k)){runs[i]=1}
      for(j in 1:(value)){
        runs[i-j]<-value
      }
      value<-1
    }
  }
  return(runs)
}
idyntify.solitons<-function(k){
  oryginal.index=seq(from=1,to=length(k),by=1)
  temp<-rbind(oryginal.index,k)
  soliton_list <- vector(mode = "list")
  runs.unique<-unique(idyntify.runs(temp[2,]))
  for(l in 1:max(runs.unique)){assign(paste('sol_list', l, sep='_'), vector(mode = "list"))}
  is.temp.empty<-1

  while(is.temp.empty==1){
    runs<-idyntify.runs(temp[2,])
    i<-runs[min(which((runs[-1]-runs[-length(runs)])>=0 & (temp[2,-1]-temp[2,-ncol(temp)])!=0))]
    index<-which(runs==i)
    x<-c(temp[1,seq(from=index[1],length.out=(2*i))])
    if(exists(paste('sol_list', i, sep='_'))==FALSE){assign(paste('sol_list', i, sep='_'), vector(mode = "list"))}
    assign(paste('sol_list', i, sep='_'),c(get(paste('sol_list', i, sep='_')),list(x)))
    soliton_list[[i]]<-get(paste('sol_list', i, sep='_'))
    if(ncol(temp)==length(x)){is.temp.empty=0}
    else {temp<-temp[,-which(temp[1,] %in% x)]}
  }
  return(soliton_list)
}

idyntify.slots<-function(k){
  solitions<-idyntify.solitons(k)
  M<-length(solitions)
  Diagram<-matrix(0,nrow=M,ncol=length(k))
  if(M==1){
    for(j in 1:length(solitions[[1]])){
      ind<-solitions[[1]][[j]]
      Diagram[M,ind[1:2]]<-1
    }
  }
  else{

    for (i in 2:M){
      for(j in 1:length(solitions[[i]])){
        ind<-solitions[[i]][[j]]
        for(k in 1:(i-1)){
          Diagram[M-(k-1),ind[c((k+1):i,(k+i+1):(2*i))]]<-k
        }
      }
    }
  }
  Diagram<-cbind(c(M:1),Diagram,c(M:1)) #adding records columns
  return(Diagram)
}

plot.slots<-function(k){
  x<-idyntify.slots(k)
  x[x==0]<-NA
  NamesRows<-vector(length=nrow(x))
  for (i in 1:nrow(x)) {
    NamesRows[nrow(x)-i+1]=paste(i,"-slots")
  }
  rownames(x)<-NamesRows

  res<-plot(x,asp=TRUE,col=brewer.pal(n = nrow(x), name = "PuRd"),axis.col=NULL, axis.row=list(cex.axis=0.7,las=1),xlab='', ylab='', yaxt = "n",main="Slot Diagram",key=NULL)

  for (i in 1:nrow(x)) {
    val<-0
    for (j in 1:ncol(x)) {
      if(!(is.na(x[i,j])) & j<ncol(x)){
        text(x=j,y=(nrow(x)-i+1),
             labels=paste(val),
             cex = 0.5, col = "black",font=2)
        val = val+1
      }
      if(!(is.na(x[i,j])) & j==ncol(x)){
        text(x=j,y=(nrow(x)-i+1),
             labels="0",
             cex = 0.5, col = "black",font=2)
      }
    }
  }
}
slot.diagram<-function(k){
  M<-idyntify.slots(k)
  M<-M[,-c(1,ncol(M))]
  if(is.vector(M)){M<-as.matrix(t(M))}
  solitions<-idyntify.solitons(k)
  for(i in (nrow(M)):2){
    assign(paste('x', (nrow(M)-i+1), sep='_'),vector(mode='numeric',length=sum(M[i,]>0)+1))
    slots.ind<-which(M[i,]>0)
    l<-1
    if(length(solitions[[nrow(M)-i+1]])>0 && length(slots.ind)>=1){
      for(k in 1:(length(slots.ind)-1)){
        solitions.attached<-0
        if(k==1){
          for (j in 1: length(solitions[[nrow(M)-i+1]]) ){
            if(max(solitions[[nrow(M)-i+1]][[j]])<slots.ind[k]){solitions.attached=solitions.attached+1}
          }
          new.vec<-get(paste('x', (nrow(M)-i+1), sep='_'))
          new.vec[0]<-solitions.attached
          assign(paste('x', (nrow(M)-i+1), sep='_'),new.vec)
          solitions.attached<-0
        }
        if(k==(length(slots.ind)-1)){
          for (j in 1: length(solitions[[nrow(M)-i+1]]) ){
            if(min(solitions[[nrow(M)-i+1]][[j]])>slots.ind[k+1]){solitions.attached=solitions.attached+1}
          }
          new.vec<-get(paste('x', (nrow(M)-i+1), sep='_'))
          new.vec[(sum(M[i,]>0)+1)]<-solitions.attached
          assign(paste('x', (nrow(M)-i+1), sep='_'),new.vec)
          solitions.attached<-0
        }
        for (j in 1: length(solitions[[nrow(M)-i+1]]) ){
          if((min(solitions[[nrow(M)-i+1]][[j]])>slots.ind[k]) && (max(solitions[[nrow(M)-i+1]][[j]])<slots.ind[k+1])){solitions.attached=solitions.attached+1}
        }
        new.vec<-get(paste('x', (nrow(M)-i+1), sep='_'))
        new.vec[l+1]<-solitions.attached
        assign(paste('x', (nrow(M)-i+1), sep='_'),new.vec)
        l<-l+1
      }
    }
  }
  assign(paste('x', (nrow(M)), sep='_'),length(solitions[[nrow(M)]]))
  diagram<- vector(mode = "list")
  for( i in 1:nrow(M)){
    x<-get(paste('x', (nrow(M)-i+1), sep='_'))
    diagram[[nrow(M)-i+1]]<-x
  }
  return (diagram)
}
phi<-function(k,q){
  stopifnot(is.numeric(q), is.vector(q),max(q)<1, min(q)>=0,length(q)>0,is.finite(sum(q)))

  diagram<-slot.diagram(k)
  M<-length(diagram)
  prod<-1
  for ( i in 1:M){
    prod<-prod*(q[i])^{sum(diagram[[i]])}*(1-q[i])^{length(diagram[[i]])}
    q[i]<-0
  }
  for( i in which(q>0)){
    prod<-prod*(1-q[i])
  }
  return(prod)
}
p_m<-function(M,q){
  stopifnot(is.numeric(q), is.vector(q),max(q)<1, min(q)>=0,length(q)>0,is.finite(sum(q)))
  if(M>0){prod<-q[M]}
  else{prod<-1} #q_{0}=1
if (M<(length(q))){
  for ( i in (M+1):length(q)){
    prod<-prod*(1-q[i])
  }
}
  return(prod)
}

p_x_m<-function(x,M,q){
  stopifnot(is.numeric(q), is.vector(q),max(q)<1, min(q)>=0,length(q)>0,is.finite(sum(q)),M<=length(q),length(x)==1,x>0)
  if(M==0){prod<-0}
  else {
    prod<-(q[M]^{x-1})*(1-q[M])
  }
  return (prod)
}

p_x_k<-function(x.k,k,q){
  stopifnot(is.numeric(q), is.vector(q),max(q)<1, min(q)>=0,length(q)>0,is.finite(sum(q)),k<=length(q),k == floor(k),k>0)
  prod<-(q[k]^{sum(x.k)})*(1-q[k])^(length(x.k))
  return (prod)
}

alpha.to.q<-function(alpha){
  stopifnot(is.numeric(alpha), is.vector(alpha),max(alpha)<1, min(alpha)>=0,length(alpha)>0,is.finite(sum(alpha)))
  q<-numeric(length=length(alpha))
  q[1]<-alpha[1]
  for(i in 2:length(alpha)){
    prod<-1
    for(j in 1:(i-1)){
      prod<-prod*((1-q[j])^(2*(i-j)))
    }
    q[i]<-alpha[i]/prod
  }

  if(is.finite(sum(q))&& max(q)<1 && min(q)>=0){return(q)}
  else cat("the resulting q vector does not satisfy the conditions, please redefinie the inpout alpha vector")
}

q.to.alpha<-function(q){
  stopifnot(is.numeric(q), is.vector(q),max(q)<1, min(q)>=0,length(q)>0,is.finite(sum(q)))
  alpha<-numeric(length=length(q))
  alpha[1]<-q[1]
  for ( i in 2:length(q)){
    prod<-1
    for (j in 1:(i-1)){
      prod<-prod*(1-q[j])^(2*(i-j))
    }
    alpha[i]<-q[i]*prod
  }
  return (alpha)
}

