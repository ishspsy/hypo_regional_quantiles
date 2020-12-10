real_main_function = function(taus, n,p,X_tr, Y_tr,t_point){
  
  if (taus == "low"){
    #fix this part for each tau
    tau = seq(0.1,0.3,length.out=11)
    tau0 = tau
    
    # length.out can be changed
    tau = seq(0.1,0.3,length.out=11)
    K=length(tau)
  }
  
  if (taus == "middle"){
    #fix this part for each tau
    tau = seq(0.4,0.6,length.out=11)
    tau0 = tau
    
    # length.out can be changed
    tau = seq(0.4,0.6,length.out=11)
    K=length(tau)
  }
  
  if (taus == "high"){
    #fix this part for each tau
    tau = seq(0.7,0.9,length.out=11)
    tau0 = tau
    
    # length.out can be changed
    tau = seq(0.7,0.9,length.out=11)
    K=length(tau)
  }
  
  plot_save=0
  
  ################################################ THIS PART HERE~
  k <- p- length(ind00)     #number of free variable

  # Applying BIC
  fun_cl = function(ttype,psp,psp2,l,tcoef_mat){
    coefm=array(0,c(nrow(psp),nrow(psp2)))
    for (j in 1:nrow(psp)){
      for (k in 1:nrow(psp2)){
        coefm[j,k] = matrix(tcoef_mat[l,], nrow=1) %*% matrix(kronecker(psp2[k,], psp[j,]), ncol=1)
      }}
    return(coefm)}    
  
  para_bic = function(tau,tau0,lout,tvas,leng_ind,t_point,ttype,X_tr,Y_tr,bsp,bsp2,p){
    bbsp=bsplineS(tau, seq(min(tau0)-0.02, max(tau0)+0.02,length.out=lout[leng_ind]), norder=2, nderiv=0, returnMatrix=FALSE)
    bbsp2=bsplineS(t_point, seq(0, 1, length.out=lout[leng_ind]), norder=2, nderiv=0, returnMatrix=FALSE)
    
    bsp=bbsp
    bsp2=bbsp2
    
    initial_saving = fn_estimation(ttype=ttype, X_tr,Y_tr,tau,bsp, bsp2, 0.01)

    coef_ini = initial_saving$coef
    
    coef_mat= matrix(coef_ini, ncol=p, byrow=T)
    tcoef_mat = t(coef_mat)
    
    coef_mat_set=mclapply(1:p, fun_cl, ttype=ttype,psp=bsp, psp2=bsp2,tcoef_mat=tcoef_mat, mc.cores=2)
    
    fit_tot=0
    for (l in 1:K){
      sumsi=0
      for (ii in 1:n){
        sums = 0
        for (jj in 1:p){
          sums= sums + X_tr[ii,jj] * coef_mat_set[[jj]][l,ii]
        }
        sumsi=sumsi+tau[l]*(Y_tr[ii] - sums)*(Y_tr[ii] - sums >=0) + (tau[l]-1)*(Y_tr[ii] - sums)*(Y_tr[ii] - sums <0)
      }
      fit_tot = fit_tot+sumsi
    }
    return(log(fit_tot) + p* (lout[leng_ind])^2 * log(n*K)/ (2*n*K))}
  
  
  lout = 6:(K-2)

  # BIC perform or not?
  BIC_index = 1
  bic_check=0
  
  if (bic_check == 1){
    BIC_set = mclapply(1:length(lout),para_bic,tau=tau,tau0=tau0,lout=lout,tvas=tvas,t_point=t_point,ttype=ttype,X_tr=X_tr,Y_tr=Y_tr,bsp=bsp,bsp2=bsp2,p=p,mc.cores=2)
    BIC_set = unlist(BIC_set)
    BIC_index = which.min(BIC_set)}
  
  opt_num = lout[BIC_index]
  
  # using BIC optimal result
  bsp=bsplineS(tau, seq(min(tau0)-0.02, max(tau0)+0.02,length.out=opt_num), norder=2, nderiv=0, returnMatrix=FALSE)
  bsp2=bsplineS(t_point, seq(0, 1, length.out=opt_num), norder=2, nderiv=0, returnMatrix=FALSE)
  
  initial_saving = fn_estimation(ttype=ttype, X_tr,Y_tr,tau,bsp, bsp2, 0.05)
  coef_ini = initial_saving$coef
  
  coef_mat= matrix(coef_ini, ncol=p, byrow=T)
  tcoef_mat = t(coef_mat)
  
  psp=bsplineS(tau0,seq(min(tau0)-0.02, max(tau0)+0.02,length.out=opt_num), norder=2, nderiv=0, returnMatrix=FALSE)
  psp2=bsplineS(seq(0,1,0.05), seq(0, 1, length.out=opt_num), norder=2, nderiv=0, returnMatrix=FALSE)
  
  coef_mat_set=mclapply(1:p, fun_cl, ttype=ttype,psp=psp, psp2=psp2,tcoef_mat=tcoef_mat, mc.cores=2)
  
  return(list(initial_saving=initial_saving,X_tr=X_tr,Y_tr=Y_tr, coef_mat_set=coef_mat_set,taus=taus,t_point=t_point,p=p))
}



fn_estimation =function(ttype,XXX_tr,Y_tr,tau,bsp, bsp2, tva)
{ 
  K <- length(tau)
  p <- dim(XXX_tr)[2]
  n=dim(XXX_tr)[1]
  
  X_tr=XXX_tr 
  X_tr= as.matrix(X_tr)

  XX_tr=NULL
  for (k in 1:K)
  {
    for (i in 1:n)
    {
      XX_tr=rbind(XX_tr,kronecker(as.vector(kronecker(as.vector(bsp2[i,]),as.vector(bsp[k,]))),as.vector(X_tr[i,])))
    }
  }
  
  v1k= XX_tr 
  ddd = fn(v1k, Y_tr, n, p= ncol(v1k), tau, K, tva)
  return(list(coef=ddd, data=v1k))  
}



fn=function(X_tr,Y_tr,n,p,tau,K,tva)
{
  I=diag(1,n)
  II=diag(1,p)
  III=diag(1,p-1)
  
  Tau=NULL
  for(i in 1:K)
  {Tau=c(Tau,c(rep(tau[i],n),rep((1-tau[i]),n)))}
  
  FF1=NULL
  for(k in 1:K)
  {
    FF=array(0,c(n,2*n*K+ 2*p))
    for(i in 1:n)
    {
      FF[i,]=c(rep(0,2*n*(k-1)),I[i,],-I[i,],rep(0,2*n*(K-k)), X_tr[(k-1)*n+i,],-X_tr[(k-1)*n+i,])
    }
    FF1=rbind(FF1,FF)
  }
  
  FF3=diag(1,2*n*K+2*p)
  
  FF=rbind(FF1,FF3)
  f.con=1*round(FF,2)
  
  f.rhs1=matrix(rep(Y_tr,K))
  f.rhs3=matrix(rep(0, 2*n*K+2*p))
  f.rhs=1*round(rbind(f.rhs1,f.rhs3),2)
  
  f.dir1=matrix(rep("=",n*K))
  f.dir3=matrix(rep(">=",2*n*K+2*p))
  f.dir=rbind(f.dir1,f.dir3)
  
  f.obj=c(round(1*Tau,2),rep(0.2,2*p))
  f.obj=c(round(1*Tau,2),rep(tva,2*p))
  
  coef=lp("min", f.obj, f.con, f.dir, f.rhs)$solution
  ccoef=coef[(2*n*K+1):(2*n*K+p)]-coef[(2*n*K+p+1):(2*n*K+2*p)]
  return(ccoef)}




dens=function(X_tr,Y_tr,tau,lenout)
{
  tau = seq(tau[1], tau[length(tau)], length.out=lenout)
  X_tr = as.matrix(X_tr)
  n=dim(X_tr)[1]
  p=dim(X_tr)[2]
  K=length(tau)
  fff2=NULL
  for(k in 1:K)
  {
    hn=min(n^(-1/6),tau[k]*(1-tau[k])/2)
    beta1=rq.fit(X_tr,Y_tr,tau=tau[k]+hn)$coefficients
    beta2=rq.fit(X_tr,Y_tr,tau=tau[k]-hn)$coefficients
    ff2=2*hn/(X_tr%*% matrix(beta1,ncol=1) -X_tr%*%matrix(beta2,ncol=1) -0.01)
    
    for(i in 1:n)
    {
      if((ff2[i]>0)*(ff2[i]<100)==1){ff2[i]=ff2[i]}
      else{ff2[i]=0}
    }
    
    fff2=rbind(fff2,ff2)
  }
  return(abs(fff2))
}

denss=function(X_tr,tau,sdd,lenout)
{
  tau = seq(tau[1], tau[length(tau)], length.out=lenout)
  n=dim(X_tr)[1]
  K=length(tau)
  fff=NULL
  for(k in 1:K)
  {
    ff=rep(dnorm(qnorm(tau[k])),n)
    fff=c(fff,ff)
  }
  return(fff)
}



