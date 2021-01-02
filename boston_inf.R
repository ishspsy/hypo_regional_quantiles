
## Boston data
data_real_ind = "boston"

## Specify the test among 1,2,3
ttest = 0

## Specify the quantile region among "low", "middle", "high"
taus = "high"   
tva = 0.01

# Specify the quantile levels 
tau =

# You can change tau levels
if (taus == "low"){
  load("resulting_files/data4_real_boston_coef_low.RData")
  tau0 = tau
  K=length(tau)}

if (taus == "middle"){
  load("resulting_files/data4_real_boston_coef_middle.RData")
  tau0 = tau
  K=length(tau)}

if (taus == "high"){
  load("resulting_files/data4_real_boston_coef_high.RData")
  tau0 = tau
  K=length(tau)}

# Recall and save the previous estimation results
data_model = data4_real_coef
coef_mat_set = data4_real_coef$coef_mat_set
initial_saving=data4_real_coef$initial_saving


# For each test, we set the corresponding test type
if (ttest ==1){
  ttype = "1"   # independent of tau and t
  type="known"}

if (ttest ==2){
  ttype = "tau" # independent of t
  type="unknown"}

if (ttest ==3){
  ttype = "t"    # independent of tau
  type="unknown"}



eps0=0; ias_sum = NULL; ias_stat = NULL; 
for (ias in 2:p){
  ind00=ias
  p=data_model$p
  n= nrow(data_model$X_tr)
  
  coef_ini = data_model$initial_saving$coef
  coef_ini = matrix(coef_ini, byrow=F, nrow=p)
  mm = ncol(coef_ini); m=sqrt(mm)
  coef_ini = coef_ini[c(setdiff(1:p,ind00),ind00),]
  
  X_tr = data_model$X_tr[,c(setdiff(1:p,ind00),ind00)]
  Y_tr = data_model$Y_tr
  
  ind00 = 7
  if (ttype == "tau"  &  type == "unknown"){
    H = diag(1,mm) -  repmat(diag(1,m),m) / m
    tH = H[, 1:(m^2-m)]
    gam = length(ind00) * (m^2-m)
  }
  
  if (ttype == "t"  &  type == "unknown"){
    gam = length(ind00) * (m^2-m)
    col_ind = m*rep(0:(m-1), m) + rep(1:m, each=m) 
    coef_ini = coef_ini[,col_ind]
    H = diag(1,mm) -  repmat(diag(1,m),m) / m
    tH = H[, 1:(m^2-m)]}
  
  if (ttype == "1"  &  type == "unknown"){
    gam = length(ind00) * (m^2-1)
    H = diag(1,mm) -  array(1, c(mm,mm)) / mm
    tH = H[,1:(mm-1)]}
  
  
  if (ttype == "1"  &  type == "known"){
    gam = length(ind00) * (m^2)
    stat_mat= coef_ini 
  }
  
  if (ttype != "1"  |  type != "known"){
    stat_mat = coef_ini %*% tH}
  
  matvec = matrix(stat_mat, ncol=1)
  matvec2 = matrix(coef_ini, ncol=1)
  
  
  ffff= dens(data_model$X_tr, data_model$Y_tr, tau, length(tau))
  
  t_point= data_model$t_point
  tau_all = seq(tau[1], tau[length(tau)], length(tau))
  

  bsp=bsplineS(tau_all, seq(min(tau)-0.02, max(tau)+0.02,length.out=6), norder=2, nderiv=0, returnMatrix=FALSE)
  bsp2=bsplineS(t_point, seq(0, 1, length.out=6), norder=2, nderiv=0, returnMatrix=FALSE)
  ffff= dens(data_model$X_tr, data_model$Y_tr, tau,length(tau))
  
  X_tr = as.matrix(X_tr)
  
  tnum = n
  En = array(0,c(p, length(ind00)))
  En[(p-length(ind00)+1):p, ] = diag(1,length(ind00))
  

  ffff2=mean(ffff)

  if (ttype == "1"  &  type == "known"){
    tH = diag(1,mm)
    
    tHEn = kronecker(tH, En)   
    
    Ttheta = t(tHEn) %*% matvec2   
    test_stat = Ttheta  
    stat_zero = (sum((1*test_stat)^2))/(sqrt(2*gam))
    pval = 2*(1-pnorm(abs(stat_zero*1)))
  }
  
  
  if (ttype != "1"  |  type != "known"){
    tHEn = kronecker(tH, En)   
    Ttheta = t(tHEn) %*% matvec2   
    test_stat = Ttheta  
    stat_zero = (sum((1*test_stat)^2))/(sqrt(2*gam))
    pval = 2*(1-pnorm(abs(stat_zero*1)))} 
   ias_sum = c(ias_sum,pval); ias_stat =c(ias_stat,stat_zero)
} 

# p-value for each covariate
round(ias_sum,3)

