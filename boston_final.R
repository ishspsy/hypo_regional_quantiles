
#### Real Boston data 

data_real_ind = "boston"


# Set "plot_save = 1" for generating coefficient functions
plot_save = 1

# In the estimation setting, must set to be "ttest = 0"
ttest = 0

# quantile intervals (choose among "low", "middle", "high")
taus = "low" 


if (taus == "low"){
  #fix this part for each tau
  tau = seq(0.1,0.3,length.out=11)
  tau0 = tau
  K=length(tau)
}

if (taus == "middle"){
  #fix this part for each tau
  tau = seq(0.4,0.6,length.out=11)
  tau0 = tau
  K=length(tau)
}

if (taus == "high"){
  #fix this part for each tau
  tau = seq(0.7,0.9,length.out=11)
  tau0 = tau
  K=length(tau)
}




# In the testing setting, must set to be "ttest != 0"
if (ttest ==1){
ttype = "tau"   # (independent of both)
type="known"}

if (ttest ==2){
  ttype = "tau"   # (independent of t)
  type="unknown"}


# load Boston data
boston = read.csv("Boston_data.csv", sep= "")  

Y = boston[,14]
U = boston[,13]
Y= Y-mean(Y)

# Choose the varaibles
X =   boston[,c(6,1,10,11,5,7)]  
U = (U - min(U)) / (max(U) - min(U))

n = dim(X)[1] 
p = dim(X)[2]

X = cbind(X, rep(1, n))

n = dim(X)[1]
p = dim(X)[2]

# Normalize data
X[,1:(p-1)] = X[,1:(p-1)] -  matrix(rep(apply(X[,1:(p-1)],2,mean),n), nrow =n, byrow=T) 

yy = sqrt(apply(X[,1:(p-1)]^2,2,mean))

# intercept, sex, smoking, height, height * smoking
X = cbind(X[,1:(p-1)] /  matrix(rep(yy,n), nrow=n,byrow=T), rep(1, n)); 
p = dim(X)[2]
U = (U - min(U)) / (max(U) - min(U))
t_point = U


## Start
# X_tr has an intercept in the last column
X_tr <- X
Y = sqrt(length(Y)) * Y / sqrt(sum(Y^2))
Y_tr <- Y

## Dimension
n = dim(X_tr)[1]
p = dim(X_tr)[2]
pp = p


ind0 <- 1:dim(X_tr)[2]
coef_names = c("Intercept","Rm","Crime","Tax","Ptratio","Nox","Age")
#1   6    % rm average number of rooms per dwelling
#2   1 % per capita crime rate by town (CRIM),
#3   10  % tax	full-value property-tax rate per USD 10,000
#4   11 % ptratio	pupil-teacher ratio by town
#5    5  % nitric oxides concentration(parts per 10 million, NOX),
#6    7  % age	proportion of owner-occupied units built prior to 1940

X_tr <- cbind(X_tr[,p], X_tr[,1:(p-1)])

# Choose variables of interest (In the estimation setting, choose all variables)
if (ttest ==0){ind00_origin = 1:7}
if (ttest == 1){ind00_origin = c(4,7)}  # test1: tax and age
if (ttest == 2){ind00_origin = 2}  # test2:rm

ind00 = ind00_origin

# Change of variables such that the last few variables are considered in the null hypothesis.
X_tr = cbind(X_tr[, -ind00],   X_tr[,ind00])
coef_names_updated = c(coef_names[-ind00], coef_names[ind00])

if (ttest ==0){ind00 = 1:7}
if (ttest == 1){ind00 = c(6,7)}  # test1: tax and age
if (ttest == 2){ind00 = c(7)}  # test2:rm

 
# Estimating for lower quantile region
if (taus == "low"){
data4_real_coef= real_main_function(taus="low", n,p,X_tr, Y_tr,t_point)
save(data4_real_coef, file="data4_real_boston_coef_low.RData")}

# Estimating for middle quantile region
if (taus == "middle"){
data4_real_coef= real_main_function(taus="middle", n,p,X_tr, Y_tr,t_point)
save(data4_real_coef, file="data4_real_boston_coef_middle.RData")}

# Estimating for upper quantile region
if (taus == "high"){
data4_real_coef= real_main_function(taus="high", n,p,X_tr, Y_tr,t_point)
save(data4_real_coef, file="data4_real_boston_coef_high.RData")}

coef_mat_set = data4_real_coef$coef_mat_set

# Generating coefficient functions
library(ggplot2)
library(dplyr)
library(tidyr)


# Arrange estimated quantile coefficient functions.
time_index= seq(0,1,0.05)

coef_spec4_set = list()
for (ind in 1:p){
  coef_spec = coef_mat_set[[ind]]
  colnames(coef_spec) = seq(0,1,length.out=21) 
  rownames(coef_spec) =  tau0 
  coef_spec2 = as.data.frame(coef_spec)
  coef_spec2 = coef_spec2[-c(1,21)]
  coef_spec2$Quantile = tau0 
  coef_spec3 = gather(coef_spec2, "Index", "Values", 1:(ncol(coef_spec2)-1))
  coef_spec4 = spread(coef_spec3, Index, Values)
  coef_spec4=coef_spec4[,-1]
  coef_spec4_set[[ind]] = coef_spec4
}


# plot lines for each coefficient and generating plots
time_indexs= seq(0.05,0.95,0.05)

for (iii in 1:length(coef_names_updated)){
  mname = coef_names_updated[iii]
  
  if (data_real_ind == "boston"){
    pdf(paste("coef_plot_boston_", taus,  "_",  mname, ".pdf", sep=""))}
  if (data_real_ind == "plasma"){
    pdf(paste("coef_plot_plasma_", taus, "_", mname, ".pdf", sep=""))
  }
  
  plot(time_indexs, coef_spec4_set[[iii]][2,],type="l", lty=1, main = mname, ylim =range(coef_spec4_set), xlab  = 'Dietary beta-carotene index (T)',  #"lstat index (T)",       #,
       ylab = expression(paste("Coefficients (", beta, "(", t, "))")),cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
  for (irow in 2:5){
    lines(time_indexs, coef_spec4_set[[iii]][2*irow,], lty=irow, cex=3)
  }
  legend("topright",legend= paste("tau", "=", tau[seq(2,10,2)],sep =" "), 
         lty=c(1:5), lwd=rep(1,5), cex=1.5)
  dev.off()
}

