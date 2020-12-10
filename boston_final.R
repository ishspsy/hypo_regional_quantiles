
#### Real Boston data 

data_real_ind = "boston"

#bic_check = 0
plot_save = 1
ttest = 0
taus = "low" #"high"     # "" # "", "high"

if (ttest ==1){
ttype = "tau"   # (independent of both)
type="known"}

if (ttest ==2){
  ttype = "tau"   # (independent of t)
  type="unknown"}



setwd("/Users/seyoungpark1/Downloads/quant_hypo")
boston = read.csv("Boston_data.csv", sep= "")  

Y = boston[,14];
U = boston[,13]

Y= Y-mean(Y)
#1   6    % rm average number of rooms per dwelling
#2   1 % per capita crime rate by town (CRIM),
#3   10  % tax	full-value property-tax rate per USD 10,000
#4   11 % ptratio	pupil-teacher ratio by town
#5    5  % nitric oxides concentration(parts per 10 million, NOX),
#6    7  % age	proportion of owner-occupied units built prior to 1940

X =   boston[,c(6,1,10,11,5,7)];   
U = (U - min(U)) / (max(U) - min(U));

n = dim(X)[1];
p = dim(X)[2];

X = cbind(X, rep(1, n));

n = dim(X)[1];
p = dim(X)[2];

X[,1:(p-1)] = X[,1:(p-1)] -  matrix(rep(apply(X[,1:(p-1)],2,mean),n), nrow =n, byrow=T)    ;

yy = sqrt(apply(X[,1:(p-1)]^2,2,mean))
# intercept, sex, smoking, height, height * smoking
X = cbind(X[,1:(p-1)] /  matrix(rep(yy,n), nrow=n,byrow=T), rep(1, n)); 
p = dim(X)[2]
U = (U - min(U)) / (max(U) - min(U)); 

dim(X)  # 506   6
length(Y)  # 506
length(U)  # 506

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

X_tr_orig = X_tr  

X_tr <-  X_tr_orig
X_tr_0 = X_tr
ind0 <- 1:dim(X_tr_0)[2]

coef_names = c("Intercept","Rm","Crime","Tax","Ptratio","Nox","Age")

#1   6    % rm average number of rooms per dwelling
#2   1 % per capita crime rate by town (CRIM),
#3   10  % tax	full-value property-tax rate per USD 10,000
#4   11 % ptratio	pupil-teacher ratio by town
#5    5  % nitric oxides concentration(parts per 10 million, NOX),
#6    7  % age	proportion of owner-occupied units built prior to 1940

X_tr <- cbind(X_tr[,p], X_tr[,1:(p-1)])
X_tr_0 <- X_tr

if (ttest == 1){ind00_origin = c(4,7)}  # test1: tax and age
if (ttest == 2){ind00_origin = 2}  # test2:rm
if (ttest ==0){ind00_origin = 1:7}

ind00 = ind00_origin

X_tr = cbind(X_tr[, -ind00],   X_tr[,ind00])
coef_names_updated = c(coef_names[-ind00], coef_names[ind00])

if (ttest ==0){ind00 = 1:7}
if (ttest == 1){ind00 = c(6,7)}  # test1: tax and age
if (ttest == 2){ind00 = c(7)}  # test2:rm



