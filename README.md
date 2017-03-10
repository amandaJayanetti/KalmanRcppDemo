# KalmanRcppDemo - An implementation of Kalman Filter with Rcpp

##### Installing KalmanRcppDemo

```R

install.packages("devtools")
library(devtools)
install_github("amandaJayanetti/KalmanRcppDemo")

```

###### Example

```R
# preparing arguments for kalmanUpdate function

# Best estimate of measurements/state at times k,k+dt,k+2dt,..
Y<-matrix(c(0.39,0.50,0.48,0.29,0.25,0.32,0.34,0.48,0.41,0.45),10,1,byrow=TRUE)
# Initial state
X<-matrix(c(0.39,0,7),3,1,byrow=TRUE)
# Prediction matrix used to obtain next state
F<-matrix(c(1,0.05,0,0,1,0.05,0,0,1),3,3,byrow = TRUE)
# Matrix to model the signal
H<-matrix(c(1,1,0),1,3,byrow = TRUE)
# Covariance matrix to model uncertainity associated with untracked influences
Q<-matrix(c(0.0025, 0.0025, 0, 0, .0025, 0, 0, 0.0025, 0),3,3,byrow=TRUE)
# Covariance matrix to model uncertainity eg. sensor noise
R<-matrix(c(3),1,1,byrow=TRUE)
# Error covariance matrix of the state estimate 
P<-matrix(c(0.05,0.05,5,0.5,10000,10,0.5,10,100),3,3,byrow = TRUE)

kalmanUpdate(Y,X,F,H,Q,R,P)

```


