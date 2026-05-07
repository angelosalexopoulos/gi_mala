rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)

ergMean<-function(x){cumsum(x)/(1:length(x))} #ergodic mean for ith dimension
set.seed(121088)

#Give name of dataset
data_name = "German"
d=25
initial_inds=c(d,1:(d-1)) #auxiliary variable for importing the data


####MCMC settings
its =10000
Burnin= 5000


###auxiliary functions
sigmoid = function(x) {exp(x)/(1+exp(x))}

log_sigmoid = function(x){
  
  y = rep(0,length(x))
  idxPos = which(x>0)
  idxNeg = which(x<=0)
  y[idxPos] = -log(1+exp(-x[idxPos]))
  y[idxNeg] = x[idxNeg] - log(1+exp(x[idxNeg]))
  
  return(y)
  
}


#target
log_target=function(z){
  
  n = length(z) 
  Xz = myX%*%z
  log_lik = sum( log_sigmoid(tY*(Xz)) )
  
  grad = t(Y -sigmoid(Xz))%*%myX
  return(list(log_lik = log_lik,grad=grad))
}


myX = as.matrix(read.table(paste("/Users/angelos/Desktop/vr_pCN/data/Logistic_",data_name,"/myX",".txt",sep=""),sep=''))
myX=myX[,initial_inds] #put the columns in the right order
Y = as.matrix(read.table(paste("/Users/angelos/Desktop/vr_pCN/data/Logistic_",data_name,"/Y",".txt",sep=""),sep=''))
tY = 2*Y-1

####obtain MLE for the proposal
fit = glm(Y~myX-1,family = binomial("logit")  )
Sprop = vcov(fit)
R = t(chol(Sprop))

mysamples = array(NA,c(its,d,10))
myESSstats = matrix(NA,3,10)

samp =0 
results = list()
resultsESS = array(NA,c(3,10,2))
likStore2= array(NA,c(10000,10,2))
for(pMALA in c(TRUE,FALSE)){
#  for(pMALA in c(FALSE)){#for variance reduction we do not need to run pMALA
  samp = samp +1
  red_fac= rep(NA,d)
  acc_store = rep(NA,100)
  means = means_cv = matrix(NA,d,100)
  for(rep in 1:10){
    
    yx = grad_store = x =matrix(NA, its,d)
    alpha = logLik=rep(NA,its)
    acceptance_rates = rep(NA,Burnin)
    xcurrent =fit$coefficients
    gamma_tune = 1
    log_post_eval = log_target(xcurrent)
    log_post = log_post_eval[[1]]
    log_grad = log_post_eval[[2]]
    
    if(pMALA){
      target_acc = 0.574
      lim_acc1 = 0.5; lim_acc2 = 0.6
    }else{
      target_acc = 0.774
      lim_acc1 = 0.7; lim_acc2 = 0.8
    }
    ###########run the chain for a burn-in period###############
    ii=0
    acc =0
    for(i in 1:(its+Burnin)){
      

        logp_grad <- log_target(xcurrent)
        
        fx <- logp_grad[[1]]
        grad_x <- logp_grad[[2]]
        
        RT_grad_x <- t(R) %*% t(grad_x)
        A_grad_x <- R %*% RT_grad_x
        
        noise <- rnorm(d, mean=0, sd=1)
        if(pMALA){
          y <- xcurrent + gamma_tune * A_grad_x + sqrt(2 * gamma_tune ) * (R %*% noise)
        }else{
          y <- xcurrent + gamma_tune * A_grad_x + sqrt(2 * gamma_tune- gamma_tune*gamma_tune) * (R %*% noise)
        }

        
        v <- y - xcurrent - (0.5 * gamma_tune) * A_grad_x
        if(pMALA){
          h_y_x <- 0.5 * sum(v * t(grad_x))
        }
        else{
          h_y_x <- (1 / (2 - gamma_tune)) * sum(v * t(grad_x))
        }

        
        logp_grad_y <- log_target(y)
        
        fy <- logp_grad_y[[1]]
        grad_y <- logp_grad_y[[2]]
        
        
        RT_grad_y <- t(R) %*% t(grad_y)
        A_grad_y <- R %*% RT_grad_y
        
        v <- xcurrent - y - (0.5 * gamma_tune) * A_grad_y
        if(pMALA){
          h_x_y <- 0.5 * sum(v * t(grad_y))
        }else{
          h_x_y <- (1 / (2 - gamma_tune)) * sum(v * t(grad_y))
        }
        
        
        delta_f <- fy - fx
        delta_q <- h_x_y - h_y_x
        logp_a <- delta_f + delta_q
        
        u <- runif(1)
        a_prob <- min(exp(logp_a), 1)
        a <- ifelse(u < a_prob, 1, 0)
        #cat(a,'\n')
        xcurrent <- a * y + (1 - a) * xcurrent
        fcurrent <- a * fy + (1-a) * fx
        grad_x_store = a*grad_y + (1-a)*grad_x 
        if(i>Burnin)acc = acc+a
      
      if(TRUE){  
        if (i <= Burnin) {
          acceptance_rates[i] <- a_prob
          
          if (i %% 200 == 0 && i > 200) {
            avg_acc <- mean(acceptance_rates[(i-199):i])
            
            
            if (!((avg_acc >= lim_acc1 && avg_acc <= lim_acc2))) {
              if (avg_acc < target_acc) {
                gamma_tune <- gamma_tune / (1 + 0.1 * (target_acc - avg_acc) / target_acc)
              } else if (avg_acc > target_acc) {
                gamma_tune <- gamma_tune * (1 + 0.1 * (avg_acc - target_acc) / target_acc)
              }
            }
          }
        }
      } 
      
      
      if(i>Burnin){
        ii =ii+1
        x[ii,] = xcurrent
        yx[ii,] = y-x[ii,]
        logLik[ii] = fcurrent
        alpha[ii] = a_prob
        grad_store[ii,] = grad_x_store
      }
      #if(i%%10000==0)cat(i,'\n')
    }
    
    acc_store[rep] =acc/its
    gradS = Sprop%*%t(grad_store)
    
    for(dd in 1:d){
      G = x[,dd]/gamma_tune
      thetaCV = matrix(NA,2,its)
    if(FALSE){     
      #b2 = sum((alpha*(yx[,dd]/gamma_tune))*(yx[,dd]/gamma_tune-gradS[dd,]))/sum((yx[,dd]/gamma_tune-gradS[dd,])^2)
      b2 = sum((alpha*(yx[,dd]/gamma_tune))*(yx[,dd]/gamma_tune-mean(yx[,dd]/gamma_tune)))/sum((yx[,dd]/gamma_tune-mean(yx[,dd]/gamma_tune))^2)
      
      #b2 = median(b2)
      ylm = alpha*yx[,dd]/gamma_tune #- gradS[dd,]
      xlm = yx[,dd]/gamma_tune #- gradS[dd,]
      #b2 = cov( ylm,xlm  )/var(xlm)
      #b2=-alpha
      b2 = 1
      #b2 = mean(ylm^2)/mean(xlm^2)
      PGhat = G + alpha*(yx[,dd])/gamma_tune -b2*(yx[,dd])/gamma_tune+b2*gradS[dd,]
      
      theta1 = ergMean(x[-1,dd]*(G[-1]+PGhat[-1]))-ergMean(x[-1,dd])*ergMean(G[-1]+PGhat[-1])
      theta2 = ergMean((x[2:its,dd]/gamma_tune - PGhat[1:(its-1)])^2)
      theta = theta1/theta2
      means[dd,rep]= mean(x[,dd]) 
      
      means_cv[dd,rep]=mymu = mean(x[-1,dd]) - theta[its-1]*mean(G -PGhat)
    }
      if(FALSE){
        #try here the second class of estimators proposed by Michalis
       g1 = alpha*(yx[,dd]/gamma_tune)
       g2 = yx[,dd]/gamma_tune - gradS[dd,]
       #yx[,dd]/gamma_tune
       K11 = ergMean(g1^2)
       K12 = K21 = ergMean(g1*g2)
       K22 = ergMean(g2^2)
       
       th1 = ergMean(x[,dd]*g1)-ergMean(x[,dd])*ergMean(g1)
       th2 = ergMean(x[,dd]*g2)-ergMean(x[,dd])*ergMean(g2)
       
       for(i in 20:its){
         Kinv = solve( rbind(c(K11[i],K12[i]),c(K21[i],K22[i])) + 10^(-6)   )
         thetaCV[,i] = Kinv%*%c(th1[i],th2[i])
       }
       means[dd,rep]= mean(x[,dd]) 
       
      means_cv[dd,rep]= mean(x[,dd]) -mean(g1)- mean(g2)
     #- thetaCV[1,its]*
      #- thetaCV[2,its]*
       
      } 
      if(FALSE){
        PGhat = G + alpha*(yx[,dd])/gamma_tune 
        
        theta1 = ergMean(x[-1,dd]*(G[-1]+PGhat[-1]))-ergMean(x[-1,dd])*ergMean(G[-1]+PGhat[-1])
        theta2 = ergMean((x[2:its,dd]/gamma_tune - PGhat[1:(its-1)])^2)
        theta = theta1/theta2
        means[dd,rep]= mean(x[,dd]) 
        
        means_cv[dd,rep]=mymu = mean(x[-1,dd]) - theta[its-1]*mean(G -PGhat)
      }
      if(FALSE){
        PGhat =x[,dd]/gamma_tune  + gradS[dd,]
        
        theta1 = ergMean(x[-1,dd]*(G[-1]+PGhat[-1]))-ergMean(x[-1,dd])*ergMean(G[-1]+PGhat[-1])
        theta2 = ergMean((x[2:its,dd]/gamma_tune - PGhat[1:(its-1)])^2)
        theta = theta1/theta2
        means[dd,rep]= mean(x[,dd]) 
        
        means_cv[dd,rep]=mymu = mean(x[-1,dd]) - theta[its-1]*mean(G -PGhat)
      }
      
      
     # cat(b2,'\n')
    }
    ############uncomment when need ESS comparisons ############
    mysamples[,,rep] = x
    myESS = effectiveSize(x)
    myESSstats[,rep] = c(min(myESS),median(myESS),max(myESS)) 
    cat(c(gamma_tune,samp,rep,avg_acc),'\n')
    likStore2[,rep,samp] = logLik
    
     cat(rep,'\n')
  }
  for(dd in 1:d)red_fac[dd]= var(means[dd,])/var(means_cv[dd,])
  
  results[[samp]] = mysamples
  resultsESS[,,samp]  =myESSstats 
  
  ################################################
}
############uncomment when need ESS comparisons ############
ESSmat2 = matrix(NA,2,3)
rownames(ESSmat2) = c("pMALA","GaussInv")
colnames(ESSmat2) = c("min","median","max")

ESSmat2[1,] = apply(resultsESS[,,1],1,mean)
ESSmat2[2,] = apply(resultsESS[,,2],1,mean)
############################################################

if(FALSE){
pdf('/Users/angelos/Desktop/LogisticToyACF.pdf')
par(mfrow=c(1,2))
acf_mala = acf(likStore[,1,1],plot=F)
acf_angelos = acf(likStore[,1,2],plot=F)
acf_mala_mean =  acf_mala$acf
acf_angelos_mean = acf_angelos$acf
for(i in 2:10){
  acf_mala = acf(likStore[,i,1],plot=F)
  acf_angelos = acf(likStore[,i,2],plot=F)
  acf_mala_mean = acf_mala_mean + acf_mala$acf
  acf_angelos_mean = acf_angelos_mean + acf_angelos$acf
}
plot(acf_mala_mean,type='l',xlab = "Lag",ylab = "Autocorrelation",
     cex.axis = 1.5,              # Change size of axis tick labels
     cex.lab = 1.5,lwd=2)
lines(acf_angelos_mean,col=2,lty=2,lwd=2)
legend("topright",c('GaussInv','pMALA'),col=c(2,1),lty=c(2,1),lwd=c(2,2))
acf_mala = acf(likStore2[,1,1],plot=F)
acf_angelos = acf(likStore2[,1,2],plot=F)
acf_mala_mean =  acf_mala$acf
acf_angelos_mean = acf_angelos$acf
for(i in 2:10){
  acf_mala = acf(likStore2[,i,1],plot=F)
  acf_angelos = acf(likStore2[,i,2],plot=F)
  acf_mala_mean = acf_mala_mean + acf_mala$acf
  acf_angelos_mean = acf_angelos_mean + acf_angelos$acf
}
plot(acf_mala_mean,type='l',xlab = "Lag",ylab = "Autocorrelation",
     cex.axis = 1.5,              # Change size of axis tick labels
     cex.lab = 1.5,lwd=2)
lines(acf_angelos_mean,col=2,lty=2,lwd=2)
legend("topright",c('GaussInv','pMALA'),col=c(2,1),lty=c(2,1),lwd=c(2,2))
dev.off()
}

