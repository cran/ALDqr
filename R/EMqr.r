################################################################################
#############       EM algorithm Quantile regression                 ###########
#                       Atualizado em 05/03/13                                 #   
################################################################################
################################################################################
###                    Iniciando o ALgoritmo EM
################################################################################
  
  
  EM.qr<-function(y,x=NULL,perc=NULL, error = 0.000001 ,iter=2000){
    ################################################################################
    ###                    MI Empirica: Veja givens
    ################################################################################
    
    MI_empirica<-function(y,x,perc,theta){
      
      p <- ncol(x)
      n <- nrow(x)
      taup2<-(2/(perc*(1-perc)))
      thep<- (1-2*perc)/(perc*(1-perc))
      
      beta<- theta[1:p]
      sigma<- theta[p+1]
      mu<- x%*%beta
      muc<- y-mu
      
      delta2<-(y-x%*%beta)^2/(taup2*sigma)
      gamma2<-(2+thep^2/taup2)/sigma
      K05P<- 2*besselK(sqrt(delta2*gamma2), 0.5)*(sqrt(delta2/gamma2)^0.5)
      K05N<- 2*besselK(sqrt(delta2*gamma2), -0.5)*(sqrt(delta2/gamma2)^(-0.5))
      K15P<- 2*besselK(sqrt(delta2*gamma2), 1.5)*(sqrt(delta2/gamma2)^(1.5))  
      
      DerG<- matrix(0,nrow=(p+1),ncol=(p+1))
      
      for (i in 1:n){
        dkibeta<- (muc[i]/(taup2*sigma))*(K05N[i])*x[i,]
        dkisigma<- sqrt(delta2[i])/(2*sigma)*K05N[i]+sqrt(gamma2)/(2*sigma)*K15P[i]
        GradBeta<- -thep/(taup2*sigma)*x[i,]+(K05P[i])^(-1)*dkibeta
        Gradsigma<- -1.5/sigma-thep*muc[i]/(taup2*sigma^2)+ (K05P[i])^(-1)*dkisigma      
        GradAux<-as.matrix(c(GradBeta,Gradsigma),p+1,1)
        DerG<-DerG+GradAux%*%t(GradAux)
      }
      
      EP<-sqrt(diag(solve(DerG)))
      
      obj.out <- list(EP = as.vector(EP))
      
      return(obj.out)    
      
    }
    ################################################################################
    
    ################################################################################
    ## Verossimilhanca da RQ: Usando a funcao de perda e usando a Bessel funcion
    ################################################################################
    
    logVQR<-function(y,x,perc,theta){
      p <- ncol(x)
      n <- nrow(x)
      
      beta<-theta[1:p]
      sigma<-theta[p+1]
      mu<-x%*%beta
      muc<- (y-mu)/sigma
      Ind<- (muc<0)+0  
      logver<-sum(-log(sigma)+log(perc*(1-perc))-muc*(perc-Ind))
      return(logver)
    }
    
    p <- ncol(x)
    n <- nrow(x)
    reg <- fastLm(y ~ x[,2:p])
    taup2<- (2/(perc*(1-perc)))
    thep<-(1-2*perc)/(perc*(1-perc))
    #Inicializa beta e sigma2 com os estimadores de m?nimos quadrados
    beta<- as.vector(coefficients(reg),mode="numeric")
    sigma <- sqrt(sum((y-x%*%beta)^2)/(n-p))
    
    teta_velho <- matrix(c(beta,sigma),ncol=1)
    cont <- 0
    criterio<-1
    
    while( criterio> error) {
      
      cont <- (cont+1)
      
      muc<-(y-x%*%beta)
      delta2<-(y-x%*%beta)^2/(taup2*sigma)
      gamma2<-(2+thep^2/taup2)/sigma
      
      vchpN<-besselK(sqrt(delta2*gamma2), 0.5-1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
      vchp1<-besselK(sqrt(delta2*gamma2), 0.5+1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
      
      xM<-c(sqrt(vchpN))*x		
      suma1<-t(xM)%*%(xM)
      suma2<-x*c(vchpN*y-thep)		
      
      sigma<-sum(vchpN*muc^2-2*muc*thep+vchp1*(thep^2+2*taup2))/(3*n*taup2)
      beta<-solve(suma1)%*%apply(suma2,2,sum)
      
      teta_novo<-matrix(c(beta,sigma),ncol=1)
      criterio<-  sqrt(sum((teta_velho-teta_novo)^2)) 
      
      if (cont==iter)
      { 
        break
      }
      
      teta_velho <- teta_novo
      
      
      
    }	
    EP<- MI_empirica(y,x,perc,teta_novo)$EP
    logver<- logVQR(y,x,perc,teta_novo)
    return(list(theta=teta_novo,EP=EP,logver=logver,iter=cont))
    
  }
  
  
