#' @export
#' @importFrom rstan stan
#' @importFrom  stats lm


GPRMortalityChild <- function( data.mortality , data.mean ,minYear,maxYear,
                               nu,rho_,n.itr=4000,n.warm=3000,verbose=FALSE){


  ################## err & warn #####################

  if( !is.data.frame(data.mortality) ) stop("data.mortality should be data.frame.")
  if( !is.data.frame(data.mean)) stop("data.mean should be data.frame.")

  if( sum(c("location","year","pop","mortality","isDR","type")%in%colnames(data.mortality))!=6  ) stop("data.mortality should contain  these variable: location, year, pop, mortality, isDR and type ")
  if( sum(c("location","year","mean")%in%colnames(data.mean))!=3  ) stop("data.mean should contain  these variable: location, year and mean ")

  data.mortality <- data.mortality[,  colnames(data.mortality)%in%c("location","year","pop","mortality","isDR","type")  ]
  data.mean <- data.mean[,  colnames(data.mean)%in%c("location","year","mean")  ]

  if(  !is.numeric(data.mortality$year) )stop("year variable in data.mortality should be numeric")
  if(  !is.numeric(data.mortality$location) )stop("location variable in data.mortality should be numeric")
  if(  !is.numeric(data.mortality$mortality) )stop("mortality variable in data.mortality should be numeric")
  if(  !is.numeric(data.mortality$isDR) )stop("isDR variable in data.mortality should be numeric")
  if(  !is.numeric(data.mortality$pop) )stop("pop variable in data.mortality should be numeric")
  if(  !is.character(data.mortality$type) )stop("type variable in data.mortality should be character")

  if( sum(  apply(data.mean,2,function(x) is.numeric(x))  ) !=3 ) stop("all variables in data.mean should be numeric")


  data.mortality$Nt <- data.mortality$pop
  data.mortality$province <- data.mortality$location
  data.mean$province <- data.mean$location

  if( dim(na.omit(data.mean))[1] != dim(data.mean)[1]  )  warning("some NA Obs removed in data.mean ")
  if( dim(na.omit(data.mortality))[1] != dim(data.mortality)[1]  )  warning("some NA Obs removed in data.mortality  ")

  data.mortality<-na.omit(data.mortality)
  data.mean<-na.omit(data.mean)

  if ( !all(data.mortality$year%in%data.mean$year) ) stop("all years in data mortality should exist in data.mean")
  if ( !all( c(minYear:maxYear)%in%data.mean$year  ) ) stop("minYear and maxYear do not match with years in data.mean")


  if( minYear>=maxYear) stop("minYear should less than maxYear")
  if ( nu>2 | nu<0.5) stop("nu should greater than 0.5 and less than two")
  if (n.warm>=n.itr) stop("n.warm should less than n.itr")
  if( sum(data.mortality$mortality<=0 | data.mortality$mortality>1  )!=0 ) stop("all mortality rate should bigger than zero and less than 1")
  if( sum(data.mortality$pop<=0  )!=0 ) stop(" population should bigger than zero ")
  if( sum(data.mortality$year<0   )!=0 ) stop(" year should bigger than zero ")
  if( sum(data.mean$mean<0 | data.mean$mean>1 )!=0 ) stop("all mean rate should bigger than zero and less than one")
  if( rho_<0 | rho_>1) stop("rho should less than 1 and bigger than zero")

  data.mortality <- data.mortality[order(data.mortality$year,data.mortality$province),]
  data.mean <- data.mean[order(data.mean$year,data.mean$province),]



  provinceID_unique  <- unique(data.mortality$province  )


  if( !all(data.mortality$year%in%c(minYear:maxYear)) ){
    data.mortality <- data.mortality[data.mortality$year%in%c(minYear:maxYear) , ]
    warning("only the years between minYear and maxyear of the time range are considered. other time mortality points are dropped.")
  }

  if( !all(data.mean$year%in%c(minYear:maxYear)) ){
    data.mean <- data.mean[data.mean$year%in%c(minYear:maxYear) , ]
    warning("only the years between minYear and maxyear of the time range are considered. other time mean points are dropped.")
  }

  i=1
  myT=maxYear-minYear+1
      for(i in provinceID_unique){
        if( nrow(data.mortality[data.mortality$province==i,])<=1 | nrow(data.mean[data.mean$province==i,])<myT ){
          stop("at least two mortality rates should exist for each location
               OR
               time span in data.mean does not match with minYear and maxyear  ")
        }
      }

 if( any(tapply(data.mortality$mortality,data.mortality$type,length)<=1)  ) stop("each data source must contain at least two data points. ")


 ########################################################## data ###############



rho_ = rho_*100

data.mortality <- data.mortality[order(data.mortality$year,data.mortality$province),]
data.mean <- data.mean[order(data.mean$year,data.mean$province),]

data.mortality$q = data.mortality$mortality
data.mortality$y = log10(data.mortality$q)

data.mean$M = log10(data.mean$mean)
  ########################################################### function ############################
  calcLambdas = function (MeanData,MortalityData)
  {
    MeanData$UID = paste0(MeanData$year,"-",MeanData$province)
    MortalityData$UID = paste0(MortalityData$year,"-",MortalityData$province)
    Mt = MeanData$M[match(MortalityData$UID,MeanData$UID)]
    ProvinceCodes = unique(MortalityData$province)
    rLambda=rep(0,length(ProvinceCodes))
    for (i in 1:length(ProvinceCodes))
    {
      S2 = (Mt[MortalityData$province==ProvinceCodes[i]] - MortalityData$y[MortalityData$province==ProvinceCodes[i]])^2
      rLambda[i]=sqrt(mean(S2))
    }
    rMat = matrix(c(ProvinceCodes,rLambda),ncol=2,byrow=F)
    colnames(rMat)=c("Province","Lambda")
    return (rMat)
  }

  calcBetaDRI = function(MortalityData)
  {
    rMat = matrix(0,ncol=4,nrow=length(unique(MortalityData$province)))
    colnames(rMat) = c("Province","Sig","BMean","BStdDev")
    rMat[,1]=unique(MortalityData$province)
    for (i in 1:nrow (rMat))
    {
      pCon = MortalityData$province==rMat[i,1] ### remove & MortalityData$type!="census90"
      DataSubset = MortalityData[pCon,]
      DataSubset$isDR = 0                      ##### need really ??
      DataSubset$isDR[DataSubset$type=="DR"]=1 #######
      if (sum(DataSubset$isDR)>0){
        rModel = lm(y ~ year + isDR,DataSubset)
        if (summary(rModel)$coefficients[3,4]<0.05)
        {
          if (rModel$coefficients[3]>0)
          {
            cat("Error! Death registration data for province #",rMat[i,1]," has negative incompleteness!")
          }else{
            rMat[i,2]=1
            rMat[i,3]=rModel$coefficients[3]
            rMat[i,4]=summary(rModel)$coefficients[3,2]
          }
        }
      }
    }
    return (rMat)

  }

  calcParameters = function (MeanData,MortalityData)
  {
    rPars = list()
    for (i in unique(MortalityData$province))
    {
      DataSubset = MortalityData[MortalityData$province==i,]
      provinceID = i
      Sources = unique(DataSubset$type)
      S=length(Sources)
      Mx = rep(0,S)
      Cx1 = rep(0,S)
      Cx2 = rep(0,S)
      for (j in 1:S)
      {
        Mx[j]=length(DataSubset$type[DataSubset$type==Sources[j]])
        Cx1[j]=floor((Mx[j]+1)/2)
        Cx2[j]=ceiling((Mx[j])/2)
      }
      Y=matrix(NA,nrow=S,ncol=max(Mx))
      YEAR=matrix(NA,nrow=S,ncol=max(Mx))
      Y2=matrix(NA,nrow=S,ncol=max(Mx))
      Nt=matrix(NA,nrow=S,ncol=max(Mx))
      ISDRi = rep(NA,S)
      for (j in 1:S)
      {
        Y[j,1:Mx[j]]=DataSubset$y[DataSubset$type==Sources[j]]
        YEAR[j,1:Mx[j]]=DataSubset$year[DataSubset$type==Sources[j]]-minYear+1
        Y2[j,1:Mx[j]]=DataSubset$q[DataSubset$type==Sources[j]]
        Nt[j,1:Mx[j]]=DataSubset$Nt[DataSubset$type==Sources[j]]
        ISDRi[j]=DataSubset$isDR[DataSubset$type==Sources[j]][1]
      }
      SIGMAy = ((1 / (Y2 * log(10)))^2) * Y2 * (1-Y2) / Nt
      Y[is.na(Y[])]=0
      YEAR[is.na(YEAR[])]=0
      subTau = ISDRi * BDRi[length(rPars)+1,4] * BDRi[length(rPars)+1,4]  + apply(SIGMAy,1,mean,na.rm=T)
      SIGMAy[is.na(SIGMAy[])]=0
      rPars[[length(rPars)+1]]=list(provinceID=provinceID,Sources=Sources,S=S,Mx=Mx,Cx1=Cx1,Cx2=Cx2,Y=Y,YEAR=YEAR,Y2=Y2,Nt=Nt,ISvri=ISDRi,SIGMAy=SIGMAy,M=MeanData$M[MeanData$province==i],subTau=subTau)
    }
    return (rPars)
  }
  calcRho = function(MeanData,MortalityData)
  {
    return (rep(rho_,length(unique(MortalityData$province))))
  }




  RunModel = function(provinceID,myT,nu,lambda,rho,BDRi,Pars,n.itr,n.warm) # remove computerID
  {
    calcSigmaFunction = function(rho,lambda,nu,size){
      Cmat = diag(size)*(lambda^2)
      for (i in 1:(size-1))
        for (j in (i+1):size)
        {
          Cmat[i,j]=	(lambda^2) *
            ((2 ^ (1-nu)) / gamma (nu)) *
            (((j-i) * (2*sqrt(nu)) / rho) ^ nu) *
            besselK((j-i) * (2*sqrt(nu)) / rho,nu)
          Cmat[j,i]=Cmat[i,j]
        }
      return (Cmat)
    }
    k = Pars[[provinceID]]
    datalist = list(Cx1=k$Cx1,Cx2=k$Cx2,myT=myT,nu=nu,lambda=lambda[provinceID,2],YEAR=k$YEAR,ISvri=k$ISvri,SIGMAy=k$SIGMAy,S=k$S,Mx=k$Mx,Y=k$Y,maxDS=max(k$Mx),subTau=k$subTau)

    #initval
    if (BDRi[provinceID,2]==1)
    {
      scode <- "
      data{
	int myT; //Number of times
      vector[myT] M; //Mean Function
      real lambda; //lambda of the covariance function
      real nu; //nu of the covariance function
      int S; //Number of sources
      int Mx[S]; //Number of data for each source
      real ISvri[S];
      real MUvri;
      real<lower=0> SIGMAvri;
      int Cx1[S];
      int Cx2[S];
      real subTau[S];
      int maxDS; //Max number of data in a source
      matrix[S,maxDS] Y;
      int YEAR[S,maxDS];
      matrix[S,maxDS] SIGMAy;
      cov_matrix[myT] C;
    }
      parameters{
      vector[myT] tmpf;
      vector[S] tmpgamma;
      real Bvri;
      }
      transformed parameters {
      matrix[myT,myT] L;
      real d;
      real fpart3;
      real MEDx;
      real MAD;
      real<lower=0> tmptau;
      vector[S] tau;

      vector[S] omega;
      vector[myT] f;

      L= cholesky_decompose(C);
      f = M + L * tmpf;

      //calculating taus
      for (i in 1:S)
      {
      vector[Mx[i]] A;
      vector[Mx[i]] B;
      vector[Mx[i]] sV;
      vector[Mx[i]] T;
      for (j in 1:Mx[i])
      {
      A[j] = (Y[i,j] - f[YEAR[i,j]] );
      T[j] =  f[YEAR[i,j]];
      }
      sV = sort_asc(A);
      MEDx = (sV[Cx1[i]] + sV[Cx2[i]])/2;
      for (j in 1:Mx[i])
      B[j] = fabs(A[j] - MEDx);
      sV = sort_asc(B);
      MAD = (sV[Cx1[i]] + sV[Cx2[i]])/2;
      tmptau = MAD*MAD - subTau[i];
      tau[i] = 1/tmptau;
      }

      for (i in 1:S)
      omega[i] = 1/tmpgamma[i];
      }
      model{
      tmpf ~ normal(0,1);

      for (i in 1:S)
      tmpgamma[i] ~ gamma(0.1 * tau[i],10);
      for (i in 1:S)
      for (j in 1:Mx[i])
      Y[i,j] ~ normal( f[YEAR[i,j]] + ISvri[i] * Bvri, sqrt(SIGMAy[i,j] + omega[i] * omega[i]));
      Bvri ~ normal(MUvri,SIGMAvri);
      }
      "
      # writeLines(scode,paste0(getwd(),"stanmodel/mymodel.stan"))
      datalist$SIGMAvri=BDRi[provinceID,4]
      datalist$MUvri=BDRi[provinceID,3]
      inlist = list(tmpf=rep(0,length(k$M)),tmpgamma=rep(100,k$S),Bvri=BDRi[provinceID,3])
    }else{
      scode <- "
      data{
	int myT; //Number of times
      vector[myT] M; //Mean Function
      real lambda; //lambda of the covariance function
      real nu; //nu of the covariance function
      int S; //Number of sources
      int Mx[S]; //Number of data for each source
      int Cx1[S];
      int Cx2[S];
      real subTau[S];
      int maxDS; //Max number of data in a source
      matrix[S,maxDS] Y;
      int YEAR[S,maxDS];
      matrix[S,maxDS] SIGMAy;
      cov_matrix[myT] C;
    }
      parameters{
      vector[myT] tmpf;
      vector[S] tmpgamma;
      }
      transformed parameters {
      matrix[myT,myT] L;
      real d;
      real fpart3;
      real MEDx;
      real MAD;
      real tmptau;
      vector[S] tau;

      vector[S] omega;
      vector[myT] f;

      L<- cholesky_decompose(C);
      f <- M + L * tmpf;

      //calculating taus
      for (i in 1:S)
      {
      vector[Mx[i]] A;
      vector[Mx[i]] B;
      vector[Mx[i]] sV;
      vector[Mx[i]] T;
      for (j in 1:Mx[i])
      {
      A[j] <- (Y[i,j] - f[YEAR[i,j]] );
      T[j] <-  f[YEAR[i,j]];
      }
      sV <- sort_asc(A);
      MEDx <- (sV[Cx1[i]] + sV[Cx2[i]])/2;
      for (j in 1:Mx[i])
      B[j] <- fabs(A[j] - MEDx);
      sV <- sort_asc(B);
      MAD <- (sV[Cx1[i]] + sV[Cx2[i]])/2;
      tmptau <- MAD*MAD - subTau[i];
      tau[i] <- 1/tmptau;
      }

      for (i in 1:S)
      omega[i] <- 1/tmpgamma[i];
      }
      model{
      tmpf ~ normal(0,1);

      for (i in 1:S)
      tmpgamma[i] ~ gamma(0.1 * tau[i],10);
      for (i in 1:S)
      for (j in 1:Mx[i])
      Y[i,j] ~ normal( f[YEAR[i,j]], sqrt(SIGMAy[i,j] + omega[i] * omega[i]));
      }
      "
      # writeLines(scode, paste0(getwd(),"stanmodel/mymodel.stan"))
      inlist = list(tmpf=rep(0,length(k$M)),tmpgamma=rep(100,k$S))
    }

    datalist$C = calcSigmaFunction(rho[provinceID],lambda[provinceID,2],nu,myT)
    datalist$M = k$M
    mcmcrun=stan(model_code=scode,data=datalist,init=list(inlist), #"mymodel.stan"
                 chains=1, cores= 1,
                 iter=n.itr,pars=c("f","tau"),warmup=n.warm)
    l=as.matrix(mcmcrun)[,1:(myT+k$S)]   # ?? 54
    # colnames(l)=c(paste0("f",1:myT),paste0("tau",1:k$S))
     colnames(l)=c(paste0("f",1:myT),paste0("tau",1:k$S))

    l=cbind(l,it=1:(n.itr-n.warm))
    out_sim = l[,1:myT]
    out_var = l[,(myT+1):(ncol(l)-1)]
    out = list(simulation=out_sim,variance = out_var)
    dim(out$simulation)
    dim(out$variance)
    out
  }

  ##################################################### calc ############################################
  #calculating parameters
  myT=maxYear-minYear+1
  #nu=2
  lambda=calcLambdas(data.mean,data.mortality)
  rho=calcRho(data.mean,data.mortality) ## should equall all????????
  BDRi = calcBetaDRI(data.mortality)
  Pars = calcParameters(data.mean,data.mortality)
  U = unique(data.mortality$province)

  OUT_sim = matrix(NA,0,myT+2)
  OUT_var = list()
  for(i in 1:length(U)){
      t = RunModel(i,myT,nu,lambda,rho,BDRi,Pars,n.itr,n.warm) ### 11 NOt province 0 1 2 ...
      u <- U[i]
      out_sim = cbind(u,1:c(n.itr - n.warm  ), 10^t$simulation)
      colnames(out_sim) = c("location","iteration",minYear:maxYear)
      OUT_sim = rbind(OUT_sim , out_sim )
      t2 = cbind( 1:(n.itr-n.warm) ,t$variance)
      colnames(t2) <- c("iteration", unique(data.mortality$type[data.mortality$province == u]) )
      OUT_var[[as.character(u)]] <- 1/t2
  }
  OUT = list(simulation = OUT_sim,variance = OUT_var)
  return(OUT)

}
