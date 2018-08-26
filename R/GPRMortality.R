#' @export
#' @importFrom stats na.omit
#' @importFrom utils head
########################################### FUN &  error & warning ##############
 GPRMortality <- function(data.mortality,data.mean,minYear,maxYear,
                  nu,rho_,product , n.itr=4000,n.warm=3000,verbose=FALSE){

   ####################### err & warning #################
   if( !is.data.frame(data.mortality) ) stop("data.mortality should be data.frame.")
   if( !is.data.frame(data.mean)) stop("data.mean should be data.frame.")

   if( sum(c("sex","age_cat","location","year","pop","mortality","completeness")%in%colnames(data.mortality))!=7  ) stop("data.mortality should contain  these variable: sex, age_cat, location, year, pop, mortality and completeness ")
   if( sum(c("sex","age_cat","location","year","mean")%in%colnames(data.mean))!=5  ) stop("data.mean should contain  these variable: sex, age_cat, location, year and mean ")

   data.mortality <- data.mortality[,  colnames(data.mortality)%in%c("sex","age_cat","location","year","pop","mortality","completeness")  ]
   data.mean <- data.mean[,  colnames(data.mean)%in%c("sex","age_cat","location","year","mean")  ]

   if( sum(  apply(data.mortality,2,function(x) is.numeric(x))  ) !=7 ) stop("all variables in data.mortality should be numeric")
   if( sum(  apply(data.mean,2,function(x) is.numeric(x))  ) !=5 ) stop("all variables in data.mean should be numeric")


 data.mortality$province = data.mortality$location
 data.mean$province = data.mean$location
 myT=maxYear-minYear+1

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
if(product<=0 | product>=1) stop("product should less than one and bigger than zero")

 data.mortality <- data.mortality[order(data.mortality$year,data.mortality$province,data.mortality$sex,data.mortality$age_cat),]
 data.mean <- data.mean[order(data.mean$year,data.mean$province,data.mean$sex,data.mean$age_cat),]


  ageID_unique <- unique(data.mortality$age_cat)
  computerID_unique <- unique(data.mortality$sex)
  provinceID_unique  <- unique(data.mortality$province  )


  if( !all(data.mortality$year%in%c(minYear:maxYear)) ){
  data.mortality <- data.mortality[data.mortality$year%in%c(minYear:maxYear) , ]
  warning("only the years between minYear and maxyear of the time range are considered. other time mortality points are dropped.")
  }

  if( !all(data.mean$year%in%c(minYear:maxYear)) ){
  data.mean <- data.mean[data.mean$year%in%c(minYear:maxYear) , ]
  warning("only the years between minYear and maxyear of the time range are considered. other time mean points are dropped.")
  }

i=1;ii=0;iii=0
for(i in ageID_unique){
for(ii in computerID_unique){
for(iii in provinceID_unique){
    if( nrow(data.mortality[data.mortality$province==iii & data.mortality$sex==ii & data.mortality$age_cat==i,])<=1 | nrow(data.mean[data.mean$province==iii & data.mean$sex==ii & data.mean$age_cat==i,])<myT ){
    stop("at least two mortality rates should exist for each age-sex-location combinatin
         OR
         time span in data.mean does not match with minYear and maxyear  ")
    }
}
}
}

	######################## prepare data ###################

rho_ = rho_*100
data.mortality$y = log10(data.mortality$mortality)
data.mean$M = log10(data.mean$mean)
data.mortality$DDM <-  data.mortality$completeness*product

	######################### FUN ##################
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

calcParameters = function (MeanData,MortalityData,IncompletnessData)
{
	rPars = list()
	for (i in unique(MortalityData$province))
	{
		DataSubset = MortalityData[MortalityData$province==i,]
		DataSubset_d = IncompletnessData[IncompletnessData$province==i,]
		provinceID = i



    		 Mx=length(DataSubset$province[DataSubset$province==i])
    		 Nx=length(DataSubset_d$province[DataSubset_d$province==i])

if(Nx>Mx)
{
DataSubset_d = na.omit(DataSubset_d)
}


  		Y=rep(NA,Mx)
  		YEAR=rep(NA,Mx)
  		Y2=rep(NA,Mx)
		Nt=rep(NA,Mx)
		ICDd=rep(NA,Mx)


    			Y=DataSubset$y
    			YEAR=(DataSubset$year)-minYear+1
    			Y2=DataSubset$mortality
    			Nt=DataSubset$pop
			ICDd=(DataSubset_d$DDM)^2


  		SIGMAy = ((1 / (Y2 * log(10)))^2) * Y2 * (1-Y2) / Nt
		SIGMAy = SIGMAy + ICDd
		Y[is.na(Y[])]=0
		YEAR[is.na(YEAR[])]=0
		SIGMAy[is.na(SIGMAy[])]=0
		rPars[[length(rPars)+1]]=list(provinceID=provinceID,Mx=Mx,Nx=Nx,Y=Y,YEAR=YEAR,Y2=Y2,Nt=Nt,ICDd=ICDd,SIGMAy=SIGMAy,M=MeanData$M[MeanData$province==i])
	}
	return (rPars)
}


calcRho = function(MeanData,MortalityData)
{
	return (rep(rho_,length(unique(MortalityData$province))))
}

RunModel = function(provinceID,myT,nu,lambda,rho,Pars,computerID,ageID,n.itr,n.warm)
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
	#library(rstan)
	#library(foreign)
	#setwd(path)


#select data
	k = Pars[[provinceID]]
	datalist = list(myT=myT,nu=nu,lambda=lambda[provinceID,2],YEAR=k$YEAR,SIGMAy=k$SIGMAy,Mx=k$Mx,Y=k$Y)

	#initval

		#fName = "stancode.txt"
		scode<-"
		data{
	int myT; //Number of times
	vector[myT] M; //Mean Function
	real lambda; //lambda of the covariance function
	real nu; //nu of the covariance function
	int Mx; //Number of data for each source
	vector[Mx] Y;
	int YEAR[Mx];
	vector[Mx] SIGMAy;
	cov_matrix[myT] C;
}

parameters{
	vector[myT] tmpf;
}

transformed parameters {
	matrix[myT,myT] L;
	real d;
	real fpart3;
	vector[myT] f;

	L =  cholesky_decompose(C);
	f = M + L * tmpf;

	}


model{
	tmpf ~ normal(0,1);

		for (i in 1:Mx)
			Y[i] ~ normal( f[YEAR[i]], sqrt(SIGMAy[i]));
}
		"

		# writeLines(scode, paste0(getwd(),"stanmodel/mymodel.stan"))
		inlist = list(tmpf=rep(0,length(k$M)))


	datalist$C = calcSigmaFunction(rho[provinceID],lambda[provinceID,2],nu,myT)
	datalist$M = k$M
	mcmcrun=stan(model_code=scode,data=datalist,init=list(inlist),chains=1,cores = 1,iter=n.itr,pars=c("f"),warmup=n.warm,verbose = verbose)
	l=as.matrix(mcmcrun)[,1:myT]
	l = 10^l
	l=cbind(it=1:(n.itr-n.warm),provinceID_unique[provinceID],computerID,ageID,l)
	colnames(l)=c( "iteration","location","sex","age_cat" ,paste0(minYear:maxYear))
  return(l)
	#write.dta(as.data.frame(l),paste0("results/sims/res_",computerID,"_",ageID,"_",k$provinceID,".dta"))

}

################################## calc ######################
 # computerID = 0 ; ageID = 3 ; provinceID = 1 ; k=1
myT=maxYear-minYear+1
out = matrix(NA, 0 , myT + 4 )

for(ageID in ageID_unique ){
for(computerID in computerID_unique ){
      data.mortality_sub = data.mortality[data.mortality$sex==computerID & data.mortality$age_cat==ageID ,]	#choose sex & age
      data.mean_sub = data.mean[data.mean$sex==computerID & data.mean$age_cat==ageID,]	#choose sex & age
      # ICD_sub = ICD[ICD$sex==computerID & ICD$age_cat==ageID ,]	#choose sex & age
      ICD_sub <- data.mortality_sub
      lambda=calcLambdas(data.mean_sub,data.mortality_sub)
      rho=calcRho(data.mean_sub,data.mortality_sub)
      Pars = calcParameters(data.mean_sub,data.mortality_sub, ICD_sub)

      for (k in  (1:length(provinceID_unique))   ){
          out <- rbind( out,
      	                RunModel(k,myT,nu,lambda,rho,Pars,computerID,ageID,n.itr,n.warm)
      			           )
      }

}
}
return(out)
 }
