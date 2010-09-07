#include           <R.h>
#include           <Rmath.h>  
#include           <Rdefines.h>
#include           <stdio.h>
#include           <math.h>
#include           <stdlib.h>
#include 		   <R_ext/Utils.h>


int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}


double ran_inv_gamma(double a, double b){
  return(1/rgamma(a, 1/b));
}

double ran_trunc_norm_upper(double mu, double sigma, double a){
  double u = runif(0,1); 
  double v = pnorm((-a+mu)/sigma,0,1,1,0);
  return(-sigma * qnorm(u*v,0,1,1,0)+mu);
}

double ran_trunc_norm_lower(double mu, double sigma, double a){
  double u = runif(0,1);
  double v = pnorm((a-mu)/sigma,0,1,1,0);
  return(sigma * qnorm(u*v,0,1,1,0)+mu);
}

double ran_trunc_norm_both(double mu, double sigma, double a, double b){
  double low = pnorm((a-mu)/sigma,0,1,1,0);
  double high = pnorm((b-mu)/sigma,0,1,1,0);
  double u = runif(low,high);
  return(sigma * qnorm(u,0,1,1,0)+mu);
}

double expXPhiA(double x,double a){
  double ret=0;
  switch(a<-5){
  case 1:
    ret=exp(x-pow(a,2)/2-log(-a)-.94/pow(a,2)-.5*log(2*M_PI));
    break;
  case 0:
     Rprintf("\n****Tried to use approximation inappropriately. x=%f, a=%f\n",x,a);
     exit(1);
    break;
  } 
  return(ret);
}

double logPhiAminusPhiB(double a, double b){
  int i;
  double c0=.2316419;
  double c[5]={.319381530,-.356563782,1.781477937,-1.821255978,1.330274429};
  double da=0,db=0,ta,tb,fa,fb,g,lza,lzb;
  if(a<0||b<0){     
    Rprintf("\n****Tried to use approximation inappropriately. a=%f, b=%f\n",a,b);
    exit(1);
  }
  ta=1/(1+c0*a);
  tb=1/(1+c0*b);  
  lza=-pow(a,2)/2;
  lzb=-pow(b,2)/2;
  for(i=0;i<5;i++){
    da+=c[i]*pow(ta,i+1);
    db+=c[i]*pow(tb,i+1);
  }
  fa=exp(lza+log(da));
  fb=exp(lzb+log(db));
  g=log(fb-fa);
  return(-.5*log(2*M_PI)+g);
}


double expXphiAminusphiB(double x, double a, double b, int returnlog){
  double ret;
  switch(abs(a)>5&&abs(b)>5){
  case 1:
    if(a>5&&b>5){
      ret=x+logPhiAminusPhiB(a,b);
    }else if(a<-5&&b<-5){
      ret=x+logPhiAminusPhiB(-b,-a);
    }else{
      ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    }
    break;
  default:
    ret=x+log(pnorm(a,0,1,1,0)-pnorm(b,0,1,1,0));
    break;
  }
  if(!returnlog){
    return(exp(ret));
  }else{
    return(ret);
  }
}




SEXP MassAtChance1P(SEXP yR, SEXP NR, SEXP nSamp, SEXP mu0R, SEXP sig20R, SEXP a0R, SEXP b0R, SEXP decorrSig2R,SEXP progress, SEXP pBar,SEXP rho)
{
     
  int ITERS,iter,i,j,k,p,m,D,sumN,regionAlpha,regionMu,J,I;  
  double a0,b0,sig20,mu0,z,sig2Decorr;

  //input and check prior parameters
  ITERS=INTEGER_VALUE(nSamp);
  mu0=NUMERIC_VALUE(mu0R);
  sig20=NUMERIC_VALUE(sig20R);
  a0=NUMERIC_VALUE(a0R);
  b0=NUMERIC_VALUE(b0R);
  sig2Decorr=NUMERIC_VALUE(decorrSig2R);
  
  SEXP R_fcall,sampCounter,chainsR,prxlt0R,returnList;

  PROTECT(sampCounter = NEW_INTEGER(1));
  PROTECT(R_fcall = lang2(pBar, R_NilValue));
  
  
  I = INTEGER_POINTER(getAttrib(yR,R_DimSymbol))[0];
  J = INTEGER_POINTER(getAttrib(yR,R_DimSymbol))[1];

  
  int y[I][J],N[I][J];

  for(i=0;i<I;i++){
	for(j=0;j<J;j++){
		y[i][j] = INTEGER_POINTER(yR)[j*I+i]; 
		N[i][j] = INTEGER_POINTER(NR)[j*I+i];
	}
  }

  PROTECT(chainsR = allocMatrix(REALSXP, I+J+2, ITERS));
  double *chains = REAL(chainsR);
  
  PROTECT(prxlt0R = allocMatrix(REALSXP,I,J));
  PROTECT(returnList = allocVector(VECSXP, 2));
  
  double sumMu,falpha[J+1],galpha[J+1],halpha[J+1],fmu[I+1],gmu[I+1],hmu[I+1],Balpha[J+1],Bmu[I+1],wSum[I][J];
  double alpha[I],mu[J],sortAlpha[I],sortMu[J],wMu=0;
  double sumBalpha=0,probBalpha=0,regionAlphaUnif,sumBmu,probBmu=0,regionMuUnif;
  double sumAlphSqr,maxBmu,maxBalpha,prxlt0[I][J],sumAlpha=0,sumAllAlpha[I],sumAllMu[J],sumAllSigma2,pDecorr,meanBDecorr=0;

  double sig2=0.1; // Starting Value

  int bDecorr;

  GetRNGstate();
  
  /*Initializing starting values*/
 
  for(i=0;i<I;i++){
    for(j=0;j<J;j++){
      wSum[i][j]=y[i][j]*.5+(N[i][j]-y[i][j])*-.5;
      prxlt0[i][j]=0;
    }
    sumAllAlpha[i]=0;
  }
  
  for(j=0;j<J;j++){
    mu[j]=runif(0,1);
    sortMu[j]=mu[j];
    sumAllMu[j]=0;
 }
  
  
  qsort(sortMu, J, sizeof (double), compare_doubles);
  
  for(iter=0;iter<ITERS;iter++){
	
	
	R_CheckUserInterrupt();
	if(INTEGER_VALUE(progress) && !((iter+1)%INTEGER_VALUE(progress))){
			INTEGER_POINTER(sampCounter)[0]=iter+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho);
	}
	
	
	
    /*Sample alpha*/
    sumAlphSqr=0;
    sumAlpha=0;
    m=0;
    for(i=0;i<I;i++){
      regionAlpha=0;
      sumBalpha=0;
      falpha[0]=sig2;
      galpha[0]=0;
      halpha[0]=0;
      Balpha[0]=.5*log(sig2)+pnorm((-sortMu[J-1]/sqrt(sig2)),0,1,1,1);
      //sumBalpha+=Balpha[0];
      maxBalpha=Balpha[0];
      for(p=1;p<(J+1);p++){
	falpha[p]=0;
	galpha[p]=0;
	halpha[p]=0;
	sumN=0;
	probBalpha=0;
	for(j=0;j<J;j++){
	  D=mu[j]>sortMu[J-p-1];
	  if(D||(p==J)){  
	    halpha[p]+=N[i][j]*pow(mu[j],2)-2*mu[j]*wSum[i][j];
	    galpha[p]+=wSum[i][j]-N[i][j]*mu[j];
	    sumN+=N[i][j];
	    }
	}
	falpha[p]=1/(1/sig2+sumN);
	galpha[p]=galpha[p]*falpha[p];
        if(p==J){  
	  Balpha[p]=-.5*(halpha[p]-pow(galpha[p],2)/falpha[p])+.5*log(falpha[p])+pnorm((sortMu[0]+galpha[p])/sqrt(falpha[p]),0,1,1,1); 

	}else{
	  Balpha[p]=expXphiAminusphiB(-.5*(halpha[p]-pow(galpha[p],2)/falpha[p])+.5*log(falpha[p]),
	  		(-sortMu[J-p-1]-galpha[p])/sqrt(falpha[p]),(-sortMu[J-p]-galpha[p])/sqrt(falpha[p]),1); 

	}
	if(Balpha[p]>maxBalpha){maxBalpha=Balpha[p];}
      }
      for(p=0;p<(J+1);p++){
	Balpha[p]=exp(Balpha[p]-maxBalpha);
	sumBalpha+=Balpha[p];
      }
      regionAlphaUnif=runif(0,1);
      for(p=0;p<(J+1);p++){
	probBalpha+=Balpha[p]/sumBalpha;
	if(regionAlphaUnif>probBalpha) regionAlpha=p+1;
      }
      if(regionAlpha==0){
	alpha[i]=ran_trunc_norm_lower(galpha[regionAlpha], sqrt(falpha[regionAlpha]), -sortMu[J-1]);
      }else if(regionAlpha==J){
	alpha[i]=ran_trunc_norm_upper(galpha[regionAlpha], sqrt(falpha[regionAlpha]), -sortMu[0]);
      }else{
	alpha[i]=ran_trunc_norm_both(galpha[regionAlpha], sqrt(falpha[regionAlpha]), -sortMu[J-regionAlpha],-sortMu[J-regionAlpha-1]);
      }
      sumAlpha+=alpha[i];
      sumAlphSqr+=pow(alpha[i],2);
      sumAllAlpha[i]+=alpha[i];
      sortAlpha[i]=alpha[i];
      //fprintf(CHAINS,"%f ",alpha[i]);
    }
    /*end alpha*/
    //for(i=0;i<I;i++){
    //  alpha[i]=alpha[i]-sumAlpha/I;
    //  sumAlphSqr+=pow(alpha[i],2);
    //  fprintf(CHAINS,"%f ",alpha[i]);
    //  sortAlpha[i]=alpha[i];
    //  sumAllAlpha[i]+=alpha[i];
    //}
    
    qsort(sortAlpha, I, sizeof (double), compare_doubles);
 

    /*Sample mu*/
    
    m=0;
    sumMu=0;
    for(j=0;j<J;j++){
      regionMu=0;
      sumBmu=0;
      fmu[0]=sig20;
      gmu[0]=mu0;
      hmu[0]=0;
      Bmu[0]=-.5*(hmu[0]-pow(gmu[0],2)/fmu[0])+.5*log(fmu[0])+pnorm((-sortAlpha[I-1]-gmu[0])/sqrt(fmu[0]),0,1,1,1);
      maxBmu=Bmu[0];
      for(p=1;p<(I+1);p++){
	fmu[p]=0;
	gmu[p]=0;
	hmu[p]=0;
	sumN=0;
	probBmu=0;
	for(i=0;i<I;i++){
	  D=alpha[i]>sortAlpha[I-p-1];
	  if(D||p==I){  
	    hmu[p]+=N[i][j]*pow(alpha[i],2)-2*alpha[i]*wSum[i][j];
	    gmu[p]+=wSum[i][j]-N[i][j]*alpha[i];
	    sumN+=N[i][j];
	  }
	}
	fmu[p]=1/(1/sig20+sumN);
	gmu[p]=(mu0/sig20+gmu[p])*fmu[p];
	if(p==I){ 
	  Bmu[p]=-.5*(hmu[p]-pow(gmu[p],2)/fmu[p])+.5*log(fmu[p])+pnorm((sortAlpha[0]+gmu[p])/sqrt(fmu[p]),0,1,1,1);
	}else{ 
          Bmu[p]=expXphiAminusphiB(-.5*(hmu[p]-pow(gmu[p],2)/fmu[p])+.5*log(fmu[p]), (-sortAlpha[I-p-1]-gmu[p])/sqrt(fmu[p]), (-sortAlpha[I-p]-gmu[p])/sqrt(fmu[p]), 1);
	}
	if(Bmu[p]>maxBmu){maxBmu=Bmu[p];}
      }
      for(p=0;p<(I+1);p++){
	Bmu[p]=exp(Bmu[p]-maxBmu);
	sumBmu+=Bmu[p];
      }
      regionMuUnif=runif(0,1);
      for(p=0;p<(I+1);p++){
	probBmu+=Bmu[p]/sumBmu;
	if(regionMuUnif>probBmu) regionMu=p+1;
      }
      if(regionMu==0){
	mu[j]=ran_trunc_norm_lower(gmu[regionMu], sqrt(fmu[regionMu]), -sortAlpha[I-1]);
      }else if(regionMu==I){
	mu[j]=ran_trunc_norm_upper(gmu[regionMu], sqrt(fmu[regionMu]), -sortAlpha[0]);
      }else{
	mu[j]=ran_trunc_norm_both(gmu[regionMu], sqrt(fmu[regionMu]), -sortAlpha[I-regionMu],-sortAlpha[I-regionMu-1]);
      }
      sumAllMu[j]+=mu[j];
      sortMu[j]=mu[j];
      sumMu+=mu[j];
    }
    qsort(sortMu, J, sizeof (double), compare_doubles);
    /*end mu*/


    /*Sample wSum*/
    for(i=0;i<I;i++){
      for(j=0;j<J;j++){
	for(k=0;k<N[i][j];k++){
	  if(k==0){
	    if(alpha[i]+mu[j]>0){
	      wMu=alpha[i]+mu[j];
	    }else{
	      wMu=0;
	      prxlt0[i][j]+=(1/(double)(ITERS));
	    }
	    if(k<y[i][j]){
	      wSum[i][j]=ran_trunc_norm_upper(wMu,1,0);
	    }else{
	      wSum[i][j]=ran_trunc_norm_lower(wMu,1,0);
	    }
	  }else{
	    if(k<y[i][j]){
	      wSum[i][j]+=ran_trunc_norm_upper(wMu,1,0);
	    }else{
	      wSum[i][j]+=ran_trunc_norm_lower(wMu,1,0);
	    }
	    
	  }
	}
      }
    }/*end wSum*/
    
    /*decorrelating step*/
    z=rnorm(0,sqrt(sig2Decorr));
   
    pDecorr=z/2*((2*sumAlpha-I*z)/sig2-(2*(sumMu-J*mu0)-J*z)/sig20);
    if(pDecorr>0){
      pDecorr=1;
    }else{
      pDecorr=exp(pDecorr);
    }
    bDecorr=rbinom(1,pDecorr);
    if(bDecorr){
      sumAlphSqr=sumAlphSqr-2*z*sumAlpha+I*pow(z,2);
      for(i=0;i<I;i++){
        alpha[i]=alpha[i]-z;
        sumAllAlpha[i]=sumAllAlpha[i]-z;
		sortAlpha[i]=sortAlpha[i]-z;
        chains[iter*(I+J+2)+i] = alpha[i];
      }
      for(j=0;j<J;j++){
		mu[j]=mu[j]+z;
		sortMu[j]=sortMu[j]+z;
		sumAllMu[j]=sumAllMu[j]+z;
		chains[iter*(I+J+2)+I+j] = mu[j];
      }
    }else{
		Memcpy(&chains[iter*(I+J+2)],alpha,I);
	    Memcpy(&chains[iter*(I+J+2)]+I,mu,J);
    }
    meanBDecorr+=bDecorr*1.0/(float)(ITERS);  
    /*Sample sig2*/
    sig2=ran_inv_gamma(a0+.5*I,b0+(.5*sumAlphSqr));
    sumAllSigma2+=sig2;
    /*end sig2*/
	chains[iter*(I+J+2)+I+J] = sig2;
    chains[iter*(I+J+2)+I+J+1] = (double)(bDecorr);


    
  }//End MCMC
  
    for(i=0;i<I;i++){
		for(j=0;j<J;j++) {
			REAL(prxlt0R)[j*I+i] = prxlt0[i][j];
		}
	}
  
   /*
  for(i=0;i<I;i++){
    printf("alpha[%d] = %f  | ",i+1,sumAllAlpha[i]/ITERS);
    for(j=0;j<J;j++) {
      printf("%f  ",prxlt0[i][j]);
      if(j==(J-1)) printf("\n");
    }
  }
  for(j=0;j<J;j++)
    printf("mu[%d] = %f\n",j+1,sumAllMu[j]/ITERS);
  
  printf("sigma^2 = %f\n",sumAllSigma2/ITERS);
  printf("Acc. rate = %f\n",meanBDecorr); 
  */

  PutRNGstate();
  
  SET_VECTOR_ELT(returnList, 0, chainsR);
  SET_VECTOR_ELT(returnList, 1, prxlt0R);
  
  UNPROTECT(5);
  
  return(returnList);
}

