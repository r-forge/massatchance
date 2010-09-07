MassAtChance <-
function(y,N,MACtype="1P",n.samp=1,mu0=0,sig20=1,a0=2,b0=1,decorr.sig2=0.3,progress=FALSE){
	
	I=dim(y)[1]
	J=dim(y)[2]
	y = as.integer(y)
	dim(y)=c(I,J)
	N = as.integer(N)
	dim(N)=c(I,J)
	n.samp = as.integer(n.samp)
	progress = as.integer(abs(progress))

## Create progress bar
	if(progress){ pb = txtProgressBar(min = 0, max = n.samp, style = 3) }else{ pb=NULL }
	pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
## Call the C code
	
	if(MACtype=="1P"){
		returnVal = .Call("MassAtChance1P",y,N,n.samp,mu0,sig20,a0,b0,decorr.sig2,progress,pbFun,new.env(),PACKAGE="MassAtChance")
	}else{
		stop("Invalid MAC model type.")
	}
	if(progress & is.list(pb)) close(pb)
	cnames = c(paste("alpha",1:I,sep=""),paste("mu",1:J,sep=""),"sig2","accDecorr")
	returnVal[[1]] = t(returnVal[[1]])
	colnames(returnVal[[1]])=cnames
	returnVal[[1]] = mcmc(returnVal[[1]])
	return(returnVal)
}
