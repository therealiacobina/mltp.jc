# a Pop MC scheme for mixtures of two Gaussian distributions (24/02/2005)
# Film of the PMC with four normal random walk proposals
# nsimu: number of simulations
# niter: maximum number of iterations
# vari: variance of normal random walk proposal number i
# ---> Model:
# xi ~ p*N(vmu1,sigma1)+(1-p)*N(vmu2,sigma2) with sigma1=1 and vmu1=0
# prior:
# vmu1 ~ N(delta,sigma1/lambda) and vmu2 ~ N(delta,sigma2/lambda)

#############################################################################################
#############################################################################################
lpmix <- function(x,mu1,mu2,p,sigma1,sigma2,delta,lambda)
{
# log(f(x1,....,xn,mu1,mu2) pi(mu1,mu2))
  sum(log(p*dnorm(x,mu1,sqrt(sigma1))+(1-p)*dnorm(x,mu2,sqrt(sigma2))))+
    log(dnorm(mu1,delta,sqrt(sigma1/lambda)))+log(dnorm(mu2,delta,sqrt(sigma2/lambda)))
}

#############################################################################################
#############################################################################################
mltp.jc <- function(newppdyn=TRUE,ncycles=1,n=100,p=.8,nsimu=1000,niter=200,genu=TRUE,noRB=TRUE,
                    vmu1=0,vmu2=2.5,sigma1=1,sigma2=1,delta=1,lambda=1,gridsiz=150,
                    var1=.5,var2=.05,var3=0.005,var4=0.0005,
                    tau=1e-2,convcri=0,howmany=1,beta=1,
                    againD=FALSE,myseedD=NULL,nfseedD="",againP=FALSE,myseedP=NULL,nfseedP="",
                    doplot=2,hor=FALSE,ndir="simres",info=FALSE,comp=FALSE)
{
  library(splancs)
  library(pastecs)

# creating dirs
  if(newppdyn){
    alldirs<-c("/noRBnew","/RBnew")
  }else{
    alldirs<-c("/noRB","/RB")
  }
  
  if(ndir!=""){
    if(noRB){
      system(paste("mkdir -p ",ndir,alldirs[1],sep=""),ignore.stderr=TRUE)
    }else{
      system(paste("mkdir -p ",ndir,alldirs[2],sep=""),ignore.stderr=TRUE)
    }
  }else{
    ndir="./"
  }
  
  ndirup<-ndir
  if(noRB){
    ndir<-paste(ndirup,alldirs[1],sep="")
  }else{
    ndir<-paste(ndirup,alldirs[2],sep="")
  }


#  set or save seeds for data distribution
  svecD<-rep(0,ncycles)
  if(!againD || (againD && nfseedD=="" && is.null(myseedD))){
    svecD<-round(runif(ncycles)*10**9)
    if(againD && nfseedD=="" && is.null(myseedD)) message("no data seed given, generated")
    if(ncycles>1){
      fsvecD<-file(paste(ndirup,"/svecD_n-",n,".p-",p,".vmu2-",vmu2,".sigma2-",sigma2,sep=""),"w+")
      cat(svecD,"\n",sep="  ",file=fsvecD,append=TRUE)
      close(fsvecD)
    }
  }else if(!is.null(myseedD)&&nfseedD==""){
      svecD<-rep(myseedD,ncycles)
  }else if(is.null(myseedD)&&nfseedD!=""){
    fsvecD<-file(paste(ndirup,nfseedD,sep="/"),"r")
    svecD<-scan(fsvecD,quiet=TRUE)
  }else{
    message("myseedD",is.null(myseedD),"nfseedD",nfseedD)
    if(ncycles<2){
      svecD<-rep(myseedD,ncycles)
    }else{
      fsvecD<-file(paste(ndirup,nfseedD,sep="/"),"r")
      svecD<-scan(fsvecD,quiet=TRUE)
    }
  }
  
#  set or save seeds for particle distribution
  svecP<-rep(0,ncycles)
  if(!againP || (againP&&nfseedP==""&&is.null(myseedP))){
    svecP<-round(runif(ncycles)*10**9)
    if(againP && nfseedP=="" && is.null(myseedP)) message("no particle seed given, generated")
    if(ncycles>1){
      fsvecP<-file(paste(ndirup,"/svecP_n-",n,".p-",p,".vmu2-",vmu2,".sigma2-",sigma2,sep=""),"w+")
      cat(svecP,"\n",sep="  ",file=fsvecP,append=TRUE)
      close(fsvecP)
    }
  }else if(!is.null(myseedP)&&nfseedP==""){
      svecP<-rep(myseedP,ncycles)
  }else if(is.null(myseedP)&&nfseedP!=""){
    fsvecP<-file(paste(ndirup,nfseedP,sep="/"),"r")
    svecP<-scan(fsvecP,quiet=TRUE)
  }else{
    message("myseedP",is.null(myseedP),"nfseedP",nfseedP)
    if(ncycles<2){
      svecP<-rep(myseedP,ncycles)
    }else{
      fsvecP<-file(paste(ndirup,nfseedP,sep="/"),"r")
      svecP<-scan(fsvecP,quiet=TRUE)
    }
  }
  
# check if there are enough seeds for the number of cycles otherwise rescale it
  if(againD||againP){
    if(length(svecD)<ncycles || length(svecP)<ncycles){
      ncycles=min(length(svecD),length(svecP))
      cat("Warning: not enough seeds for the number of cycles...\n",
          "--> ncycles is now:",ncycles,"\n",sep=" ")
    }
  }

# various inputs checks  
  if(sum(doplot==c(0,1,2))==0){
    message("doplot not in {0,1,2}... setting doplot=2 (plot on screen)")
    doplot=2
  }
  if(sum(convcri==c(0,1,2,3))==0){ 
    convcri=1
    beta=1e-05
    message("convcri not in {0,1,2,3}... setting convcri= ",convcri,"  with beta=",beta)
  }
  if(convcri==1&&sum(howmany==c(1,2,3,4))==0){
    howmany=1
    message("howmany not in {1,2,3,4}... setting howmany=1 (weakest condition)")
  }
  
  for(icycle in 1:ncycles){
    myseedD=svecD[icycle]
    myseedP=svecP[icycle]
    if(info) fiosum<-file(paste(paste(ndir,"iosum_n",sep="/"),n,"c",icycle,"sD",myseedD,sep=""),"w")
    message(" myseedD=",myseedD,"   myseedP=",myseedP)

#  writing parameters in output file
    if(info){
      if(convcri==0){
        criterio="no convergence criterion"
      }else if(convcri==1){
        criterio="proposal weights variation rate"
      }else if(convcri==2){
        criterio="entropy variation rate"
      }else if(convcri==3){
        criterio="ess variation rate"
      }
      cat("# parameters ","\n",sep="",
          "#============","\n",
          "# newppdyn = ",newppdyn,"\n",
          "#        n = ",n,"\n",
          "#    nsimu = ",nsimu,"\n",
          "#  convcri = ",criterio,"\n", # choice of the conv. criterion
          "#  howmany = ",howmany,"\n",  # number of proposals to be cons. by the conv. criterion
          "#     beta = ",beta,"\n",     # value of parameter for the improved conv. equation
          "#      tau = ",tau,"\n",      # convergence threshold
          "#     genu = ",genu,"\n",
          "#     noRB = ",noRB,"\n",
          "#     vmu1 = ",vmu1,"\n",
          "#     vmu2 = ",vmu2,"\n",
          "#        p = ",p,"\n",
          "#   sigma1 = ",sigma1,"\n",
          "#   sigma2 = ",sigma2,"\n",
          "#    delta = ",delta,"\n",
          "#   lambda = ",lambda,"\n",
          "#     ndir = ",ndir,"\n",
          "#   againD = ",againD,"\n",
          "#  myseedD = ",myseedD,"\n",
          "#  nfseedD = ",nfseedD,"\n",
          "#   againP = ",againP,"\n",
          "#  myseedP = ",myseedP,"\n",
          "#  nfseedP = ",nfseedP,"\n",
          "#    niter = ",niter,"\n",
          "# (var1,var2,var3,var4) = ","(",var1,", ",var2,", ",var3,", ",var4,")","\n",
          file=fiosum,append=FALSE)
    }    

    jcenv<-jc(newppdyn,n,nsimu,niter,genu,noRB,vmu1,vmu2,p,sigma1,sigma2,delta,lambda,var1,var2,
              var3,var4,doplot,hor,icycle,againD,myseedD,againP,myseedP,info,gridsiz,ndir,
              ndirup,fiosum,tau,convcri,howmany,beta)
    
#   upload jc() outputs
    environment(jcenv)
    nu1<-get("nu1",envir=jcenv,mode="any")
    nu2<-get("nu2",envir=jcenv,mode="any")
    z<-get("z",envir=jcenv,mode="any")
    basin<-get("basin",envir=jcenv,mode="any")
    sppts<-get("sppts",envir=jcenv,mode="any")
    mapust<-get("mapust",envir=jcenv,mode="any")
    gridsiz<-get("gridsiz",envir=jcenv,mode="any")
    w<-get("w",envir=jcenv,mode="any")
    wmodeev<-get("wmodeev",envir=jcenv,mode="any")
    wmode<-get("wmode",envir=jcenv,mode="any")
    
    if(sppts$npeaks>0){
      if(sppts$npeaks>1){
# when the algorithm does not find the highest mode:    
        if(wmode[1,2]!=1){ # check for relevant asymmetry
          deltaz<-dist(sppts$coordpeaks[3,])
          for(im in 1:sppts$npeaks){
            if(wmode[im,2]==1){
              if(sum(deltaz[wmode[1:(im-1),2]-1]>(1e-4*abs(range(z)[1])))>0){
                asymeas<-c(abs(sppts$coordpeaks[3,2:(1+(im-1))]/range(z)[1]))
                if(info){
                  cat("\n","Algorithm failure...",sep=" ",file=fiosum,append=TRUE)
                  cat("\n","asymmetry index : ",asymeas,sep=" ",file=fiosum,append=TRUE)
                }else{
                  message("Algorithm failure...")
                  message("asymmetry index : ",asymeas)
                }     
              }
              break
            }
          }
        }
      }
      if(info){
        if(sppts$npeaks==1){
          cat("\n\n","total weight on modes :",wmode[,1],"\n\n",sep=" ",file=fiosum,append=TRUE)
        }else{
          cat("\n\n","total weight on modes :",wmode[order(wmode[,2],decreasing=FALSE),1],"\n\n",sep=" ",
              file=fiosum,append=TRUE)
        }
      }else{
        if(sppts$npeaks==1){
          cat("\n\n","total weight on modes :",wmode[,1],"\n\n",sep=" ")
        }else{
          cat("\n\n","total weight on modes :",wmode[order(wmode[,2],decreasing=FALSE),1],"\n\n",sep=" ")
        }
      }
    }
    
    if(doplot!=0){
      if(icycle!=ncycles){ graphics.off()}
    }
    
# before uncommenting one or more of the following lines, the relative "assign(...)" lines in function jc()  
# must also be uncommented or created.
#    
#  save(nu1,f1,file=paste(paste(ndir,"f1_n",sep="/"),n,"c",icycle,"sD",myseedD,sep=""),ascii=TRUE)
#  save(nu2,f2,file=paste(paste(ndir,"f2_n",sep="/"),n,sep=""),"c",icycle,sep=""),"sD",myseedD,sep=""),ascii=TRUE)
#  save(mapust,file=paste(paste(ndir,"mapust_n",sep="/"),n,sep=""),"c",icycle,sep=""),"sD",myseedD,sep=""),ascii=TRUE)
    
  if(sppts$npeaks>0){
    if(info){
      for(i in 1:niter) cat("# iter",i,"wm --> ",wmodeev[i,],"\n",sep=" ",file=fiosum,append=TRUE)
    }else{
      for(i in 1:niter) cat("# iter",i,"wm --> ",wmodeev[i,],"\n",sep=" ")      
    }
  }
 
    if(doplot==0) comp=FALSE
    if(comp){      
# possibility of compressing files in case of big jobs
      system(paste("tar cvf ",paste(ndir,"allres_n",sep="/"),n,"c",icycle,"sD",myseedD,
                   ".tar ",paste(ndir,"sim*.ps",sep="/"),sep=""))
# comment/uncomment/add/modify the appropriate lines (depending on the saved quantities)
#      system("gzip -f *.tar; rm f1* f2* z* mapust* basin* sppts* wmode* output* sim*.ps")
      system(paste("gzip -f ",paste(ndir,"*.tar",sep="/"),"; rm ",
                   paste(ndir,"sim*.ps",sep="/"),sep=""))
    }
    if(info) close(fiosum)
  }
}
       
#############################################################################################
#############################################################################################
jc <- function(newppdyn=TRUE,n=100,nsimu=1000,niter=200,genu=TRUE,noRB=TRUE,vmu1=0,vmu2=1,
               p=0.8,sigma1=1,sigma2=1,delta=1,lambda=1,var1=.5,var2=.05,var3=.005,
               var4=.0005,doplot=2,hor=FALSE,icycle=0,againD=FALSE,myseedD=NULL,
               againP=FALSE,myseedP=NULL,info=FALSE,gridsiz=150,ndir="simres",ndirup="./",
               fiosum="iosum000",tau=1e-2,convcri=0,howmany=1,beta=1)
{  
# simulation of the sample (x1,....,xn) from true or false distribution (true mean=vmu2)  
# data generation
  set.seed(myseedD)
  uu=runif(n)
  if(genu){
    x=rnorm(n,sd=1)+(uu<p)*vmu2
  }else{
#    x=rnorm(n,sd=1.)+(uu>.25&uu<.5)*vmu2-
#      (uu>.5&uu<.75)*2*vmu2+(uu>.75)*4*vmu2 #special sample with 3 modes
#    x=rnorm(n,sd=1.)+(uu>.25&uu<.5)*vmu2-(uu>.5&uu<.75)*2*vmu2+(uu>.75)*4*vmu2
    x=rnorm(n,sd=.1)-(uu>.2&uu<.4)*2*vmu2+(uu>.4&uu<.6)*vmu2-
      (uu>.6&uu<.8)*vmu2+(uu>.8)*2*vmu2 #special sample with many modes
  }

#-------------------------------- info/graphical stuff -------------------------------------
  esteso=FALSE
  titolo<-paste("n = ",n,",  p = ",p,",  vmu2 = ",vmu2,",  sigma2 = ",sigma2,",  nsimu = ",nsimu,
                ",  genu = ",genu,",  myseedD = ",myseedD,",  noRB = ",noRB,sep="")
  if(doplot==2){hor=FALSE}
  if(doplot!=0){
    if(doplot==1){
      nfplots<-paste(ndir,"/sim_n",n,"c",icycle,"sD",myseedD,"_1.ps",sep="")
      postscript(file=nfplots,horizontal=hor,paper="a4")
      mipos=dev.cur()
      if(hor){
        split.screen(c(1,2)) # one row two columns
        split.screen(c(2,1),screen=1)
      }else{
#      par(mfrow=c(3,1),pin=c(7,2),mar=c(2,2,2,0))  # if horizontal=F
        split.screen(c(2,1))
        split.screen(c(2,1),screen=1)
      }
    }else if(doplot==2){
      x11(height=10,width=12)
      mihist=dev.cur()
    }    
# data histogram  
    hist(x,br=80,xlab="",ylab="",col="gold3",asp=.5,main="")
  }
#--------------------------------------------------------------------------------------------
  
# Simulated mu's and logweights
  mu <- matrix(1,nsimu,2)
  lw <- 1:(nsimu)	   # logweights
  llw <- 1:(nsimu)	   # logweights
  la=matrix(0,ncol=4,nrow=nsimu)
  lala=lw
  res <-rep(0,nsimu) # component indices
  v <- c(var1,var2,var3,var4)*var(x) # properly scaled
  stdd <- sqrt(v)
  
# Previous mus for RB purposes
# First iteration : grid over potential space
  mupast <- expand.grid(seq(mean(x)-2*sd(x),mean(x)+2*sd(x),length=sqrt(nsimu)),
                        seq(mean(x)-2*sd(x),mean(x)+2*sd(x),length=sqrt(nsimu)))
  mupast=mupast[sample(1:dim(mupast)[1],nsimu),] # makes sure we get the right size
  mapust <- mupast # = x in paper notations
  
# Initial weights, entropy and effective sample size
  w <- rep(1,nsimu)/nsimu
  for (j in 1:nsimu)
    w[j]=lpmix(x,mapust[j,1],mapust[j,2],p,sigma1,sigma2,delta,lambda) # - 0 for uniform
  w <- exp(w-max(w))
  w <- w/sum(w)
  if(convcri==2){
    snew=0
    for(j in 1:nsimu){
      if(w[j]!=0){
        snew=snew-w[j]*log(w[j])
      }else{
        next
      }
    }
    s<-snew
  }else if(convcri==3){
    if(mean(w)!=0){
      ess<-1/(1+(var(w)/(mean(w)*mean(w))))
    }else{
      ess<-0
    }
  }

# Remembrance of past things (= x tilde in paper notations)
  mupast <- mapust[sample(1:nsimu,size=nsimu,rep=TRUE,prob=w),]

  pp <- matrix(1/4,1,4)
  ppnew <- matrix(0,1,4)
  pprate <- matrix(0,1,4)
  ppratenew <- matrix(0,1,4)
  
  nu1 <- seq(mean(x)-3*sd(x),mean(x)+3*sd(x),length=gridsiz)
  nu2 <- nu1
  
# Calculate or upload log-posterior surface
  namefzeta<-paste(ndirup,"/z_n-",n,".p-",p,".vmu2-",vmu2,".sigma2-",sigma2,".sD-",myseedD,sep="")
# watch out: the following flags are set to 'FALSE' if the file exists.
  nofzeta<-as.logical(system(paste("ls ",namefzeta),ignore.stderr=TRUE))
  if(nofzeta){
    z <- matrix(1,gridsiz,gridsiz)
    znodiv <- matrix(1,gridsiz,gridsiz)  
    for (i in 1:gridsiz){
      for (j in 1:gridsiz){
        z[i,j] <- lpmix(x,nu1[i],nu2[j],p,sigma1,sigma2,delta,lambda)
        znodiv[i,j]=z[i,j]
        if(z[i,j]==-Inf){znodiv[i,j]=-z[i,j]}
      }
    }
# Getting rid of possible divergences
    for (i in 1:gridsiz){
      for (j in 1:gridsiz){
        if(z[i,j]==-Inf){z[i,j]=min(znodiv)}
      }
    }
    zmax0<-z-max(z)
    fzeta<-file(namefzeta,"wb") 
    save(zmax0,file=fzeta)
    close(fzeta)
  }else{
    load(file=namefzeta)
  }    
 
# Fast&gridy approximation to the density
  f1=apply(exp(zmax0),1,sum) #*diff(range(nu1))/gridsiz
  f1=f1/(sum(f1)) #/diff(range(nu1))/gridsiz)
  f2=apply(exp(zmax0),2,sum) #*diff(range(nu1))/gridsiz
  f2=f2/(sum(f2)) #/diff(range(nu1))/gridsiz)

# Calculate or upload sppts' positions and their basins of influence
  namefsppts<-paste(ndirup,"/sppts_n-",n,".p-",p,".vmu2-",vmu2,".sigma2-",sigma2,".sD-",myseedD,sep="")
# watch out: the following flags are set to 'FALSE' if the file exists.
  nofsppts<-as.logical(system(paste("ls ",namefsppts),ignore.stderr=TRUE))                                  
  if(nofsppts){
    sppts<-findspecialpoints(nu1,nu2,zmax0,gridsiz)
    fsppts<-file(namefsppts,"wb")
    save(sppts,file=fsppts)
    close(fsppts)
  }else{
    load(file=namefsppts)
  }
  
  namefbasin<-paste(ndirup,"/basin_n-",n,".p-",p,".vmu2-",vmu2,".sigma2-",sigma2,".sD-",myseedD,sep="")
# watch out: the following flags are set to 'FALSE' if the file exists.
  nofbasin<-as.logical(system(paste("ls ",namefbasin),ignore.stderr=TRUE))
  if(nofbasin){
    basin<-findcontours(nu1,nu2,zmax0,sppts)
    fbasin<-file(namefbasin,"wb")
    save(basin,file=fbasin)
    close(fbasin)
  }else{
    load(file=namefbasin)
  }
  
#-------------------------------- info/graphical stuff -------------------------------------
  if(info){
    cat("# init   pp --> ",pp[1,],"\n\n",sep="      ",file=fiosum,append=TRUE)
    cat("z-range :",abs(range(zmax0)[1]),"\n",sep=" ",file=fiosum,append=TRUE)
    cat("modes :",sppts$npeaks,"\n",sep=" ",file=fiosum,append=TRUE)
    if(sppts$npeaks>0){
      for (im in 1:sppts$npeaks){
        cat(label=as.character(cat("mode",im,sep=" ",file=fiosum,append=TRUE)),
            sppts$coordpeaks[,im],sep="   ",fill=TRUE,file=fiosum,append=TRUE)
      }
    }
    if(sppts$npeaks>1){
      cat("\n",sep="",file=fiosum,append=TRUE)
      cat("pits :",sppts$npits,"\n",sep=" ",file=fiosum,append=TRUE)
      if(sppts$npits>0){
        for (im in 1:sppts$npits)
          cat(label=as.character(cat("pit",im,sep=" ",file=fiosum,append=TRUE)),
              sppts$coordpits[,im],sep="   ",fill=TRUE,file=fiosum,append=TRUE)
      }
      cat("\n",sep="",file=fiosum,append=TRUE)
      cat("saddle points:",sppts$nspts,"\n",sep=" ",file=fiosum,append=TRUE)
      if(sppts$nspts>0){
        for (im in 1:sppts$nspts)
          cat(label=as.character(cat("sadpt",im,sep=" ",file=fiosum,append=TRUE)),
              sppts$coordspts[,im],sep="   ",fill=TRUE,file=fiosum,append=TRUE)
      }
    }
    cat("\n",sep="",file=fiosum,append=TRUE)
  }else{
    cat("# init   pp --> ",pp[1,],"\n\n",sep="      ")
    cat("z-range :",abs(range(zmax0)[1]),"\n",sep=" ")
    cat("modes :",sppts$npeaks,"\n",sep=" ")
    if(sppts$npeaks>0){
      for (im in 1:sppts$npeaks)
        cat(label=as.character(cat("mode",im,sep=" ")),
            sppts$coordpeaks[,im],sep="   ",fill=TRUE)
    }
    if(sppts$npeaks>1){
      cat("pits :",sppts$npits,"\n",sep=" ")
      if(sppts$npits>0){
        for (im in 1:sppts$npits)
          cat(label=as.character(cat("pit",im,sep=" ")),
              sppts$coordpits[,im],sep="   ",fill=TRUE)
      }
      cat("saddle points:",sppts$nspts,"\n",sep=" ")
      if(sppts$nspts>0){
        for (im in 1:sppts$nspts)
          cat(label=as.character(cat("sadpt",im,sep=" ")),
              sppts$coordspts[,im],sep="   ",fill=TRUE)
      }
    }
  }
#-------------------
  if(doplot!=0){ 
    colo <- c("skyblue","blue","slategray4","black")
    if(doplot==1){
      if(hor){
        screen(4)
      }else{
        par(pin=c(7,2),mar=c(2,2,2,0))
#        screen(4)
      }
    }else if(doplot==2){    
      x11(height=10,width=12)
      midense=dev.cur()
    }    
# data densities
    plot(nu1,f1,type="l",lwd=2,xlab="",ylab="")
    lines(nu2,f2,type="l",col="gray",lwd=2)
    lines(density(mupast[,1]),lty=3,lwd=1)
    lines(density(mupast[,2]),col="gray",lty=3,lwd=1)

    if(doplot==1){
      if(hor){
        screen(2)
        par(pin=c(5,5))
      }else{
        screen(2)
        par(pin=c(5,5),mar=c(2,2,2,0))     
      }
    }else if(doplot==2){
      title(main=titolo,cex.main=1.25,font.main=2,col.main="black")
      x11(height=12,width=12)
      mimaj=dev.cur()
    }    
# z surface    
    image(nu1,nu2,zmax0,xlab=expression(mu[1]),ylab=expression(mu[2]),col=heat.colors(120));
    points(mupast[,1],mupast[,2],col=colo[1],pch=20)
    if(sppts$npeaks>0){
      for(ipp in 1:sppts$npeaks) points(sppts$coordpeaks[1,ipp],sppts$coordpeaks[2,ipp],
                                      pch=as.character(ipp),cex=2,col="red")
      for(im in 1:sppts$npeaks){lines(basin[im,2][[1]],basin[im,3][[1]],'l',lwd=2,col="green")}
    }
      
    if(doplot==1){
      close.screen(all=TRUE)
    }
    pagina=2 # per salvare una pag per file.ps.forse utile conservare opzione contraria...
  }
#-----------------------------------------------------------------------------------------
  
  oneseed=FALSE
  if(!oneseed){set.seed(myseedP)}


#=========================================================================================
#=====>            STARTING ITERATIONS  
#=========================================================================================
  i=0
  conv=FALSE  # label convergenza
  while(!conv){
    i=i+1
    res=sample(1:4,nsimu,rep=TRUE,prob=pp[i,]) # component index
    mu[,1] = mupast[,1]+stdd[res]*rnorm(nsimu)
    mu[,2] = mupast[,2]+stdd[res]*rnorm(nsimu)
      
    for (j in 1:nsimu){
#      mu[j,]=c(rnorm(1,mupast[j,1],stdd[res[j]],rnorm(1,mupast[j,2],stdd[res[j]]))
      if (noRB){ # No add'al RB
        for(ll in 1:4){la[j,ll]=dnorm(mu[j,1],mupast[j,1],stdd[ll])*
                                dnorm(mu[j,2],mupast[j,2],stdd[ll])}
      }else{ # Double RB
        for (ll in 1:4){la[j,ll]= sum(w*dnorm(mu[j,1],mapust[,1],stdd[ll])*
                                        dnorm(mu[j,2],mapust[,2],stdd[ll]))}
      }
      lala[j]=sum(pp[i,]*la[j,])
      lw[j]=lpmix(x,mu[j,1],mu[j,2],p,sigma1,sigma2,delta,lambda)-log(lala[j])
    }
    
    w <- exp(lw-max(lw))
    w <- w/sum(w)

    if(i==1){
      wmodeev<- countpart(nsimu,sppts,basin,mapust,w)
    }else{
      wmodeev <- rbind(wmodeev,countpart(nsimu,sppts,basin,mapust,w))
    }

# updating of the proposal weights: new pp dyn or old pp dyn
    if(newppdyn){
      for(ip in 1:4){
        if(i<=3){
          ppnew[ip] <- sum(w*pp[i,ip]*la[,ip]/
                           (pp[i,1]*la[,1]+pp[i,2]*la[,2]+pp[i,3]*la[,3]+pp[i,4]*la[,4]))
        }else{
          ppnew[ip] <- (1-beta)*pp[i,ip]+beta*
            sum(w*pp[i,ip]*la[,ip]/
                (pp[i,1]*la[,1]+pp[i,2]*la[,2]+pp[i,3]*la[,3]+pp[i,4]*la[,4]))
        }
      }
    }else{ # old pp dyn
      if (noRB){
        for(ip in 1:4){
          if(i<=3){
            ppnew[ip] <- sum(w[res==ip]) #length(res[res==1])/nsimu
          }else{
            ppnew[ip] <- (1-beta)*pp[i,ip]+beta*sum(w[res==ip])
          }
        }
      }else{
        for (ip in 1:4){
          if(i<=3){        
            ppnew[ip]=pp[i,ip]*sum(w*(sum(la[,ip])/lala))
          }else{
            ppnew[ip] <- (1-beta)*pp[i,ip]+beta*pp[i,ip]*sum(w*(sum(la[,ip])/lala))
          }
        }
      } 
    }
    
    ppnew <- ppnew/sum(ppnew)
    pp<-rbind(pp,ppnew)

    if(convcri==2){
#Shannon's entropy    
      snew=0
      for(j in 1:nsimu){
        if(w[j]!=0){
          snew=snew-w[j]*log(w[j])
        }else{
          next
        }
      }
      s<-c(s,snew)
    }else if(convcri==3){
# Effective Sample Size
      if(mean(w)!=0){
        ess<-c(ess,1/(1+(var(w)/(mean(w)*mean(w)))))
      }else{
        ess<-c(ess,0)
      }
    }

#-------------------------------- info/graphical stuff ------------------------------------    
    if(info){
      cat(paste("# iter ",i,sep=""),"pp --> ",pp[i,],"\n",sep=" ",
          file=fiosum,append=TRUE)
    }else{
      cat(paste("# iter ",i,sep=""),"pp --> ",pp[i,],"\n",sep=" ")
    }
      
#===============================================================
#          Analysis of convergence    
#===============================================================
    
    if(convcri==0){
      if(i==1) message("No convergence criterium imposed.")
    }else{
      dp=1e-6 # Maximum value for p[i]-p[i-1]
      iterm=3 # "thermalizing" time... do we need it???
    }
    if(convcri==1){      
# compute proposals weights rates
      if(i==1){
	message("Convergence criterion: proposal weights")
        pprate<-rep(0,4)
      }else{
        for(j in 1:4){
          if(pp[i,j]==0 && pp[i-1,j]==0){ # pbs for (0-0)/0=NaN
            ppratenew[j]=0
          }else{ 
            ppratenew[j]<-abs(1-pp[i-1,j]/pp[i,j])  # ==> dt=1
          }
        } 
        pprate<-rbind(pprate,ppratenew)
      }      
# checking convergence condition on pp()
      if(i>iterm&&(sum(pprate[i,]<tau)>=howmany&&sum(pprate[i-1,]<tau)>=howmany&&sum(pprate[i-2,]<tau)>=howmany&&
                sum(pprate[i-3,]<tau)>=howmany&&sum(pprate[i-4,]<tau)>=howmany)
         ||i>iterm&&(sum(abs(pp[i,]-pp[i-1,])<dp)>=howmany&&sum(abs(pp[i-1,]-pp[i-2,])<dp)>=howmany&&
                sum(abs(pp[i-2,]-pp[i-3,])<dp)>=howmany)){
        conv=TRUE
        break
      }
    }else if(convcri==2){
# compute entropy rate s()
      if(s[i]==0&&s[i+1]==0){
        if(i==1){
	  message("Convergence criterion: entropy")
          srate<-0
        }else{
          srate<-c(srate,0)
        }
      }else{
        if(i==1){
	  message("Convergence criterion: entropy")
          srate<-abs(1-s[i+1]/s[i])
        }else{
          srate<-c(srate,abs(1-s[i+1]/s[i]))
        }
      }
# checking convergence condition on s()      
      if(i>iterm&&srate[i]<tau&&srate[i-1]<tau&&srate[i-2]<tau&&srate[i-3]<tau&&srate[i-4]<tau){
        conv=TRUE
        break
      }
    }else if(convcri==3){
# compute ess rate    
      if(ess[i]==0&&ess[i+1]==0){
        if(i==1){
          message("Convergence criterion: ESS")
          essrate<-0
        }else{
          essrate<-c(essrate,0)
        }
      }else{
        if(i==1){
          message("Convergence criterion: ESS")
          essrate<-abs(1-ess[i+1]/ess[i])
        }else{
          essrate<-c(essrate,abs(1-ess[i+1]/ess[i]))
        }
      }
# checking convergence confition on ess()
      if(i>iterm&&essrate[i]<tau&&essrate[i-1]<tau&&essrate[i-2]<tau&&essrate[i-3]<tau&&essrate[i-4]<tau){
        conv=TRUE
        break
      }
    }
      
#-------------------------------- info/graphical stuff ------------------------------------    
    if(doplot!=0){
      if(doplot==1){     
        if(as.logical(i==1) || as.logical((i-1)/12==i%/%12)){
          dev.off()
          nfplots<-paste(paste(ndir,"sim_n",sep="/"),n,"c",icycle,"sD",myseedD,"_",pagina,".ps",sep="")
          postscript(file=nfplots,horizontal=hor,paper="a4")
          mipos2=dev.cur()
          if(hor){
            par(mfrow=c(3,4),mai=c(.8,.6,0,0),pty="s")
          }else{
            par(mfrow=c(4,3),mai=c(.8,.6,0,0),pty="s")
          }
          pagina=pagina+1
        }
        image(nu1,nu2,zmax0,xlab=expression(mu[1]),ylab=expression(mu[2]),col=heat.colors(120))
        for (ll in 1:4){points(mu[res[]==ll,1],mu[res[]==ll,2],col=colo[ll],pch=20)}
        if(sppts$npeaks>0){
          for(ipp in 1:sppts$npeaks){
            points(sppts$coordpeaks[1,ipp],sppts$coordpeaks[2,ipp],pch=as.character(ipp),cex=2,col="red")}
          for(im in 1:sppts$npeaks){
            lines(basin[im,2][[1]],basin[im,3][[1]],'l',lwd=2,col="green")}
        }
      }else if(doplot==2){
        dev.set(mimaj)
        title(main=titolo,cex.main=1.25,font.main=2,col.main="black")
        if(i==niter){image(nu1,nu2,zmax0,xlab=expression(mu[1]),
                           ylab=expression(mu[2]),col=heat.colors(120))}
        for (ll in 1:4)
          points(mu[res[]==ll,1],mu[res[]==ll,2],col=colo[ll],pch=20)
        if(sppts$npeaks>0){
          for(ipp in 1:sppts$npeaks){
            points(sppts$coordpeaks[1,ipp],sppts$coordpeaks[2,ipp],pch=as.character(ipp),cex=2,col="red")}
          for(im in 1:sppts$npeaks){
            lines(basin[im,2][[1]],basin[im,3][[1]],'l',lwd=2,col="green")}
        }
      }
    }
#----------------------------------------------------------------------------------
    
    samp <- sample(1:nsimu,nsimu,prob=w,rep=TRUE)
    mapust <- mu
    mupast <- mu[samp,]

#-------------------------------- info/graphical stuff ------------------------------------        
    if(doplot!=0){
      if(esteso){
# Partial plot of the output linear density and the resampled particles
        if(i!=niter){
          if(doplot==2){dev.set(midense)}          
          plot(nu1,f1,type="l",lwd=1)
          lines(nu2,f2,type="l",col="darkgray",lwd=1)
          lines(density(mupast[,1]),lty=3,lwd=1)
          lines(density(mupast[,2]),col="darkgray",lty=3,lwd=1)
          
          if(doplot==2){dev.set(mimaj)}    
          image(nu1,nu2,zmax0,xlab=expression(mu[1]),
                ylab=expression(mu[2]),col=heat.colors(120))
          points(mupast[,1],mupast[,2],col="sienna4",pch=20)
          if(sppts$npeaks>0){
            for(ipp in 1:sppts$npeaks){
              points(sppts$coordpeaks[1,ipp],sppts$coordpeaks[2,ipp],pch=as.character(ipp),cex=2,col="red")}
            for(im in 1:sppts$npeaks){
              lines(basin[im,2][[1]],basin[im,3][[1]],'l',lwd=2,col="green")}
          }
        }
      }
    }
#-------------------------------------------------------------------------------------------    

    if(i==niter){
      if(convcri!=0){ 
        message(criterio,"Maximum number of iteration reached... ")
      }
      break
    }
  }   # end loop on iterations

  wmode<-as.matrix(cbind(wmodeev[i,],seq(1:sppts$npeaks)))
  if(sppts$npeaks>1){wmode<-wmode[order(wmode[,1],decreasing=TRUE),]}

#-------------------------------- info/graphical stuff ------------------------------------        
# Final plot of the output linear density and the resampled particles
  if(doplot!=0){
    if(esteso){   
      if(doplot==2){dev.set(mimaj)}    
      image(nu1,nu2,zmax0,xlab=expression(mu[1]),ylab=expression(mu[2]),col=heat.colors(120))
      points(mupast[,1],mupast[,2],col="sienna4",pch=20)
      if(sppts$npeaks>0){
        for(ipp in 1:sppts$npeaks){
          points(sppts$coordpeaks[1,ipp],sppts$coordpeaks[2,ipp],pch=as.character(ipp),cex=2,col="red")}
        for(im in 1:sppts$npeaks){
          lines(basin[im,2][[1]],basin[im,3][[1]],'l',lwd=2,col="green")}
      }
      
      if(doplot==2){dev.set(midense)}      
      plot(nu1,f1,type="l",lwd=2)
      lines(nu2,f2,type="l",col="darkgray",lwd=2)
      lines(density(mupast[,1]),lty=3,lwd=1)
      lines(density(mupast[,2]),col="darkgray",lty=3,lwd=1)
    }

    if(doplot==1){
      dev.off()
      graphics.off()
    }
  }  
#-------------------------------------------------------------------------------------------    
  
#  pass various variables to mltp.jc() as an environment
  jcenv <- new.env(parent=globalenv())
#  assign("x",x,envir=jcenv)
  assign("nu1",nu1,envir=jcenv)
#  assign("f1",f1,envir=jcenv)
  assign("nu2",nu2,envir=jcenv)
#  assign("f2",f2,envir=jcenv)
  assign("z",zmax0,envir=jcenv)
  assign("sppts",sppts,envir=jcenv)
  assign("basin",basin,envir=jcenv)
  assign("gridsiz",gridsiz,envir=jcenv)
  assign("mapust",mapust,envir=jcenv)
  assign("w",w,envir=jcenv)
  assign("wmodeev",wmodeev,envir=jcenv)
  assign("wmode",wmode,envir=jcenv)
  return(jcenv)
}

#############################################################################################
#############################################################################################
findspecialpoints <- function(nu1,nu2,z,gridsiz)
{
#   compute the coordinates of the logposterior surface peaks

#  initialize the vector of peaks positions: the maximum number of peaks a priori is npeaksf1*npeaksf2
  totpeaks<-matrix(0,3,2*round(gridsiz/2))
  totpits<-matrix(0,3,2*round(gridsiz/2))
  totspts<-matrix(0,3,2*round(gridsiz/2))
  oucaturn<-matrix(0,gridsiz,1)             # vector of unique positions for orthogonal sounding

#  count the total number of peaks and save their coordinates
  npeaks=0
  npits=0
  nspts=0
  for(itp in 1:gridsiz){
    tpz0<-turnpoints(z[itp,])   # search for turnpoints in vertical directions
#  count the peaks of current z row
    if(tpz0$nturns==0) next  # itp
    newturns<-matrix(0,2,tpz0$nturns)
    for(inpz in 1:tpz0$nturns){
      newturns[1,inpz]=itp
      newturns[2,inpz]=tpz0$tppos[inpz]
    }
    if(itp==1){
      turns<-newturns
    }else{
      turns<-cbind(turns,newturns)
    } 
  }  #itp
# Maxima and saddle points
  if(length(turns)>0){
    oucaturn<-as.matrix(unique(turns[2,]))
    for (itp in 1:as.numeric(dim(oucaturn)[1])){
      tpz<-turnpoints(z[,oucaturn[itp,1]])      
      for(inpz in 1:tpz$nturns){
        tpzaux<-turnpoints(z[tpz$tppos[inpz],])
        if(tpz$peak[tpz$tppos[inpz]]){
          if(tpzaux$peak[oucaturn[itp,1]]){ # maximum
            npeaks=npeaks+1
            totpeaks[1,npeaks]=nu1[tpz$tppos[inpz]]
            totpeaks[2,npeaks]=nu2[oucaturn[itp,1]]
            totpeaks[3,npeaks]=z[tpz$tppos[inpz],oucaturn[itp,1]]
          }else if(tpzaux$pit[oucaturn[itp,1]]){
            nspts=nspts+1
            totspts[1,nspts]=nu1[tpz$tppos[inpz]]
            totspts[2,nspts]=nu2[oucaturn[itp,1]]
            totspts[3,nspts]=z[tpz$tppos[inpz],oucaturn[itp,1]]
          }
        }else if(tpz$pit[tpz$tppos[inpz]]){
          if(tpzaux$pit[oucaturn[itp,1]]){  # minimum
            npits=npits+1
            totpits[1,npits]=nu1[tpz$tppos[inpz]]
            totpits[2,npits]=nu2[oucaturn[itp,1]]
            totpits[3,npits]=z[tpz$tppos[inpz],oucaturn[itp,1]]
          }else if(tpzaux$peak[oucaturn[itp,1]]){
            nspts=nspts+1
            totspts[1,nspts]=nu1[tpz$tppos[inpz]]
            totspts[2,nspts]=nu2[oucaturn[itp,1]]
            totspts[3,nspts]=z[tpz$tppos[inpz],oucaturn[itp,1]]
          }
        }
      }
    }
    if(npeaks>0){
      peaks<-totpeaks[,1:npeaks]
      if(length(peaks)<=3){   # meaning only one peak was found
        totpeaks<-as.matrix(peaks)
      }else{  # order modes in decreasing z values
        totpeaks<-peaks[,order(peaks[3,],decreasing=TRUE)]
      }
    }else{
      totpeaks=0
    }
    if(npits>0){    
      pits<-totpits[,1:npits]
      if(length(pits)<=3){   # meaning only one pit was found
        totpits<-as.matrix(pits)
      }else{  # order sppts in decreasing z values
        totpits<-pits[,order(pits[3,],decreasing=TRUE)]
      }
    }else{
      totpits=0
    }
    if(nspts>0){
      spts<-totspts[,1:nspts]
      if(length(spts)<=3){
        totspts<-as.matrix(spts)
      }else{
        if(sum(duplicated(t(spts)))>0){
          uniquespts<-spts[,!duplicated(t(spts))]
          if(length(uniquespts)>=3){
            totspts<-as.matrix(uniquespts)
          }else{
            totspts<-uniquespts[,order(uniquespts[3,],decreasing=TRUE)]
          }
        }else{
          totspts<-spts[,order(spts[3,],decreasing=TRUE)]
        }
      }
    }else{
      totspts=0
    }
  }else{
    totpeaks=totpits=totspts=0
  }  
  sppts<-list(npeaks=npeaks,coordpeaks=totpeaks,npits=npits,coordpits=totpits,nspts=nspts,coordspts=totspts)
  return(sppts)
}

#############################################################################################
#############################################################################################
findcontours <- function(nu1,nu2,z,sppts)
{
# compute the largest contourlines encircling each mode.
  basin<-matrix(list(0),sppts$npeaks,3)
  lavora<-as.matrix(sppts$coordpeaks)
  livello<-numeric(0)
  imm<-seq(1:sppts$npeaks)
  modeOK<-rep(FALSE,sppts$npeaks)
  checkmode<-rep(FALSE,sppts$npeaks)
# starting level choice (there is maybe a more clever choice)
  if(sppts$nspts==0&&sppts$npits==0){ # c'e' un solo modo (eventualmente "multiplo")
    livello=min(sppts$coordpeaks[3,])+(min(sppts$coordpeaks[3,])-abs(min(z)))/89.1 # a default value fixed by hand
                                                                                   # on the base of many trials
  }else{  # more modes
    if(sppts$nspts>0){
      livello<-c(livello,sppts$coordspts[3,1:sppts$nspts])
    }
    if(sppts$npits>0){
      livello<-c(livello,sppts$coordpits[3,1:sppts$npits])
    }
    livello<-livello[order(livello,decreasing=TRUE)]
  }

  il=0
  while(sum(modeOK)<sppts$npeaks){
    il=il+1
    if(il>length(livello))
      stop("findcontours not working as expected :-(  .... stopping...")
    pappero<-contourLines(nu1,nu2,z,levels=livello[il]*.9)
    for(im1 in imm){
      imm2=imm[imm!=im1]
      for(icont in 1:length(pappero)){
        checkmode[im1]=inout(t(sppts$coordpeaks[,im1]),as.points(pappero[[icont]]$x,pappero[[icont]]$y),bound=TRUE)
        if(checkmode[im1]){
          for(im2 in imm2){
            checkmode[im2]=inout(t(sppts$coordpeaks[,im2]),as.points(pappero[[icont]]$x,pappero[[icont]]$y),
                       bound=TRUE)
            if(checkmode[im2]){
              if(!modeOK[im2]){
                modeOK[im2]=TRUE
              }else{
                modeOK[im1]=TRUE
                break
              }
            } # if(!modeOK[im2])
          } # for im2
          if(!modeOK[im1]){
            basin[im1,]=pappero[[icont]]
            modeOK[im1]=TRUE
            break  # next im1
          }else{ # if(!modeOK[im1])
            break
          }
        } # if(checkmode[im1]) 
      } # for icont
      if(icont==length(pappero)&&il==length(livello)) modeOK[im1]=TRUE 
    } # for im1
  } # while
  return(basin)
}

#############################################################################################
#############################################################################################
countpart <- function(nsimu,sppts,basin,mapust,w)
{
#  count the total weigth associated to each mode. it uses the function inpip() 
#  of library "splancs" to select particles inside the basin of attraction of each mode    

  wmode<-rep(0,sppts$npeaks)
  checkpar<-rep(FALSE,nsimu)
  for(im in 1:sppts$npeaks){
    pup<-as.points(basin[im,2][[1]],basin[im,3][[1]])
    dove<-inpip(mapust,pup,bound=TRUE,quiet=TRUE)
    if(length(dove)>0){
      if(sum(checkpar[dove])!=0)
        message("Some particles are counted twice... maybe basins bounds overlap?")
      wmode[im]=sum(w[dove])
      checkpar[dove]=TRUE
    }
    if(sum(checkpar)==nsimu) break
  }
  return(wmode)
}
