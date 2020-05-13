# Warning!
# This software comes as it is. No guarantee whatsoever is given. 
# The authors cannot be held responsible for any misunderstanding, 
# incorrect use, false scientific conclusions or other problems using 
# these programs. 
# Please give proper reference to the corresponding paper when used. 
# Charkhi, A. and Claeskens, G. (2018). Asymptotic post-selection inference 
# for Akaike's information criterion. Biometrika.
# This code is written by A. Charkhi and further adjusted by A. Garcia Angulo and G. Claeskens
# May 2018, allowing for factor variables and an offset.
# May 2020, fixing a bug with empty factor levels

library('tmg')

PostAIC=function(y,X,offset=NULL,model.set='user',quant=0.95, common=c(),intercept=TRUE, 
			family=gaussian, Pi=NULL, linearcomb=FALSE,xnew=X[1,]){

# Required arguments:  
# y: a nx1 vector of response values
# X: a nxp matrix containing the covariates, it is not required to have a column with ones for  
#   the intercept in this matrix.   
#   Inclusion of an intercept can be set via 'intercept=TRUE'. 

# Options:
# intercept=TRUE: if there is no column of 1s in the X matrix, but an intercept should
#   be included to fit each model. 
#   The intercept is not assumed to be part of every model with allsubsets.
#   Use model.set='partsubsets' with common=c() to have all subsets with common intercept
# intercept=FALSE: No automatic intercept is added in the model fits. Use this option if
#   a column of 1s is already contained in the X matrix or no intercept is needed.

# model.set='user'  this is for a user-specified set of models, Pi is required
# model.set='allsubsets' all possible models for all parameters
# model.set='partsubsets' all possible models for a subset of the parameters  
# model.set='nested' a sequence of nested models, starting with including the variable
#          in column 1 of X, then adding the variable in column 2 etc.

# family: for glm model fitting. Extensions to general likelihood are straightforward, 
#   by using the Fisher information matrix Jfull corresponding to the largest fitted model, 
#   but this is not (yet) included in this function.
  
# quant: confidence level, default is 0.95.
  
# Pi: indicates which variables are in each model, one row for each model. This is required
#  for model.set="user" to know which models are used.
#  This matrix is automatically computed for model.set equal to allsubsets, partsubsets and nested.
#  The option intercept=TRUE requires to have the first column of Pi equal to 1s;
#  if X already contains a column of 1s for the intercept, Pi has p columns, otherwise, in that case 
#  Pi has p+1 columns.

# common: a vector of length at most p, consisting of those column numbers of X for which the
#   variables are forced to be present in every model (e.g. c(2,5) for variables 2 and 5). 
#   If provided, it is used in model.set="partsubsets" to construct the Pi matrix, and in the model selection.

# linearcomb: When TRUE, xnew is required;
#    confidence intervals for linear combinations of the form xnew%*%beta, xnew needs to be 
#    specified as a matrix with the same number of columns as X, an intercept can be added via intercept=T.
#    When linearcomb=FALSE, confidence intervals for the components of beta.

# Output:
# A list containing the "Model selected by AIC and inference ignoring model selection"
#  "PostAIC", consisting of the separate quantiles of the specified level for each component in the 
#  selected model, the simultaneous quantile, the constraints from the model comparison, and the 
#  "PostAIC intervals" either for each parameter appearing in the selected model or for the linear combinations.


    ####################### Pi matrix ###################################
  
  if (is.null(Pi)){
    if(model.set=='user'){stop('Pi matrix is required with a user-specified set of models')}    

    if (model.set=='allsubsets'){   ####### All possible models for all parameters
      if (intercept==T){   
        Pi1=combinations(ncol(X)+1)
      }else {
        Pi1=combinations(ncol(X))
      }
      Pi = Pi1[rowSums(Pi1 != 0) > 0,]  # no 'empty' model
    }
    
    else if(model.set=='partsubsets'){  ######## All Possible for a subset
      if(intercept==T){ 
        if(is.null(common)){
          Pi1=combinations(ncol(X))  
          Pi=cbind(1,Pi1)
        } else {
          Pi1=combinations(ncol(X)-length(common))
          common1<-matrix(1,nrow(Pi1),length(common))
          Pi<-cbind(1,common1,Pi1)
        }}
      if(intercept==F){
        if(!is.null(common)){
          Pi1=combinations(ncol(X)-length(common))
          common1<-matrix(1,nrow(Pi1),length(common))
          Pi<-cbind(common1,Pi1)
        }
        if(is.null(common)) stop('Either the intercept or another parameter should be common to all models in model.set=partsubsets')
        
      }
    }
    
    else if(model.set=='nested'){  ######## Nested models
      m2 <- matrix(1, ncol(X), ncol(X))
      m2[upper.tri(m2, diag = FALSE)]<-0
      Pi1<-m2
      if (intercept==T){   
        Pi=cbind(1,Pi1)
      }else {
        Pi=Pi1
      }
    }
    
  } else {
    
    Pi=Pi
  }
  
  
  
  ####################### Model selection ###################################
  
  Mn=nrow(Pi)
  aic=rep(NA,Mn)
  
  if(!is.null(common)){
    Xcommon=X[,common]
    Xother= X[,-common]
    Xcommon2<-cbind(Xcommon,Xother)
    
    if (intercept==T){    
      Xfull=cbind(Intercept=1,Xcommon2)
    }
    else{
      Xfull=Xcommon2}}
  
  if(is.null(common)) {
    if (intercept==T){
      Xfull=cbind(Intercept=1,X)
    }
    else{
      Xfull=X
    }
  }
  
  
  #if(class(Xfull)!='matrix'){
  #  Xfull=matrix(unlist(Xfull),ncol=ncol(Pi))}
  
  for (k in 1:Mn){
    data = as.data.frame(cbind(Xfull[,which(Pi[k,]!=0)],y=y))#sweep(Xfull,MARGIN=2,Pi[k,],'*') 
    model = glm(y~.-1,offset=offset,data=data,family=family) 
    aic[k] = AIC(model) 
  }
  
  pselected= which.min(aic)
  X.aic = as.data.frame(Xfull[,which(Pi[pselected,]!=0)]) #sweep(Xfull,MARGIN=2,Pi[pselected,],'*') 
  X.aic = X.aic[, colSums(X.aic != 0) > 0]
  newdata = as.data.frame(cbind(y=y,X.aic))
  
  if(!is.null(common)){
    
    if (is.null(names(X))){
      namesX=c(paste('X',common,sep=''),paste('X', (1:dim(X)[2])[-common], sep=''))
    }else { namesX=c(names(X)[common],names(X)[-common])}
    
  }
  
  if(is.null(common)){
    if (is.null(names(X))){
      namesX=paste('X', 1:ncol(X), sep='')
    }else { namesX=names(X)}}
  
  
  if(intercept==T){
    namesXfull<-c('Intercept',namesX)
    names.newdata<-c('y',namesXfull[which(Pi[pselected,]!=0)])
    names(newdata)<-names.newdata
  } else {  namesXfull<-c(namesX)
  names.newdata<-c('y',namesXfull[-which(Pi[pselected,]==0)])
  names(newdata)<-names.newdata}
  
  
  fit.aic= glm(y~.-1, offset=offset, data = newdata, family=family)
  summary.fit=summary(fit.aic)$coefficients
  
  
  ########################## Jfull, Fisher information matrix ####################
  
  data.full=as.data.frame(cbind(y,Xfull))
  fit.full = glm(y~.-1,offset=offset,family=family, data=data.full)
  
  if (length(y)<30){
    Jfull=solve((vcov(fit.full,complete = F)*(length(y)-length(fit.full$coefficients))/length(y))*length(y))  
    
  }else {
    Jfull=solve(length(y)*vcov(fit.full,complete = F)) 
  }


  #factors
  
  splits<-list()
  
  for (i in 1:ncol(Pi)){
    
    if(is.factor(Xfull[,i])){
      factori = droplevels(Xfull[,i])
      lev<-length(levels(factori))-1
      splits[[i]]<-matrix(rep(Pi[,i],each=lev), ncol=lev, byrow=TRUE)
    }else{
      splits[[i]]<-Pi[,i]
    }
    
  } 
  
  Pi.new2<-do.call(cbind,splits)
  
   
if(linearcomb==FALSE){
  
  results1<-Aicpost(Pi=Pi.new2,p=pselected,npar=ncol(Pi.new2),model.set=model.set,Information=Jfull,quant=quant,nsamp=200000,simplify.constraints=TRUE,morecons=FALSE,unsure=c(),show.constraints=FALSE)
  
  quantiles = unlist(results1[1])
  if(sum(is.na(fit.aic$coefficients))>0){
    coefs<-fit.aic$coefficients[-which(is.na(fit.aic$coefficients))]
  } else {
    coefs<-fit.aic$coefficients
  }
 
  intervals=matrix(NA,length(coef),2)
  
  intervals=cbind(coefs-quantiles/sqrt(length(y)), coefs+quantiles/sqrt(length(y)))
  rownames(intervals)<- names(coefs[is.na(as.vector(coefs))==F])

  outlist=list(summary.fit,results1,intervals)
  names(outlist)=c('Model selected by AIC and inference ignoring model selection','PostAIC','PostAIC intervals')

} # if linearcomb false

if(linearcomb==TRUE){

  if(!is.null(common)){
    xnewcommon=xnew[,common]
    xnewother= xnew[,-common]
    xnewcommon2<-cbind(xnewcommon,xnewother)
    
    if (intercept==T){    
      xnewfull=cbind(Intercept=1,xnewcommon2)
    }
    else{
      xnewfull=xnewcommon2}}
  
  if(is.null(common)) {
    if (intercept==T){
      xnewfull=cbind(Intercept=1,xnew)
    }
    else{
      xnewfull=xnew
    }
  }
  
  if(class(xnewfull)!='matrix'){
    xnewfull=matrix(unlist(xnewfull),ncol=ncol(Pi))}
  
  results2=Aicpostpred(xp=xnewfull,Pi=Pi,p=pselected,npar=ncol(Pi), model.set=model.set,Information=Jfull,quant=quant,
		nsamp=200000,simplify.constraints=TRUE,morecons=FALSE,show.constraints=FALSE,intercept=intercept)

  quantilespred = unlist(results2[1])

  intervalspred = cbind(xnewfull[,which(Pi[pselected,]>0)]%*%fit.aic$coefficients - 
                         rep(quantilespred,nrow(xnewfull))/sqrt(length(y)), 
				xnewfull[,which(Pi[pselected,]>0)]%*%fit.aic$coefficients + 
				 rep(quantilespred,nrow(xnewfull))/sqrt(length(y)))
   
  outlist=list(summary.fit,results2,intervalspred)
  names(outlist)=c('AIC selected model and inference ignoring model selection','PostAIC',
                    'PostAIC intervals for xnew')

} # if linearcomb true


  
 return(outlist)
  
  
}

Aicpost=function(Pi,p,npar,model.set='user',Information=diag(npar),quant=0.95,nsamp=200000,simplify.constraints=TRUE,morecons=FALSE,unsure=c(),show.constraints=FALSE){
  
  selected=Pi[p,]
  ind=which(selected>0)
  
  if (length(ind)==0){
    
    print('No variable has been selected')
    
  } else {
    
    rtmginput=constraint(Pi=Pi,p=p,modelset=model.set,m.cons.=morecons,simplify.const=simplify.constraints,show.const=show.constraints,Unsure=unsure)
    const=rtmginput[[1]]
    initials=rtmginput[[2]][[1]]
    checkcov=unique(rtmginput[[2]][[2]])
    
    n=nsamp
    M=diag(npar)
    rlin=c(rep(0,npar))
    eachinitial=ceiling(n/nrow(initials))
    samples=c()
    
    if (nrow(initials)>1){
      
      for(i in 1:nrow(initials)){
        
        initial=initials[i,]
        
        seed=npar
        
        while(TRUE){
          
          set.seed(seed)
          sample = rtmg(eachinitial, M, r=rlin, initial, q=const)
          initialsign=sign(initial[checkcov])
          covsign=sign(sample[,checkcov])
          if(all(covsign==matrix(initialsign, ncol=length(initialsign),nrow=eachinitial,byrow=T))==TRUE) break()
          else seed=seed+npar
          
        }
        
        samples=rbind(samples,sample)
        
      }
      
    }else{
      
      seed=npar
      
      set.seed(seed)
      initial=as.vector(initials)
      sample = rtmg(eachinitial, M, r=rlin, initial, q=const)
      posnum=colSums(sample>0)
      negnum=colSums(sample<0)
      
      while(TRUE){
        
        posnum=colSums(sample>0)
        negnum=colSums(sample<0)
        if(all(abs(posnum-negnum)<1000)==TRUE) break()
        
        else{ 
          
          posnegdif=abs(posnum-negnum)
          indp=which(posnegdif>1000)
          indm=which.max(posnegdif[indp])
          ind=indp[indm]
          aa=posnum[ind] ; bb=negnum[ind]
          indm1=which.max(c(aa,bb)[c(1,2)])
          
          if(indm1==1){
            
            pcol=sample[,ind]
            samind=which(pcol>0)[1:posnegdif[ind]]
            samplen=rtmg(posnegdif[ind], M, r=rlin, initial, q=const)
            sample[samind,]=samplen
            
          }else{
            
            pcol=sample[,ind]
            samind=which(pcol<0)[1:posnegdif[ind]]
            samplen=rtmg(posnegdif[ind], M, r=rlin, initial, q=const)
            sample[samind,]=samplen  
            
          }
        }
      }
      samples=sample 
    }
    
    
    
    Jp=Information[ind,ind]
    Jpinv=solve(Jp)
    
    e <- eigen(Jpinv)
    V <- e$vectors
    
    mm=diag(length(ind))
    diag(mm)=sqrt(e$values)
    Jinvhalfun=V %*% mm %*% t(V)
    samplestheta=t(apply(samples[,ind],1,function(x) Jinvhalfun%*%t(t(x))))
    qquant=quant/2+1/2
    quantiles=apply(samplestheta,2,function(x) quantile(x,qquant)[[1]])
    
    simulsamp=t(apply(samples[,ind],1,function(x) t(x)%*%t(t(x))))
    simulquantile=quantile(simulsamp,quant)[[1]]
     
    outlist=list(quantiles,simulquantile,rtmginput[3])
    names(outlist)=c(paste(quant,'quantile',sep=''),paste(quant,'simultaneous quantile',sep=''),'constraints')
    
    names(outlist[[1]])=ind
    
    return(outlist)
    
  }
  
}

Aicpostpred=function(xp,Pi,p,npar,model.set='user',Information=diag(npar),quant=0.95,nsamp=200000,
            simplify.constraints=TRUE,recons=FALSE,morecons=FALSE,show.constraints=FALSE,intercept){
  
  selected=Pi[p,]
  ind=which(selected>0)
    
  if (length(ind)==0){
    
    print('No variable has been selected')
    
  } else {
    
    if(intercept==T & ncol(as.matrix(xp))==(ncol(Pi)-1)) { xp=matrix(unlist(cbind(Intercept=1,xp)),ncol=ncol(Pi)) }

    rtmginput=constraint(Pi=Pi,p=p,modelset=model.set,m.cons.=FALSE,simplify.const=simplify.constraints,show.const=show.constraints,Unsure=c())
    const=rtmginput[[1]]
    initials=rtmginput[[2]][[1]]
    checkcov=unique(rtmginput[[2]][[2]])
    
    n=nsamp
    M=diag(npar)
    rlin=c(rep(0,npar))
    eachinitial=ceiling(n/nrow(initials))
    samples=c()
    
    if (nrow(initials)>1){
      
      for(i in 1:nrow(initials)){
        
        initial=initials[i,]
        
        seed=npar
        
        while(TRUE){
          
          set.seed(seed)
          sample = rtmg(eachinitial, M, r=rlin, initial, q=const)
          initialsign=sign(initial[checkcov])
          covsign=sign(sample[,checkcov])
          if(all(covsign==matrix(initialsign, ncol=length(initialsign),nrow=eachinitial,byrow=T))==TRUE) break()
          else seed=seed+npar
          
        }
        
        samples=rbind(samples,sample)
        
      }
      
    }else{
      
      seed=npar
      
      set.seed(seed)
      initial=as.vector(initials)
      sample = rtmg(eachinitial, M, r=rlin, initial, q=const)
      posnum=colSums(sample>0)
      negnum=colSums(sample<0)
      
      while(TRUE){
        
        posnum=colSums(sample>0)
        negnum=colSums(sample<0)
        if(all(abs(posnum-negnum)<1000)==TRUE) break()
        
        else{ 
          
          posnegdif=abs(posnum-negnum)
          indp=which(posnegdif>1000)
          indm=which.max(posnegdif[indp])
          ind=indp[indm]
          aa=posnum[ind] ; bb=negnum[ind]
          indm1=which.max(c(aa,bb)[c(1,2)])
          
          if(indm1==1){
            
            pcol=sample[,ind]
            samind=which(pcol>0)[1:posnegdif[ind]]
            samplen=rtmg(posnegdif[ind], M, r=rlin, initial, q=const)
            sample[samind,]=samplen
            
          }else{
            
            pcol=sample[,ind]
            samind=which(pcol<0)[1:posnegdif[ind]]
            samplen=rtmg(posnegdif[ind], M, r=rlin, initial, q=const)
            sample[samind,]=samplen  
            
          }
          
        }
        
        samples=sample 
        
      }
      
    }
    
    if (morecons==FALSE){
      
      selected=Pi[p,]
      ind=which(selected>0)
      Jp=Information[ind,ind]
      Jpinv=solve(Jp)
      
      e <- eigen(Jpinv)
      V <- e$vectors
      
      mm=diag(length(ind))
      diag(mm)=sqrt(e$values)
      Jinvhalfun=V %*% mm %*% t(V)
      
      samplepred=xp[,ind]%*%Jinvhalfun%*%t(samples[,ind])
      
    }else{
      
      selected=Pi[which.max(rowSums(Pi)),]
      ind=which(selected>0)
      Jp=Information[ind,ind]
      Jpinv=solve(Jp)
      
      e <- eigen(Jpinv)
      V <- e$vectors
      
      mm=diag(length(ind))
      diag(mm)=sqrt(e$values)
      Jinvhalfun=V %*% mm %*% t(V)
      
      samplepred=xp[,ind]%*%Jinvhalfun%*%t(samples[,ind])
      
    }
    
    qquant=quant/2+1/2
    quantiles=quantile(samplepred,qquant)[[1]]
    
    
    outlist=list(quantiles,rtmginput[3])
    names(outlist)=c(paste(quant,'quantile',sep=''),'constraints')
    
    
    return(outlist)
  }
  
  
}

constraint=function(Pi,p,modelset='user',m.cons.=FALSE,simplify.const=TRUE,show.const=TRUE,Unsure=c()){  
  
  if (modelset=='allsubsets'){   ####### All possible models for all parameters
    
    selected=Pi[p,]
    diffline=ifelse(selected>0, selected, -1)
    diff=diag(diffline)
    constants=rowSums(diff)*(-2)-0.000001
    pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)    
    
  }else if(modelset=='partsubsets'){     ####### All Possible for a subset
    
    if (m.cons.==FALSE){
      
      common=which((colSums(Pi)==nrow(Pi))==TRUE)
      selected=Pi[p,]
      diffline=ifelse(selected>0, selected, -1)
      diff=diag(diffline)[-c(common),]
      constants=rowSums(diff)*(-2)-0.000001
      pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
      
    }else{   #### if Unsure is Null, then more conservative for all parr. else, more conservative for the Unsure parameters
      
      if (is.null(Unsure)==TRUE){
        
        selected=Pi[p,]
        diffline=ifelse(selected>0, selected, -1)
        diff=diag(diffline)
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
        
      }else{
        
        common=which((colSums(Pi)==nrow(Pi))==TRUE)
        selected=Pi[p,]
        diffline=ifelse(selected>0, selected, -1)
        diff=diag(diffline)[-c(common),]
        constants=rowSums(diff)*(-2)-0.000001
        
        for (i in 1:length(Unsure)){
          
          if(selected[Unsure[i]]==1){
            
            attachdiff=rep(0,ncol(Pi))
            attachdiff[Unsure[i]]=1
            diff=rbind(diff,attachdiff)
            
          }
        }
        
        diff=unique(diff,MARGIN=1)
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
        
      }
      
    }
    
  }else if(modelset=='nested'){       ###### for conditional remove lines with negative values
    
    if(m.cons.==FALSE){        ##### inference for nested set
      
      selected=Pi[p,]
      diff=matrix(NA,nrow=nrow(Pi),ncol=ncol(Pi))   ##### difference between selected and all models
      
      diff=t(apply(Pi,1,function(x)selected-x))
      
      diff=matrix(diff[-p,],ncol=ncol(Pi))
      constants=rowSums(diff)*(-2)-0.000001
      pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
      
    }else {     
      
      if (is.null(Unsure)==TRUE){  ##### more conservative inference for non-common parameters
        
        common=which((colSums(Pi)==nrow(Pi))==TRUE)
        selected=Pi[p,]
        
        diff=matrix(NA,nrow=nrow(Pi),ncol=ncol(Pi))   ##### difference between selected and all models
        
        diff=t(apply(Pi,1,function(x)selected-x))
        
        diff=matrix(diff[-p,],ncol=ncol(Pi))
        
        for(i in 1:length(common)){
          
          difflinei=rep(0,ncol(Pi))
          difflinei[common[i]]<-1
          diff=rbind(diff,difflinei)
          rownames(diff)=c()
          
        }
        
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
        
      }else{  ##### more conservative inference for Unsure parameters
        
        selected=Pi[p,]
        diff=matrix(NA,nrow=nrow(Pi),ncol=ncol(Pi))   ##### difference between selected and all models
        
        diff=t(apply(Pi,1,function(x)selected-x))
        
        diff=matrix(diff[-p,],ncol=ncol(Pi))
        common=which((colSums(Pi)==nrow(Pi))==TRUE)
        
        for (i in 1:length(Unsure)){
          
          if(selected[Unsure[i]]==1){
            
            if (Unsure[i]%in% common){
              
              attachdiff=rep(0,ncol(Pi))
              attachdiff[Unsure[i]]=1
              diff=rbind(attachdiff,diff)            
              rownames(diff)=c()
              
            }else{
              
              attachdiff=rep(0,ncol(Pi))
              attachdiff[Unsure[i]]=1
              diff=rbind(attachdiff,diff)            
              rownames(diff)=c()
              print('The unfixed unsure parameter(s) may lead to false confidence for the other parameter(s)')
              
            }
            
          }
          
        }
        
        diff=unique(diff,MARGIN=1)
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=simplify.const,silent=!show.const)
        
        if (simplify.const==TRUE){
          
          diff=diff[pconst[[3]],]
          constants=constants[pconst[[3]]]
          
        }
        
      } 
      
    }
    
  }else {
    
    if(m.cons.==FALSE){        ##### inference for nested set
      
      selected=Pi[p,]
      diff=matrix(NA,nrow=nrow(Pi),ncol=ncol(Pi))   ##### difference between selected and all models
      
      diff=t(apply(Pi,1,function(x)selected-x))
      
      diff=matrix(diff[-p,],ncol=ncol(Pi))
      constants=rowSums(diff)*(-2)-0.000001  
      pconst=print.orig.const04(diff,constants,simplify.const=simplify.const,silent=!show.const)
      
      if (simplify.const==TRUE){
        
        diff=diff[pconst[[3]],]
        constants=constants[pconst[[3]]]
        
      }
      
    }else{
      
      if (is.null(Unsure)==TRUE){
        
        common=which((colSums(Pi)==nrow(Pi))==TRUE)
        selected=Pi[p,]
        diffline=ifelse(selected>0, selected, -1)
        
        if (length(common)>0){
          diff=diag(diffline)[-c(common),]
        }else{
          diff=diag(diffline)
        }
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=FALSE,silent=!show.const)
        
      }else {
        
        selected=Pi[p,]
        diff=matrix(NA,nrow=nrow(Pi),ncol=ncol(Pi))   ##### difference between selected and all models
        
        diff=t(apply(Pi,1,function(x)selected-x))
        
        diff=matrix(diff[-p,],ncol=ncol(Pi))
        
        for (i in 1:length(Unsure)){
          
          if(selected[Unsure[i]]==1){
            
            attachdiff=rep(0,ncol(Pi))
            attachdiff[Unsure[i]]=1
            diff=rbind(attachdiff,diff)            
            
          }
          
        }
        
        diff=unique(diff,MARGIN=1)
        constants=rowSums(diff)*(-2)-0.000001
        pconst=print.orig.const04(diff,constants,simplify.const=simplify.const,silent=!show.const)
        
        if (simplify.const==TRUE){
          
          diff=diff[pconst[[3]],]
          constants=constants[pconst[[3]]]
          
        }
        
      }
      
    }
    
  }
  
  Blinear=c(rep(0,ncol(Pi)))                    ##### No linear constraint  
  diff=matrix(diff,ncol=ncol(Pi))
  totconstlist=vector('list', nrow(diff))       ##### M matrices for constraints
  
  for(i in (1:nrow(diff))){
    
    totconstlist[[i]]=diag(diff[i,])
    
  }
  
  compconstlist=vector('list', nrow(diff))      ##### list of all elements of constraints
  
  for(j in 1:nrow(diff)){
    
    compconstlist[[j]]= list(totconstlist[[j]],Blinear,constants[j])
    
  }
  
  initialcov=initialvalues(diff,sel=selected)
  
  
  tmginputs=list(compconstlist,initialcov,pconst)
  return(tmginputs)
  
}

filtercombn=function(com,orderc,testleve){
  
  firstsum=com%*%orderc
  indexpos=which(firstsum%in%testleve)
  return(indexpos)
  
}

possiblecombn = function(adeq,orderb,uniqueleveltest){
  
  if (length(adeq)>1){
    
    vecs <- mapply(seq, 0, adeq,SIMPLIFY = FALSE)
    tmp0 <- as.matrix(do.call(expand.grid, vecs))
    index <- filtercombn(com=tmp0,orderc=orderb,testleve=uniqueleveltest)
    poscombn=matrix(tmp0[index,],ncol=nrow(adeq))
    colnames(poscombn)=orderb
    return(poscombn)
    
  }else{
    
    vecs <- mapply(seq, 0, adeq)
    index <- filtercombn(com=vecs,orderc=orderb,testleve=uniqueleveltest)
    poscombn=matrix(vecs[index,],ncol=nrow(adeq))
    colnames(poscombn)=orderb
    return(poscombn)
    
  }
  
}

initialvalues=function(diff,sel){
  
  absdiff=abs(diff)
  
  strictconstraints=which(rowSums(absdiff)==1)     ###### models with strict constraints (boundries for covariates)
  positivestrconst=which(diff[strictconstraints,]==1,arr.ind=T)#[,2] ###### Z^2_i > 2
  
  if(length(positivestrconst)>1){
    
    posstr=positivestrconst[,2]
    
  }else{
    
    if (length(positivestrconst)==1){
      
      posstr=positivestrconst
      
    }else
      
      posstr=c(NULL)
    
  }
  
  posnegoneeach= which(rowSums(diff)==0 & rowSums(absdiff)==2)
  posnegind=which(diff[which(rowSums(diff)==0 & rowSums(absdiff)==2),]==1,arr.ind=T)#[,2]  
          ###### Z^2_i>Z^2_j
  #dublenegpos=which(diff[posnegoneeach[posnegind],]==1,arr.ind=T)
  if (length(posnegind)>1){
    
    posneg=posnegind[,2]
    
  }else{
    
    if (length(posnegind)==1){
      
      posneg=posnegind
      
    }else{
      
      posneg=c(NULL)
    }
    
  }
  
  posduble=c(posstr,posneg)
  
  ###### haminja duble ro ba posneg combine kon
  
  if (length(posduble)>0){
    
    incombinations=rep(c(2,-2),length(unique(posduble)))        
    initialcombinations=unique(combn(incombinations,length(unique(posduble))),MARGIN=2)   ##### initial combinations
    
    initialmatrix=matrix(0,nrow=ncol(initialcombinations),ncol=ncol(diff))              ##### matrix of initials
    selectedcov=sel
    initialmatrix[,selectedcov==1]<-2
    
    for(i in 1:ncol(initialcombinations)){
      
      initialmatrix[i,unique(posduble)]=initialcombinations[,i]
    }
    
    
  }else{
    
    selectedcov=sel  #### find covariates that are common to all models 
    initialmatrix=matrix(0,nrow=1,ncol=ncol(diff))
    initialmatrix[,selectedcov==1]<-2
    
  }
  
  if (length(posduble)>0){
    
    initialcov=list(initialmatrix,sort(unique(posduble)))
    
  }else{
    
    initialcov=list(initialmatrix,c()) 
    
  }
  
  return(initialcov)
  
}  

simplify.const04=function(diff){
  
  rownames(diff)=c(1:nrow(diff))
  base1=c()
  base=c()
  
  if (nrow(diff)>2){
    
    leveldiff=rowSums(abs(diff))  ####### number of nonzero values
    
    firstminleveldiff=names(which((leveldiff==min(leveldiff))==T))  ###### smallest numbers of nonzero 
    
    if (length(firstminleveldiff)>1){
      
      base=diff[firstminleveldiff,]
      minleveldiff=firstminleveldiff
      ncommon=length(which((colSums(abs(diff))==0)==TRUE))
      
      if(length(firstminleveldiff)==(ncol(diff)-ncommon) & all(rowSums(abs(base[firstminleveldiff,]))==1)){
        
        base1=base
        return(rownames(base1))
        
      }
      
    }else{
      
      levelremdiff=matrix(abs(diff[rownames(diff)%in%firstminleveldiff==F,]),ncol=ncol(diff))
      levelremdiff=sort(rowSums(abs(diff[rownames(diff)%in%firstminleveldiff==F,])))
      secondminleveldiff=names(which((levelremdiff==min(levelremdiff))==T))   ##### second smallest number of nonzero
      minleveldiff=c(firstminleveldiff,secondminleveldiff)
      base=diff[minleveldiff,]
      
    }
    
  }else{
    
    base1=diff
    return(rownames(base1))
    
  }
  
  remdiffloopnames=rownames(diff)[rownames(diff)%in%minleveldiff==F]
  #remdiffloopnames= rownames(matrix(diff[rownames(diff)%in%minleveldiff==F,],ncol=ncol(diff)))
  remaindiffloop=matrix(diff[rownames(diff)%in%minleveldiff==F,],ncol=ncol(diff))   ##### all const. for checking
  rownames(remaindiffloop)=remdiffloopnames
  levelremdiff=sort(rowSums(abs(remaindiffloop)))           ##### level of models for checking
  
  while (nrow(remaindiffloop)!= 0){
    
    if(nrow(remaindiffloop)>1){
      
      leveltest=sort(rowSums(abs(remaindiffloop)))
      
    }else{
      
      leveltest=sum(abs(remaindiffloop))
      
    }
    
    baseorder=rowSums(abs(base))
    
    if (sum(baseorder) < min(leveltest)){
      
      attnames=names(which((rowSums(abs(remaindiffloop))==min(rowSums(abs(remaindiffloop))))==T))
      remainname1=rownames(remaindiffloop[!rownames(remaindiffloop)%in%attnames,])
      baseatt=remaindiffloop[attnames,]
      base=rbind(base,baseatt)
      rownames(base)[(nrow(base)- length(attnames)+1):nrow(base)]=attnames
      remaindiffloop=matrix(remaindiffloop[remainname1,],ncol=ncol(diff))
      rownames(remaindiffloop)=remainname1
      
    }
    
    maxleveltest=max(leveltest)
    uniqleveltest=unique(leveltest)
    justbaseorder=matrix(as.numeric(names(table(baseorder))),nrow=length(names(table(baseorder))))
    quantbaseorder=matrix(table(baseorder),ncol=length(names(table(baseorder))))
    
    adeqcombn=floor(pmin(maxleveltest/justbaseorder,quantbaseorder))
    
    posblcombn=as.matrix(possiblecombn(adeq=adeqcombn,orderb=justbaseorder,uniqueleveltest=uniqleveltest))
    
    baseordernames=lapply(justbaseorder, function(x) names(baseorder[(baseorder==x)==T]))
    names(baseordernames)=justbaseorder
    
    if (nrow(posblcombn)!=0){
      
      for (i in 1:nrow(posblcombn)){
        
        combin=matrix(posblcombn[i,],ncol=ncol(posblcombn))
        colnames(combin)=colnames(posblcombn)
        
        makecombin=vector('list',ncol(combin))
        names(makecombin)=colnames(combin)
        
        for (j in colnames(combin)){
          
          makecombin[names(makecombin)%in%j][[1]]=combn(baseordernames[names(baseordernames)%in%j][[1]],combin[,j])
          
        }
        
        if(length(colnames(combin))==1){
          
          potentialcomb0=matrix(unlist(makecombin),ncol=ncol(makecombin[[1]]))
          
        }else if(length(colnames(combin))==2){
          
          potentialcomb0=matrix(unlist(apply(makecombin[[1]],2,function(x) apply(makecombin[[2]],2, function(y) c(x,y)))),ncol=ncol(makecombin[[1]])*ncol(makecombin[[2]]))
          
        }else{
          
          potentialcomb0=matrix(apply(makecombin[[1]],2,function(x) apply(makecombin[[2]],2, function(y) c(x,y))),ncol=ncol(makecombin[[1]])*ncol(makecombin[[2]]))
          
          for (k in 3:length(names(makecombin))){
            
            potentialcomb0=matrix(unlist(apply(potentialcomb0,2,function(x) apply(makecombin[[k]],2, function(y) c(x,y)))),ncol=ncol(potentialcomb0)*ncol(makecombin[[k]]))      
            
          }
          
          
          
        }
        
        if (ncol(potentialcomb0)>1){
          
          sumofpotentials=t(apply(potentialcomb0,2,function(x) colSums(base[x,])))
          
          redund=which(sumofpotentials>1 | sumofpotentials<(-1),arr.ind=T)[,1]
          
          if (length(redund)>0){
            
            superpotentials=matrix(sumofpotentials[-redund,],ncol=ncol(diff))
            
          }else{
            
            superpotentials=sumofpotentials
            
          }
          
          matchedindex=c()
          
          if (nrow(superpotentials)>1){
            
            for (l in 1:nrow(remaindiffloop)){
              
              matching=t(apply(superpotentials,1,function(x) x-remaindiffloop[l,]))
              matched=any(apply(matching,1,function(x) all(x==0))==T)
              
              if (matched==TRUE){
                
                matchedindex=c(matchedindex,l)
                #  baseattach=rbind(baseattach,remaindiffloop[l,],make.row.names = TRUE)
              }
              
            }
            
          }else
            
            if(nrow(superpotentials)==1){
              
              for (l in 1:nrow(remaindiffloop)){
                
                matching=superpotentials-remaindiffloop[l,]
                matched=all(matching==0)
                
                if (matched==TRUE){
                  
                  matchedindex=c(matchedindex,l)
                  
                }
                
              }
              
            }else{
              
              matchedindex=c()
              
            }
          
          if (length(matchedindex)>0){
            
            rmnames=rownames(remaindiffloop)          
            matchedname=rmnames[matchedindex]
            remainname=rmnames[!rmnames%in%matchedname]
            remaindiffloop=matrix(remaindiffloop[remainname,],ncol=ncol(diff))
            rownames(remaindiffloop)=remainname
            
          }
          
        }else if(ncol(potentialcomb0)==1){
          
          superpotentials=colSums(base[potentialcomb0,])
          
          matchedindex=c()
          
          for (ll in 1:nrow(remaindiffloop)){
            
            matching=superpotentials-remaindiffloop[ll,] 
            matched=all(matching==0)==T
            
            if (matched==TRUE){
              
              matchedindex=c(matchedindex,ll)
              
            }
            
          }
          
          if (length(matchedindex)>0){
            
            rmnames=rownames(remaindiffloop)          
            matchedname=rmnames[matchedindex]
            remainname=rmnames[!rmnames%in%matchedname]
            remaindiffloop=matrix(remaindiffloop[remainname,],ncol=ncol(diff))
            rownames(remaindiffloop)=remainname
            
          }
          
        }
        
        if ((nrow(remaindiffloop)==0)==TRUE) break
        
      }
      
      
      
      if (nrow(remaindiffloop)>1){
        
        attachindex=names(which((rowSums(abs(remaindiffloop))==min(rowSums(abs(remaindiffloop))))==T))
        diffnames=rownames(remaindiffloop)
        remaindiffindex=diffnames[!diffnames%in%attachindex]
        baseattach=remaindiffloop[attachindex,]
        base=rbind(base,baseattach)
        rownames(base)[(nrow(base)- length(attachindex)+1):nrow(base)]=attachindex
        remaindiffloop=matrix(remaindiffloop[remaindiffindex,],ncol=ncol(diff))
        rownames(remaindiffloop)=remaindiffindex
        
      }else{
        
        base=rbind(base,remaindiffloop)
        
        if(nrow(remaindiffloop)==1){
          
          rownames(base)[nrow(base)]=rownames(remaindiffloop)
          remaindiffloop=matrix(remaindiffloop[-1,],ncol=ncol(diff))
          
        }
        
      }
      
    }else{
      
      if (nrow(remaindiffloop)>1){
        
        attachindex=names(which((rowSums(abs(remaindiffloop))==min(rowSums(abs(remaindiffloop))))==T))
        diffnames=rownames(remaindiffloop)
        remaindiffindex=diffnames[!diffnames%in%attachindex]
        baseattach=remaindiffloop[attachindex,]
        base=rbind(base,baseattach)
        rownames(base)[(nrow(base)- length(attachindex)+1):nrow(base)]=attachindex
        remaindiffloop=matrix(remaindiffloop[remaindiffindex,],ncol=ncol(diff))
        rownames(remaindiffloop)=remaindiffindex
        
      }else{
        
        base=rbind(base,remaindiffloop)
        
        if(nrow(remaindiffloop)==1){
          
          rownames(base)[nrow(base)]=rownames(remaindiffloop)
          remaindiffloop=matrix(remaindiffloop[-1,],ncol=ncol(diff))
          
        }
        
      }
      
    }
    
  }
  
  base1=base
  
  
  
  return(rownames(base1)) 
  
}

print.orig.const04=function(diff,constants,simplify.const=TRUE,silent=TRUE){
  
  obj=matrix(NA,nrow=nrow(diff))
  for(i in 1:nrow(diff)){
    obj[i]=round(constants[i],1)
  }
  obj=replace(obj,obj==0,values=c(' ') )
  for(i in 1:nrow(diff)){
    for(j in 1:length(diff[i,])) {
      if (diff[i,j]>0) obj[i]=c(paste(obj[i],paste('+(z',j,sep='',collapse=''),')^2',sep='',collapse=''))
      if (diff[i,j]<0) obj[i]=c(paste(obj[i],paste('-(z',j,sep='',collapse=''),')^2',sep='',collapse=''))
    }
  }
  
  for(i in 1:nrow(obj)){
    obj[i]=paste(obj[i,],'>' ,'0',sep='',collapse='')
  }
  
  if(simplify.const==TRUE){
    
    reducedconst=simplify.const04(diff)
    reducedobj=matrix(obj[as.numeric(reducedconst),])
    
  }
  
  if(silent==FALSE){
    
    print('original constraints')
    for(i in 1:nrow(obj)){
      cat(paste('(',i,')',obj[i,],sep='',collapse=''),sep='\n')
      
    }
    
    cat(c(),sep='\n')
    
    if(simplify.const==TRUE){
      
      print('reduced constraints')
      for(i in 1:nrow(reducedobj)){
        cat(paste('(',i,')',reducedobj[i,],sep=' ',collapse=''),sep='\n')
        
      }
      
      cat(c(),sep='\n')
      print('Indices of necessary and sufficient constraints from originals')
      cat(sort(as.numeric((reducedconst))))
      
    }
    
    #returns=c()
    
  }
  
  if(simplify.const==TRUE){
    
    allconstraints=list(obj[order(nchar(obj))],reducedobj[order(nchar(reducedobj))],sort(as.numeric(reducedconst)))
    names(allconstraints)=c('all constraints','reduced constraints','sufficeint indicies')
    
    return(allconstraints)
    
  }else{
    
    allconstraints=list(obj[order(nchar(obj))])  
    names(allconstraints)='all constraints'
    
    return(allconstraints)
    
  }
  
}

combinations= function(n){
  comb = NULL
  for(i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
  return(comb)
}
