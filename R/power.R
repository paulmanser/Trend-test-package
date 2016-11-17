#' Calculate power for multiple contrast tests
#'
#' Calculate power for a multiple contrast test for a specified trend type.
#' @param mu, altModel mu is mean values in all doses groups including placebo. It could be a vector of arbitrary values or an object from CPL. Need to be a
#'           matrix with rows equal to number of doses. altModel is an object of class "Mods", see Mods in DoseFinding package for detail. It is used to define mean values.
#'           Either mu or altModel needs to be specified. Specify mu if you want to an arbitrary mu or an object from CPL, or specify altModel if you want to use an object from MCP Mods.
#' @param n, sigma, S n is a numeric vector of sample sizes in all doses groups including placebo. Either a vector n and sigma or S need to be specified it is assumed computation are made
#'          for a normal homoscedastic model with group sample sizes given by n and residual standard deviation sigma, i.e. the covariance matrix used for estiamtes is thus sigma^2*diag(1/n)
#'          and the degrees of freedom are caculated as sum(n)- length(mu). When S is specified this will be used as convariance matrix for the estimates.
#' @param df degrees of freedom n, need to specify when S is given.
#' @param doses a numeric vector of specified doses. When this argument is missing, doses are the object in the specified model if model is given for MCP Mod trend test, or doses are 0 to k by 1
#'        for all other trend test. You can assign less levels of doses. mu will be adjusted correspondingly.
#' @param trend a single character string, determining the trend for the multiple contrast trend test, one of "CPL", "MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear".
#' Note that when "MCP_Mod" is specified for trend, model needs to be specified. When "Arbitrary" is specified, muMat needs to be specified.
#' @param model specify models if you choose trend to be MCP_Mod trend test
#' @param alternative a single character string, specifying the alternative for the multiple contrast trend test. If this argument is missing, alternative="one.sided".
#' @param alpha significance level to use. If this argument is missing, alpha=0.05.
#' @return true mean, optimal contrast, mean matrix, power.
#' @export
#' @examples N<-150
#' n0<-30
#' n1<-30
#' k<-5
#' n<-c(n0, rep(n1,k))
#' transdose<-function(F, k, T){
#' dose<-F^(seq(0, k-1, by=1)/(k-1)	)
#' log(1+T*dose)
#' }
#'
#'transDose<-transdose(F=50, k=k, T=2)
#'
#'
#' models<-Mods(emax=c(log(1+2*0.01),log(1+2*0.25),log(1+2*4)),
#'           linear=NULL,
#'           exponential=c(log(1+2*0.25), log(1+2*0.75), log(1+2*3.5)),
#'           sigEmax=rbind(c(log(1+2*1),25), c(log(1+2*5),25), c(log(1+2*15),25),c(log(1+2*1),5), c(log(1+2*5),5), c(log(1+2*15),5) ),
#'           doses=c(0, transDose), placEff=0, maxEff=0.7)
#'
#'
#'MCP_models<-Mods(emax=c(log(1+2*0.01),log(1+2*0.25),log(1+2*4)),
#'                 linear=NULL,
#'                 exponential=c(log(1+2*0.25), log(1+2*0.75), log(1+2*3.5)),
#'                 sigEmax=rbind( c(log(1+2*1),5), c(log(1+2*5),5), c(log(1+2*15),5) ),
#'                 doses=c(0,transDose), placEff=0, maxEff=0.7)
#'
#'
#'trendpow<-trendpowfn(altModel=models,n=rep(30,6), sigma=1,doses=c(0,transDose), alpha=0.05, trend="MCP Mod", model=MCP_models)
#'
#'mufrmCPL<-t(mu.piecewise(5,doses=transDose))
#'
#'trendpow<-trendpowfn(mu=mufrmCPL,n=rep(30,6), sigma=1,doses=c(0,transDose), alpha=0.05, trend="MCP Mod", model=MCP_models)


trendpowfn<-function( mu=NULL, altModel=NULL,  n=NULL, sigma=NULL, S=NULL, df=NULL,doses=NULL,trend=c("CPL","MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear"), model=NULL, alternative="one.sided",alpha=0.05){

  #mu: used for arbitrary mu or object from MCP Mod
  #altModel: an object from MCP Mod
  if(trend[1]=="MCP Mod"){
    if (is.null(model))
      stop("model needs to specified")
  } else {
    if(!is.null(model))
      stop ("when model is specified, trend needs to be MCP Mod")
  }

  ## extract covariance matrix
  if(is.null(S)){
    if(is.null(n) | is.null(sigma))
      stop("Either S or n and sigma need to be specified")
    #    if(length(n) != nD)
    #      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2*diag(1/n)
    df<-as.integer(sum(n)-ncol(S))
  } else {
    if(!is.null(n)|!is.null(sigma))
      stop("Need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
    if(nrow(S) != ncol(S))
      stop("S needs to be a square matrix")
    # if(nrow(S) != nD)
    #   stop("S needs to have as many rows&cols as there are doses")
    if(missing(df))
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }


  #### when use an object from CPL with less doses, need doses levels

  if(is.null(mu)){
    if(is.null(altModel)){
      stop("Either mu or altModel needs to be specified")
    }else{
      #if altModel is given
      if(!is.null(doses)){
        #doses is given
        if(length(which(attributes(altModel)$doses%in%as.character(doses)))!=length(attributes(altModel)$doses))
          stop("The specified doses is not the same as the object in given altModel ")
      }
      doses<-attributes(altModel)$doses
      print("Use doses in the specified model")
      mu<-getResp(altModel)
      nD<-nrow(mu)
      k<-nD-1
    }
  }else{
    # mu is given
    if (!is.null(altModel)){
      stop("Need to specify only one of \"mu\" or \"altModel\"")
    }else{
      # print("An arbitrary mu or object from CPL is given")
      # print("mu need to be a matrix with nrow= number of doses")
      if(is.null(doses)){
        ### dose is not given
        nD<-nrow(mu)  # of doses
        k<-nD-1
        doses<-seq(0,k)
        print("Use default doses from 0 to k by 1")
      }else{
        nD<-length(doses)
        k<-nD-1
        clnm<-rownames(mu)
        idx<-which(clnm%in%as.character(doses))
        mu<-matrix(mu[idx,], nrow=length(idx))
        S<-S[idx, idx]
       # df<-as.integer(sum(n[idx])-nD)
        if(is.null(n) & is.null(sigma)){
          df<-df
        }else{
          df<-as.integer(sum(n[idx])-nD)
        }
      }
    }
  }







  ####prepare muMat and contMat for williams
  cont.Mat.Will<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat.Will<-matrix(0, nrow=k, ncol=k+1)
  for(j in 1:k){

    Sin<-solve(S)[c(1,((k-j+2):(k+1))),c(1,((k-j+2):(k+1)))]

    #Sin<-diag(1/(diag(S)[c(1,((k-j+2):(k+1)))]))
    mu.Mat.Will[j,]<-c(0, rep(NA, k-j), rep(1,j))
    mu1<-c(0, rep(1,j))
    cont.Mat.Will[j,c(1,((k-j+2):(k+1)))]<-optC(mu1, Sin)
  }
  #####prepare muMat and contMat for Aug-Williams
  cont.Mat1<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat1<-matrix(0, nrow=k, ncol=k+1)
  for(j in 1:k){
    Sin<-solve(S)[c(1:j,k+1),c(1:j,k+1)]
    #	Sin<-diag(1/(diag(S)[c(1:j,k+1)]))
    mu1<-c(rep(0,j),1)
    mu.Mat1[j,]<-c(rep(0,j), rep(NA, k-j), 1)
    cont.Mat1[j,c(1:j,k+1)]<-optC(mu1, Sin)
  }
  cont.Mat2<-matrix(0, nrow=k-1, ncol=k+1)
  mu.Mat2<-matrix(0, nrow=k-1, ncol=k+1)
  for(j in 1:(k-1)){

    Sin<-solve(S)[c(1,((k-j+1):(k+1))), c(1,((k-j+1):(k+1)))]
    #Sin<-diag(1/(diag(S)[c(1,((k-j+1):(k+1)))]))

    mu1<-c(0, rep(1, j+1))
    mu.Mat2[j,]<-c(0, rep(NA, k-j-1),rep(1, j+1))
    cont.Mat2[j,c(1,((k-j+1):(k+1)))]<-optC(mu1, Sin)
  }
  mu.Mat.AugWill<-rbind(mu.Mat1, mu.Mat2)
  cont.Mat.AugWill<-rbind(cont.Mat1, cont.Mat2)
  ####prepare muMat and contMat for Dunnett
  cont.Mat.Dun<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat.Dun<-matrix(NA, nrow=k, ncol=k+1 )
  for (j in 1:k){
    #Sin<-solve(S)[c(1, j+1),c(1, j+1)]
    Sin<-diag(1/(diag(S)[c(1, j+1)]))
    mu1<-c(0,1)
    mu.Mat.Dun[j, c(1, (j+1))]<-c(0,1)
    cont.Mat.Dun[j,c(1, j+1)]<-optC(mu1, Sin)
  }


  muMat <- switch(trend[1],
                  "MCP Mod"=as.matrix(t(getResp(model))),
                  "Linear"=mu.linear(k,doses),
                  "Williams"=mu.Mat.Will,
                  "Aug Williams"=mu.Mat.AugWill,
                  "Dunnett"=mu.Mat.Dun,
                  "CPL"=mu.piecewise(k,doses)
  )

  contMat <- switch(trend[1],
                    "MCP Mod"=optContr(model, doses=doses, S=S)$contMat,
                    "Linear"=apply(muMat,1, optC,Sin=solve(S)),
                    "Williams"=t(cont.Mat.Will),
                    "Aug Williams"=t(cont.Mat.AugWill),
                    "Dunnett"=t(cont.Mat.Dun),
                    "CPL"=apply(muMat,1, optC, Sin=solve(S))
  )





  cont<-contMat
  nC<-ncol(cont)  # of contrasts


  #calculate non-centrality parameter
  deltaMat<-t(cont)%*%mu
  covMat<-t(cont)%*%S%*%cont
  den<-sqrt(diag(covMat))
  deltaMat<-deltaMat/den
  if(alternative == "two.sided"){
    deltaMat <- abs(deltaMat)
  }

  diag(covMat)<-diag(covMat)+min(diag(covMat))*0.0001
  #corMat<-cov2cor(covMat+min(diag(covMat))*0.000001)
  corMat<-cov2cor(covMat)

  tail<-ifelse(alternative == "two.sided",
               "both.tails", "lower.tail")
  critV<-qmvt(1-alpha, tail=tail, df=df,delta=0, corr=corMat,algorithm=GenzBretz() )$quantile
  if(alternative == "two.sided"){
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper<-rep(critV, nC)

  pow<-rep(NA,ncol(deltaMat))
  for (i in 1: ncol(deltaMat)){
    pow[i]<-1-pmvt(lower=lower, upper=upper,df=df,corr=corMat, delta=deltaMat[,i], algorithm=GenzBretz())
  }
  mylist<-list(mu, contMat, muMat, pow)
  names(mylist)<-c("true_mu", "contrast_matrix", "mu_matrix", "power")
  mylist
}



