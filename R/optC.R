

optC<-function(mu, Sin){
  optcont<-Sin%*%(mu-sum(mu*rowSums(Sin))/sum(Sin))
  contMat<-(optcont-sum(optcont))/sqrt(sum(optcont^2))
}
