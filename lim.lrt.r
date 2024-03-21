### Left-Inflated Mixture model Likelihood Ratio Test (LIM-LRT)
#
# by Andreas Gleiss (2014)
#
# andreas.gleiss@meduniwien.ac.at
#

plot.lod <- function(row=1, xmax=25,ymax=1, xlabel="log2-intensity",
                      data,group,massat=-8, dl=0) {
  
  #
  # Plots histograms of both groups including bar for pointmass values
  #
  # row       if a single number then interpreted as row number of data matrix
  #           otherwise vector containing the data for both groups
  # xmax      maximum value at x axis (intensity)
  # ymax      maximum value at y axis (rel. frequency)
  # xlabel    label of x axis
  # data      data matrix (only used if row = single number)
  # group     vector of length of row (if >1) or of row length in data matrix
  # massat    numerical code indicating pointmass values in the data
  # dl        detection limit where the pointmass value bar is to be placed
  #
  
  breaks=seq(dl,xmax,1) # works only for breaks in unit steps!!! fault of hist() function
  gr0<-(group==0)
  gr1<-(group==1)
  if(length(row)==1) {
    test<-data[row,]  
    test[test==massat]<-dl
    tit0<-paste("Prot",row,"(group 0)")
    tit1<-paste("Prot",row,"(group 1)")
    brs<-c(dl-(breaks[2]-breaks[1]),dl+breaks)
    }
  else { # row contains data
    test<-row
    test[test==massat]<-dl
    tit0<-paste("(group 0)")
    tit1<-paste("(group 1)")
    brs<-c(dl-(breaks[2]-breaks[1]),dl+breaks)
    }
  # use myhist() instead of hist() to show correct y-axis scaling!!!
  par(mfrow=c(2,1))
  hist(test[gr0],xlim=c(dl-(breaks[2]-breaks[1]),xmax),ylim=c(0,ymax),breaks=brs, freq=F,
        main=tit0, xlab=xlabel,ylab="rel.frequency",
        col=c("grey",rep(0,length(breaks))))
  hist(test[gr1],xlim=c(dl-(breaks[2]-breaks[1]),xmax),ylim=c(0,ymax),breaks=brs, freq=F,
        main=tit1, xlab=xlabel,ylab="rel.frequency",
        col=c("grey",rep(0,length(breaks))))
  par(mfrow=c(1,1))
  }

lim0 <- function (arg,z,d,w,h,lod=0,quickret=100000,print=F) {
  #
  # internal function used by lim.lrt
  # delivers value of the likelihood under null
  #
  
  mu<-arg[1]
  sd1<-arg[2]
  sd2<-arg[3]
  p<-arg[4]
  Phi1<-pnorm(lod,mu,sd1)
  Phi2<-pnorm(lod,mu,sd2)

  if(print) cat(arg,Phi1,Phi2,"\n")
  if(sd1<=0 | sd2<=0) {
    return(quickret)
    }
  if(Phi1>0.999 | Phi2>0.999 | p+Phi1>0.999 | p+Phi2>0.999 | p+Phi1<=0 | p+Phi2<=0) { # account for numeric cases of "near miss" (despite marginal conditions)
    return(quickret)
    }

  # general case
  erg<- -((length(z)-sum(d))*log(p+(1-p)*Phi1)+
      sum(d)*log(1-p)-
      log(sd1)*sum(d)-
      sum(d*(z-mu)^2/(2*sd1^2)))-
      ((length(w)-sum(h))*log(p+(1-p)*Phi2)+
      sum(h)*log(1-p)-
      log(sd2)*sum(h)-
      sum(h*(w-mu)^2/(2*sd2^2)))
  if (is.na(erg) | erg==Inf | erg==-Inf) {
    cat("ERROR (H0): p =",p,", Phi1 =",Phi1,", Phi2 =",Phi2,"\n")
    print(z)
    print(d)
    print(w)
    print(h)
    }
  return(erg)
  }
 
lim1 <- function (arg,z,d,w,h,lod=0,print=F, quickret=100000) {
  #
  # internal function used by lim.lrt
  # delivers value of the likelihood under alternative
  #
  
    mu1<-arg[1] # group1
    sd1<-arg[2]
    mu2<-arg[3] # group0
    sd2<-arg[4]
    p1<-arg[5]
    p2<-arg[6]
    Phi1<-pnorm(lod,mu1,sd1)
    Phi2<-pnorm(lod,mu2,sd2)

  if(print) cat(arg,Phi1,Phi2,"\n")
  if(sd1<=0 | sd2<=0) {
    return(quickret)
    }
  if(Phi1>0.999 | Phi2>0.999 | p1+Phi1>0.999 | p2+Phi2>0.999 | p1+Phi1<=0 | p2+Phi2<=0) { # account for numeric cases of "near miss" (despite marginal conditions)
    return(quickret)
    }

  # general case
  erg<- -((length(z)-sum(d))*log(p1+(1-p1)*Phi1)+
      sum(d)*log(1-p1)-
      log(sd1)*sum(d)-
      sum(d*(z-mu1)^2/(2*sd1^2)))-
      ((length(w)-sum(h))*log(p2+(1-p2)*+Phi2)+
      sum(h)*log(1-p2)-
      log(sd2)*sum(h)-
      sum(h*(w-mu2)^2/(2*sd2^2)))
  if (is.na(erg) | erg==Inf | erg==-Inf) {
    cat("ERROR (H1): p1 =",p1,", p2 =",p2,", Phi1 =",Phi1,", Phi2 =",Phi2,"\n")
    print(z)
    print(d)
    print(w)
    print(h)
    }

  return(erg)
  }

tobit.single<-function(data, point.mass=0, dl=0) {
  #
  # internal function used by lim.lrt
  # estimation of Tobit model
  #
  library(survival)
  data[data==point.mass]<-dl
  survdat<-Surv(as.vector(data),as.vector(data!=dl),type='left')
  survfit<-survreg(survdat ~ 1, dist='gaussian')
  result<-list(mean=survfit$coefficients[1],sd=survfit$scale)
  }

lim.lrt <- function (data, group, massat=-8, quickret=100000, pprint=F,
                            mu_low=-Inf,sd_low=.1,print=F, eps=0.1,
                            mu_start=1,sd_start=1,preest="none") {
  #
  # R-function for Left-inflated mixture model likelihood ratio test (LIM-LRT)
  #
  # ARGUMENTS:
  # data      vector containing data of both groups
  # group     vector of same length containing group codes (1 and 0)
  # massat    numerical code indicating pointmass values in the data vector 
  # quickret  value of likelihood to be returned by internal sub-functions in case of 'illegal' parameters
  # pprint    if TRUE then details about iterations are output
  # mu_low    lower limit for point estimate of continuous subdistribution 
  # sd_low    lower limit for standard deviation of continuous subdistribution
  # print     if TRUE then results are output
  # eps       margin of limit of detection below data minimum
  # mu_start  starting value for point estimates in case of preest="none" and degenarated cases of other preest options  
  # sd_start  starting value for standard deviation in case of preest="none" and degenarated cases of other preest options
  # preest    "none"  => use given starting values mu_start and sd_start
  #           "2part" => use TwoPart estimates as starting values
  #           "Tobit" => use Tobit model estimates as starting values
  #
  # DETAILS:
  # see http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/limlrt/
  #
  # VALUE:
  # pvalue    p-value of the LIM-LRT
  # diffest   estimated difference of continuous subdistributions' point estimates
  # techx     estimated percentage of technical pointmass values in group 1
  # techy     estimated percentage of technical pointmass values in group 0
  # biolx     estimated percentage of biological pointmass values in group 1
  # bioly     estimated percentage of biological pointmass values in group 0
  #
  # Author:
  # Andreas Gleiss (2014)
  #
  
  z<-data[group==1] # group1 = "x"
  d<-(z>massat)
  w<-data[group==0] # group0 = "y"
  h<-(w>massat)

  lod<-min(c(w[h],z[d]))-eps

  if(preest=="2part") {
    if(sum(d)==0) { # only PMVs => p_ML=1, mu and sd undetermined => optimization gives starting values
      mu_start1<-mu_start
      sd_start1<-sd_start
      p_start1<-1
      }
    else {
      mu_start1<-median(z[d])
      sd_start1<-mad(z[d])
      p_start1<-max(0.05,1-sum(d)/(length(d)*(1-pnorm(lod,mu_start1,sd_start1))))
      }
    if(sum(h)==0) { # all PMVs => p_ML=1, mu and sd undetermined => optimization gives starting values
      mu_start2<-mu_start
      sd_start2<-sd_start
      p_start2<-1
      }
    else {
      mu_start2<-median(w[h])
      sd_start2<-mad(w[h])
      p_start2<-max(0.05,1-sum(h)/(length(h)*(1-pnorm(lod,mu_start2,sd_start2))))
      }
    if(sum(d)==0 & sum(h)==0) { # all PMVs => p_ML=1, mu and sd undetermined => optimization gives starting values
      mu_start0<-mu_start
      sd_start0<-sd_start
      p_start0<-1
      }
    else {
      mu_start0<-median(c(z[d],w[h]))
      sd_start0<-mad(c(z[d],w[h]))
      p_start0<-max(0.05,1-sum(c(d,h))/(length(c(d,h))*(1-pnorm(lod,mu_start0,sd_start0))))
      }
    }
  else if(preest=="Tobit") { 
    if(sum(d)==0) {
      mu_start1<-mu_start
      sd_start1<-sd_start
      p_start1<-1
      }
    else {
      tob1<-tobit.single(z,point.mass=massat,dl=lod)
      mu_start1<-tob1$mean
      sd_start1<-tob1$sd
      p_start1<-max(0.05,1-sum(d)/(length(d)*(1-pnorm(lod,mu_start1,sd_start1))))
      }
    if(sum(h)==0) {
      mu_start2<-mu_start
      sd_start2<-sd_start
      p_start2<-1
      }
    else {
      tob2<-tobit.single(w,point.mass=massat,dl=lod)
      mu_start2<-tob2$mean
      sd_start2<-tob2$sd
      p_start2<-max(0.05,1-sum(h)/(length(h)*(1-pnorm(lod,mu_start2,sd_start2))))
      }
    if(sum(d)==0 & sum(h)==0) { # all PMVs => p_ML=1, mu and sd undetermined => optimization gives starting values
      mu_start0<-mu_start
      sd_start0<-sd_start
      p_start0<-1
      }
    else {
      tob0<-tobit.single(c(z,w),point.mass=massat,dl=lod)
      mu_start0<-tob0$mean
      sd_start0<-tob0$sd
      p_start0<-max(0.05,1-sum(c(d,h))/(length(c(d,h))*(1-pnorm(lod,mu_start0,sd_start0))))
      }
    }
  else {
    mu_start1<-mu_start2<-mu_start0<-mu_start
    sd_start1<-sd_start2<-sd_start0<-sd_start
    p_start1<-p_start2<-p_start0<-0.5
    }

  if(print) {
    cat("\nlod =",lod)
    cat("\nPre-estimators group0:",c(mu_start2,sd_start2,p_start2),"\n")
    cat("Pre-estimators group1:",c(mu_start1,sd_start1,p_start1),"\n")
    cat("Pre-estimators null:  ",c(mu_start0,sd_start0,p_start0),"\n")
    }

  optxy<-optim(par=c(mu_start1,sd_start,mu_start2,sd_start,p_start1,p_start2),fn=lim1,
                z=z,d=d,w=w,h=h,method="L-BFGS-B",
                lower=c(mu_low,sd_low,mu_low,sd_low,0,0),
                upper=c(30,30,30,30,1,1), lod=lod,
                print=pprint,quickret=quickret)
  if(print & optxy$convergence!=0)
    cat("\nConvergence problems under H1:",optxy$message)
  pmvy<-100*(1-sum(h)/length(h))
  bioly<-100*optxy$par[6]
  techy<-100*(1-optxy$par[6])*pnorm(lod,optxy$par[3],optxy$par[4])
  if(print)
    cat("\ngroup0 (",round(pmvy,1),"% obs.PMVs,",round(techy,1),"% tech.PMVs,",round(bioly,1),"% biol.PMVs):",
        round(optxy$par[3:4],2),"(lod =",lod,")")
  pmvx<-100*(1-sum(d)/length(d))
  biolx<-100*optxy$par[5]
  techx<-100*(1-optxy$par[5])*pnorm(lod,optxy$par[1],optxy$par[2])
  if(print)
    cat("\ngroup1 (",round(pmvx,1),"% obs.PMVs,",round(techx,1),"% tech.PMVs,",round(biolx,1),"% biol.PMVs):",
        round(optxy$par[1:2],2)," -> -logL =",optxy$value,"\n")

  opt0<-optim(par=c(mu_start0,sd_start0,sd_start0,p_start0),fn=lim0,
              z=z,d=d,w=w,h=h,method="L-BFGS-B",
              lower=c(mu_low,sd_low,sd_low,0),
              upper=c(30,30,30,1), lod=lod,
              print=pprint,quickret=quickret)
  if(print & opt0$convergence!=0)
    cat("\nConvergence problems under H0:",opt0$message)
  if(print)  cat("null:")
  bioly0<-100*opt0$par[4]
  techy0<-100*(1-opt0$par[4])*pnorm(lod,opt0$par[1],opt0$par[3])
  if(print)
    cat("\ngroup0 (",round(techy0,1),"% tech.PMVs,",round(bioly0,1),"% biol.PMVs):",
        round(opt0$par[c(1,3)],2),"(lod =",lod,")")
  biolx0<-100*opt0$par[4] 
  techx0<-100*(1-opt0$par[4])*pnorm(lod,opt0$par[1],opt0$par[2]) 
  if(print)
    cat("\ngroup1 (",round(techx0,1),"% tech.PMVs,",round(biolx0,1),"% biol.PMVs):",
       round(opt0$par[1:2],2)," -> -logL =",opt0$value,"\n\n")

  erg<-2*opt0$value-2*optxy$value 
  
  pvalue<-1-pchisq(erg,2)  
  diffest=optxy$par[1]-optxy$par[3]
  
  if(print)
    cat("statistic =",erg,", p-value =",pvalue,"\n")

  return(list(pvalue=pvalue, diffest=diffest, techx=techx, techy=techy,biolx=biolx, bioly=bioly))
  }

v<-c(   0.69247375,   -0.38525328,   -0.09195440,    0.51218987,   -0.17191259, -999.00000000,   -0.23233105,   -0.39769759,
       -0.13172148, -999.00000000,    0.42270656,   -0.39672140,    0.05429913,   -0.07517673, -999.00000000, -999.00000000,
     -999.00000000,    0.01362811, -999.00000000, -999.00000000,   -0.01328998,    1.14730704,   -0.00702700, -999.00000000,
     -999.00000000,    2.15011347,    0.28037684,    0.91463191,   -0.43666377,    2.27688526,    2.53985537, -999.00000000,
        0.65461774, -999.00000000,    0.66475825,    2.26210668,    1.71774215,   -0.02573570,    0.75957574,    0.52330104,
        2.63435338,    1.90234394,    1.38587776,    1.65364339,    3.07758557,    1.63079534,    2.04622760, -999.00000000,
        1.79002209,    1.46890535)

plot.lod(row=v,group=c(rep(0,25),rep(1,25)),dl=-10,xmax=20,massat=-999)
lim.lrt(data=v, group=c(rep(0,25),rep(1,25)), massat=-999, pprint=F,print=T,preest="2part")
 