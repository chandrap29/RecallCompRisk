rm(list=ls(all=TRUE))
#Install Required Packages
#library("VGAM");library("zipfR");library("stats");
#library("stats4");library("coda")

# Define Code Inside Loop
ff=function(){
  #---------------------------------Data Generation------------------------#
  #Initialize Parameters
  n           #sample size
  alp1=alp11; #shape 
  bet1=bet11; #scale
  alp2=alp22; #shape
  bet2=bet22; #scale
  lmd1=lmd11; #non-recall parameter
  lmd2=lmd22; #non-recall parameter
  alpha=0.05  #Type-I error
  
  n3=0
  while(n3<n*0.10){ #to ensure minimum observations in no-recall
    
    x=rweibull(n,shape=alp1,scale=bet1^(-1/alp1));x
    y=rweibull(n,shape=alp2,scale=bet2^(-1/alp2));y 
    tt=pmin(x,y);tt                               #failure time
    cc1=ifelse(x==tt,1,2);cc1
    #S=rexp(n,.);S                                #Generate Monitoring 
    #S=runif(n,.,.)                              
    cc=ifelse(tt<S,cc1,0);cc                      #censoring Indicator                   
    cs=S[which(cc==0)];cs                         #right-censored observations
    ncs=data.frame(tt[which(cc!=0)],cc[which(cc!=0)],S[which(cc!=0)]);ncs #non-censored observations
    
    y11=which(ncs[,2]==1);y11       #Cause 1
    y12=which(ncs[,2]==2);y12       #Cause 2
    
    g1=ncs[y11,1];g1                #Cause 1 Failure time
    g2=ncs[y11,3];g2                #monitoring time
    
    g3=ncs[y12,1];g3                #cause 2 Failure time
    g4=ncs[y12,3];g4                #monitoring time
    
    pp1=exp(-(g2-g1)*lmd1);pp1        # non-recall probability for cause 1   
    pp2=exp(-(g4-g3)*lmd2);pp2        # non-recall probability for cause 2
    z11=rbinom(length(y11),1,pp1);z11 #binomial variate
    z22=rbinom(length(y12),1,pp2);z22 #binomial variate
    
    r1=which(z11==1);r1
    r2=which(z22==1);r2
    nr1=which(z11==0);nr1
    nr2=which(z22==0);nr2
    nr=c(nr1,nr2);nr  
    
    t1=g1[r1];t1
    s1=g2[r1];s1
    t2=g3[r2];t2
    s2=g4[r2];s2
    
    s3=c(g2[nr1],g4[nr2]);s3
    s4=cs;s4
    
    z1=t1;S1=s1;z2=t2;S2=s2;S3=s3;S4=s4
    n1=length(z1);n2=length(z2);n3=length(S3);n4=length(S4)
  }
  z1;S1;z2;S2;S3;S4
  n1;n2;n3;n4 #observations under different categories
  
  
  #==================================================================#
  #
  #                        MLE Calculation
  #==================================================================#
  m.thh1=m.thh2=m.thh3=m.thh4=m.thh5=m.thh6=rep()
  #initial values
  m.th1=0.94;m.th2=1.12;m.th3=0.18;m.th4=0.95;m.th5=1.11;m.th6=0.20
  it=100
  for(j in 2:it){
    #E-step of EM
    i1=function(S,u){m.th1*m.th2*(u^(m.th1-1))*exp(-m.th2*(u^m.th1)-m.th5*(u^m.th4))*(1-exp(-(S-u)*m.th3))}
    i11=sapply(1:length(S3),function(ii){integrate(function(u){i1(S3[ii],u)},0,S3[ii],subdivisions=100L,stop.on.error=FALSE)$value});i11
    i2=function(S,u){m.th4*m.th5*(u^(m.th4-1))*exp(-m.th2*(u^m.th1)-m.th5*(u^m.th4))*(1-exp(-(S-u)*m.th6))}
    i12=sapply(1:length(S3),function(ii){integrate(function(u){i2(S3[ii],u)},0,S3[ii],subdivisions=100L,stop.on.error=FALSE)$value});i12
    p_i<<-i11/(i11+i12);p_i
    P<<-sum(p_i);P
    
    aa=m.th2*(S3^m.th1);aa
    D1=1-exp(-aa);D1
    
    bb=m.th5*(S3^m.th4);bb
    D2=1-exp(-bb);D2
    
    # Expectations
    int1=function(x){log(x)*exp(-x)}
    int2=function(x){x*log(x)*exp(-x)}
    int3=function(x){(x^(m.th4/m.th1))*log(x)*exp(-x)}
    int4=function(x){x*(log(x)^2)*exp(-x)}
    int5=function(x){(x^(m.th4/m.th1))*(log(x)^2)*exp(-x)}
    int6=function(x){(x^(m.th1/m.th4))*log(x)*exp(-x)}
    int7=function(x){(x^(m.th1/m.th4))*(log(x)^2)*exp(-x)}
    
    intgd1=array(0,length(aa))
    for (i in 1:length(aa)){
      l=0
      u=aa[i]
      intgd1[i]=integrate(int1,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K1=intgd1;K1   
    
    intgd2=array(0,length(aa))
    for (i in 1:length(aa)){
      l=0
      u=aa[i]
      intgd2[i]=integrate(int2,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K2=intgd2;K2
    
    intgd3=array(0,length(aa))
    for (i in 1:length(aa)){
      l=0
      u=aa[i]
      intgd3[i]=integrate(int3,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K3=intgd3;K3    
    
    intgd4=array(0,length(aa))
    for (i in 1:length(aa)){
      l=0
      u=aa[i]
      intgd4[i]=integrate(int4,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K4=intgd4;K4
    
    intgd5=array(0,length(aa))
    for (i in 1:length(aa)) {
      l=0
      u=aa[i]
      intgd5[i]=integrate(int5,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K5=intgd5;K5
    
    intgd6=array(0,length(bb))
    for (i in 1:length(bb)){
      l=0
      u=bb[i]
      intgd6[i]=integrate(int1,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K6=intgd6;K6
    
    intgd7=array(0,length(bb))
    for (i in 1:length(bb)){
      l=0
      u=bb[i]
      intgd7[i]=integrate(int6,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K7=intgd7;K7
    
    intgd8=array(0,length(bb))
    for (i in 1:length(bb)){
      l=0
      u=bb[i]
      intgd8[i]=integrate(int2,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K8=intgd8;K8
    
    intgd9=array(0,length(bb))
    for (i in 1:length(bb)){
      l=0
      u=bb[i]
      intgd9[i]=integrate(int7,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K9=intgd9;K9
    
    intgd10=array(0,length(bb))
    for (i in 1:length(bb)){
      l=0
      u=bb[i]
      intgd10[i]=integrate(int4,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
    }
    K10=intgd10;K10
    
    #Based on Component 1
    xi01=((m.th2^(-1/m.th1))*(gamma(1/m.th1+1)*pgamma(aa,1/m.th1+1,lower.tail = TRUE)))/D1;xi01
    xi02=((m.th5^(-1/m.th4))*(gamma(1/m.th4+1)*pgamma(bb,1/m.th4+1,lower.tail = TRUE)))/D2;xi02
    xi1=(1-(1+aa)*exp(-aa))*(1/m.th2)*(1/D1);xi1
    xi2=(gamma((m.th1+m.th4)/m.th1)*(1-pgamma(aa,(m.th1+m.th4)/m.th1,lower.tail=FALSE)))*(1/D1)*(m.th2^(-m.th4/m.th1));xi2
    xi3=(K1-log(m.th2)*D1)*(1/m.th1)*(1/D1);xi3
    xi4=(K2-log(m.th2)*(1-(1+aa)*exp(-aa)))*(1/m.th2)*(1/m.th1)*(1/D1);xi4
    xi5=(K3-log(m.th2)*(gamma((m.th1+m.th4)/m.th1)*(1-pgamma(aa,(m.th1+m.th4)/m.th1,lower.tail=FALSE))))*(m.th2^(-m.th4/m.th1))*(1/m.th1)*(1/D1);xi5
    DDD=m.th3*(S3-xi01);DDD; EEE=m.th6*(S3-xi02);EEE
    xi6=(1/m.th3)*(1-(1+DDD)*exp(-DDD))*(1/(1-exp(-DDD)));xi6
    xi7=(K4+(log(m.th2)^2)*(1-(1+aa)*exp(-aa))-2*log(m.th2)*K2)*(1/m.th2)*(1/(m.th1^2))*(1/D1);xi7
    xi8=(K5+(log(m.th2)^2)*(gamma((m.th1+m.th1)/m.th1)*(1-pgamma(aa,(m.th1+m.th4)/m.th1,lower.tail=FALSE)))-2*log(m.th2)*K3)*(m.th2^(-m.th4/m.th1))*(1/m.th1^2)*(1/D1);xi8
    
    #Based on Component 2
    xi9=(gamma((m.th1+m.th4)/m.th4)*(1-pgamma(bb,(m.th1+m.th4)/m.th4,lower.tail=FALSE)))*(1/D2)*(m.th5^(-m.th1/m.th4));xi9
    xi10=(1-(1+bb)*exp(-bb))*(1/m.th5)*(1/D2);xi10
    xi11=(K6-log(m.th5)*D2)*(1/m.th1)*(1/D2);xi11
    xi12=(K7-log(m.th5)*(gamma((m.th1+m.th4)/m.th4)*(1-pgamma(bb,(m.th1+m.th4)/m.th4,lower.tail=FALSE))))*(m.th5^(-m.th1/m.th4))*(1/m.th4)*(1/D2);xi12
    xi13=(K8-log(m.th5)*(1-(1+bb)*exp(-bb)))*(1/m.th5)*(1/m.th4)*(1/D2);xi13
    xi14=(1/m.th6)*(1-(1+EEE)*exp(-EEE))*(1/(1-exp(-EEE)));xi14
    xi15=(K9+(log(m.th5)^2)*(gamma((m.th1+m.th4)/m.th4)*(1-pgamma(bb,(m.th1+m.th4)/m.th4,lower.tail=FALSE)))-2*log(m.th5)*K7)*(m.th5^(-m.th1/m.th4))*(1/m.th4^2)*(1/D2);xi15
    xi16=(K10+(log(m.th5)^2)*(1-(1+bb)*exp(-bb))-2*log(m.th5)*K8)*(1/m.th5)*(1/(m.th4^2))*(1/D2);xi16
    
    #M Step of EM
    m.thh1[j]=m.th1;m.thh2[j]=m.th2;m.thh3[j]=m.th3;m.thh4[j]=m.th4;m.thh5[j]=m.th5;m.thh6[j]=m.th6
    m.thh1[j+1]=(n1+P)/(-sum(log(z1))-sum(p_i*xi3)+m.thh2[j]*sum((z1^m.thh1[j])*log(z1))+m.thh2[j]*sum((z2^m.thh1[j])*log(z2))+m.thh2[j]*sum(p_i*xi4)+m.thh2[j]*sum((1-p_i)*xi12)+m.thh2[j]*sum((S4^m.thh1[j])*log(S4)))
    m.thh2[j+1]=(n1+P)/(sum(z1^m.thh1[j+1])+sum(z2^m.thh1[j+1])+sum(p_i*xi1)+sum((1-p_i)*xi9)+sum(S4^m.thh1[j+1]))
    m.thh3[j+1]=P/(sum(S1-z1)+sum(p_i*xi6))
    m.thh4[j+1]=(n2+n3-P)/(-sum(log(z2))-sum((1-p_i)*xi11)+m.thh5[j]*sum((z1^m.thh4[j])*log(z1))+m.thh5[j]*sum((z2^m.thh4[j])*log(z2))+m.thh5[j]*sum(p_i*xi5)+m.thh5[j]*sum((1-p_i)*xi13)+m.thh5[j]*sum((S4^m.thh4[j])*log(S4)))
    m.thh5[j+1]=(n2+n3-P)/(sum(z1^m.thh4[j+1])+sum(z2^m.thh4[j+1])+sum(p_i*xi2)+sum((1-p_i)*xi10)+sum(S4^m.thh4[j+1]))
    m.thh6[j+1]=(n3-P)/(sum(S2-z2)+sum((1-p_i)*xi14))
    m.th1=m.thh1[j+1];m.th2=m.thh2[j+1];m.th3=m.thh3[j+1];m.th4=m.thh4[j+1];m.th5=m.thh5[j+1];m.th6=m.thh6[j+1]
  }
  m.thh1;m.thh2;m.thh3;m.thh4;m.thh5;m.thh6
  #par(mfrow=c(2,3))
  #plot(m.thh1,type='l');plot(m.thh2,type='l');plot(m.thh3,type='l');plot(m.thh4,type='l');plot(m.thh5,type='l');plot(m.thh6,type='l')
  
  #Estimates at Final Iteration
  ml=c(m.thh1[it],m.thh2[it],m.thh3[it],m.thh4[it],m.thh5[it],m.thh6[it]);ml
  ML_MSE=c((ml[1]-alp11)^2,(ml[2]-bet11)^2,(ml[3]-lmd11)^2,(ml[4]-alp22)^2,(ml[5]-bet22)^2,(ml[6]-lmd22)^2);ML_MSE
  ML_AB=c(abs(ml[1]-alp11),abs(ml[2]-bet11),abs(ml[3]-lmd11),abs(ml[4]-alp22),abs(ml[5]-bet22),abs(ml[6]-lmd22));ML_AB  
  
  #====================================================================#
  #
  #                ACI Calculations
  #=====================================================================#
  aa1=ml[2]*(S3^ml[1]);aa1
  DD1=1-exp(-aa1);DD1
  
  bb1=ml[5]*(S3^ml[4]);bb1
  DD2=1-exp(-bb1);DD2
  
  # Expectations
  intt1=function(x){log(x)*exp(-x)}
  intt2=function(x){x*log(x)*exp(-x)}
  intt3=function(x){(x^(ml[4]/ml[1]))*log(x)*exp(-x)}
  intt4=function(x){x*(log(x)^2)*exp(-x)}
  intt5=function(x){(x^(ml[4]/ml[1]))*(log(x)^2)*exp(-x)}
  intt6=function(x){(x^(ml[1]/ml[4]))*log(x)*exp(-x)}
  intt7=function(x){(x^(ml[1]/ml[4]))*(log(x)^2)*exp(-x)}
  
  #Based on Component 1
  inttgd1=array(0,length(aa1))
  for (i in 1:length(aa1)){
    l=0
    u=aa1[i]
    inttgd1[i]=integrate(intt1,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK1=inttgd1;KK1   
  
  inttgd2=array(0,length(aa1))
  for (i in 1:length(aa1)){
    l=0
    u=aa1[i]
    inttgd2[i]=integrate(intt2,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK2=inttgd2;KK2
  
  inttgd3=array(0,length(aa1))
  for (i in 1:length(aa1)){
    l=0
    u=aa1[i]
    inttgd3[i]=integrate(intt3,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK3=inttgd3;KK3    
  
  inttgd4=array(0,length(aa1))
  for (i in 1:length(aa1)){
    l=0
    u=aa1[i]
    inttgd4[i]=integrate(intt4,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK4=inttgd4;KK4
  
  inttgd5=array(0,length(aa1))
  for (i in 1:length(aa1)) {
    l=0
    u=aa1[i]
    inttgd5[i]=integrate(intt5,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK5=inttgd5;KK5
  
  #Based on Component 2
  inttgd6=array(0,length(bb1))
  for (i in 1:length(bb1)) {
    l=0
    u=bb1[i]
    inttgd6[i]=integrate(intt1, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK6=inttgd6;KK6
  
  inttgd7=array(0,length(bb1))
  for (i in 1:length(bb1)) {
    l=0
    u=bb1[i]
    inttgd7[i]=integrate(intt6, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK7=inttgd7;KK7
  
  inttgd8=array(0,length(bb1))
  for (i in 1:length(bb1)) {
    l=0
    u=bb1[i]
    inttgd8[i]=integrate(intt2, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK8=inttgd8;KK8
  
  inttgd9=array(0,length(bb1))
  for (i in 1:length(bb1)) {
    l=0
    u=bb1[i]
    inttgd9[i]=integrate(intt7, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK9=inttgd9;KK9
  
  inttgd10=array(0,length(bb1))
  for (i in 1:length(bb1)) {
    l=0
    u=bb1[i]
    inttgd10[i]=integrate(intt4, l, u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  KK10=inttgd10;KK10
  
  xxi01=((ml[2]^(-1/ml[1]))*(gamma(1/ml[1]+1)*pgamma(aa1,1/ml[1]+1,lower.tail = TRUE)))/DD1;xxi01
  xxi02=((ml[5]^(-1/ml[4]))*(gamma(1/ml[4]+1)*pgamma(bb1,1/ml[4]+1,lower.tail = TRUE)))/DD2;xxi02
  xxi1=(1-(1+aa1)*exp(-aa1))*(1/ml[2])*(1/DD1);xxi1
  xxi2=(gamma((ml[1]+ml[4])/ml[1])*(1-pgamma(aa1,(ml[1]+ml[4])/ml[1],lower.tail=FALSE)))*(1/DD1)*(ml[2]^(-ml[4]/ml[1]));xxi2
  xxi3=(KK1-log(ml[2])*DD1)*(1/ml[1])*(1/DD1);xxi3
  xxi4=(KK2-log(ml[2])*(1-(1+aa1)*exp(-aa1)))*(1/ml[2])*(1/ml[1])*(1/DD1);xxi4
  xxi5=(KK3-log(ml[2])*(gamma((ml[1]+ml[4])/ml[1])*(1-pgamma(aa1,(ml[1]+ml[4])/ml[1],lower.tail=FALSE))))*(ml[2]^(-ml[4]/ml[1]))*(1/ml[1])*(1/DD1);xxi5
  
  DDD1=ml[3]*(S3-xxi01);DDD1; EEE1=ml[6]*(S3-xxi02);EEE1
  xxi6=(1/ml[3])*(1-(1+DDD1)*exp(-DDD1))*(1/(1-exp(-DDD1)));xxi6
  xxi7=(KK4+(log(ml[2])^2)*(1-(1+aa1)*exp(-aa1))-2*log(ml[2])*KK2)*(1/ml[2])*(1/(ml[1]^2))*(1/DD1);xxi7
  xxi8=(KK5+(log(ml[2])^2)*(gamma((ml[1]+ml[4])/ml[1])*(1-pgamma(aa1,(ml[1]+ml[4])/ml[1],lower.tail=FALSE)))-2*log(ml[2])*KK3)*(ml[2]^(-ml[4]/ml[1]))*(1/ml[1]^2)*(1/DD1);xxi8
  
  xxi9=(gamma((ml[1]+ml[4])/ml[4])*(1-pgamma(bb1,(ml[1]+ml[4])/ml[4],lower.tail=FALSE)))*(1/DD2)*(ml[5]^(-ml[1]/ml[4]));xxi9
  xxi10=(1-(1+bb1)*exp(-bb1))*(1/ml[5])*(1/DD2);xxi10
  xxi11=(KK6-log(ml[5])*DD2)*(1/ml[1])*(1/DD2);xxi11
  xxi12=(KK7-log(ml[5])*(gamma((ml[1]+ml[4])/ml[4])*(1-pgamma(bb1,(ml[1]+ml[4])/ml[4],lower.tail=FALSE))))*(ml[5]^(-ml[1]/ml[4]))*(1/ml[4])*(1/DD2);xxi12
  xxi13=(KK8-log(ml[5])*(1-(1+bb1)*exp(-bb1)))*(1/ml[5])*(1/ml[4])*(1/DD2);xxi13
  xxi14=(1/ml[6])*(1-(1+EEE1)*exp(-EEE1))*(1/(1-exp(-EEE1)));xxi14
  xxi15=(KK9+(log(ml[5])^2)*(gamma((ml[1]+ml[4])/ml[4])*(1-pgamma(bb1,(ml[1]+ml[4])/ml[4],lower.tail=FALSE)))-2*log(ml[5])*KK7)*(ml[5]^(-ml[1]/ml[4]))*(1/ml[4]^2)*(1/DD2);xxi15
  xxi16=(KK10+(log(ml[5])^2)*(1-(1+bb1)*exp(-bb1))-2*log(ml[5])*KK8)*(1/ml[5])*(1/(ml[4]^2))*(1/DD2);xxi16
  
  #Complete Information Matrix
  #alpha1
  I11=-(n1+P)/(ml[1]^2)-ml[2]*sum((z1^ml[1])*(log(z1)^2))-ml[2]*sum((z2^ml[1])*(log(z2)^2))-ml[2]*sum(p_i*xxi7)-ml[2]*sum((1-p_i)*xxi15)-ml[2]*sum((S4^ml[1])*(log(S4)^2));I11
  I12=-sum((z1^ml[1])*log(z1))-sum((z2^ml[1])*log(z2))-sum(p_i*xxi4)-sum((1-p_i)*xxi12)-sum((S4^ml[1])*log(S4));I12
  I13=0;I13
  I14=0;I14
  I15=0;I15
  I16=0;I16
  
  #beta1
  I21=I12;I21
  I22=-(n1+P)/(ml[2]^2);I22
  I23=0;I23
  I24=0;I24
  I25=0;I25
  I26=0;I26
  
  #lambda1  
  I31=0;I31
  I32=0;I32
  I33=-P/(ml[3]^2);I33
  I34=0;I34
  I35=0;I35
  I36=0;I36
  
  #alpha2
  I41=0;I41
  I42=0;I42
  I43=0;I43
  I44=-(n2+n3-P)/(ml[4]^2)-ml[5]*sum((z1^ml[4])*(log(z1)^2))-ml[5]*sum((z2^ml[4])*(log(z2)^2))-ml[5]*sum(p_i*xxi8)-ml[5]*sum((1-p_i)*xxi16)-ml[5]*sum((S4^ml[4])*(log(S4)^2));I44
  I45=-sum((z1^ml[4])*log(z1))-sum((z2^ml[4])*log(z2))-sum(p_i*xxi5)-sum((1-p_i)*xxi13)-sum((S4^ml[4])*log(S4));I45
  I46=0;I46
  
  #beta2
  I51=0;I51 
  I52=0;I52
  I53=0;I53
  I54=I45;I54
  I55=-(n2+n3-P)/(ml[5]^2);I55
  I56=0;I56
  
  #lambda2
  I61=0;I61
  I62=0;I62
  I63=0;I63
  I64=0;I64
  I65=0;I65
  I66=-(n3-P)/(ml[6]^2);I66
  
  Mat1=-matrix(c(I11,I12,I13,I14,I15,I16,I21,I22,I23,I24,I25,I26,I31,I32,I33,I34,I35,I36,I41,I42,I43,I44,I45,I46,
                 I51,I52,I53,I54,I55,I56,I61,I62,I63,I64,I65,I66),6,6,byrow=TRUE);Mat1
  
  
  #Incomplete Information Matrix
  #grad
  GD1=(n1+P)/ml[1]+sum(log(z1))-ml[2]*sum((z1^ml[1])*log(z1))-ml[2]*sum((z2^ml[1])*log(z2))+sum(p_i*xxi3)-ml[2]*sum(p_i*xxi4)-ml[2]*sum((1-p_i)*xxi12)-ml[2]*sum((S4^ml[1])*log(S4));GD1
  GD2=(n1+P)/ml[2]-sum(z1^ml[1])-sum(z2^ml[1])-sum(p_i*xxi1)-sum((1-p_i)*xxi9)-sum(S4^ml[1]);GD2
  GD3=P/ml[3]-sum(S1-z1)-sum(p_i*xxi6);GD3
  GD4=(n2+n3-P)/ml[4]+sum(log(z2))-ml[5]*sum((z1^ml[4])*log(z1))-ml[5]*sum((z2^ml[4])*log(z2))-ml[5]*sum(p_i*xxi5)+sum((1-p_i)*xxi11)-ml[5]*sum((1-p_i)*xxi13)-ml[5]*sum((S4^ml[4])*log(S4));GD4
  GD5=(n2+n3-P)/ml[5]-sum(z1^ml[4])-sum(z2^ml[4])-sum(p_i*xxi2)-sum((1-p_i)*xxi10)-sum(S4^ml[4]);GD5
  GD6=(n3-P)/ml[6]-sum(S2-z2)-sum((1-p_i)*xxi14);GD6
  
  GD=c(GD1,GD2,GD3,GD4,GD5,GD6);GD
  
  Mat2=GD%*%t(GD);Mat2
  Var_Cov=solve(Mat1-Mat2);Var_Cov
  
  #-------------100(1-alpha)% CIs and CP--------#  
  ci_l=ci_u=rep();alpha=0.05
  for (i in 1:6) {
    ci_l[i]=ml[i]-1.96*sqrt((diag(Var_Cov)[i]))
    ci_u[i]=ml[i]+1.96*sqrt((diag(Var_Cov)[i]))
  }
  
  cl1=c(ci_l[1],ci_u[1]);cl1
  cl2=c(ci_l[2],ci_u[2]);cl2
  cl3=c(ci_l[3],ci_u[3]);cl3
  cl4=c(ci_l[4],ci_u[4]);cl4
  cl5=c(ci_l[5],ci_u[5]);cl5
  cl6=c(ci_l[6],ci_u[6]);cl6
  #ACI_lim=c(cl1,cl2,cl3,cl4,cl5,cl6);ACI_lim
  
  l1=c(ci_u[1]-ci_l[1]);l1
  l2=c(ci_u[2]-ci_l[2]);l2
  l3=c(ci_u[3]-ci_l[3]);l3
  l4=c(ci_u[4]-ci_l[4]);l4
  l5=c(ci_u[5]-ci_l[5]);l5
  l6=c(ci_u[6]-ci_l[6]);l6
  ACI_L=c(l1,l2,l3,l4,l5,l6);ACI_L
  
  cp1=ifelse(ci_l[1]<alp11 && alp11<ci_u[1],1,0);cp1
  cp2=ifelse(ci_l[2]<bet11 && bet11<ci_u[2],1,0);cp2
  cp3=ifelse(ci_l[3]<lmd11 && lmd11<ci_u[3],1,0);cp3
  cp4=ifelse(ci_l[4]<alp22 && alp22<ci_u[4],1,0);cp4
  cp5=ifelse(ci_l[5]<bet22 && bet22<ci_u[5],1,0);cp5
  cp6=ifelse(ci_l[6]<lmd22 && lmd22<ci_u[6],1,0);cp6
  
  coverage_ACI=c(cp1,cp2,cp3,cp4,cp5,cp6);coverage_ACI
  ML=c(ml,ML_MSE,ML_AB,ACI_L,1,1,1,1,1,1,coverage_ACI);ML
  
  #===========================================================================#
  #
  #               Bayesian Estimation
  #
  #===========================================================================#
  
  #Define Hyper-parameters
  mu1=mu2=mu3=mu4=mu5=mu6=2;
  nu1=mu1/alp1;nu1
  nu2=mu2/bet1;nu2
  nu3=mu3/lmd1;nu3
  nu4=mu4/alp2;nu4
  nu5=mu5/bet2;nu5
  nu6=mu6/lmd2;nu6
  b.thh1=b.thh2=b.thh3=b.thh4=b.thh5=b.thh6=rep()
  it1=35000 #iterations
  for(j in 1:it1){
    b.thh1[1]=ml[1];b.thh2[1]=ml[2];b.thh3[1]=ml[3];b.thh4[1]=ml[4];b.thh5[1]=ml[5];b.thh6[1]=ml[6]
    j1=function(S,u){b.thh1[j]*b.thh2[j]*(u^(b.thh1[j]-1))*exp(-b.thh2[j]*(u^b.thh1[j])-b.thh5[j]*(u^b.thh4[j]))*(1-exp(-(S-u)*b.thh3[j]))}
    l11=sapply(1:length(S3),function(ii){integrate(function(u){j1(S3[ii],u)},0,S3[ii],subdivisions=100L,stop.on.error=FALSE)$value});l11
    j2=function(S,u){b.thh4[j]*b.thh5[j]*(u^(b.thh4[j]-1))*exp(-b.thh2[j]*(u^b.thh1[j])-b.thh5[j]*(u^b.thh4[j]))*(1-exp(-(S-u)*b.thh6[j]))}
    l12=sapply(1:length(S3),function(ii){integrate(function(u){j2(S3[ii],u)},0,S3[ii],subdivisions=100L,stop.on.error=FALSE)$value});l12
    q_i=l11/(l11+l12);q_i
    
    z_i=rbinom(length(S3),1,q_i);z_i
    Z=sum(z_i);Z
    
    tt1=(-(1/b.thh2[j])*log(1-runif(1)*(1-exp(-b.thh2[j]*(S3^b.thh1[j])))))^(1/b.thh1[j]);tt1  
    tt2=(-(1/b.thh5[j])*log(1-runif(1)*(1-exp(-b.thh5[j]*(S3^b.thh4[j])))))^(1/b.thh4[j]);tt2  
    
    u_i=-(1/b.thh3[j])*log(1-runif(1)*(1-exp(-b.thh3[j]*(S3-tt1))));u_i 
    v_i=-(1/b.thh6[j])*log(1-runif(1)*(1-exp(-b.thh6[j]*(S3-tt2))));v_i
    
    post1=function(al1,be1){
      zz1=(al1^(n1+Z+mu1-1))*(prod(z1^(al1-1)))*(prod((tt1^z_i)^(al1-1)))
      zz2=exp(-be1*sum(z1^al1)-be1*sum(z2^al1)-be1*sum(z_i*(tt1^al1))-be1*sum((1-z_i)*(tt2^al1))-be1*sum(S4^al1)-al1*nu1)
      zzz=zz1*zz2
      return(zzz)}
    
    last1=b.thh1[j]
    cand1=(rnorm(1,alp1,Var_Cov[1,1]))
    r1=ifelse(is.finite(post1(cand1,b.thh2[j])/post1(last1,b.thh2[j]))==TRUE,post1(cand1,b.thh2[j])/post1(last1,b.thh2[j]),0.01)
    if(runif(1)<min(r1,1)) last1<-cand1
    b.thh1[j+1]=last1
    
    b.thh2[j+1]<-rgamma(1,shape=n1+Z+mu2,rate=sum(z1^b.thh1[j+1])+sum(z2^b.thh1[j+1])+sum(z_i*(tt1^b.thh1[j+1]))+sum((1-z_i)*(tt2^b.thh1[j+1]))+sum(S4^b.thh1[j+1])+nu2)
    b.thh3[j+1]<-rgamma(1,shape=Z+mu3,rate=sum(S1-z1)+sum(z_i*u_i)+nu3)
    
    post2=function(al2,be2){
      ww1=(al2^(n2+n3-Z+mu4-1))*(prod(z2^(al2-1)))*(prod((tt2^(1-z_i))^(al2-1)))
      ww2=exp(-be2*sum(z1^al2)-be2*sum(z2^al2)-be2*sum(z_i*(tt1^al2))-be2*sum((1-z_i)*(tt2^al2))-be2*sum(S4^al2)-al2*nu4)
      www=ww1*ww2
      return(www)}
    
    last2=b.thh4[j]
    cand2=(rnorm(1,alp2,Var_Cov[4,4]))
    r2=ifelse(is.finite(post2(cand2,b.thh5[j])/post2(last2,b.thh5[j]))==TRUE,post2(cand2,b.thh5[j])/post2(last2,b.thh5[j]),0.01)
    if(runif(1)<min(r2,1)) last2<-cand2
    b.thh4[j+1]=last2
    
    b.thh5[j+1]<-rgamma(1,shape=n2+n3-Z+mu5,rate=sum(z1^b.thh4[j+1])+sum(z2^b.thh4[j+1])+sum(z_i*(tt1^b.thh4[j+1]))+sum((1-z_i)*(tt2^b.thh4[j+1]))+sum(S4^b.thh4[j+1])+nu5)
    b.thh6[j+1]<-rgamma(1,shape=n3-Z+mu6,rate=sum(S2-z2)+sum((1-z_i)*v_i)+nu6)
  }
  b.thh1;b.thh2;b.thh3;b.thh4;b.thh5;b.thh6
  length(b.thh1);length(b.thh2);length(b.thh3);length(b.thh4);length(b.thh5);length(b.thh6)
  #remove first 5K observations from chain 
  ch1=b.thh1[5001:35000];ch2=b.thh2[5001:35000];ch3=b.thh3[5001:35000];ch4=b.thh4[5001:35000];ch5=b.thh5[5001:35000];ch6=b.thh6[5001:35000]
  #par(mfrow=c(2,3))
  #acf(ch1);acf(ch2);acf(ch3);acf(ch4);acf(ch5);acf(ch6)
  
  #Remove Lag from Chains
  zz1=seq(10,length(ch1),10);zz2=seq(10,length(ch2),10);zz3=seq(10,length(ch3),10);zz4=seq(10,length(ch4),10);zz5=seq(10,length(ch5),10);zz6=seq(10,length(ch6),10)
  ch111=ch1[zz1];ch222=ch2[zz2];ch333=ch3[zz3];ch444=ch4[zz4];ch555=ch5[zz5];ch666=ch6[zz6]
  length(ch111);length(ch222);length(ch333);length(ch444);length(ch555);length(ch666)
  #Final Chains corresponding to each parameters
  ch11=ch111[1:3000];ch22=ch222[1:3000];ch33=ch333[1:3000];ch44=ch444[1:3000];ch55=ch555[1:3000];ch66=ch666[1:3000];
  
  #Bayes estimates under SELF
  b.th1=b.th2=b.th3=b.th4=b.th5=b.th6=rep()
  b.th1=mean(ch11);b.th2=mean(ch22);b.th3=mean(ch33);b.th4=mean(ch44);b.th5=mean(ch55);b.th6=mean(ch66)
  b.th=c(b.th1,b.th2,b.th3,b.th4,b.th5,b.th6);b.th  
  b.mse=c((b.th1-alp1)^2,(b.th2-bet1)^2,(b.th3-lmd1)^2,(b.th4-alp2)^2,(b.th5-bet2)^2,(b.th6-lmd2)^2);b.mse
  b.bias=c(abs(b.th1-alp1),abs(b.th2-bet1),abs(b.th3-lmd1),abs(b.th4-alp2),abs(b.th5-bet2),abs(b.th6-lmd2));b.bias
  
  #HPD Intervals
  h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33));h.th4=HPDinterval(mcmc(ch44));h.th5=HPDinterval(mcmc(ch55));h.th6=HPDinterval(mcmc(ch66))
  #HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2],h.th4[,1],h.th4[,2],h.th5[,1],h.th5[,2],h.th6[,1],h.th6[,2]);HPD_lim
  
  #HPD Lengths
  h.th1_l=h.th1[,2]-h.th1[,1];h.th2_l=h.th2[,2]-h.th2[,1];h.th3_l=h.th3[,2]-h.th3[,1];h.th4_l=h.th4[,2]-h.th4[,1];h.th5_l=h.th5[,2]-h.th5[,1];h.th6_l=h.th6[,2]-h.th6[,1]
  HPD_L=c(h.th1_l,h.th2_l,h.th3_l,h.th4_l,h.th5_l,h.th6_l);HPD_L                                #length of HPD
  
  #Shape
  shp1=(h.th1[,2]-b.th1)/(b.th1-h.th1[,1])
  shp2=(h.th2[,2]-b.th2)/(b.th2-h.th2[,1])
  shp3=(h.th3[,2]-b.th3)/(b.th3-h.th3[,1])
  shp4=(h.th4[,2]-b.th4)/(b.th4-h.th4[,1])
  shp5=(h.th5[,2]-b.th5)/(b.th5-h.th5[,1])
  shp6=(h.th6[,2]-b.th6)/(b.th6-h.th6[,1])
  
  H_shape=c(shp1,shp2,shp3,shp4,shp5,shp6);H_shape
  
  #Coverage Probability
  coverage_HPD1=ifelse(alp1>h.th1[,1] && alp1<h.th1[,2],1,0)                      #HPD interval Estimate 
  coverage_HPD2=ifelse(bet1>h.th2[,1] && bet1<h.th2[,2],1,0)
  coverage_HPD3=ifelse(lmd1>h.th3[,1] && lmd1<h.th3[,2],1,0)
  coverage_HPD4=ifelse(alp2>h.th4[,1] && alp2<h.th4[,2],1,0)                      #HPD interval Estimate 
  coverage_HPD5=ifelse(bet2>h.th5[,1] && bet2<h.th5[,2],1,0)
  coverage_HPD6=ifelse(lmd2>h.th6[,1] && lmd2<h.th6[,2],1,0)
  coverage_HPD=c(coverage_HPD1,coverage_HPD2,coverage_HPD3,coverage_HPD4,coverage_HPD5,coverage_HPD6);coverage_HPD              #coverage HPD
  
  Bayes=c(b.th,b.mse,b.bias,HPD_L,H_shape,coverage_HPD);Bayes
  return(c(ML,Bayes))
}
ff()
start=Sys.time()
#k=500 (1000) define desired iterations 
rr=replicate(k,ff());rr
rm=rr[,colSums(is.na(rr))!= nrow(rr)];rm
rm1=rm[ ,!apply(rm,2,function(x) any(is.na(x)))];rm1 #remove if any NaN or NA cases in any iterations
rm2=rowMeans(rm1);rm2      #take row-means to get average results based on k iterations
#write.csv(rm2,"file.csv") #save file
end=Sys.time()
end-start