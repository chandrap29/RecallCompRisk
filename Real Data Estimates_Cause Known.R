#--------------------------------Data Analysis--------------------------------#
rm(list=ls(all=TRUE))
#library("VGAM");library("zipfR");library("coda")
#Importing and Filtering Data
Data_2<-read.csv("file.name",header=T)
#View(Data_2)
dim(Data_2)

#categorise data as per recall, non-recall and right-censored
z1;S1
z2;S2
S3;S4

n1<-length(z1);n1      #recall observation associated with cause 1
n2<-length(S2);n2      #recall observation associated with cause 2
n31<-length(S31);n31   #non-recall observation associated with cause 1
n32<-length(S32);n32   #non-recall observation associated with cause 2
n4<-length(S4);n4      #right-censored observations
n<-n1+n2+n31+n32+n4;n  #total number of observations

aa=60 #scale data
z1=z1/aa;S1=S1/aa;z2=z2/aa;S2=S2/aa;S31=S31/aa;S32=S32/aa;S4=S4/aa  # Scaling of data

summary(z1);summary(z2)
summary(S31);summary(S31);
summary(S4)
#---------------------------------------------------#
#Initialization of Parameters
#---------------------------------------------------#
#install.packages("fitdistrplus")
library(fitdistrplus)
est=fitdist(z1,"weibull")
sc=est$estimate;sc
a1=sc[1];b1=sc[2]
a2=sc[1]+0.1;a2;b2=sc[2]+0.2;b2

lower1=min(S31)/10;lower1
lower2=min(S32)/10;lower2

#=====================================================================#
#                           MLE CALCULATION
#=====================================================================#
m.thh1=m.thh2=m.thh3=m.thh4=m.thh5=m.thh6=rep()
m.th1=a1;m.th2=b1;m.th3=0.5;m.th4=a2;m.th5=b2;m.th6=0.25
it=150
for(j in 1:it){
  #E-Step of EM
  dd1=exp(-m.th2*(lower1^m.th1))-exp(-m.th2*(S31^m.th1));dd1  # Denominator 1
  dd2=exp(-m.th5*(lower2^m.th4))-exp(-m.th5*(S32^m.th4));dd2  # Denominator 2
  
  int01=function(u){m.th1*m.th2*(u^m.th1)*exp(-m.th2*(u^m.th1))}
  intt01=array(0,length(S31))
  for (i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt01[i]=integrate(int01,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt01
  xi01=intt01/dd1;xi01
  
  int02=function(u){m.th4*m.th5*(u^m.th4)*exp(-m.th5*(u^m.th4))}
  intt02=array(0,length(S32))
  for (i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt02[i]=integrate(int02,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt02
  xi02=intt02/dd2;xi02
  
  int1=function(u){(u^(2*m.th1-1))*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt1=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt1[i]=integrate(int1,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt1
  xi1=intt1/dd1;xi1
  
  int2=function(u){(u^(m.th1+m.th4-1))*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt2=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt2[i]=integrate(int2,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt2
  xi2=intt2/dd1;xi2
  
  int3=function(u){log(u)*m.th1*m.th2*(u^(m.th1-1))*exp(-m.th2*(u^m.th1))}
  intt3=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt3[i]=integrate(int3,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt3
  xi3=intt3/dd1;xi3
  
  int4=function(u){(u^(2*m.th1-1))*log(u)*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt4=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt4[i]=integrate(int4,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt4
  xi4=intt4/dd1;xi4
  
  int5=function(u){(u^(m.th1+m.th4-1))*log(u)*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt5=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt5[i]=integrate(int5,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt5
  xi5=intt5/dd1;xi5
  
  dd3=exp(-m.th3*lower1)-exp(-m.th3*(S31-xi01));dd3 
  ee3=exp(-m.th3*lower2)-exp(-m.th3*(S32-xi02));ee3
  
  SS1=S31-xi01;SS1
  int6=function(u){u*m.th3*exp(-m.th3*u)}
  intt6=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=SS1[i]
    intt6[i]=integrate(int6,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt6
  xi6=intt6/dd3;xi6
  
  int7=function(u){(u^(2*m.th1-1))*(log(u)^2)*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt7=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt7[i]=integrate(int7,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt7
  xi7=intt7/dd1;xi7
  
  int8=function(u){(u^(m.th1+m.th4-1))*(log(u)^2)*m.th1*m.th2*exp(-m.th2*(u^m.th1))}
  intt8=array(0,length(S31))
  for(i in 1:length(S31)){
    l=lower1
    u=S31[i]
    intt8[i]=integrate(int8,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt8
  xi8=intt8/dd1;xi8
  
  #Based on Component 2
  int9=function(u){(u^(m.th1+m.th4-1))*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt9=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt9[i]=integrate(int9,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt9
  xi9=intt9/dd2;xi9
  
  int10=function(u){(u^(2*m.th4-1))*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt10=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt10[i]=integrate(int10,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt10
  xi10=intt10/dd2;xi10    
  
  int11=function(u){(u^(m.th4-1))*log(u)*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt11=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt11[i]=integrate(int11,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt11
  xi11=intt11/dd2;xi11
  
  int12=function(u){(u^(m.th1+m.th4-1))*log(u)*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt12=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt12[i]=integrate(int12,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt12
  xi12=intt12/dd2;xi12
  
  int13=function(u){(u^(2*m.th4-1))*log(u)*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt13=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt13[i]=integrate(int13,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt13
  xi13=intt13/dd2;xi13
  
  SS2=S32-xi02;SS2
  int14=function(u){u*m.th6*exp(-m.th6*u)}
  intt14=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=SS2[i]
    intt14[i]=integrate(int14,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt14
  xi14=intt14/ee3;xi14
  
  int15=function(u){(u^(m.th1+m.th4-1))*(log(u)^2)*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt15=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt15[i]=integrate(int15,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt15
  xi15=intt15/dd2;xi15
  
  int16=function(u){(u^(2*m.th4-1))*(log(u)^2)*m.th4*m.th5*exp(-m.th5*(u^m.th4))}
  intt16=array(0,length(S32))
  for(i in 1:length(S32)){
    l=lower2
    u=S32[i]
    intt16[i]=integrate(int16,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
  }
  intt16
  xi16=intt16/dd2;xi16
  
  #M Step of E-M
  m.thh1[j]=m.th1;m.thh2[j]=m.th2;m.thh3[j]=m.th3;m.thh4[j]=m.th4;m.thh5[j]=m.th5;m.thh6[j]=m.th6
  m.thh1[j+1]=(n1+n31)/(-sum(log(z1))-sum(xi3)+m.thh2[j]*sum((z1^m.thh1[j])*log(z1))
                      +m.thh2[j]*sum((z2^m.thh1[j])*log(z2))+m.thh2[j]*sum(xi4)
                      +m.thh2[j]*sum(xi12)+m.thh2[j]*sum((S4^m.thh1[j])*log(S4)))
  m.thh2[j+1]=(n1+n31)/(sum(z1^m.thh1[j+1])+sum(z2^m.thh1[j+1])+sum(xi1)+sum(xi9)+sum(S4^m.thh1[j+1]))
  m.thh3[j+1]=n31/(sum(S1-z1)+sum(xi6))
  m.thh4[j+1]=(n2+n32)/(-sum(log(z2))-sum(xi11)+m.thh5[j]*sum((z1^m.thh4[j])*log(z1))+m.thh5[j]*sum((z2^m.thh4[j])*log(z2))+m.thh5[j]*sum(xi5)+m.thh5[j]*sum(xi13)+m.thh5[j]*sum((S4^m.thh4[j])*log(S4)))
  m.thh5[j+1]=(n2+n32)/(sum(z1^m.thh4[j+1])+sum(z2^m.thh4[j+1])+sum(xi2)+sum(xi10)+sum(S4^m.thh4[j+1]))
  m.thh6[j+1]=(n32)/(sum(S2-z2)+sum(xi14))
  m.th1=m.thh1[j+1];m.th2=m.thh2[j+1];m.th3=m.thh3[j+1];m.th4=m.thh4[j+1];m.th5=m.thh5[j+1];m.th6=m.thh6[j+1]
}
m.thh1;m.thh2;m.thh3;m.thh4;m.thh5;m.thh6
par(mfrow=c(2,3))
plot(m.thh1,type='l');plot(m.thh2,type='l');plot(m.thh3,type='l');plot(m.thh4,type='l');plot(m.thh5,type='l');plot(m.thh6,type='l')

#ML Estimates after final iterations
ml=c(m.thh1[it],m.thh2[it],m.thh3[it],m.thh4[it],2+m.thh5[it],m.thh6[it]);ml

# Mean duration associated with causes
mn1=aa*ml[2]^(-1/ml[1])*gamma(1+1/ml[1]);mn1 
mn2=aa*ml[5]^(-1/ml[4])*gamma(1+1/ml[4]);mn2 

# Median duration associated with causes
md1=aa*ml[2]^(-1/ml[1])*(log(2))^(1/ml[1]);md1 
md2=aa*ml[5]^(-1/ml[4])*(log(2))^(1/ml[4]);md2 

#===========================================================================#
#
#                                  ACI Calculation
#============================================================================#
#Expected Values
ddd1=exp(-ml[2]*(lower1^ml[1]))-exp(-ml[2]*(S31^ml[1]));ddd1  # Denominator 1
ddd2=exp(-ml[5]*(lower2^ml[4]))-exp(-ml[5]*(S32^ml[4]));ddd2  # Denominator 2

int011=function(u){ml[1]*ml[2]*(u^ml[1])*exp(-ml[2]*(u^ml[1]))}
intt011=array(0,length(S31))
for (i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt011[i]=integrate(int011,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt011
xxi01=intt011/ddd1;xxi01

int022=function(u){ml[4]*ml[5]*(u^ml[4])*exp(-ml[5]*(u^ml[4]))}
intt022=array(0,length(S32))
for (i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt022[i]=integrate(int022,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt022
xxi02=intt022/ddd2;xxi02

int11=function(u){(u^(2*ml[1]-1))*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt11=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt11[i]=integrate(int11,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt11
xxi1=intt11/ddd1;xxi1

int22=function(u){(u^(ml[1]+ml[4]-1))*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt22=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt22[i]=integrate(int22,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt22
xxi2=intt22/ddd1;xxi2

int33=function(u){log(u)*ml[1]*ml[2]*(u^(ml[1]-1))*exp(-ml[2]*(u^ml[1]))}
intt33=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt33[i]=integrate(int33,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt33
xxi3=intt33/ddd1;xxi3

int44=function(u){(u^(2*ml[1]-1))*log(u)*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt44=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt44[i]=integrate(int44,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt44
xxi4=intt44/ddd1;xxi4

int55=function(u){(u^(ml[1]+ml[4]-1))*log(u)*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt55=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt55[i]=integrate(int55,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt55
xxi5=intt55/ddd1;xxi5

ddd3=exp(-ml[3]*lower1)-exp(-ml[3]*(S31-xxi01));ddd3 
eee3=exp(-ml[3]*lower2)-exp(-ml[3]*(S32-xxi02));eee3

SS11=S31-xxi01;SS11
int66=function(u){u*ml[3]*exp(-ml[3]*u)}
intt66=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=SS11[i]
  intt66[i]=integrate(int66,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt66
xxi6=intt66/ddd3;xxi6

int77=function(u){(u^(2*ml[1]-1))*(log(u)^2)*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt77=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt77[i]=integrate(int77,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt77
xxi7=intt77/ddd1;xxi7

int88=function(u){(u^(ml[1]+ml[4]-1))*(log(u)^2)*ml[1]*ml[2]*exp(-ml[2]*(u^ml[1]))}
intt88=array(0,length(S31))
for(i in 1:length(S31)){
  l=lower1
  u=S31[i]
  intt88[i]=integrate(int88,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt88
xxi8=intt88/ddd1;xxi8

#Based on Component 2
int99=function(u){(u^(ml[1]+ml[4]-1))*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt99=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt99[i]=integrate(int99,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt99
xxi9=intt99/ddd2;xxi9

int101=function(u){(u^(2*ml[4]-1))*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt101=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt101[i]=integrate(int101,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt101
xxi10=intt101/ddd2;xxi10    

int111=function(u){(u^(ml[4]-1))*log(u)*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt111=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt111[i]=integrate(int111,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt111
xxi11=intt111/ddd2;xxi11

int121=function(u){(u^(ml[1]+ml[4]-1))*log(u)*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt121=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt121[i]=integrate(int121,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt121
xxi12=intt121/ddd2;xxi12

int131=function(u){(u^(2*ml[4]-1))*log(u)*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt131=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt131[i]=integrate(int131,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt131
xxi13=intt131/ddd2;xxi13

SS22=S32-xxi02;SS2
int141=function(u){u*ml[6]*exp(-ml[6]*u)}
intt141=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=SS2[i]
  intt141[i]=integrate(int141,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt141
xxi14=intt141/eee3;xxi14

int151=function(u){(u^(ml[1]+ml[4]-1))*(log(u)^2)*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt151=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt151[i]=integrate(int151,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt151
xxi15=intt151/ddd2;xxi15

int161=function(u){(u^(2*ml[4]-1))*(log(u)^2)*ml[4]*ml[5]*exp(-ml[5]*(u^ml[4]))}
intt161=array(0,length(S32))
for(i in 1:length(S32)){
  l=lower2
  u=S32[i]
  intt161[i]=integrate(int161,l,u,subdivisions=100L,stop.on.error = FALSE,abs.tol = 1E-4)$value
}
intt161
xxi16=intt161/ddd2;xxi16

#Complete Information Matrix
#alpha1
I11=-(n1+n31)/(ml[1]^2)-ml[2]*sum((z1^ml[1])*(log(z1)^2))-ml[2]*sum((z2^ml[1])*(log(z2)^2))-ml[2]*sum(xxi7)-ml[2]*sum(xxi15)-ml[2]*sum((S4^ml[1])*(log(S4)^2));I11
I12=-sum((z1^ml[1])*log(z1))-sum((z2^ml[1])*log(z2))-sum(xxi4)-sum(xxi12)-sum((S4^ml[1])*log(S4));I12
I13=0;I13
I14=0;I14
I15=0;I15
I16=0;I16

#beta1
I21=I12;I21
I22=-(n1+n31)/(ml[2]^2);I22
I23=0;I23
I24=0;I24
I25=0;I25
I26=0;I26

#lambda1  
I31=0;I31
I32=0;I32
I33=-n31/(ml[3]^2);I33
I34=0;I34
I35=0;I35
I36=0;I36

#alpha2
I41=0;I41
I42=0;I42
I43=0;I43
I44=-(n2+n32)/(ml[4]^2)-ml[5]*sum((z1^ml[4])*(log(z1)^2))-ml[5]*sum((z2^ml[4])*(log(z2)^2))-ml[5]*sum(xxi8)-ml[5]*sum(xxi16)-ml[5]*sum((S4^ml[4])*(log(S4)^2));I44
I45=-sum((z1^ml[4])*log(z1))-sum((z2^ml[4])*log(z2))-sum(xxi5)-sum(xxi13)-sum((S4^ml[4])*log(S4));I45
I46=0;I46

#beta2
I51=0;I51 
I52=0;I52
I53=0;I53
I54=I45;I54
I55=-(n2+n32)/(ml[5]^2);I55
I56=0;I56

#lambda2
I61=0;I61
I62=0;I62
I63=0;I63
I64=0;I64
I65=0;I65
I66=-(n32)/(ml[6]^2);I66

Mat1=-matrix(c(I11,I12,I13,I14,I15,I16,I21,I22,I23,I24,I25,I26,I31,I32,I33,I34,I35,I36,I41,I42,I43,I44,I45,I46,
               I51,I52,I53,I54,I55,I56,I61,I62,I63,I64,I65,I66),6,6,byrow=TRUE);Mat1


#Incomplete Information Matrix
#grad
GD1=(n1+n31)/ml[1]+sum(log(z1))-ml[2]*sum((z1^ml[1])*log(z1))-ml[2]*sum((z2^ml[1])*log(z2))+sum(xxi3)-ml[2]*sum(xxi4)-ml[2]*sum(xxi12)-ml[2]*sum((S4^ml[1])*log(S4));GD1
GD2=(n1+n31)/ml[2]-sum(z1^ml[1])-sum(z2^ml[1])-sum(xxi1)-sum(xxi9)-sum(S4^ml[1]);GD2
GD3=n31/ml[3]-sum(S1-z1)-sum(xxi6);GD3
GD4=(n2+n32)/ml[4]+sum(log(z2))-ml[5]*sum((z1^ml[4])*log(z1))-ml[5]*sum((z2^ml[4])*log(z2))-ml[5]*sum(xxi5)+sum(xxi11)-ml[5]*sum(xxi13)-ml[5]*sum((S4^ml[4])*log(S4));GD4
GD5=(n2+n32)/ml[5]-sum(z1^ml[4])-sum(z2^ml[4])-sum(xxi2)-sum(xxi10)-sum(S4^ml[4]);GD5
GD6=(n32)/ml[6]-sum(S2-z2)-sum(xxi14);GD6

GD=c(GD1,GD2,GD3,GD4,GD5,GD6);GD

Mat2=GD%*%t(GD);Mat2
Var_Cov=solve(Mat1-Mat2);Var_Cov

SE=c(sqrt(Var_Cov[1,1]),sqrt(Var_Cov[2,2]),sqrt(Var_Cov[3,3]),
     sqrt(Var_Cov[4,4]),sqrt(Var_Cov[5,5]),sqrt(Var_Cov[6,6]));SE

#-------------100(1-alpha)% CIs--------#  
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
ACI_lim=c(cl1,cl2,cl3,cl4,cl5,cl6);round(ACI_lim,4)

l1=c(ci_u[1]-ci_l[1]);l1
l2=c(ci_u[2]-ci_l[2]);l2
l3=c(ci_u[3]-ci_l[3]);l3
l4=c(ci_u[4]-ci_l[4]);l4
l5=c(ci_u[5]-ci_l[5]);l5
l6=c(ci_u[6]-ci_l[6]);l6
ACI_L=c(l1,l2,l3,l4,l5,l6);ACI_L

se=c(sqrt(diag(Var_Cov)[1]),sqrt(diag(Var_Cov)[2]),sqrt(diag(Var_Cov)[3]),sqrt(diag(Var_Cov)[4]),sqrt(diag(Var_Cov)[5]),sqrt(diag(Var_Cov)[6]));se

#===================================================================================#
#
#               Bayesian Estimations
#====================================================================================#
#Define hyper-parameters
mu1=mu2=mu3=mu4=mu5=mu6=0.11;nu1=1/mu1;nu2=1/mu2;
nu3=1/mu3;nu4=1/mu4;nu5=1/mu5;nu6=1/mu6
b.thh1=b.thh2=b.thh3=b.thh4=b.thh5=b.thh6=rep()
b.th1=b.th2=b.th3=b.th4=b.th5=b.th6=rep()
b.th1=ml[1];b.th2=ml[2];b.th3=ml[3];b.th4=ml[4];b.th5=ml[5];b.th6=ml[6]
it1=250000 #iterations
for(j in 1:it1){

  ww1=runif(length(S31));ww1
  ww2=runif(length(S32));ww2
  
  tt1=(-(1/b.th2)*log(exp(-b.th2*(lower1^b.th1))-ww1*(exp(-b.th2*(lower1^b.th1))-exp(-b.th2*(S31^b.th1)))))^(1/b.th1);tt1  
  tt2=(-(1/b.th5)*log(exp(-b.th5*(lower2^b.th4))-ww2*(exp(-b.th5*(lower2^b.th4))-exp(-b.th5*(S32^b.th4)))))^(1/b.th4);tt2  
  
  u_i=-(1/b.th3)*log(exp(-b.th3*lower1)-ww1*(exp(-b.th3*lower1)-exp(-b.th3*(S31-tt1))));u_i 
  v_i=-(1/b.th6)*log(exp(-b.th6*lower2)-ww2*(exp(-b.th6*lower2)-exp(-b.th6*(S32-tt2))));v_i
  
  post1=function(al1,be1){
    zz1=(al1^(n1+n31+mu1-1))*(prod(z1^(al1-1)))*(prod((tt1^n31)^(al1-1)))
    zz2=exp(-be1*sum(z1^al1)-be1*sum(z2^al1)-be1*sum((tt1^al1))-be1*sum((tt2^al1))-be1*sum(S4^al1)-al1*nu1)
    zzz=zz1*zz2
    return(zzz)}
  
  last1=b.th1
  cand1=abs(rnorm(1,ml[1],Var_Cov[1,1]))
  #cand1=abs(rnorm(1,ml[1],0.5))
  r1=ifelse(is.finite(post1(cand1,b.th2)/post1(last1,b.th2))==TRUE,post1(cand1,b.th2)/post1(last1,b.th2),0.01)
  if(runif(1)<min(r1,1)) last1<-cand1
  b.thh1[j+1]=last1
  
  b.thh2[j+1]<-rgamma(1,shape=n1+n31+mu2,rate=sum(z1^b.thh1[j+1])+sum(z2^b.thh1[j+1])+sum((tt1^b.thh1[j+1]))+sum((tt2^b.thh1[j+1]))+sum(S4^b.thh1[j+1])+nu2)
  b.thh3[j+1]<-rgamma(1,shape=n31+mu3,rate=sum(S1-z1)+sum(u_i)+nu3)
  
  post2=function(al2,be2){
    ww11=(al2^(n2+n32+mu4-1))*(prod(z2^(al2-1)))*(prod((tt2^(n32))^(al2-1)))
    ww22=exp(-be2*sum(z1^al2)-be2*sum(z2^al2)-be2*sum((tt1^al2))-be2*sum((tt2^al2))-be2*sum(S4^al2)-al2*nu4)
    www=ww11*ww22
    return(www)}
  
  last2=b.th4
  cand2=abs(rnorm(1,ml[4],Var_Cov[4,4]))
  r2=ifelse(is.finite(post2(cand2,b.th5)/post2(last2,b.th5))==TRUE,post2(cand2,b.th5)/post2(last2,b.th5),0.01)
  if(runif(1)<min(r2,1)) last2<-cand2
  b.thh4[j+1]=last2
  
  b.thh5[j+1]<-rgamma(1,shape=n2+n32+mu5,rate=sum(z1^b.thh4[j+1])+sum(z2^b.thh4[j+1])+sum((tt1^b.thh4[j+1]))+sum((tt2^b.thh4[j+1]))+sum(S4^b.thh4[j+1])+nu5)
  b.thh6[j+1]<-rgamma(1,shape=n32+mu6,rate=sum(S2-z2)+sum(v_i)+nu6)
  
  b.thh1[j]=b.th1;b.thh2[j]=b.th2;b.thh3[j]=b.th3;b.thh4[j]=2+b.th4;b.thh5[j]=b.th5;b.thh6[j]=b.th6
  b.th1=b.thh1[j+1];b.th2=b.thh2[j+1];b.th3=b.thh3[j+1];b.th4=b.thh4[j+1];b.th5=b.thh5[j+1];b.th6=b.thh6[j+1] 
  
}
b.thh1;b.thh2;b.thh3;b.thh4;b.thh5;b.thh6
length(b.thh1);length(b.thh2);length(b.thh3);length(b.thh4);length(b.thh5);length(b.thh6)

#discard first few values 
ch1=b.thh1[10001:it1];ch2=b.thh2[10001:it1];ch3=b.thh3[10001:it1];ch4=b.thh4[10001:it1];ch5=2+b.thh5[10001:it1];ch6=b.thh6[10001:it1]

#checking ACF
acf(ch1);acf(ch2);acf(ch3);acf(ch4);acf(ch5);acf(ch6)

#Remove Lag in MCMC chains
zz1=seq(150,length(ch1),150);zz2=seq(90,length(ch2),90)
ch111=ch1[zz1];ch222=ch2[zz2];ch333=ch3[zz2];ch444=ch4[zz2];ch555=ch5[zz2];ch666=ch6[zz2]
length(ch111);length(ch222);length(ch333);length(ch444);length(ch555);length(ch666)
ch11=ch111[1:1500];ch22=ch222[1:1500];ch33=ch333[1:1500];ch44=ch444[1:1500];ch55=ch555[1:1500];ch66=ch666[1:1500]

#quantile plots
gx1=na.omit(as.vector(ch11));gx1
ii=1:length(gx1)
q1.25=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.25)));q1.25
q1.50=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.5)));q1.50
q1.75=sapply(ii,function(ii) quantile((gx1[1:ii]), probs = c(.75)));q1.75

gx2=na.omit(as.vector(ch22));gx2
ii=1:length(gx2)
q2.25=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.25),na.rm=TRUE));q2.25
q2.50=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.5),na.rm=TRUE));q2.50
q2.75=sapply(ii,function(ii) quantile((gx2[1:ii]), probs = c(.75),na.rm=TRUE));q2.75

gx3=as.vector(ch33)
ii=1:length(gx3)
q3.25=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.25)));q3.25
q3.50=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.5)));q3.50
q3.75=sapply(ii,function(ii) quantile((gx3[1:ii]), probs = c(.75)));q3.75

gx4=as.vector(ch44);gx4
ii=1:length(gx4)
q4.25=sapply(ii,function(ii) quantile((gx4[1:ii]), probs = c(.25)));q4.25
q4.50=sapply(ii,function(ii) quantile((gx4[1:ii]), probs = c(.5)));q4.50
q4.75=sapply(ii,function(ii) quantile((gx4[1:ii]), probs = c(.75)));q4.75

gx5=na.omit(as.vector(ch55));gx5
ii=1:length(gx5)
q5.25=sapply(ii,function(ii) quantile((gx5[1:ii]), probs = c(.25),na.rm=TRUE));q5.25
q5.50=sapply(ii,function(ii) quantile((gx5[1:ii]), probs = c(.5),na.rm=TRUE));q5.50
q5.75=sapply(ii,function(ii) quantile((gx5[1:ii]), probs = c(.75),na.rm=TRUE));q5.75

gx6=as.vector(ch66)
ii=1:length(gx6)
q6.25=sapply(ii,function(ii) quantile((gx6[1:ii]), probs = c(.25)));q6.25
q6.50=sapply(ii,function(ii) quantile((gx6[1:ii]), probs = c(.5)));q6.50
q6.75=sapply(ii,function(ii) quantile((gx6[1:ii]), probs = c(.75)));q6.75

par(mfrow=c(3,2))

#Trace Plot
matplot(ii,cbind(q1.25,q1.50,q1.75),main="",type="l",col=1,xlab="iteration",ylab=expression(alpha[1]))
matplot(ii,cbind(q2.25,q2.50,q2.75),main="",type="l",col=1,xlab="iteration",ylab=expression(beta[1]))
matplot(ii,cbind(q3.25,q3.50,q3.75),main="",type="l",col=1,xlab="iteration",ylab=expression(lambda[1]))
matplot(ii,cbind(q4.25,q4.50,q4.75),main="",type="l",col=1,xlab="iteration",ylab=expression(alpha[2]))
matplot(ii,cbind(q5.25,q5.50,q5.75),main="",type="l",col=1,xlab="iteration",ylab=expression(beta[2]))
matplot(ii,cbind(q6.25,q6.50,q6.75),main="",type="l",col=1,xlab="iteration",ylab=expression(lambda[2]))

#acf plots
acf(ch11,main="",col="black",xlab="Lag",ylab=expression("ACF"(alpha[1])))
acf(ch22,main="",col="black",xlab="Lag",ylab=expression("ACF"(beta[1])))
acf(ch33,main="",col="black",xlab="Lag",ylab=expression("ACF"(lambda[1])))
acf(ch44,main="",col="black",xlab="Lag",ylab=expression("ACF"(alpha[2])))
acf(ch55,main="",col="black",xlab="Lag",ylab=expression("ACF"(beta[2])))
acf(ch66,main="",col="black",xlab="Lag",ylab=expression("ACF"(lambda[2])))

#matplots
matplot(ch11,type="l",col="black",xlab="iterations",ylab=expression(alpha[1]),main="")
matplot(ch22,type="l",col="black",xlab="iterations",ylab=expression(beta[1]),main="")
matplot(ch33,type="l",col="black",xlab="iterations",ylab=expression(lambda[1]),main="")
matplot(ch44,type="l",col="black",xlab="iterations",ylab=expression(alpha[2]),main="")
matplot(ch55,type="l",col="black",xlab="iterations",ylab=expression(beta[2]),main="")
matplot(ch66,type="l",col="black",xlab="iterations",ylab=expression(lambda[2]),main="")

#density plots
plot(density(ch11),col="black",type="l",ylab=expression("Density"(alpha[1])),main="")
plot(density(ch22),col="black",type="l",ylab=expression("Density"(beta[1])),main="")
plot(density(ch33),col="black",type="l",ylab=expression("Density"(lambda[1])),main="")
plot(density(ch44),col="black",type="l",ylab=expression("Density"(alpha[2])),main="")
plot(density(ch55),col="black",type="l",ylab=expression("Density"(beta[2])),main="")
plot(density(ch66),col="black",type="l",ylab=expression("Density"(lambda[2])),main="")

#Bayes estimates
b.th1=b.th2=b.th3=b.th4=b.th5=b.th6=rep()
b.th1=mean(ch11);b.th2=mean(ch22);b.th3=mean(ch33);b.th4=mean(ch44);b.th5=mean(ch55);b.th6=mean(ch66)
b.th=c(b.th1,b.th2,b.th3,b.th4,b.th5,b.th6);b.th  

#hpd intervals
h.th1=HPDinterval(mcmc(ch11));h.th2=HPDinterval(mcmc(ch22));h.th3=HPDinterval(mcmc(ch33));h.th4=HPDinterval(mcmc(ch44));h.th5=HPDinterval(mcmc(ch55));h.th6=HPDinterval(mcmc(ch66))
HPD_lim=c(h.th1[,1],h.th1[,2],h.th2[,1],h.th2[,2],h.th3[,1],h.th3[,2],h.th4[,1],h.th4[,2],h.th5[,1],h.th5[,2],h.th6[,1],h.th6[,2]);HPD_lim

h.th1_l=h.th1[,2]-h.th1[,1];h.th2_l=h.th2[,2]-h.th2[,1];h.th3_l=h.th3[,2]-h.th3[,1];h.th4_l=h.th4[,2]-h.th4[,1];h.th5_l=h.th5[,2]-h.th5[,1];h.th6_l=h.th6[,2]-h.th6[,1]
HPD_L=c(h.th1_l,h.th2_l,h.th3_l,h.th4_l,h.th5_l,h.th6_l);HPD_L                                #length of HPD

#Mean and Median Durations 
Mean11=mean(aa*(ch22^(-1/ch11))*gamma(1+1/ch11));Mean11  
Mean22=mean(aa*(ch55^(-1/ch44))*gamma(1+1/ch44));Mean22  
Med11=mean(aa*(ch22^(-1/ch11))*(log(2)^(1/ch11)));Med11  
Med22=mean(aa*(ch55^(-1/ch44))*(log(2)^(1/ch44)));Med22  
