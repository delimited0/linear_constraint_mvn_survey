import numpy as np
from scipy.interpolate import CubicSpline


#If we are generating our own set of x values:
#Number of data points:


###################
#If we are importing the x values from somwhere else:
xval=np.loadtxt("Q2.txt")

#Number of files with noisy data you want to generate, each file is 9 fold since we have 9 generators:

files=1


#Size of the Gaussian Error:

sig=0.005

#Size of the percentage error:
#(0.5 means 0.5% error)
per=0.5


###################
#Physical constants

#4*Proton Mass2 in fm2 inverse
mp=90.4478

#Converter from GeV2 to fm-2

co=25.689

###################

#Following the "Robust" Paper we have 9 generators:

#Dipole
def f1(x,p1):
  return (1+x*1.0/p1)**(-2)

#Monopole
def f2(x,p1):
  return (1+x*1.0/p1)**(-1)

#Gaussian
def f3(x,p1):
  return np.exp(-x/p1)


#Kelly2004
#a1= -0.24, bi = 10.98. b2 = 12.8, b3 = 21.97, rp= 0.863
#for 0.85 -0.111, 10.781, 12.82, 21.97
def f4(x,a1,b1,b2,b3):
  return (1+a1*x/(mp))/(1+1.0*b1*x/mp+b2*(x*1.0/mp)**2+b3*(x/mp)**3)

#Arrington2004
#Parameters in Gev:
#p2=3.226,p4=1.508, p6=-0.3773,p8=0.611,p10=-0.1853,p12=0.01596.
def f5(x,p2,p4,p6,p8,p10,p12):
  return 1.0/(1+p2*x+p4*x**2+p6*x**3+p8*x**4+p10*x**5+p12*x**6)


#Arrignton2007
#p1=3.440, p2= -0.178, p3= -1.212 . p4 = 1.176 , p5 = -0.284
def f6(x,p1,p2,p3,p4,p5):
  return 1.0/(1+ 1.0*x*p1/(1+ 1.0*x*p2/(1+ 1.0*x*p3/(1+ 1.0*x*p4/(1+ 1.0*x*p5  )  )  )  )  )


#Venkat2011
def f7(x,a1,a2,a3,b1,b2,b3,b4,b5):
  return 1.0*(1+a1*x/mp+a2*(x/mp)**2+a3*(x/mp)**3)/(1+b1*x/mp+b2*(x/mp)**2+b3*(x/mp)**3+b4*(x/mp)**4+b5*(x/mp)**5)


#Bernauer2014
#p1=-3.3686,p2=14.5606,p3=-88.1912,p4=453.6244,p5=-1638.7911,p6=3980.7114,p7=-6312.6333,p8=6222.3646,p9=-3443.2251,p10=+814.4112;
def f8(x,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10):
  return 1+x*p1+p2*x**2+p3*x**3+p4*x**4+p5*x**5+p6*x**6+p7*x**7+p8*x**8+p9*x**9+p10*x**10


#Alarcon2017 (WE ARE NOT GOING TO USE ALARCON FOR THIS TESTS):
#xal=np.loadtxt("AlarconXValues.txt")
#yal=np.loadtxt("AlarconYValues.txt")
#alcs=CubicSpline(xal,yal,bc_type='not-a-knot')


#Ye-2018

t0=-17.9823
tcut=2.001
def z(xc):
  return (np.sqrt(tcut+xc)-np.sqrt(tcut-t0))/(np.sqrt(tcut+xc)+np.sqrt(tcut-t0)) 
def f10(x,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12):
  return p0+z(x)*p1+p2*z(x)**2+p3*z(x)**3+p4*z(x)**4+p5*z(x)**5+p6*z(x)**6+p7*z(x)**7+p8*z(x)**8+p9*z(x)**9+p10*z(x)**10+p11*z(x)**11+p12*z(x)**12


#We define a master function F[x,i], where i is between 0 and 8 and selects which function to use. We set all the free parameters to match a radius of 0.85.


def F(x,i):
  if(i==0):
    return f1(x,12.0/(0.85**2))
  if(i==1):
    return f2(x,6.0/(0.85**2))
  if(i==2):
    return f3(x,6.0/(0.85**2))
  if(i==3):
    return f4(x,-0.111, 10.781, 12.82, 21.97)
  if(i==4):
    return f5(x,0.120416,1.508/(co**2),-0.3773/(co**3),0.611/(co**4),-0.1853/(co**5),0.01596/(co**6))
  if(i==5):
    return f6(x,3.09338/co,-0.178/co,-1.212/co,1.176/co,-0.284/co)
  if(i==6):
    return f7(x,3.10858,-1.1154222,0.03866171,14,40.8833,99.99999,0.00004579,10.35804)
  if(i==7):
    return f8(x,-3.09338/co,14.5606/co**2,-88.1912/co**3,453.6244/co**4,-1638.7911/co**5,3980.7114/co**6,-6312.6333/co**7,6222.3646/co**8,-3443.2251/co**9,814.4112/co**10)
  if(i==8):
    return 0
  if(i==9):
    return f10(x,0.27767280155272844, -1.0454893852325884, 1.4344338057502042,0.49035313941611075, -2.294998339931553, 1.132143141088303,1.2470712638274697, -3.628821019050907,4.080886445171562,0.504887446241421, -5.085582522974911, 3.9680089231149216,-0.9816813830273542)



yvals=[[],[],[],[],[],[],[],[],[],[],[]]
for j in range(len(xval)):
  yvals[0].append(xval[j])
for k in range(10):
    for j in range(len(xval)):
      yvals[k+1].append(F(xval[j],k))
    np.savetxt("YF"+str(k)+".txt",yvals[k+1])

np.savetxt("yvalsAll.txt",np.transpose(yvals))
delt=0.00001

for j in range(10):
  print(j)
  print(1.0*np.sqrt(-6*(F(delt,j)-F(0,j))/delt))

np.savetxt("xvalsGenerated.txt",yvals[0])

