import matplotlib.pyplot as plt
import numpy as np

def Euler_V1(T=10,dt=0.001,m=1,k=1):

  x,v,t,Ep,Ek,Ec,i = [0],[1],[0],[0],[1],[0.5],0

  while t[-1]<=T:
    F = -k*x[i]
    xn = x[i]+v[i]*dt+F/(2*m)*dt**2 
    vn = v[i] + F/m*dt
    nEp = (m*(vn)**2)/2
    nEk = (k*(xn)**2)/2
    nEc = nEp+nEk
    x.append(xn)
    v.append(vn)
    t.append(t[i]+dt)
    Ep.append(nEp)
    Ek.append(nEk)
    Ec.append(nEc)
    i+=1

  return x,v,t,Ep,Ek,Ec


def Euler_V2(T=10,dt=0.001,m=1,k=1):

  x,v,t,Ep,Ek,Ec,i = [0],[1],[0],[0],[1],[0.5],0

  while t[-1]<=T:
    F = -(1/6)*x[i]**3+x[i]
    xn = x[i]+v[i]*dt+F/(2*m)*dt**2 
    vn = v[i]+F/m*dt
    nEp = k*(1/24*xn**4-0.5*xn**2)
    nEk = (m*(vn)**2)/2
    nEc = nEp+nEk
    x.append(xn)
    v.append(vn)
    t.append(t[i]+dt)
    Ep.append(nEp)
    Ek.append(nEk)
    Ec.append(nEc)
    i+=1
  
  return x,v,t,Ep,Ek,Ec


def Verlet_V1(T=10,dt=0.001,m=1,k=1):

  x,v,t,Ep,Ek,Ec,i = [0],[1],[0],[0],[1],[0.5],0

  while t[-1]<=T:
    F = -k*x[i]
    if i == 0:
      vn=v[i]+(F/m)*dt
      xn=x[i]+vn*dt
      nEp = (m*(vn)**2)/2
      nEk = (k*(xn)**2)/2
      nEc = nEp+nEk
    else:
      xn=2*x[i]-x[i-1]+F/m*dt**2
      vn=(x[i]-x[i-2])/(m*2*dt)
      nEp = (m*(vn)**2)/2
      nEk = (k*(xn)**2)/2
      nEc = nEp+nEk
    x.append(xn)
    v.append(vn)
    t.append(t[i]+dt)
    Ep.append(nEp)
    Ek.append(nEk)
    Ec.append(nEc)
    i+=1

  return x,v,t,Ep,Ek,Ec


def Verlet_V2(T=10,dt=0.001,m=1,k=1):

  x,v,t,Ep,Ek,Ec,i = [0],[1],[0],[0],[1],[0.5],0

  while t[-1]<=T:
    F = -(1/6)*x[i]**3+x[i]
    if i == 0:
      vn=v[i]+(F/m)*dt
      xn=x[i]+vn*dt
      nEp = k*(1/24*xn**4-0.5*xn**2)
      nEk = (m*(vn)**2)/2
      nEc = nEp+nEk
    else:
      xn=2*x[i]-x[i-1]+F/m*dt**2
      vn=(x[i]-x[i-2])/(m*2*dt)
      nEp = k*(1/24*xn**4-0.5*xn**2)
      nEk = (m*(vn)**2)/2
      nEc = nEp+nEk
    x.append(xn)
    v.append(vn)
    t.append(t[i]+dt)
    Ep.append(nEp)
    Ek.append(nEk)
    Ec.append(nEc)
    i+=1
    
  return x,v,t,Ep,Ek,Ec

def Plot(x,v,t,Ep,Ek,Ec):
  
  plt.plot(t, x)
  plt.grid(True)
  plt.title('położenie w funkcji czasu', fontweight = 'light', fontsize = 16)
  plt.xlabel('t', fontweight = 'light', fontsize = 14)
  plt.ylabel('x', fontweight = 'light', fontsize = 14)
  plt.show()

  plt.plot(t,v)
  plt.grid(True)
  plt.title('pręskość w funkcji czasu', fontweight = 'light', fontsize = 16)
  plt.xlabel('t', fontweight = 'light', fontsize = 14)
  plt.ylabel('v', fontweight = 'light', fontsize = 14)
  plt.show()

  plt.plot(x,v)
  plt.grid(True)
  plt.title('położenie w funkcji pręskości', fontweight = 'light', fontsize = 16)
  plt.xlabel('v', fontweight = 'light', fontsize = 14)
  plt.ylabel('x', fontweight = 'light', fontsize = 14)
  plt.show()

  plt.plot(t,Ep)
  plt.plot(t,Ek)
  plt.plot(t,Ec)
  plt.grid(True)
  plt.title('energia w funkcji czasu', fontweight = 'light', fontsize = 16)
  plt.xlabel('t', fontweight = 'light', fontsize = 14)
  plt.ylabel('Ep, Ek & Ec', fontweight = 'light', fontsize = 14)
  plt.ylim([-2, 2.5])
  plt.show()


#V1 = 0.5*x**2
#V2 = (1/24)*x**4-0.5*x**2
#F1 = -x
#F2 = -(1/6)*x**3+x

Plot(*Euler_V1())
Plot(*Euler_V2())
Plot(*Verlet_V1())
Plot(*Verlet_V2())
