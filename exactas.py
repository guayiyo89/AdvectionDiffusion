import numpy as np

a0 = 1
sigma = 2
lamb = 0.1
cx = 1.5
cy = 1.5
xc = -2
yc = -4
s = 1.0

x0 = 2.5
y0 = 2.5
x1 = 0
y1 = 0
tau0 = 7.0
tau1 = 0.001
alpha = 0.4

d = 0.001
kappa = 0.01

def vx(x,y,t,T0):
    vel = cx + lamb*(x-y) * np.sin(t/T0)
    return vel

def vy(x,y,t,T0):
    vel = cy + lamb*(x+y) * np.sin(t/T0)
    return vel

def delta(x):
    delt = 1 / (alpha * np.sqrt(np.pi)) * np.exp(-(x**2) / (alpha**2))
    return delt

def exacta(xi,yi,t,T0):
    exac = a0 + np.exp( -sigma * ( (xi - vx(xi,yi,t,T0)*t - xc)**2 + (yi - vy(xi,yi,t,T0)*t - yc)**2 ) ) +  s * t* (delta(xi-t-x0) * delta(yi-t-y0) * delta( t - tau0 ) + delta(xi-t-x1)*  delta(yi-t-y1) * delta( t - tau1 )) 
    return exac

def source(x,y,t,T0):
    src = -d*((sigma**2)*(2*(1-t*np.sin(t/T0)*lamb)*(-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)+2*t*np.sin(t/T0)*lamb*(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x))**2 * np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2)) + (sigma**2)*(2*(1-t*np.sin(t/T0)*lamb)*(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)-2*t*np.sin(t/T0)*lamb*(-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y))**2 * np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2 + (-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))-2*sigma*(2*((1-t*np.sin(t/T0)*lamb)**2)+2 * (t**2) * (np.sin(t/T0))**2 * lamb**2 )*np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))+(4*s*t*((-y1+y-t)**2)*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**7)*np.pi**(3/2))+(4*s*t*((-x1+x-t)**2)*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**7)*np.pi**(3/2))-(4*s*t*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2))+(4*s*t*((-y0+y-t)**2)*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**7)*np.pi**(3/2))+(4*s*t*(-x0+x-t)**2*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**7)*np.pi**(3/2))-(4*s*t*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2)))+(np.sin(t/T0)*(y+x)*lamb+cy)*(-sigma*(2*(1-t*np.sin(t/T0)*lamb)*(-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)+2*t*np.sin(t/T0)*lamb*(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x))*np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))-(2*s*t*(-y1+y-t)*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2))-(2*s*t*(-y0+y-t)*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2)))+(np.sin(t/T0)*(x-y)*lamb+cx)*(-sigma*(2*(1-t*np.sin(t/T0)*lamb)*(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)-2*t*np.sin(t/T0)*lamb*(-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y))*np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))-(2*s*t*(-x1+x-t)*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2))-(2*s*t*(-x0+x-t)*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**5)*np.pi**(3/2)))+kappa*(np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))+(s*t*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))+(s*t*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/alpha**2-(t-tau0)**2/alpha**2))/((alpha**3)*np.pi**(3/2))+a0)-sigma*(2*(-np.sin(t/T0)*(y+x)*lamb-(t*np.cos(t/T0)*(y+x)*lamb)/T0-cy)*(-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)+2*(-np.sin(t/T0)*(x-y)*lamb-(t*np.cos(t/T0)*(x-y)*lamb)/T0-cx)*(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x))*np.exp(-sigma*((-t*(np.sin(t/T0)*(y+x)*lamb+cy)-yc+y)**2+(-t*(np.sin(t/T0)*(x-y)*lamb+cx)-xc+x)**2))+(s*t*((2*(-y1+y-t))/(alpha**2)+(2*(-x1+x-t))/(alpha**2)-(2*(t-tau1))/(alpha**2))*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))+(s*np.exp(-(-y1+y-t)**2/(alpha**2)-(-x1+x-t)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))-(s*np.exp(-(y-y1)**2/(alpha**2)-(x-x1)**2/(alpha**2)-(t-tau1)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))+(s*t*((2*(-y0+y-t))/(alpha**2)+(2*(-x0+x-t))/(alpha**2)-(2*(t-tau0))/(alpha**2))*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))+(s*np.exp(-(-y0+y-t)**2/(alpha**2)-(-x0+x-t)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))-(s*np.exp(-(y-y0)**2/(alpha**2)-(x-x0)**2/(alpha**2)-(t-tau0)**2/(alpha**2)))/((alpha**3)*np.pi**(3/2))
    return src

vall = source(2,3,45,2)