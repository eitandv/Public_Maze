import numpy as np
from math import pi
from scipy.integrate import quad
from scipy.special import jv

a = 1/137
c = 1
m1 = 1
m2 = m1*(np.sqrt(0.48**2 + 17.3**2)*4*pi*10**(-7) + 1)
e1 = 1
eaux = 8.85418782*10**(-12)
maux = 1.25663706*10**(-6)
n1 = 1
n2 = 1.65
theta1 = 0
theta2 = pi
delta = theta1 - theta2
skip = 15
scale = (600/(2*pi*c))**(-1)
w = 1
e2 = 1.65**2/(c**2*m2)

def f1(theta):
    return m1*m2*(e2*m1 + e1*m2 + n1*n2*(np.sqrt((1-np.sin(theta)**2)/(1-(n1**2/n2**2)*np.sin(theta)**2)) + np.sqrt((1-(n1**2/n2**2)*np.sin(theta)**2)/(1-np.sin(theta)**2))))

def f2(theta):
    return m1*m2*(e2*m1 - e1*m2 + n1*n2*(np.sqrt((1-np.sin(theta)**2)/(1-(n1**2/n2**2)*np.sin(theta)**2)) - np.sqrt((1-(n1**2/n2**2)*np.sin(theta)**2)/(1-np.sin(theta)**2))))

def f3(theta):
    return m1*m2*(e1*m2 - e2*m1 + n1*n2*(np.sqrt((1-np.sin(theta)**2)/(1-(n1**2/n2**2)*np.sin(theta)**2)) - np.sqrt((1-(n1**2/n2**2)*np.sin(theta)**2)/(1-np.sin(theta)**2))))

def R_tmte(theta): 
    return (-2*m2*n1*delta)/(f1(theta) + delta**2)

def R_tmtm(theta): 
    return (f2(theta)+delta**2)/(f1(theta)+delta**2)

def R_tete(theta):
    return (f3(theta)+delta**2)/(f1(theta)+delta**2)

def R_tetm(theta):
    return (-2*m2*n1*delta)/(f1(theta) + delta**2)

def rho_z_integrand(theta, z): 
    return 3*(n1**3/4)*(np.sin(theta)**3)*((-2*m2*n1*delta)/(f1(theta) + delta**2)**2 + (np.abs(1+np.e**(2j*n1*w*np.cos(theta)*z/c)*(f2(theta)+delta**2)/(f1(theta)+delta**2)))**2)

def rho_x_integrand(theta, z): 
    return 3*(n1**3/4)*(np.sin(theta)**3)*((-2*m2*n1*delta)/(f1(theta) + delta**2)**2 + (np.abs(1-np.e**(2j*n1*w*np.cos(theta)*z/c)*(f2(theta)+delta**2)/(f1(theta)+delta**2)))**2)

def rho_y_integrand(theta, z): 
    return 3*(n1**3/4)*(np.sin(theta)**3)*((-2*m2*n1*delta)/(f1(theta) + delta**2)**2 + (np.abs(1+np.e**(2j*n1*w*np.cos(theta)*z/c)*(f3(theta)-delta**2)/(f1(theta)+delta**2)))**2)

def rho_z(z):
    return quad(rho_z_integrand, 0, pi/2, args=(z,))

def rho_x(z):
    return quad(rho_x_integrand, 0, pi/2, args=(z,))

def rho_y(z):
    return quad(rho_y_integrand, 0, pi/2, args=(z,))




def integrand_rho_z(theta, z):
    return ((np.sin(theta))**3)*(R_tmte(theta)**2 + np.abs((np.e**(-1j*n1*w*np.cos(theta)*z/c) + R_tmtm(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)) - ((R_tmte(theta)*(R_tmtm(theta) + R_tete(theta) + (2*np.sin(theta)**2)*np.e**(2j*n1*w*np.cos(theta)*z/c)))/(R_tmte(theta)**2+((1+np.cos(2*n1*w*np.cos(theta)*z/c))*R_tete(theta) + 1)**2))*(R_tmte(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)))**2)

def integrand_rho_x(theta, z):
    return ((np.sin(theta))*np.cos(theta)**2)*(R_tmte(theta)**2 + np.abs((np.e**(-1j*n1*w*np.cos(theta)*z/c) - R_tmtm(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)) - ((R_tmte(theta)*(R_tmtm(theta) + R_tete(theta) + (2*np.sin(theta)**2)*np.e**(2j*n1*w*np.cos(theta)*z/c)))/(R_tmte(theta)**2+((1+np.cos(2*n1*w*np.cos(theta)*z/c))*R_tete(theta) + 1)**2))*(-1*R_tmte(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)))**2)

def integrand_rho_y(theta, z):
    return (np.abs(1+R_tete(theta)*np.e**(-2j*n1*w*np.cos(theta)*z/c))**2 + np.abs(((-1)*R_tetm(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)) + ((R_tmte(theta)*(R_tmtm(theta) + R_tete(theta) + (2*np.sin(theta)**2)*np.e**(2j*n1*w*np.cos(theta)*z/c)))/(R_tmte(theta)**2+((1+np.cos(2*n1*w*np.cos(theta)*z/c))*R_tete(theta) + 1)**2))*(np.e**(-1j*n1*w*np.cos(theta)*z/c) + R_tete(theta)*np.e**(1j*n1*w*np.cos(theta)*z/c)))**2)

def Rtmtm(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c):
  thingy = ( (e2*np.sqrt(k1**2-kp**2)-e1*np.sqrt(k2**2-kp**2))*(m1*m2*(m2*np.sqrt(k1**2-kp**2)+m1*np.sqrt(k2**2-kp**2))) + np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta**2)/()
  return thingy

def Rtetm(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c):
  thingy = -2*m2*n1*np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta/( (e2*np.sqrt(k1**2-kp**2) + e1*np.sqrt(k2**2-kp**2) )*(m1*m2*(m2*np.sqrt(k1**2-kp**2)+m1*np.sqrt(k2**2-kp**2))) + np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta**2)
  return thingy

def Rtete(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c):
  thingy = ( (m2*np.sqrt(k1**2-kp**2)-m1*np.sqrt(k2**2-kp**2) )*(m1*m2*(e2*np.sqrt(k1**2-kp**2)+e1*np.sqrt(k2**2-kp**2))) + np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta**2)/( (m2*np.sqrt(k1**2-kp**2) + m1*np.sqrt(k2**2-kp**2) )*(m1*m2*(e2*np.sqrt(k1**2-kp**2)+e1*np.sqrt(k2**2-kp**2))) + np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta**2)
  return thingy

def Rtmte(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c):
  thingy = -2*m2*n1*np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta/( (m2*np.sqrt(k1**2-kp**2) + m1*np.sqrt(k2**2-kp**2) )*(m1*m2*(e2*np.sqrt(k1**2-kp**2)+e1*np.sqrt(k2**2-kp**2))) + np.sqrt(k1**2-kp**2)*np.sqrt(k2**2-kp**2)*delta**2)
  return thingy

def G_z_integrand(kp, k, x, y, z): 
    return (1j/(4*pi))*(np.e**(1j*np.sqrt(k**2-kp**2)*(z + 1.5)))*jv(0,kp*np.sqrt(x**2+y**2))*(kp**3/(k**2*np.sqrt(k**2-kp**2)))*Rtmtm(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c)

def G_x_integrand(theta, z): 
    return None

def G_y_integrand(theta, z): 
    return None

def G_z(k, x, y, z):
    return quad(rho_z_integrand, 0, 1, args=(k,x,y,z,))

def G_x(z):
    return quad(rho_x_integrand, 0, 1, args=(z,))

def G_y(z):
    return quad(rho_y_integrand, 0, 1, args=(z,))

def GRxz_integrand(kp, x, y, z, zp, data=data1):
    k1 = data["k1"]
    return np.imag((1/(4*pi))*(np.e**(1j*np.sqrt(k1**2-kp**2)*(z + zp)))*jv(1,kp*np.sqrt(x**2+y**2))*(kp**2/(k1*np.sqrt(k1**2-kp**2)))*Rtetm(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"]))

def GRxx_integrand(kp, x, y, z, zp, data=data1):
    rtmtm = Rtmtm(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"])
    rtete = Rtete(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"])
    k1 = data["k1"]
    return np.imag((1j/(8*pi))*(np.e**(1j*np.sqrt(k1**2-kp**2)*(z + zp)))*((jv(0,kp*np.sqrt(x**2+y**2)) + jv(2,kp*np.sqrt(x**2+y**2)))*(kp/(np.sqrt(k1**2-kp**2)))*rtete - (jv(0,kp*np.sqrt(x**2+y**2)) - jv(2,kp*np.sqrt(x**2+y**2)))*(kp*np.sqrt(k1**2-kp**2)/k1**2)*rtmtm))

def GRyx_integrand(kp, x, y, z, zp, data=data1):
    rtetm = Rtetm(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"])
    rtmte = Rtmte(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"])
    k1 = data["k1"]
    return np.imag((-1j/(8*pi))*(np.e**(1j*np.sqrt(k1**2-kp**2)*(z + zp)))*((jv(0,kp*np.sqrt(x**2+y**2)) - jv(2,kp*np.sqrt(x**2+y**2)))*(kp/k1)*rtetm + (jv(0,kp*np.sqrt(x**2+y**2)) + jv(2,kp*np.sqrt(x**2+y**2)))*(kp/k1)*rtmte))

def GRzx_integrand(kp, x, y, z, zp, data=data1):
    rtmtm = Rtmtm(kp, data["e1"], data["e2"], data["k1"], data["k2"], data["m1"], data["m2"], data["n1"], data["n2"], data["delta"], data["c"])
    k1 = data["k1"]
    return np.imag((-1/(4*pi))*(np.e**(1j*np.sqrt(k1**2-kp**2)*(z + zp)))*jv(1,kp*np.sqrt(x**2+y**2))*(kp**2/k1**2)*rtmtm)

def ReflectiveMatrix(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c):
    # print(e1, e2, k1, k2, m1, m2, n1, n2)
    return np.array([[Rtete(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c), Rtetm(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c)], [Rtmte(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c), Rtmtm(kp, e1, e2, k1, k2, m1, m2, n1, n2, delta, c)]])