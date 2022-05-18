from Functions import *
from timeit import default_timer as timer

if __name__ == '__main__':
    start = timer()
    r = 3
    p = 30
    par = 2*r/p
    zc = 0.5*(600*10**(-9))
    zp = 0.1*(600*10**(-9))
    xs = np.arange(-r, r+par, par)*(600*10**(-9))
    ys = np.arange(-r, r+par, par)*(600*10**(-9))
    x,y = np.meshgrid(xs,ys)

    s0_values = S0v_dipole(x, y, zc, zp=zp, data=data1)
    s_values = Sv_dipole(x, y, zc, zp=zp, data=data1)
    s_values_notopo = Sv_dipole(x, y, zc, zp=zp, data=data0) 

    load_data('test.npy', [s0_values, s_values, s_values_notopo])

    end = timer()
    print('elapsed time: {} min'.format((end - start)/60))

