from Functions import *
from timeit import default_timer as timer
from threading import Thread

class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None
    def run(self):
        print(type(self._target))
        if self._target is not None:
            self._return = self._target(*self._args,
                                                **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return

if __name__ == '__main__':
    # start = timer()
    # r = 3
    # p = 30
    # par = 2*r/p
    # zc = 0.5*(600*10**(-9))
    # zp = 0.1*(600*10**(-9))
    # xs = np.arange(-r, r+par, par)*(600*10**(-9))
    # ys = np.arange(-r, r+par, par)*(600*10**(-9))
    # x,y = np.meshgrid(xs,ys)

    # s0_values = S0v_dipole(x, y, zc, zp=zp, data=data1)
    # s_values = Sv_dipole(x, y, zc, zp=zp, data=data1)
    # s_values_notopo = Sv_dipole(x, y, zc, zp=zp, data=data0) 

    # load_data('test.npy', [s0_values, s_values, s_values_notopo])

    # end = timer()
    # print('elapsed time: {} min'.format((end - start)/60))
    
    start = timer()
    zcs = np.arange(0, 3.05, 0.05)*(600*10**(-9))
    zps = np.arange(0, 3.05, 0.05)*(600*10**(-9))
    zc,zp = np.meshgrid(zcs,zps)
    y_range = 'Place Holder'

    t1 = ThreadWithReturnValue(target=Sv_dipole_max, args=(0, y_range,zc,), kwargs={'zp': zp, 'data': data1})
    t2 = ThreadWithReturnValue(target=Sv_dipole_max, args=(0, y_range,zc,), kwargs={'zp': zp, 'data': data0})

    # starting thread 1
    t1.start()
    # starting thread 2
    t2.start()
  
    # wait until both threads are completely executed
    S_values = t1.join() 
    S_notopo_values = t2.join() 

    load_data('test2.npy', [S_values, S_notopo_values])

    end = timer()
    print('elapsed time: {} min'.format((end - start)/60))

