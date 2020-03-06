import ppmt
import ppmt_interface
import numpy as np

np.random.seed(1634120)

ndata,nvars = 1000,2

data = np.random.random(size=(ndata,nvars))
data = np.asfortranarray(data,dtype=np.float64)
w = np.ones(ndata)
w = np.asfortranarray(w,dtype=np.float64)

lo = 10

if False:
    dir = np.array([1.0,0.0])
    dir = np.asfortranarray(dir,dtype=np.float64)

    projection = data@dir
    r = ppmt.r_function(projection)

    for a,b in zip(projection,r):
        print(a,b,sep=',')

    print("LEGENDRE:")
    lg = ppmt.legendre_polinomial(r,lo)
    for row in lg:
        print(*row,sep=',')
    print("LEGENDRE:")

    er = np.empty(lo+1)
    for k in range(0,lo+1):
        er[k] = np.mean(lg[:,k])**2

    ret = 0.0
    for k in range(1,lo+1):
        ix = ((2*(k+1)+1)/2)*er[k]
        ret += ix
        print(er[k],ix,ret,sep=',')

    print("FORTRAN:")

    pp = ppmt.projection_index(projection,lo)
    pf = ppmt_interface.ppmt_ext.pp_index(dir,data,w,lo)

    print("PYTHON",pp,"FORTRAN:",pf)

    quit()

min_index = 1e-4 #None

if min_index is None:
    min_index = ppmt.bootstrap_ppmt(nvars,nboot=100,percentile=50,lo=10)
    print('min_index:',min_index)

print('ppmt.fit...')
iterations,y,snorm,z,sph,sph_inv,U,projections,zd2 = ppmt.fit(data,min_index=min_index)
yd2 = y.copy()
print('snorm min/max:',np.min(snorm),np.max(snorm))
print('ppmt.fit...DONE iterations=',iterations)

x_back = ppmt.back_transform(y,iterations,data,snorm,z,sph_inv,U,projections,zd2,yd2,trace=0)
new_samples = np.random.normal(size=(10000,nvars))
x2_back = ppmt.back_transform(new_samples,iterations,data,snorm,z,sph_inv,U,projections,zd2,yd2,trace=0)
print('x2_back.shape:',x2_back.shape)


error = (data-x_back)**2
print(np.mean(error))

import matplotlib.pyplot as plt

plt.subplot(2, 2, 1)
plt.scatter(data[:,0],data[:,1])
plt.subplot(2, 2, 2)
plt.scatter(y[:,0],y[:,1])
plt.subplot(2, 2, 3)
plt.scatter(x_back[:,0],x_back[:,1])
plt.subplot(2, 2, 4)
plt.scatter(x2_back[:,0],x2_back[:,1])
plt.show()
