import ppmt_interface
import numpy as np
import scipy
import scipy.stats

from sklearn.base import BaseEstimator, TransformerMixin

def generate_directions(d,n=100):
    v = np.random.normal(size=(n,d))
    s = np.sqrt(np.sum(v**2,axis=1))

    for i in range(d):
        v[:,i] /= s
        #v[:,i] /= np.linalg.norm(v[:,i])

    for i in range(n):
        v[i,:] /= np.linalg.norm(v[i,:])

    return v

def sphere(X,centered=True,trace=False):
    #This function will sphere the data. This means the data will now have a
    #mean of 0 and a covariance matrix of 1.
    n,p = X.shape
    muhat = np.mean(X,axis=0)

    if trace: print ("muhat",muhat)

    #covariance matrix
    cov = np.cov(X.T)
    if trace: print(cov)
    #eigenvalue/eigenvector decomposition
    D,V = np.linalg.eig(cov)
    if trace: print("D",D,np.diag(D))
    if trace: print("V",V)
    Dinv = np.linalg.inv(np.diag(D))
    if trace: print("Dinv",Dinv)
    sq = scipy.linalg.sqrtm(Dinv)

    S = V@(sq@V.T)

    if trace: print("S",S)


    Xc = X - muhat

    Z = Xc@S

    return Z,S

def bootstrap_ppmt(nvars,ndirs=10,nboot=100,nsamples=100000,percentile=50,lo=20):
    #the idea is to generate samples from mutigaussian and calculate the projectino index onto the main axis
    np.random.seed(1634120)

    #random samples from mutivariate Gaussian
    mean = np.zeros(nvars)
    cov = np.diag(np.ones(nvars))

    #random direction for projection
    directions = generate_directions(nvars,n=ndirs)
    #directions[:,0] =
    #directions = np.diag(np.ones(nvars))
    #directions = np.asfortranarray(directions,dtype=np.float64)

    w = np.ones(nsamples)
    w = np.asfortranarray(w,dtype=np.float64)

    #compute pp index
    indices = np.empty((nboot,ndirs))
    #indices2 = np.empty((nboot,ndirs))

    for i in range(nboot):
        samples = np.random.multivariate_normal(mean,cov,size=nsamples)
        #samples,S = sphere(samples)

        samples = np.asfortranarray(samples,dtype=np.float64)

        for j in range(ndirs):
            #dir = np.zeros(nvars)
            dir = directions[j].copy()
            dir = np.asfortranarray(dir,dtype=np.float64)

            #project
            projection = samples@dir

            #import matplotlib.pyplot as plt
            #plt.hist(projection)
            #plt.show()

            pi = ppmt_interface.ppmt_ext.pp_index(dir,samples,w,lo)
            #p2 = projection_index(projection,lo)
            #print(i,j,pi,p2,dir,np.mean(projection),np.var(projection),np.min(projection),np.max(projection))

            #assert p1==p2,"dir=%s INDEX: %f vs %f"%(dir,p1,p2)

            ##quit()
            #indices2[i,j] = p2
            indices[i,j] = pi

    #calculate the percentil of the distribution of generated indices
    #per1 = np.percentile(indices1,percentile,axis=0)
    per = np.percentile(indices,percentile,axis=0)
    #print(per1,per2)
    return np.mean(per)

class PPMT(TransformerMixin, BaseEstimator):
    def __init__(self,n_vars=None,min_index=1e-4,maxiters=10,seed=1634129,trace=0):
        self.n_vars = n_vars
        self.min_index = min_index
        self.maxiters = maxiters
        self.seed = seed
        self.trace = trace

        assert 0 <= maxiters <= 1000

    def fit(self, X, y=None):
        ndata,nvars = X.shape

        if self.n_vars is not None:
            assert nvars == self.n_vars
        else:
            self.n_vars = nvars


        #make them Fortranable
        data = np.asfortranarray(X,dtype=np.float64)
        self._yd2 = np.asfortranarray(np.empty_like(data),dtype=np.float64)
        self._snorm = np.asfortranarray(np.empty(ndata),dtype=np.float64)
        self._zd = np.asfortranarray(np.empty_like(data),dtype=np.float64)
        self._sph = np.asfortranarray(np.empty((nvars,nvars)),dtype=np.float64)
        self._sph_inv = np.asfortranarray(np.empty((nvars,nvars)),dtype=np.float64)
        self._Us = np.asfortranarray(np.empty((nvars,nvars,self.maxiters)),dtype=np.float64)
        self._projections = np.asfortranarray(np.empty((ndata,self.maxiters)),dtype=np.float64)
        self._zd2 = np.asfortranarray(np.empty((ndata,nvars)),dtype=np.float64)

        #print(ppmt_interface.ppmt_ext.ppmt_fit.__doc__)

        self._iterations = ppmt_interface.ppmt_ext.ppmt_fit(data,self.min_index,self.maxiters,self._yd2,self._snorm,self._zd,self._sph,self._sph_inv,self._Us,self._projections,self._zd2,self.seed,self.trace)

        #reshape objects to hold just actual iterations
        self._Us = self._Us[:,:,:self._iterations]
        self_projections = self._projections[:,:self._iterations]

        return self

    def fit_transform(self,X, y=None):
        self.fit(X, y=None)
        return self._yd2

    def transform(self, X, copy=None):
        ndata,nvars = X.shape
        assert nvars == self.n_vars


        #make them Fortanable
        data = np.asfortranarray(X,dtype=np.float64)
        y = np.asfortranarray(np.empty_like(data),dtype=np.float64)

        st = ppmt_interface.ppmt_ext.ppmt_transform(data,self._iterations,self._snorm,self._zd,self._sph,self._projections,self._Us,y,self.trace)

        return y

    def inverse_transform(self, X, copy=None):
        ndata,nvars = X.shape
        assert nvars == self.n_vars

        #make them Fortanable
        X = np.asfortranarray(X,dtype=np.float64)
        ret = np.asfortranarray(np.empty_like(X),dtype=np.float64)
        ppmt_interface.ppmt_ext.ppmt_back_transform(X,ret,self._iterations,self._snorm,self._zd,self._sph_inv,self._Us,self._projections,self._zd2,self._yd2,self.trace)
        return ret
