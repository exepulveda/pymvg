import ppmt_interface
import numpy as np
import scipy
import scipy.stats
import sklearn
import sklearn.decomposition
import os.path

from scipy.stats import norm

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

def max_pi(data,ndirs=100,lo=20):
    nsamples,dmin = data.shape
    #random direction for projection
    directions = generate_directions(dmin,n=ndirs)
    #directions[:,0] =
    #directions = np.diag(np.ones(nvars))
    #directions = np.asfortranarray(directions,dtype=np.float64)

    w = np.ones(nsamples)
    w = np.asfortranarray(w,dtype=np.float64)

    data = np.asfortranarray(data,dtype=np.float64)

    #compute pp index
    indices = np.empty(ndirs)

    for j in range(ndirs):
        #dir = np.zeros(nvars)
        direction = directions[j].copy()
        direction = np.asfortranarray(direction,dtype=np.float64)

        pi = ppmt_interface.ppmt_ext.pp_index(direction,data,w,lo)
        indices[j] = pi

    return np.max(indices)


def bootstrap_ppmt(nvars,ndirs=10,nboot=100,nsamples=100000,percentile=50,lo=20):
    #the idea is to generate samples from mutigaussian and calculate the projection index onto random directions
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
    def __init__(self,n_vars=None,min_index=1e-4,maxiters=500,seed=1634129,trace=0):
        self.n_vars = n_vars
        self.min_index = min_index
        self.maxiters = maxiters
        self.seed = seed
        self.trace = trace

        assert 0 <= maxiters <= 1000

    def load_from_file(self,filename):
        if not os.path.exists(filename):
            raise Exception("not found")

        #get the dimensions stored in the binary file
        nvars, ndata, iters = ppmt_interface.ppmt_ext.ppmt_extract_param_dim(filename)
        print(ndata, nvars, iters)

        if ndata == 0 and ndata == 0 and iters == 0:
            raise Exception("poblem to read binary file")


        #make them Fortranable
        self._yd2 = np.asfortranarray(np.empty((ndata,nvars)),dtype=np.float64)
        self._snorm = np.asfortranarray(np.empty(ndata),dtype=np.float64)
        self._zd = np.asfortranarray(np.empty((ndata,nvars)),dtype=np.float64)
        self._sph_inv = np.asfortranarray(np.empty((nvars,nvars)),dtype=np.float64)
        self._Us = np.asfortranarray(np.empty((nvars,nvars,iters)),dtype=np.float64)
        self._projections = np.asfortranarray(np.empty((ndata,iters)),dtype=np.float64)
        self._zd2 = np.asfortranarray(np.empty((ndata,nvars)),dtype=np.float64)

        #print('calling ppmt_interface.ppmt_ext.ppmt_extract_param...')
        self._iterations = ppmt_interface.ppmt_ext.ppmt_extract_param(filename,nvars,ndata,iters,self._snorm,self._zd,self._sph_inv,self._Us,self._projections,self._zd2,self._yd2,self.trace)
        self.n_vars = nvars
        self.ndata = ndata
        #print('calling ppmt_interface.ppmt_ext.ppmt_extract_param...DONE')
        #print(self._iterations)
        #print(self._sph_inv)
        #reshape objects to hold just actual iterations
        self._Us = np.asfortranarray(self._Us[:,:,:self._iterations],dtype=np.float64)
        self_projections = np.asfortranarray(self._projections[:,:self._iterations],dtype=np.float64)
        self._sph = np.asfortranarray(np.linalg.inv(self._sph_inv),dtype=np.float64)


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
        self._projections = self._projections[:,:self._iterations]

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

        st = ppmt_interface.ppmt_ext.ppmt_transform(
                data,
                self._iterations,
                self._snorm,
                self._zd,
                self._sph,
                self._projections,
                self._Us,
                y,
                self.trace)

        return y

    def inverse_transform(self, X, copy=None):
        ndata,nvars = X.shape
        assert nvars == self.n_vars

        #make them Fortanable
        X = np.asfortranarray(X,dtype=np.float64)
        ret = np.asfortranarray(np.empty_like(X),dtype=np.float64)
        ppmt_interface.ppmt_ext.ppmt_back_transform(
                X,
                ret,
                self._iterations,
                self._snorm,
                self._zd,
                self._sph_inv,
                self._Us,
                self._projections,
                self._zd2,
                self._yd2,
                self.trace)
        return ret


def departure_from_normality(X):
    V = np.cov(X,rowvar=False)
    var = np.var(X,axis=0)
    D = np.eye(X.shape[1])*var
    D_inv = np.linalg.inv(D)
    R = D_inv@V@D_inv

    E=np.diag(R)
    EE=R-np.diag(E)

    return np.sum((E-1.0)**2) + np.sum(EE**2)

#==================== Normal Score =============================================
def normal_scores(o_data,o_weights=None,indices=None,na_value=np.nan):
    epsilon = 1.0e-7

    m = len(o_data)

    if o_weights is None:
        o_weights = np.ones(m)

    if indices is not None:
        data = o_data[indices]
        weights = o_weights[indices]
    else:
        data = o_data
        weights = o_weights

    n = len(data)

    data = data + np.random.random(n) * epsilon

    #sort data and weights
    table = np.empty((n,2))
    table[:,0] = data
    table[:,1] = weights

    table = table[table[:,0].argsort()]

    sorted_data = table[:,0]
    sorted_weights = table[:,1]

    #normalize weights
    wsum = np.sum(sorted_weights)
    nweights = sorted_weights/wsum #normalized weights

    #Cummulative distribution
    cumsum = np.cumsum(nweights) - 0.5/wsum #centroids
    weights = norm.ppf(cumsum)

    #transformation table
    table[:,0] = sorted_data
    table[:,1] = weights

    #Interpolate
    transformed_data = np.interp(data, sorted_data, weights)


    if indices is not None:
        tmp = np.empty(m)
        tmp[:] = na_value
        tmp[indices] = transformed_data
        transformed_data = tmp

    return transformed_data, table


def interpolator(xval,xlow,xhigh,ylow,yhigh,pow=1.0):
    if (xhigh - xlow) < EPSILON:
        ret = (yhigh + ylow) / 2.0
    else:
        ret = ylow + (yhigh-ylow) * (((xval-xlow)/(xhigh-xlow))**pow)

    return ret

def back_normal_scores(o_data,table,indices=None,na_value=np.nan,min_value=0.0,max_value=6.0,ltail=1.0,utail=2.0):
    '''exponential extrapolation'''
    epsilon = 1.0e-7

    m = len(o_data)

    if indices is not None:
        data = o_data[indices]
    else:
        data = o_data

    n = len(data)

    #apply linear interpolation to all data
    ntable = len(table)
    new_table = np.empty((ntable+2,2))
    #add extremmes
    new_table[0,0] = min_value
    new_table[0,1] = np.min(data)
    new_table[-1,0] = max_value
    new_table[-1,1] = np.max(data)
    new_table[1:-1,:] = table

    backtransformed_data = np.interp(data, new_table[:,1], new_table[:,0])

    #take care of tails

    #find tails
    lower_tails = np.where(data < table[0,1])[0]
    if len(lower_tails) > 0:
        if ltail != 1:
            cpow = 1.0 / ltpar
        else:
            cpow = 1.0

        z1 = table[0,0]
        y1 = table[0,1]
        b1 = (z1-min_value)*np.exp(-ltail*y1)
        for i in lower_tails:
            backtransformed_data[i] = min_value + b1*np.exp(ltail*data[i])
            #print('LB',i,data[i],table[0,1],z1,y1,b1,backtransformed_data[i])


    upper_tails = np.where(data > table[-1,1])[0]
    if len(upper_tails) > 0:
        zn = table[-1,0]
        yn = table[-1,1]

        bn = (zn-max_value)*np.exp(utail*yn)

        for i in upper_tails:
            backtransformed_data[i] = max_value + bn*np.exp(-utail*data[i])
            #print('UB',i,data[i],table[-1,1],zn,yn,bn,backtransformed_data[i])


    return backtransformed_data

def marginal_gaussian(X):
    n,nd = X.shape

    G = np.empty_like(X)
    tables = []
    for i in range(nd):
        transformed_data, table = normal_scores(X[:,i])
        G[:,i] = transformed_data
        tables += [table]

    return G,tables

def apply_marginal_gaussian(X,tables):
    n,nd = X.shape

    G = np.copy(X)
    for i in range(nd):
        table = tables[i]
        backtransformed_data = np.interp(G[:,i], table[:,0], table[:,1])
        G[:,i] = backtransformed_data

    return G

class MarginalMT():
    def __init__(self,max_iter=100,tol=1e-3,rot_max_iter=200,rot_tol=1e-6,seed=1634120,method='pca'):
        self.max_iter = max_iter
        self.tol = tol
        self.rot_max_iter = rot_max_iter
        self.rot_tol = rot_tol
        self.ica = []
        self.ica_algorithm = 'parallel'#'deflation'
        self.ica_fun = 'logcosh'
        self.seed = seed
        self.method = method

    def fit(self,X,trace=0):
        self.fit_transform(X,trace=trace)

    def fit_transform(self,X,trace=0):
        n,nd = X.shape
        self.ica = []
        self.tables = []

        Z = X.copy()
        #Z, tables_initial = marginal_gaussian(X)

        #whitening
        #Z = whiten(X)

        if trace>0:
            print(np.mean(Z,axis=0))
            print(np.std(Z,axis=0))

            pi = max_pi(Z)
            print('initial projection index:',pi)

        np.random.seed(self.seed)
        seeds = np.int32(np.random.random(self.max_iter)* 100000)

        for i in range(self.max_iter):
            if self.method == 'ica':
                ica = sklearn.decomposition.FastICA(n_components=None,max_iter=self.rot_max_iter,tol=self.rot_tol,algorithm=self.ica_algorithm,fun=self.ica_fun,random_state=seeds[i],whiten=True)
            else:
                ica = sklearn.decomposition.PCA(n_components=None,random_state=seeds[i],whiten=True)
            #ica.fit(Z)
            G = ica.fit_transform(Z)

            if trace>0:
                pi = max_pi(G)
                print('projection index at iteration %d after ICA:'%(i+1),pi)

            newG, tables = marginal_gaussian(G)
            self.ica += [ica]
            self.tables += [tables]

            Z = newG

            pi = max_pi(Z)
            if trace>0:
                print('projection index at iteration %d after marginal gaussian transform:'%(i+1),pi)

            if pi < self.tol:
                break

        self._iterations = i

        if trace>0:
            print(np.mean(Z,axis=0))
            print(np.std(Z,axis=0))

        self.Y = Z

        return self.Y

    def transform(self,X,trace=0):
        n,nd = X.shape

        Z = X.copy()

        if trace>0:
            print(np.mean(Z,axis=0))
            print(np.std(Z,axis=0))

        if trace>0: 
            pi = max_pi(Z)
            print('initial pi:',pi)

        for i in range(self._iterations):
            G = self.ica[i].transform(Z)

            Z = apply_marginal_gaussian(G,self.tables[i])

            if trace>0: 
                pi = max_pi(Z)
                print('pi at iteration %d after ICA:'%(i+1),pi)

        if trace>0:
            print(np.mean(Z,axis=0))
            print(np.std(Z,axis=0))

        return Z

    def inverse_transform(self,Y,trace=0):
        n,nd = Y.shape

        Z = Y.copy()
        for i in range(self._iterations-1,-1,-1):
            #back transform normal scores
            G = np.empty_like(Z)
            for k in range(nd):
                G[:,k] = np.interp(Z[:,k],self.tables[i][k][:,1],self.tables[i][k][:,0])

            #back transform ICA
            Z = self.ica[i].inverse_transform(G)

        return Z
