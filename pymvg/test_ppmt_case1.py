import ppmt
import ppmt_interface
import numpy as np
from sklearn.model_selection import train_test_split



np.random.seed(1634120)

ndata,nvars = 10000,2

data = np.random.random(size=(ndata,nvars))
dummy = data[:,0].copy()

X_train, X_test, _, _ = train_test_split(data, dummy, test_size=0.5, random_state=1634120)


X_train = np.asfortranarray(X_train,dtype=np.float64)
w = np.ones(ndata)
w = np.asfortranarray(w,dtype=np.float64)

lo = 10
trace=0

min_index = 1e-4 #None

model = ppmt.PPMT(n_vars=nvars,min_index=min_index,maxiters=40,seed=1634129,trace=0)

print('ppmt.fit...')
yd2 = model.fit_transform(X_train)
pi = ppmt.test_ppmt(yd2)
print(np.min(pi),np.mean(pi),np.max(pi))
print('ppmt.fit...DONE')

print('ppmt.back_transform...')
X_train_back = model.inverse_transform(yd2)
print('ppmt.back_transform...DONE')

print('ppmt.transform Xtest...')
Y_test = model.transform(X_test)

pi = ppmt.test_ppmt(Y_test)
print(np.min(pi),np.mean(pi),np.max(pi))

X_test_back = model.inverse_transform(Y_test)

import matplotlib.pyplot as plt

#plot X_train
plt.subplot(3, 2, 1)
plt.scatter(X_train[:,0],X_train[:,1])
plt.subplot(3, 2, 2)
plt.scatter(yd2[:,0],yd2[:,1])
#plot X_test
plt.subplot(3, 2, 3)
plt.scatter(X_test[:,0],X_test[:,1])
plt.subplot(3, 2, 4)
plt.scatter(Y_test[:,0],Y_test[:,1])
#plot backs
plt.subplot(3, 2, 5)
plt.scatter(X_train_back[:,0],X_train_back[:,1])
plt.subplot(3, 2, 6)
plt.scatter(X_test_back[:,0],X_test_back[:,1])
plt.show()
