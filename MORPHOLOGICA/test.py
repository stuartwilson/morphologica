import numpy as np
import pylab as pl


A = np.zeros(30000)

tswitch = 20000
alpha = 1.

for i in range(30000):
    alpha -= 1./(1.*tswitch)
    alpha = np.fmax(alpha,0.)
    A[i] = alpha

F = pl.figure()
f = F.add_subplot(111)
f.plot(A)
pl.show()
