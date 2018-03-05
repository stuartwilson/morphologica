import numpy as np
import pylab as pl
import sys

fs = 20

Fn = 'data.npz'

D = np.load(Fn)

Q=np.array(D['Q'],dtype=float)
Ta=np.array(D['Ta'],dtype=float)
TaScale = 50
Ta *= TaScale



correl = Q[:,38,:]


fs = 20
ms = 8
F = pl.figure(1,figsize=(12,6))
#F.subplots_adjust(wspace=.3,hspace=0.1,left=0.1,right=0.95,top=0.95,bottom=0.05)
F.subplots_adjust(wspace=.35,hspace=0.35)


n = correl.shape[0]
N = 1.*n
f1 = F.add_subplot(111)
for i in range(n):
    f1.plot(correl[i,:].T,color=(1.*i/N,0,1.-(1.*i/N)))
f1.axis([0,27000,-1.5,0.5])
f1.set_xlabel('time',fontsize=fs)
f1.set_ylabel(r'correlation$(w_L-w_R,w_M)$',fontsize=fs)
F.savefig('coolTing.pdf',dpi=300)
#pl.show()

