import numpy as np
import pylab as pl
import sys

fs = 20

Fn = 'data.npz'

D = np.load(Fn)

Q=np.array(D['Q'],dtype=float)
Ta=np.array(D['Ta'],dtype=float)

# Use these for early
s1 = np.arange(500,1500)
s2 = np.arange(2000,3000)
s3 = np.arange(3500,4500)
s4 = np.arange(5000,6000)
s5 = np.arange(6500,7500)
s6 = np.arange(8000,9000)
s7 = np.arange(9500,10500)
s8 = np.arange(11000,12000)
s9 = np.arange(12500,13500)
s10 = np.arange(14000,15000)
#slice1 = np.hstack([s1,s2,s3,s4,s5])
slice1 = np.hstack([s1,s2,s3,s4,s5,s6,s7,s8,s9,s10])

# Use these for late
s1 = 20000 + np.arange(1000)
s2 = 21500 + np.arange(1000)
s3 = 23000 + np.arange(1000)
s4 = 24500 + np.arange(1000)
s5 = 26000 + np.arange(1000)
slice2 = np.hstack([s1,s2,s3,s4,s5])


slice = slice1

Tbody = Q[:,0:12,slice]
ctact = Q[:,12:24,slice]
learn = Q[:,24:36,slice]
macro = Q[:,36,slice]
micro = Q[:,37,slice]


TbodyM = np.mean(Tbody,2) # Average over time
ctactM = 1.-np.mean(ctact,2) # Average over time
learnM = np.mean(learn,2) # Average over time

macroMean = np.mean(macro,1) # Average over time
microMean = np.mean(micro,1) # Average over time

macroStd = np.std(macro,1) # Standard Deviation over time
microStd = np.std(micro,1) # Standard Deviation over time

TbodyMean = np.mean(TbodyM,1) # Average over pups
ctactMean = np.mean(ctactM,1) # Average over pups
learnMean = np.mean(learnM,1) # Average over pups

TbodyStd = np.std(TbodyM,1) # Standard deviation over pups
ctactStd = np.std(ctactM,1) # Standard deviation over pups
learnStd = np.std(learnM,1) # Standard deviation over pups


F = pl.figure(1, figsize=(12,4))
f = F.add_subplot(131)
pl.plot(Ta,TbodyMean,'.-')
for i,j in enumerate(TbodyStd):
    pl.plot([Ta[i],Ta[i]],[TbodyMean[i]-j,TbodyMean[i]+j],color=(0,0,0))
pl.plot([0,1],[0,1],color=(1,0,0))
#f.axis([min(Ta),max(Ta),min(Ta),max(Ta)])
f.axis([0,1.0,0,1.0])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel(r'$T_b$',fontsize=fs)
#f.set_title('time '+str(To)+'--'+str(Te))


f = F.add_subplot(132)
pl.plot(Ta,ctactMean,'.-')
for i,j in enumerate(ctactStd):
    pl.plot([Ta[i],Ta[i]],[ctactMean[i]-j,ctactMean[i]+j],color=(0,0,0))
f.axis([0,1.0,f.get_ylim()[0],f.get_ylim()[1]])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel(r'$1-A$',fontsize=fs)


f = F.add_subplot(133)
pl.plot(Ta,microMean,'.-')
for i,j in enumerate(microStd):
    pl.plot([Ta[i],Ta[i]],[microMean[i]-j,microMean[i]+j],color=(0,0,0))
f.axis([0,1.0,f.get_ylim()[0],f.get_ylim()[1]])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel('# subgroups',fontsize=fs)

'''
f = F.add_subplot(224)
pl.plot(Ta,learnMean,'.-')
for i,j in enumerate(learnStd):
    pl.plot([Ta[i],Ta[i]],[learnMean[i]-j,learnMean[i]+j],color=(0,0,0))
f.axis([0,1.0,f.get_ylim()[0],f.get_ylim()[1]])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$')
f.set_ylabel('learning iterations')
'''
F.tight_layout()
F.savefig(Fn+'.pdf',dpi=300)
pl.show()

