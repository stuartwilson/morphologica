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
slice1 = np.hstack([s1,s2,s3,s4,s5])

# Use these for late
s1 = 20000 + np.arange(1000)
s2 = 21500 + np.arange(1000)
s3 = 23000 + np.arange(1000)
s4 = 24500 + np.arange(1000)
s5 = 26000 + np.arange(1000)
slice2 = np.hstack([s1,s2,s3,s4,s5])


slice = slice1



macro = Q[:,36,slice]
micro = Q[:,37,slice]

macroMean = np.mean(macro,1) # Average over time
microMean = np.mean(micro,1) # Average over time
macroStd = np.std(macro,1) # Standard Deviation over time
microStd = np.std(micro,1) # Standard Deviation over time

#A

Male = np.array([0,2,4,6,8,10])
Female = np.array([1,3,5,7,9,11])

Tbody = Q[:,0:12,slice]
ctact = Q[:,12:24,slice]
learn = Q[:,24:36,slice]

Tbody_M = Tbody[:,Male,:] #Q[:,0+Male,slice]
ctact_M = ctact[:,Male,:] #Q[:,12+Male,slice]
learn_M = learn[:,Male,:] #Q[:,24+Male,slice]
TbodyM_M = np.mean(Tbody_M,2) # Average over time
ctactM_M = 1.-np.mean(ctact_M,2) # Average over time
learnM_M = np.mean(learn_M,2) # Average over time
TbodyMean_M = np.mean(TbodyM_M,1) # Average over pups
ctactMean_M = np.mean(ctactM_M,1) # Average over pups
learnMean_M = np.mean(learnM_M,1) # Average over pups
TbodyStd_M = np.std(TbodyM_M,1) # Standard deviation over pups
ctactStd_M = np.std(ctactM_M,1) # Standard deviation over pups
learnStd_M = np.std(learnM_M,1) # Standard deviation over pups


Tbody_F = Tbody[:,Female,:] #Q[:,0+Male,slice]
ctact_F = ctact[:,Female,:] #Q[:,12+Male,slice]
learn_F = learn[:,Female,:] #Q[:,24+Male,slice]
TbodyM_F = np.mean(Tbody_F,2) # Average over time
ctactM_F = 1.-np.mean(ctact_F,2) # Average over time
learnM_F = np.mean(learn_F,2) # Average over time
TbodyMean_F = np.mean(TbodyM_F,1) # Average over pups
ctactMean_F = np.mean(ctactM_F,1) # Average over pups
learnMean_F = np.mean(learnM_F,1) # Average over pups
TbodyStd_F = np.std(TbodyM_F,1) # Standard deviation over pups
ctactStd_F = np.std(ctactM_F,1) # Standard deviation over pups
learnStd_F = np.std(learnM_F,1) # Standard deviation over pups


F = pl.figure(1, figsize=(12,4))
f = F.add_subplot(131)
pl.plot(Ta,TbodyMean_M,'.-',color=(0,0,1))
for i,j in enumerate(TbodyStd_M):
    pl.plot([Ta[i],Ta[i]],[TbodyMean_M[i]-j,TbodyMean_M[i]+j],color=(0,0,1))

pl.plot(Ta,TbodyMean_F,'.-',color=(1,0,0))
for i,j in enumerate(TbodyStd_F):
    pl.plot([Ta[i],Ta[i]],[TbodyMean_F[i]-j,TbodyMean_F[i]+j],color=(1,0,0))

pl.plot([0,1],[0,1],'--',color=(0,0,0))

#f.axis([min(Ta),max(Ta),min(Ta),max(Ta)])
f.axis([0,1.0,0,1.0])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel(r'$T_b$',fontsize=fs)
#f.set_title('time '+str(To)+'--'+str(Te))


# "Males" (lower G) in red
f = F.add_subplot(132)
pl.plot(Ta,ctactMean_M,'.-',color=(0,0,1))
for i,j in enumerate(ctactStd_M):
    pl.plot([Ta[i],Ta[i]],[ctactMean_M[i]-j,ctactMean_M[i]+j],color=(0,0,1))

pl.plot(Ta,ctactMean_F,'.-',color=(1,0,0))
for i,j in enumerate(ctactStd_F):
    pl.plot([Ta[i],Ta[i]],[ctactMean_F[i]-j,ctactMean_F[i]+j],color=(1,0,0))

f.axis([0,1.0,f.get_ylim()[0],f.get_ylim()[1]])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel(r'$1-A$',fontsize=fs)


f = F.add_subplot(133)
pl.plot(Ta,microMean,'.-',color=(0,0,0))
for i,j in enumerate(microStd):
    pl.plot([Ta[i],Ta[i]],[microMean[i]-j,microMean[i]+j],color=(0,0,0))
f.axis([0,1.0,f.get_ylim()[0],f.get_ylim()[1]])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_axisbelow(True)
f.set_xlabel(r'$T_a$',fontsize=fs)
f.set_ylabel('# subgroups',fontsize=fs)

F.tight_layout()
F.savefig(Fn+'.pdf',dpi=300)
pl.show()
