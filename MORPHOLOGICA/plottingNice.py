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

slice = slice2

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

TbodyMean *= TaScale
TbodyStd *= TaScale

flow = Q[:,12:24,slice]
#flowS = np.mean(np.mean(np.abs(np.diff(1.-flow,2)),2),1) # Average over time

#flowS = np.mean(np.mean(np.abs(np.diff(flow,2)),2),1) # Average over time

flowS = np.mean(np.abs(np.diff(np.mean(flow,1),2)),1) # Average over time
flowStd = np.std(np.abs(np.diff(np.mean(flow,1),2)),1)/np.sqrt(12) # Average over time
print np.shape(flow)
print np.shape(np.mean(flow,1))
print np.shape(np.abs(np.diff(np.mean(flow,1),2)))

#print np.shape(np.mean(np.diff(1.-flow,2),2))


#flowS = np.mean(np.std(1.-flow,2),1) # Average over time
#print np.shape(flow)
#print np.shape(np.diff(1.-flow,2))
#print np.shape(np.mean(np.diff(1.-flow,2),2))


fs = 20
ms = 8
F = pl.figure(1,figsize=(10,10))
#F.subplots_adjust(wspace=.3,hspace=0.1,left=0.1,right=0.95,top=0.95,bottom=0.05)
F.subplots_adjust(wspace=.35,hspace=0.35)

f1 = F.add_subplot(221)
f1.plot(Ta,ctactMean,'ko',color=(0,0,0))#,markersize=ms)
for i,j in enumerate(ctactStd):
    f1.plot([Ta[i],Ta[i]],[ctactMean[i]-j,ctactMean[i]+j],color=(0,0,0))
#f1.axis([0,1.0,f1.get_ylim()[0],f1.get_ylim()[1]])
f1.axis([0,TaScale,0,0.6])
f1.set_aspect(np.diff(f1.get_xlim())/np.diff(f1.get_ylim()))
f1.set_axisbelow(True)
f1.set_xlabel(r'$T_A (^{\circ}{\rm C})$',fontsize=fs)
#f1.set_ylabel(r'$1-A$',fontsize=fs)
f1.set_ylabel(r'huddling',fontsize=fs)


f2 = F.add_subplot(222)
f2.plot(Ta,microMean,'ko',color=(0,0,0))#,markersize=ms)
for i,j in enumerate(microStd):
    f2.plot([Ta[i],Ta[i]],[microMean[i]-j,microMean[i]+j],color=(0,0,0))
f2.axis([0,TaScale,1,9])#f2.get_ylim()[0],f2.get_ylim()[1]])
f2.set_aspect(np.diff(f2.get_xlim())/np.diff(f2.get_ylim()))
f2.set_axisbelow(True)
f2.set_xlabel(r'$T_A (^{\circ}{\rm C})$',fontsize=fs)
f2.set_ylabel('groups',fontsize=fs)

f3 = F.add_subplot(223)
f3.plot(Ta,flowS,'ko',color=(0,0,0))#,markersize=ms)
for i,j in enumerate(flowStd):
    f3.plot([Ta[i],Ta[i]],[flowS[i]-j,flowS[i]+j],color=(0,0,0))
#f3.axis([0,1.0,f3.get_ylim()[0],f3.get_ylim()[1]])
f3.axis([0,TaScale,0.005,0.035])
f3.set_aspect(np.diff(f3.get_xlim())/np.diff(f3.get_ylim()))
f3.set_axisbelow(True)
f3.set_xlabel(r'$T_A (^{\circ}{\rm C})$',fontsize=fs)
f3.set_ylabel('pup flow',fontsize=fs)

f4 = F.add_subplot(224)
f4.plot(Ta,TbodyMean,'ko',color=(0,0,0))#,markersize=ms)
for i,j in enumerate(TbodyStd):
    f4.plot([Ta[i],Ta[i]],[TbodyMean[i]-j,TbodyMean[i]+j],color=(0,0,0))
f4.plot([0,TaScale],[0.6*TaScale,0.6*TaScale],'--',color=(0,0,0))
#f.axis([min(Ta),max(Ta),min(Ta),max(Ta)])
f4.axis([0,TaScale,0,60])
f4.set_aspect(np.diff(f4.get_xlim())/np.diff(f4.get_ylim()))
f4.set_axisbelow(True)
f4.set_xlabel(r'$T_A (^{\circ}{\rm C})$',fontsize=fs)
f4.set_ylabel(r'$T_B (^{\circ}{\rm C})$',fontsize=fs)
#f.set_title('time '+str(To)+'--'+str(Te))

# ANOTATION
f1.annotate('A',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-55, 25), va='top',xycoords='axes fraction', textcoords='offset points')

f2.annotate('B',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-55, 25), va='top',xycoords='axes fraction', textcoords='offset points')

f3.annotate('C',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-55, 25), va='top',xycoords='axes fraction', textcoords='offset points')

f4.annotate('D',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-55, 25), va='top',xycoords='axes fraction', textcoords='offset points')

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
#F.tight_layout()
#F.savefig(Fn+'.pdf',dpi=300)
F.savefig('Post1.pdf',dpi=600)
pl.show()

