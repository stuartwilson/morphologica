import pylab as pl
import numpy as np

# model

tmin = -1.5
tmax = 60
t = np.linspace(tmin,tmax,100)

k = 8.31
c1 = 19.
c2 = 3.
c3 = 1./40.
Tp = 35.7

P = np.exp(-t/k)
M = k*(1.-P)
S = -k*P*np.log(P)
Smax = k*np.exp(-1)
G = k*(1.+S)
T1 = c2*P
N = c1*np.exp(-k*T1)
T2 = c3*N*G
W = M+N

#Mscale = 3.5


lagerX = np.array([-1.32673611,-0.33840278,1.64395833,3.67326389,5.67034722,9.60895833,11.60895833,13.603125,15.61770833,17.64993056,33.60590278])
lagerY = np.array([1.02707883,1.17001608,0.91548949,0.62592582,0.50506952,0.46018333,0.45676491,0.49481427,0.49566718,0.45590523,0.3602843])

eedyX = np.array([1,3,5,7,14,21,35,49.])
eedyY = np.array([3.01899547,4.58730481,6.06548142,6.38995921,5.89422925,3.99242885,2.3970797,2.05457536,2.05457536])

eedyY *= 0.001      # convert cm^3 O2 / g/h to liters 02 /g/h
eedyY *= 5000.      # convert liters 02 /g/h to cal /g/h


eedyXa = np.array([1,3,5,7,10,14,21,28,35,49.])
eedyYa = np.array([0.79733333,1.96666667,2.472,4.07111111,4.82888889,6.73688889,10.97377778,15.91244444,21.74888889,25.62088889,32.09288889])
eedyYb = np.array([38.10860234,37.85330412,35.53773585,36.73535724,35.2981758,35.28225879,34.91500492,32.7752392,31.41800054,31.91133864,32.40798533])


# PLOTTING
fs = 15
ms = 8
F = pl.figure(1,figsize=(12,12))
#F.subplots_adjust(wspace=.3,hspace=0.1,left=0.1,right=0.95,top=0.95,bottom=0.05)

F.subplots_adjust(wspace=.3,hspace=0.1,left=0.1,right=0.95,top=0.95,bottom=0.05)


Ymax = lagerY[0]#np.max(Y)
Ymin = np.min(lagerY)
Yfit = P*(Ymax-Ymin)+Ymin

f1 = F.add_subplot(221)
f1.plot(t,Yfit,'-',color=(0,0,0))
f1.plot(lagerX,lagerY,'ko',color=(0,0,0),markersize=ms)
f1.axis([tmin,tmax,0.3,1.2])

f1.legend([r'$P=e^{-t/k}$'],frameon=False,loc='upper right',fontsize=fs)

f1.set_xlabel('$t$ (days)',fontsize=fs,labelpad=20)
f1.set_ylabel(r'brown fat / 100 g$\,$',fontsize=fs)
f1.set_aspect(np.diff(f1.get_xlim())/np.diff(f1.get_ylim()),adjustable='box-forced')
f1b = f1.twinx()
f1b.spines['top'].set_visible(False)
f1b.spines['bottom'].set_visible(False)
f1b.spines['left'].set_visible(False)
f1b.spines['right'].set_visible(False)
f1b.set_aspect(np.diff(f1b.get_xlim())/np.diff(f1b.get_ylim()),adjustable='box-forced')
ylim = f1.get_ylim()
y0 = ylim[0]*1.
y1 = ylim[1]*1.
f1b.set_yticks([(Ymax-y0)/(y1-y0),(Ymin-y0)/(y1-y0)])
f1b.set_yticklabels([r'$P=1$',r'$P=0$'],fontsize=fs)

f2 = F.add_subplot(222)
f2.plot(t,G,'-',color=(0,0,0))
f2.plot(0,0,'.',color=(1,1,1))
f2.plot(eedyX,eedyY[:-1],'ko',color=(0,0,0),markersize=ms)

'''
offset = pm[1]#-2.5
Sscale = 60
f2.plot(t,offset+(Smax-S)*Sscale,color=(0,0,0))
for i in np.arange(len(pm)):
    f2.plot(tB[i],pm[i],'ko',color=(0,0,0),markersize=ms)
    f2.plot([tB[i],tB[i]],[pm[i]-ps[i],pm[i]+ps[i]],'-',color=(0,0,0))
f2.axis([2,tmax,0,30])
'''

f2.legend([r'$G=k(1+S),$',r'$S=-kP\,\ln P$'],frameon=False,loc='upper right',fontsize=fs)
f2.set_xlabel('$t$ (days)',fontsize=fs,labelpad=20)
f2.set_ylabel(r'metabolism (cal/g/h)$\,$',fontsize=fs)
f2.axis([tmin,tmax,0,35])
f2.set_aspect(np.diff(f2.get_xlim())/np.diff(f2.get_ylim()),adjustable='box-forced')
f2b = f2.twinx()
f2b.spines['top'].set_visible(False)
f2b.spines['bottom'].set_visible(False)
f2b.spines['left'].set_visible(False)
f2b.spines['right'].set_visible(False)
f2b.set_aspect(np.diff(f2b.get_xlim())/np.diff(f2b.get_ylim()),adjustable='box-forced')
ylim = f2.get_ylim()
y0 = ylim[0]*1.
y1 = ylim[1]*1.
f2b.set_yticks([(k-y0)/(y1-y0)])
f2b.set_yticklabels([r'$k$'],fontsize=fs)


############

f3 = F.add_subplot(223)
f3.plot(t,M,'--',color=(0,0,0))
f3.plot(t,W,'-',color=(0,0,0))
f3.plot(eedyXa,eedyYa[:-1],'ko',color=(0,0,0),markersize=ms)
f3.axis([tmin,tmax,0,30])

f3.legend([r'$M=k(1-P)$',r'$W=M+ce^{-kT_1}$'],frameon=False,loc='upper left',fontsize=fs)

f3.set_xlabel('$t$ (days)',fontsize=fs,labelpad=20)
f3.set_ylabel(r'weight (g)$\,$',fontsize=fs)
f3.set_aspect(np.diff(f3.get_xlim())/np.diff(f3.get_ylim()),adjustable='box-forced')
f3b = f3.twinx()
f3b.spines['top'].set_visible(False)
f3b.spines['bottom'].set_visible(False)
f3b.spines['left'].set_visible(False)
f3b.spines['right'].set_visible(False)
f3b.set_aspect(np.diff(f3b.get_xlim())/np.diff(f3b.get_ylim()),adjustable='box-forced')
ylim = f3.get_ylim()
y0 = ylim[0]*1.
y1 = ylim[1]*1.
f3b.set_yticks([(k-y0)/(y1-y0)])
f3b.set_yticklabels([r'$k$'],fontsize=fs)

f4 = F.add_subplot(224)
f4.plot(t,Tp+T1,'--',color=(0,0,0))
f4.plot(t,Tp-T2,'-',color=(0,0,0))
f4.plot(eedyXa,eedyYb[:-1],'ko',color=(0,0,0),markersize=ms)


f4.legend([r'$T_p+\left[T_1\propto P\,\right]$',r'$T_p-\left[T_2\propto G(W-M)\right]$'],frameon=False,loc='upper right',fontsize=fs)
f4.set_xlabel('$t$ (days)',fontsize=fs,labelpad=20)
f4.set_ylabel(r'preference ($^{\circ}$C)',fontsize=fs)
f4.axis([tmin,tmax,31,39])
f4.set_aspect(np.diff(f4.get_xlim())/np.diff(f4.get_ylim()),adjustable='box-forced')
f4b = f4.twinx()
f4b.spines['top'].set_visible(False)
f4b.spines['bottom'].set_visible(False)
f4b.spines['left'].set_visible(False)
f4b.spines['right'].set_visible(False)
f4b.set_aspect(np.diff(f4b.get_xlim())/np.diff(f4b.get_ylim()),adjustable='box-forced')
ylim = f4.get_ylim()
y0 = ylim[0]*1.
y1 = ylim[1]*1.
f4b.set_yticks([(Tp-y0)/(y1-y0)])
f4b.set_yticklabels([r'$T_p$'],fontsize=fs)



# ANOTATION
f1.annotate('A',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-60, 27), va='top',xycoords='axes fraction', textcoords='offset points')

f2.annotate('B',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-60, 27), va='top',xycoords='axes fraction', textcoords='offset points')

f3.annotate('C',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-60, 27), va='top',xycoords='axes fraction', textcoords='offset points')

f4.annotate('D',fontsize=20, fontweight='bold', xy=(0, 1), xytext=(-60, 27), va='top',xycoords='axes fraction', textcoords='offset points')

pl.savefig('Fig1_2.pdf',dpi=600)
#pl.show()

