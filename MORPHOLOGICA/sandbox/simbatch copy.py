import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000

s=11*3

worlds = [P.pTemp('processes/sim','world00','logs/log00',s+0,basePort+0),
          P.pTemp('processes/sim','world01','logs/log01',s+1,basePort+1),
          P.pTemp('processes/sim','world02','logs/log02',s+2,basePort+2),
          P.pTemp('processes/sim','world03','logs/log03',s+3,basePort+3),
          P.pTemp('processes/sim','world04','logs/log04',s+4,basePort+4),
          P.pTemp('processes/sim','world05','logs/log05',s+5,basePort+5),
          P.pTemp('processes/sim','world06','logs/log06',s+6,basePort+6),
          P.pTemp('processes/sim','world07','logs/log07',s+7,basePort+7),
          P.pTemp('processes/sim','world08','logs/log08',s+8,basePort+8),
          P.pTemp('processes/sim','world09','logs/log09',s+9,basePort+9),
          P.pTemp('processes/sim','world10','logs/log10',s+10,basePort+10),
          P.pTemp('processes/sim','world11','logs/log11',s+11,basePort+11),
          P.pTemp('processes/sim','world12','logs/log12',s+12,basePort+12),
          P.pTemp('processes/sim','world13','logs/log13',s+13,basePort+13),
          P.pTemp('processes/sim','world14','logs/log14',s+14,basePort+14),
          P.pTemp('processes/sim','world15','logs/log15',s+15,basePort+15),
          P.pTemp('processes/sim','world16','logs/log16',s+16,basePort+16),
          P.pTemp('processes/sim','world17','logs/log17',s+17,basePort+17),
          P.pTemp('processes/sim','world18','logs/log18',s+18,basePort+18),
          P.pTemp('processes/sim','world19','logs/log19',s+19,basePort+19),
          P.pTemp('processes/sim','world20','logs/log20',s+20,basePort+20),
          P.pTemp('processes/sim','world21','logs/log21',s+21,basePort+21),
          P.pTemp('processes/sim','world22','logs/log22',s+22,basePort+22),
          P.pTemp('processes/sim','world23','logs/log23',s+23,basePort+23)]


Npup = 12
Nsim = len(worlds)

Ta = np.linspace(0.0,1.0,Nsim)                          # Define Ambient temperatures

k1 = np.tile(np.linspace(0.9,0.9,Nsim),[Npup,1]).T      # Define k1 (decay to ambient)
k2 = np.tile(np.linspace(1.0,1.0,Nsim),[Npup,1]).T      # Define k2 (decay to contact)

G = np.tile(np.linspace(0.055,0.055,Nsim),[Npup,1]).T   # Define G (thermogenesis)
G[:,np.array([0,2,4,6,8,10])] = 0.045
G[:,np.array([1,3,5,7,9,11])] = 0.065


# SET INITIAL CONDITIONS
for w, world in enumerate(worlds):

    world.stream('3,0,0,1,0,0,'+str(pi*0.5)+',16,16') # Delete first animat
    world.stream('1,0,')
    
    # Set environment params
    world.stream('5,0,'+str(Ta[w])+',')
    world.stream('5,1,2.0,')


    # Set pup params
    for i in range(Npup):
        world.stream('4,0,'+str(i)+','+str(k1[w,i])+',')
        world.stream('4,1,'+str(i)+','+str(k2[w,i])+',')
        world.stream('4,2,'+str(i)+','+str(G[w,i])+',')

    # Set initial SOM parameters
    world.stream('6,0,0.025,')                             # set learning rate INPUT
    world.stream('6,1,0.025,')                             # set learning rate OUTPUT
    world.stream('6,2,0.5,')                               # set smoothing rate of noise
    world.stream('6,3,0.5,')                               # set smoothing of threshold
    world.stream('6,4,0.5,')                               # set neighbourhood sigma
    world.stream('6,5,0.75,')                              # set cutoff
    world.stream('6,6,1000,')                              # set timescale for decay rates ???


# Timing
trialLength = 1500
numTrials = 13
numTests = 5
tSwitch = trialLength * numTrials
tTotal = tSwitch + trialLength * numTests + 1

# Storage
K = np.zeros([Nsim,38])
Q = K
alpha = 1.

''' 
# THIS FOR LOADING
for w, world in enumerate(worlds):
    worlds[w].stream('12,dataweights'+str(w)+'.bin')
'''

# MAIN SIMULATION LOOP
for t in range(tTotal):

    alpha -= 1./(1.*tSwitch)
    alpha = np.fmax(alpha,0.)
    for w, world in enumerate(worlds):
        world.stream('1,'+str(alpha)+',')
  
    for w, world in enumerate(worlds):
        K[w,:] = world.out()
    Q = np.dstack([Q,K])

    if(t%100==0):
        for w, world in enumerate(worlds):
            world.stream('8,')
            world.stream('9,')
        print 't='+str(t)


    if(t%trialLength==0):
        np.savez('data',Q=Q,Ta=Ta)
        print 'saved'


    if(t%trialLength==0):
        for w, world in enumerate(worlds):
            world.stream('7,')

    if(t==tSwitch):
        print "Networks operating without supervision; alpha = "+str(alpha)
        for w, world in enumerate(worlds):
            worlds[w].stream('11,data/weights'+str(w)+'.bin')
        print "Saved weights to file"

for world in worlds:
    world.quit()



