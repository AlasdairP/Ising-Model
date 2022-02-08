"""
Checkpoint 1 MVP

Ising model on 2D square lattice with periodic BCs.
User can choose dynamics (Glauber or Kawasaki), size, temperature and number of sweeps.
Animation via matplotlib.
Measurements are taken every 10 sweeps, after an initial period of 100 sweeps to reach equilibrium.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import math
import sys

def Glauber(state):
    """
    First a random spin is chosen with location (i,j).
    Then consider the energy change of flipping the spin.
    If deltaE is negative, exp(-deltaE/kT) is always greater than 1,
       so accept flip with probability 1 (always)
    If deltaE is positive, exp(-deltaE/kT) is always less than 1,
       so accept flip with probability exp(-deltaE/kT),
       done by seeing if a random number (0-1) is greater or less than exp(-deltaE/kT).
    """
    i = random.randint(0,(N-1))
    j = random.randint(0,(N-1))
    
    deltaE = 2*state[i,j]*(state[(i+1)%N,j] + state[i,(j+1)%N] + state[(i-1)%N,j] + state[i,(j-1)%N])
    
    r = random.random()
    
    if deltaE <= 0:
        state[i,j] = -state[i,j]
            
    elif r < math.exp(-deltaE/(T)):
        state[i,j] = -state[i,j]
        
    return state
    
    
def Kawasaki(state):
    """
    First, two random spins are chosen, with locations (i1,j1) and (i2,j2)
    If both are +1 or both are -1 then skip because swapping does nothing
    Else, calculate energy if the spins are swapped (equivalent to flipping both).
    Then if the chosen pair are nearest neighbours, add correction of 4J.
    Then accept the swap with probability min(1, exp(-deltaE/kT))
    If deltaE is negative, exp(-deltaE/kT) is always greater than 1,
       so accept flip with probability 1 (always)
    If deltaE is positive, exp(-deltaE/kT) is always less than 1,
       so accept flip with probability exp(-deltaE/kT)
       done by seeing if a random number (0-1) is greater or less than exp(-deltaE/kT).
    """
    i1 = random.randint(0,N-1)
    j1 = random.randint(0,N-1)
    i2 = random.randint(0,N-1)
    j2 = random.randint(0,N-1)
    
    if state[i1,j1] == state[i2,j2]:
        return state
                                                           
    else:
        deltaE1 = 2*state[i1,j1]*(state[(i1+1)%N,j1] + state[i1,(j1+1)%N] + state[(i1-1)%N,j1] + state[i1,(j1-1)%N])
        deltaE2 = 2*state[i2,j2]*(state[(i2+1)%N,j2] + state[i2,(j2+1)%N] + state[(i2-1)%N,j2] + state[i2,(j2-1)%N])
    
    # If nearest neighbours...
    if (i1,j1) == ((i2+1)%N,j1) or (i1,j1) == ((i2-1)%N,j2) or (i1,j1) == (i2,(j2+1)%N) or (i1,j1) == (i2,(j2-1)%N):
        # Add correction of 4(J)
        deltaE = deltaE1 + deltaE2 + 4    
    
    else:
        deltaE = deltaE1 + deltaE2
        
    # "flip towards order"
    if deltaE <= 0:
        
        state[i1,j1] = -state[i1,j1]
        state[i2,j2] = -state[i2,j2]
    
    # "flip towards disorder" - more likely if T is big, then Boltzmann ~ 1    
    elif r < math.exp(-deltaE/(T)):
        state[i1,j1] = -state[i1,j1]
        state[i2,j2] = -state[i2,j2]
        
    return state

def totalE(state):
    
    E = 0
    for i in range(N):
        for j in range(N):
            E_ij = -state[i,j]*(state[(i+1)%N,j] + state[i,(j+1)%N])
            E += E_ij
            
    return E
    
if len(sys.argv)!=5:              
    print("Wrong number of arguments.")
    print("Usage: " + sys.argv[0] + " <Dynamics> <Size> <T> <Numsweeps>")
    # In the end may want to remove numsteps input
    quit()
    
else:
    dynamics = sys.argv[1]
    N = int(sys.argv[2])
    T= float(sys.argv[3])
    nsweeps = int(sys.argv[4])
    sweep_time = N**2
    nsteps = nsweeps * sweep_time 
    
########## main #########

# Initialise state randomly

state = np.zeros((N,N), dtype=int)
for i in range(N):    
    for j in range(N):
        r = random.random()
        if (r<0.5): state[i,j] = 1
        if (r>0.5): state[i,j] = -1
    
fig = plt.figure()
im = plt.imshow(state,animated=True)

if dynamics == "Glauber":
    
    for n in range(nsteps):
        
        state = Glauber(state)
        if n%(5*sweep_time) == 0:
            plt.cla()
            im = plt.imshow(state, animated=True)
            plt.draw()
            plt.pause(0.00001)
            E = totalE(state)
            M = abs(sum(sum(state)))
            print("Energy = ",E)
    
elif dynamics == "Kawasaki":
    
    for n in range(nsteps):
        
        state = Kawasaki(state)
        
        if n%(5*sweep_time) == 0:
            plt.cla()
            im = plt.imshow(state, animated=True)
            plt.draw()
            plt.pause(0.00001)
            E = totalE(state)
            print("Energy = ",E)
        
else:
    print("Dynamics must be Glauber or Kawasaki")
    quit()


