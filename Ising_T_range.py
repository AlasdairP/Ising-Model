"""
Checkpoint 1 MVP

Ising model on 2D square lattice with periodic BCs.
User can choose dynamics (Glauber or Kawasaki), size, temperature range and number of sweeps.
Measurements are taken every 10 sweeps, after an initial period of 100 sweeps to reach equilibrium.
J = 1.0 and k = 1.0 so they are not included in the code for simplicity.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import math
import sys
    
def Glauber(state,T):
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


def Kawasaki(state,T):
    """
    First, two random spins are chosen, with locations (i1,j1) and (i2,j2)
    If both are +1 or both are -1 then skip because swapping does nothing
    Else, consider energy if the spins are swapped (equivalent to flipping both)
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
    
    # If chosen pair are nearest neighbours, add correction of 4(J)
    if abs(i1-i2) + abs(j1-j2) == 1 :
        deltaE = deltaE1 + deltaE2 + 4    
    
    else:
        deltaE = deltaE1 + deltaE2
    
    r = random.random()
    
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
    """
    Calculate the total energy of the state    
    """
    E = 0
    for i in range(N):
        for j in range(N):
            E_ij = -state[i,j]*(state[(i+1)%N,j] + state[i,(j+1)%N])
            E += E_ij
            
    return E
    
############################## input #################################
    
if len(sys.argv)!=7:              
    print("Wrong number of arguments.")
    print("Usage: " + sys.argv[0] + " <Dynamics> <Size> <Tstart (included)> <Tend(not included)> <Tinterval> <Numsweeps>")
    quit()
    
else:
    dynamics = sys.argv[1]
    N = int(sys.argv[2])
    Tstart = float(sys.argv[3])
    Tend = float(sys.argv[4])
    Tinterval = float(sys.argv[5])
    nsweeps = int(sys.argv[6])
    sweep_time = N**2
    nsteps = nsweeps * sweep_time
    T_array = np.arange(Tstart,Tend,Tinterval)  
    
################################### main ######################################

if dynamics == "Glauber":
    
    # Initialise all up for first temperature
    state = np.ones((N,N),dtype=int)
    num_meas = int((nsweeps-100)/10)
    
    for T in T_array:
        
        # Needed for measurements
        i=0                                  # Index of the measurement arrays
        E_array = np.zeros(num_meas)
        M_array = np.zeros(num_meas)
        
        for n in range(nsteps):
            
            # GLAUBER UPDATE - ATTEMPT FLIP ON EVERY "STEP"
            state = Glauber(state,T)
            
            # Measurements every 10 sweeps, after 100 sweeps equilib time
            
            equilib_time = 100*sweep_time
            if n >= equilib_time:
                if n%(10*sweep_time) == 0:
                    print("Progress: T = " + str(T) + ": " + str(n/sweep_time) + "/" + str(nsweeps) + " sweeps")
                    M_array[i] = abs(sum(sum(state)))
                    E_array[i] = totalE(state)
                    i += 1
                    
        avg_M = np.mean(M_array)
        avg_M_error = np.std(M_array) / np.sqrt(len(M_array))
        chi = 1/(N**2*T) * (np.mean(M_array**2) - np.mean(M_array)**2)
        avg_E = np.mean(E_array)
        avg_E_error = np.std(E_array) / np.sqrt(len(E_array))
        c = 1/(N**2*T**2) * (np.mean(E_array**2) - np.mean(E_array)**2)
            
        # Bootstrap for errors on heat capacity (c) and susceptiility (chi)
        
        num_resamplings = 100
        c_array = np.zeros(num_resamplings)
        chi_array = np.zeros(num_resamplings)
        for k in range(num_resamplings):
            E_resampled = np.random.choice(E_array,num_meas)
            M_resampled = np.random.choice(E_array,num_meas)
            c_array[k] = 1/(N**2*T**2) * (np.mean(E_resampled**2) - np.mean(E_resampled)**2)
            chi_array[k] = 1/(N**2*T) * (np.mean(M_resampled**2) - np.mean(M_resampled)**2)
                
        c_error = np.sqrt(np.mean(c_array**2) - np.mean(c_array)**2)
        chi_error = np.sqrt(np.mean(chi_array**2) - np.mean(chi_array)**2)
        
        with open("Glauber_data.txt","a+") as f:
            f.write(str(T)+" "+str(avg_E)+" "+str(avg_E_error)+" "+str(c)+" "+str(c_error)
            +" "+str(avg_M)+" "+str(avg_M_error)+" "+str(chi)+" "+str(chi_error)+"\n")
    
###################################### Kawasaki #########################################
    
elif dynamics == "Kawasaki":
    
    # Initialise half up, half down for first temperature
    state = np.ones((N,N),dtype=int)
    state[int(N/2): ,:] = -1
    
    for T in T_array:
        
        # Needed for measurements
        i=0
        num_meas = int((nsweeps-100)/10)
        E_array = np.zeros(num_meas)
        
        for n in range(nsteps):
            
            # KAWASAKI UPDATE - ATTEMPT SWAP ON EVERY "STEP"
            state = Kawasaki(state,T)
            
            # Measurements every 10 sweeps, after 100 sweeps equilib time
            
            equilib_time = 100*sweep_time
            if n >= equilib_time:
                if n%(10*sweep_time) == 0:
                    print("Progress: T = " + str(T) + ": " + str(n/sweep_time) + "/" + str(nsweeps) + " sweeps")
                    E_array[i] = totalE(state)
                    i += 1
    
        avg_E = np.mean(E_array)
        avg_E_error = np.std(E_array) / len(E_array)
        c = 1/(N**2*T**2) * (np.mean(E_array**2) - np.mean(E_array)**2)
        
        # Bootstrap for errors on heat capacity (c) 
        
        num_resamplings = 100
        c_array = np.zeros(num_resamplings)
        for k in range(num_resamplings):
            E_resampled = np.random.choice(E_array,num_meas)
            c_array[k] = 1/(N**2*T**2) * (np.mean(E_resampled**2) - np.mean(E_resampled)**2)
        
        c_error = np.sqrt(np.mean(c_array**2) - np.mean(c_array)**2)
        
        with open("Kawasaki_data.txt","a+") as f:
            f.write(str(T)+" "+str(avg_E)+" "+str(avg_E_error)+" "+str(c)+" "+str(c_error)+"\n")
    
else:
    print("Dynamics must be Glauber or Kawasaki")
    quit()







"""
def Initialise_random():
    
    state = np.zeros((N,N), dtype=int)
    for i in range(N):    
        for j in range(N):
            r = random.random()
            if (r<0.5): state[i,j] = 1
            if (r>0.5): state[i,j] = -1
        
    return state
"""


