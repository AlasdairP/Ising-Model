# Modelling-and-Visualisation
Spring 2022 course

There are three Python scripts:

1. Ising.py

Usage: Ising.py <Glauber/Kawasaki> <System Size N (total spins=N**2)> <Temperature> <Number of sweeps>

Use this for a single temperature. It will show an animation and print the evolution of the energy.
It initialises the state randomly. J and k are set to 1 implicitly (i.e. not included in the code).
It is not used for measurements of the average energy <E>, magnetisation <M>, heat capacity c or susceptibility chi.

2. Ising_T_range.py

Usage: Ising_T_range.py <Glauber/Kawasaki> <System Size N (total spins=N**2)> <T_start> <T_end> <T_interval> <Number of sweeps>

This is used to generate the data for plotting - it loops through an array of temperature values.
J and k are set to 1 implicitly (i.e. not included in the code).
The state is initialised "all up" for the first temperature of the array (so a good choice is a low temperature like 1.0),
then subsequent temperatures are not re-initialised, they use the final state of the previous temperature.
It also takes measurements of E and M, and calculates the heat capacity and susceptibility and associated errors (via Bootstrap).
Measurements are taken every 10 sweeps, after an equilibrium period of 100 sweeps.
The Glauber and Kawasaki updates are done in exactly the same way as Ising.py.

3. Plotting.py

This is used to plot the data from the text files. It reads them in to a pandas dataframe and makes a plot.
It only plots one thing at a time so you have to edit it for plotting different things (dataG = Glauber, dataK = Kawasaki).
Usage: Plotting.py
