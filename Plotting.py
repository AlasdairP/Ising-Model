"""
Plotting for MVP Checkpoint 1 - Ising model

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

col_names_G = ["T","E","E_error","c", "c_error","M","M_error","chi","chi_error"] 
dataG = pd.read_csv("Glauber_data.txt",sep="\s+",dtype=np.float64,names=col_names_G)

col_names_K = ["T","E","E_error","c","c_error"]
dataK = pd.read_csv("Kawasaki_data.txt",sep="\s+",dtype=np.float64,names=col_names_K)

plt.errorbar(dataG['T'],dataG['chi'],dataG['chi_error'], ecolor='k', capsize=2)
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.title("Glauber, susceptibility")
plt.show()

