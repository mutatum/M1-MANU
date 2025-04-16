# %%

import numpy as np
import matplotlib.pyplot as plt


i = open('GAUSS_PULSE_initial.txt')
f = open('GAUSS_PULSE_final.txt')

i_file_content = i.readlines()
Xi,Yi = np.array([list(map(float, l.split())) for l in i_file_content[2:-2]]).T
plt.figure(figsize=(6,6))
plt.plot(Xi, Yi, label="init")
f_file_content = f.readlines()
Xf,Yf = np.array([list(map(float, l.split())) for l in f_file_content[2:-2]]).T
plt.plot(Xf, Yf, label="final")
plt.legend()