import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

H_Q2_dict = {} #python dictionary that will store each data file
Q2_arr = np.array([2, 5, 10, 15, 20, 50, 70, 100]) #array for all the Q2 values
N = Q2_arr.size #number of data files

for i in np.arange(0, N, 1):
    Q2_val = str(int(Q2_arr[i]))
    file_name = "/Users/tosh/Documents/Codes/Femtonet/Others/gpd-calc/gpd-data/Qsq%s/GPD_H000.dat" % Q2_val #need to specify the path for your own files
    H_Q2_dict[Q2_arr[i]] = np.loadtxt(file_name, unpack=True) #loading all the data files into the dictionary

x_values = H_Q2_dict[2][0] #an array for all the x values

n = 7 # column number for Hg; n = 1 for Hu, n = 4 for Hd

limited_x_indices = np.array([1, 85, 167, 250, 288, 345, 433])
for x_index in limited_x_indices:
    Hg_x_func_of_Q2 = np.array([H_Q2_dict[2][n][x_index], H_Q2_dict[5][n][x_index], H_Q2_dict[10][n][x_index], H_Q2_dict[15][n][x_index], H_Q2_dict[20][n][x_index], H_Q2_dict[50][n][x_index], H_Q2_dict[70][n][x_index], H_Q2_dict[100][n][x_index]])
    x_val_str = str(x_values[x_index])
    plt.scatter(Q2_arr, Hg_x_func_of_Q2, label='$X = %s $' % x_val_str)


plt.xlabel(r'$Q^{2}$', fontsize = 20)
plt.ylabel(r'$H^{g}$', fontsize = 20)
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()
