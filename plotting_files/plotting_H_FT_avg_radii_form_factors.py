import pandas as pd
import numpy as np
from numpy import sqrt, pi, exp, log
import os
import re
import matplotlib.pyplot as plt
from scipy import fft, interpolate

#The way this is written currently, it has to be run within the Q2equalsWhatever directories. 
#You can do this anywhere on your machine if you just fully specify the path, using the method for obtaining
#the data file information given in H_versus_Q2_fixed_xt.py

#function for all the x values
def xvals_func(zeta): 
    xval = []
    i_range = np.arange(0, 490, 1)

    for i in i_range: 
        if(i <= 290): 
            yyy = log(1.e+4)*(330.-(i)+1.)/330.
            xv = exp(-yyy)
            xval.append(xv)
        else: 
            xstart = exp(-log(1.e+4)*41./330.)
            xv = xstart + ((i)-290.)*(1.-xstart)/201.
            xval.append(xv)
    
    return np.array(xval)


#first, store all the files in df_tval kind of dataframes
#then, create 490 df_Xval kind of data frames
#finally, store the nth row of the df_tval dataframes into the nth df_Xval dataframe. add a column to each df_Xval for the t values

df_tval = {} #creating dictionary for dataframes

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

path = os.getcwd()
files_in_dir = sorted(os.listdir(path))

for file_name in files_in_dir: #loop that creates all the tval dataframes
        if file_name.endswith('.dat') and file_name.startswith('GPD_H'):
            tval = get_numbers_from_filename(file_name) #of the form 000, 005, 010, ..., 200
            tval = float(tval)/100
            df_tval[tval] = pd.read_csv(file_name, header=None, delim_whitespace=True)
            df_tval[tval] = df_tval[tval][:-1]

#the above way of obtaining the data files is not as versatile as the way done in H_versus_Q2_fixed_xt.py
#I will probably end up rewriting the above in terms of what's given in H_versus_Q2_fixed_xt.py, but it works for now

df_Xval = {} #keys are numbers 0, 1, 2, ..., 489

for x in np.arange(490): #create all the dataframes organized by x values
    df_Xval[x] = pd.DataFrame(columns = np.arange(12))

for l in df_Xval.keys():
    for k in df_tval.keys():
        df_Xval[l] = pd.concat([df_Xval[l], df_tval[k].loc[[l]]])
        df_Xval[l] = df_Xval[l].reset_index(drop=True)

def f_2D(Delta_x, Delta_y, func): #create a 2D function given a 1D function; seems superfluous, but works with the way we've been doing the Fourier transforms with the python package
    minus_t = Delta_x**2 + Delta_y**2
    output_arr = np.zeros_like(minus_t)
    for i in np.arange(0, len(minus_t), 1): 
        for j in np.arange(0, len(minus_t[i]), 1):
            if minus_t[i][j] <= 2:
                output_arr[i][j] = func(minus_t[i][j])
            else: 
                output_arr[i][j] = 0
    return output_arr

#now for the actual Fourier transformation! 
def fourier(x, n): #x is the ~index~ of the x value you're interested in, n is the column of the data file
    H_for_some_x = df_Xval[x].iloc[:, n].to_numpy()
    minus_t_arr = np.arange(0, 2.05, 0.05)
    f = interpolate.interp1d(minus_t_arr, H_for_some_x) #interpolating H as a function of -t

    Nx = 100
    Ny = 100

    Delta_x = np.linspace(-10, 10, Nx) #change back to 10 #1000 for Green data
    Delta_y = np.linspace(-10, 10, Ny)

    Delta_x_v, Delta_y_v = np.meshgrid(fft.ifftshift(Delta_x), fft.ifftshift(Delta_y))
    fk = f_2D(Delta_x_v, Delta_y_v, f) #creating a 2D function, which we need for the Fourier transform 

    fk_FFT = fft.fftshift(np.diff(Delta_x)[0]*np.diff(Delta_y)[0]*fft.fft2(fk)/((2*np.pi)**2))
    xnew = fft.fftshift(fft.fftfreq(len(Delta_x), np.diff(Delta_x)[0] / (2*np.pi)))
    ynew = fft.fftshift(fft.fftfreq(len(Delta_y), np.diff(Delta_y)[0] / (2*np.pi)))

    xvnew, yvnew = np.meshgrid(xnew, ynew)

    H_FT = np.real(fk_FFT)
    bx = xnew/5.068
    by = ynew/5.068
    return [bx, by, H_FT, xvnew, yvnew]

#A few of the x values specificed as (index of x, x value)
##(2, 0.0001) (85, 0.001) (167, 0.01) (250, 0.1) (288, 0.3) (345, 0.5) (433, 0.8)

#A few of the columns of the data files
#n = 1 is Hu, n = 4 is Hd, and n = 7 is Hg

def plot_FT(x, n, num_rows, num_cols, index_k): #need to write plt.show() after using the function
    plt.subplot(num_rows, num_cols, index_k)
    plt.contourf(fourier(x, n)[0], fourier(x, n)[1], fourier(x, n)[2], cmap = 'Greens')
    plt.colorbar()
    plt.xlim(-0.6, 0.6)
    plt.ylim(-0.6, 0.6)
    plt.xlabel(r'$b_{x} [fm]$', size = 20)
    plt.xticks(size = 10)
    plt.ylabel(r'$b_{y} [fm]$', size = 20)
    plt.yticks(size = 10)

#An inefficient way to make subplots; you could just edit this to make individual plots
num_rows, num_cols = 2, 2
plot_FT(2, 7, num_rows, num_cols, 1)
plt.title(r'X = 0.0001')
plot_FT(85, 7, num_rows, num_cols, 2)
plt.title(r'X = 0.001')
plot_FT(167, 7, num_rows, num_cols, 3)
plt.title(r'X = 0.01')
plot_FT(250, 7, num_rows, num_cols, 4)
plt.title(r'X = 0.1')
plt.suptitle(r'$H^{g}, $' r'$Q^{2} = 5$' " GeV" r'$^{2}$', fontsize = 15)
plt.subplots_adjust(top= 0.92)
plt.show()

def avg_radius_func(x, n): #function for computing the average radius given the index of the x value and the column n

    bx = fourier(x, n)[3]
    FT_H_grid = fourier(x, n)[2]

    bT = []
    FT_H = []

    for i in np.arange(0, bx[50].size, 1): 
        if bx[50][i] >= 0:
            bT.append(bx[50][i])
            FT_H.append(FT_H_grid[50][i])

    bT = np.array(bT)
    FT_H = np.array(FT_H)

    int = 0
    int_num = 0
    int_den = 0
    delta_bT = (1.)*0.444

    for j in np.arange(bT.size): 
        int_temp_num = delta_bT*(bT[j]**3)*FT_H[j]
        int_num = int_temp_num + int_num
        int_temp_den = delta_bT*bT[j]*FT_H[j]
        int_den = int_temp_den + int_den
        
    int = int_num/int_den
    avg_radius = sqrt(int)/5.068
    return avg_radius

#The next few functions are just for checking that we obtain the right F1 and Ag
def integrate_data(X, Y): 
    int_arr = 0
    for i in np.arange(0, len(X) - 1, 1): 
        width = X[i+1] - X[i]
        height = (Y[i+1] + Y[i])/2
        area = width*height
        int_arr = int_arr + area
    return int_arr

def F_1_check(n, minus_t):
    H_func_of_X = np.array(df_tval[minus_t].iloc[:, n])
    X_vals = np.array(df_tval[minus_t].iloc[:, 0])
    return integrate_data(X_vals, H_func_of_X)

def A_g_check(n, minus_t): 
    H_func_of_X = np.array(df_tval[minus_t].iloc[:, n])
    X_vals = np.array(df_tval[minus_t].iloc[:, 0])
    return integrate_data(X_vals, H_func_of_X/X_vals)
#end functions meant for checking form factors----------------------------------------------

#You could ignore these next commented out lines - they're just for sending the calculations
#of the average radii to data files

#obtaining the average radii ----------------------------------------------------------------------------------------------------------------------------------
#avg_radius_arr_u = []
#n = 1 #n = 7 Hg, n = 1 Hu, n = 4 Hd
#for x in np.arange(490): 
#    avg_radius_arr_u.append(avg_radius_func(x, n))


#avg_radius_arr_d = []
#n = 4 #n = 7 Hg, n = 1 Hu, n = 4 Hd
#for x in np.arange(490): 
#    avg_radius_arr_d.append(avg_radius_func(x, n))


#avg_radius_arr_g = []
#n = 7 #n = 7 Hg, n = 1 Hu, n = 4 Hd
#for x in np.arange(490): 
#    avg_radius_arr_g.append(avg_radius_func(x, n))

#----------------------------------------------------------------------------------------------------------------------------------

#plotting the average radii: ----------------------------------------------------------------------------------------------------------------------------------
#plt.plot(xvals_func(0), avg_radius_arr_u, c = 'Orange', label = 'Hu')
#plt.plot(xvals_func(0), avg_radius_arr_d, c = 'Blue', label='Hd')
#plt.plot(xvals_func(0), avg_radius_arr_g, c = 'Green', label='Hg')

#plt.xlabel('X [0.0001, 0.5]', fontsize=20)
#plt.xlim(right = 0.5)
#plt.xlim(left = 0.0001)
#plt.ylim(bottom = 0.3)
#plt.xticks(fontsize = 18)
#plt.ylabel(r'$<b_{T}^{2}>^{1/2}$' ' [fm]', fontsize=20)
#plt.yticks(fontsize = 18)
#plt.xscale("log")
#plt.title('Average radii of gluons and quarks', fontsize=20)
#plt.text(0.1, 0.7,  r'$Q^{2} = 5$' " " r'GeV$^{2}$', fontsize=25)
#plt.legend(prop={"size":20})
#plt.show()
#----------------------------------------------------------------------------------------------------------------------------------