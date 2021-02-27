#-------------------------------------------------------------
# Script for plotting upper limits on the 21 cm power spectrum
# Author: Pablo Villanueva Domingo
#-------------------------------------------------------------

import glob, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

# Reduced Hubble rate (for units)
hlittle = 0.7
# Modeling error for sensitivity (see 1501.06576)
moderror = 0.1

# Load the power spectrum from 21cmFAST, at a given scale k (indicated by the index kindex)
def powerspectrum(path, kindex):

    powk = []
    for nz, zz in enumerate(z):
        file = glob.glob(path+"/Deldel_T_power_spec/ps_z{:06.2f}_nf*".format(zz))[0]
        tabpow = np.loadtxt(file, unpack=True)
        k = tabpow[0][kindex]
        powkk = tabpow[1][kindex]
        powk.append(powkk)

    return k, np.array(powk)

# Load sensitivity file (from 21cmSense) for a given experiment
def sensitivity(exp, moderror):

    table = np.loadtxt("data/Sensitivity_"+exp+".dat", unpack="True")
    zz, kk = np.unique(table[0]), np.unique(table[1])*hlittle   # NoiseVals_SKA in h*Mpc^-1, but 21cmFAST in Mpc^-1
    pow = table[2].reshape(len(zz), len(kk))
    logsenseSKA = interp2d(zz,np.log(kk),np.log(np.transpose(pow)),kind="linear")
    sense = np.exp(logsenseSKA( z, np.log(k) ))
    err = np.sqrt( moderror**2.*powk**2. + sense**2. )
    return err

# Experiments with sensitivity from 21cmSense
experiments = ["HERA", "SKA"]
alpha_exp = [0.3, 0.7]


fig, (axPs) = plt.subplots( figsize=(10,6) )


#--- FIDUCIAL MODEL ---#
# Fiducial theoretical model generated from 21cmFast

# Astrophysical parameters
# Mturn
mturn_vals = [1.e8]
# Number of ionizing photons per baryon
ngamma_vals = [4.e3]
# X-ray luminosity
lx_vals = [40.]

# Element index of the k wavenumber at which the power spectrum is plotted
k_index_list = [5]

for mturn in mturn_vals:
    for ngamma in ngamma_vals:
        for lx in lx_vals:

            #print("Mturn_{:.3e}_lumx_{:.3f}_ngamma_{:.3e}".format(mturn,lx,ngamma))
            path = "21cmFAST_files/Outputs_Mturn_{:.3e}_lumx_{:.3f}_ngamma_{:.3e}".format(mturn,lx,ngamma)

            if (os.path.exists(path)==0):
                print("Path "+path+" does not exist.")
                continue

            # Import the global file
            globfile = glob.glob(path+"/Ts_outs/glob*")[0]
            tab = np.loadtxt( globfile, unpack=True )

            z, xH, Tk, xe, Ts, Tcmb = tab[0], tab[1], tab[2], tab[3], tab[4], tab[5]
            xHtot = xH*(1.-xe)
            dTb = 27.*xHtot*(1.-(Tcmb/Ts))*np.sqrt((1.+z)/10.)

            # Plot the power spectrum
            for i, kindex in enumerate(k_index_list):

                k, powk = powerspectrum(path, kindex)
                #print(k)

                axPs.plot(z, powk, linestyle="-", color="black", label="Fiducial model")

                for j, exp in enumerate(experiments):
                    err = sensitivity(exp, moderror)
                    axPs.fill_between(z, powk + err, powk - err, color="grey", alpha=alpha_exp[j], label=exp+" sensitivity")



#--- EXPERIMENTS ---#
# Plot upper limits from different interferometers
# See README for the references of the experiments, the range of k and other details

# Load experiments list
list_data_file = "list_data.txt"
listdata = np.loadtxt(list_data_file, dtype=str, delimiter=";", unpack=True)

# Remove tabs and add 'r' prefix for LaTeX formulae
for i, data in enumerate(listdata):
    for j, text in enumerate(data):
        listdata[i,j]=r""+text.replace("\t","")

datanames = listdata[0]
labels = listdata[1]
colors = listdata[2]
markers = listdata[3]

for i, data in enumerate(datanames):

    data = np.loadtxt("data/"+data+".txt", unpack=True)
    axPs.scatter(data[0], data[1], marker=markers[i], color=colors[i], label=labels[i])
    #axPs.arrow(data[0], data[1], 0., -data[1]/10., color=colors[i])

    # Plot bars for extended ranges of redsfhit
    if len(data)>2:
        axPs.hlines(y=data[1], xmin=data[2], xmax=data[3], color=colors[i], zorder=1.e-1)


axPs.grid(True, linestyle=":", zorder=1.e-2)
axPs.set_yscale("log")
axPs.legend(bbox_to_anchor=(1.01, 0.9), borderaxespad=0., fontsize=10)
axPs.set_xlim([6.,29.])
axPs.set_ylabel(r"$\overline{\delta T_{\rm b}}^2 \Delta_{21}^2 \; [mK^2]$")
axPs.set_xlabel(r"$z$")

plt.savefig("plot_21ps_constraints.pdf", bbox_inches='tight')
plt.show()
plt.gcf().clear()
