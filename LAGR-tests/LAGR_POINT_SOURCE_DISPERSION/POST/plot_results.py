
#####################################################################################
#                                                                                   #
#                   Program to plot the results of the Lagrangian module            #
#                   in the case of point-source dispersion                          #
#                   with precribed values of T_L, tau_p and B_x                     #
#                                                                                   #
#####################################################################################

# Load libraries
import numpy as np
import matplotlib.pyplot as plt

#####################################################################################
#                                                                                   #
#           Parameters to be filled by Users                                        #
#                                                                                   #
#####################################################################################

# Simulation
# Folder with results
resu_folder = "%s" % ("../SCHEME_ORDER_1/RESU/")
simu_folder = "%s" % ("20200318-0941/")

# Name of the file to be browsed
file_name = "results_to_plot"           # File extension for results
path = "%s%s%s" % (resu_folder, simu_folder, file_name)

#####################################################################################
#                                                                                   #
#           Main program                                                            #
#           (DO NOT MODIFY)                                                         #
#                                                                                   #
#####################################################################################

# Load and analyze results
resu = np.loadtxt(path, skiprows = 1)

N_rows = resu.shape[0]
N_colums = resu.shape[1]

# Plot results
# 4 figures (for Xp_Xp, Up_Up, Us_Us and Up_Us correlations)

# Figure 1: Xp_Xp correlations
for iplot in range(0, 4):
    plt.figure(1)
    plt.plot(resu[:,0],
             resu[:,iplot+1],
             linestyle = '-',
             linewidth = 2,
             label = 'Analytical solution' )
    plt.plot(resu[:,0],
             resu[:,iplot+5],
             linestyle = '--',
             linewidth = 2,
             label = 'Scheme order 1' )
    # Plot environment
    plt.xlabel('Time (in s)')
    if (iplot == 0):
        plt.ylabel('$<X_p^2>$ in $m^2/s^2$')
        plt.title('Particle position')
        plt.legend(loc=2, borderaxespad=0.)
    elif (iplot == 1):
        plt.ylabel('$<U_p^2>$ in $m^2/s^2$')
        plt.title('Fluid velocity')
        plt.legend(loc=2, borderaxespad=0.)
    elif (iplot == 2):
        plt.ylabel('$<U_s^2>$ in $m^2/s^2$')
        plt.title('Fluid velocity seen')
        plt.legend(loc=2, borderaxespad=0.)
    elif (iplot == 3):
        plt.ylabel('$<U_p U_s>$ in $m^2/s^2$')
        plt.title('Velocity covariance')
        plt.legend(loc=2, borderaxespad=0.)
    plt.grid(True)
    plt.show()
