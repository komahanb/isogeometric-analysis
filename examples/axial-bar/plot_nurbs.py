import matplotlib.pyplot as plt
import numpy as np
from basis import N, Nprime

###################################################################
# plot results
###################################################################

# Configure 
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing
# font errors
params = {
  'axes.labelsize': 20,
  'legend.fontsize': 14,
  'xtick.labelsize': 20,
  'ytick.labelsize': 20,
  'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']
#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# Set marker size
markerSize = 7.0 #11.0
mew = 2.0

# bar chart settings
lalpha    = 0.9
rev_alpha = 0.9/1.5

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)


xi1 = np.linspace(0.0, 0.5, num=10)
N1 = N(xi1,1,1)
N2 = N(xi1,2,1)
N3 = N(xi1,3,1)

xi2 = np.linspace(0.5, 1.0, num=10)
N4 = N(xi2,1,2)
N5 = N(xi2,2,2)
N6 = N(xi2,3,2)

plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(xi1, N1 , '-o', label=r'$N_1^1$' , mew=mew, ms=markerSize, color=tableau20[2], mec='black')
plt.plot(xi1, N2 , '-o', label=r'$N_2^1$' , mew=mew, ms=markerSize, color=tableau20[4], mec='black')
plt.plot(xi1, N3 , '-o', label=r'$N_3^1$' , mew=mew, ms=markerSize, color=tableau20[6], mec='black')
plt.plot(xi2, N4 , '-s', label=r'$N_1^2$' , mew=mew, ms=markerSize, color=tableau20[2], mec='black')
plt.plot(xi2, N5 , '-s', label=r'$N_2^2$' , mew=mew, ms=markerSize, color=tableau20[4], mec='black')
plt.plot(xi2, N6 , '-s', label=r'$N_3^2$' , mew=mew, ms=markerSize, color=tableau20[6], mec='black')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$N(\xi)$')
plt.legend(loc='best', ncol=2, frameon=False)
plt.savefig('bar-basis-functions.pdf', bbox_inches='tight', pad_inches=0.05)

xi1 = np.linspace(0.0, 0.5, num=10)
N1 = Nprime(xi1,1,1)
N2 = Nprime(xi1,2,1)
N3 = Nprime(xi1,3,1)

xi2 = np.linspace(0.5, 1.0, num=10)
N4 = Nprime(xi2,1,2)
N5 = Nprime(xi2,2,2)
N6 = Nprime(xi2,3,2)

plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(xi1, N1 , '-o', label=r'$N_1^{1,\xi}$' , mew=mew, ms=markerSize, color=tableau20[2], mec='black')
plt.plot(xi1, N2 , '-o', label=r'$N_2^{1,\xi}$' , mew=mew, ms=markerSize, color=tableau20[4], mec='black')
plt.plot(xi1, N3 , '-o', label=r'$N_3^{1,\xi}$' , mew=mew, ms=markerSize, color=tableau20[6], mec='black')
plt.plot(xi2, N4 , '-s', label=r'$N_1^{2,\xi}$' , mew=mew, ms=markerSize, color=tableau20[2], mec='black')
plt.plot(xi2, N5 , '-s', label=r'$N_2^{2\xi}$' , mew=mew, ms=markerSize, color=tableau20[4], mec='black')
plt.plot(xi2, N6 , '-s', label=r'$N_3^{2,\xi}$' , mew=mew, ms=markerSize, color=tableau20[6], mec='black')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\partial{N(\xi)}/\partial {\xi}$')
plt.legend(loc='upper center', ncol=2, frameon=False)
plt.savefig('bar-basis-function-derivatives.pdf', bbox_inches='tight', pad_inches=0.05)





