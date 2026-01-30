## Compare CDM spectra from CLASS binary and Classy wrapper

import numpy as np
from classy import Class
import matplotlib.pyplot as plt
import pandas as pd ## probably easier to just use np.loadtxt 

## Default parameters for LambdaCDM we want to use. 
params_default = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'z_pk':99,
    'h':0.6781,
    'A_s':2.1e-09, ## Check precision settings on this? 
    'omega_b':0.02238280,
    'T_cmb':2.7255,
    'omega_cdm':0.1201075, # This will set the total mass density to the correct value
    'k_per_decade_for_pk': 50, 
    'P_k_max_1/Mpc': 500  ## WARNING: MUST BE SET HIGHER THAN kk(1/Mpc)*h!!!
    ## Careful with P_k_max_h/Mpc vs P_k_max_1/Mpc.
}

## Generate the CLASS binary file using the exact same parameters. 
## Should be in params_default.ini
## Two files generated are class_default_00_pk.dat and class_default_00_tk.dat

h=0.6781 ## Should be identitcal for all runs. 

## Binary CLASS output
default_binary = np.loadtxt("class_default_00_pk.dat", skiprows=4)

k_bin_h = default_binary[:,0]  ## in h/Mpc
Pk_bin_h = default_binary[:,1]  ## in (Mpc/h)^3


## Use classy to compute
default = Class()
default.set(params_default)
default.compute()


k_bin = k_bin_h / h  ## convert to 1/Mpc for classy
## Compute P(k) using classy
## We want to use the same k values as the binary output (some will be too large so just truncate)
## Let's just do it one value at a time instead of using get_pk_all.

wrapper_pk = []
for k_val in k_bin:
    try:
        pk = default.pk(k_val, 99) ## should return in Mpc^3
        wrapper_pk.append(pk)
    except:
            wrapper_pk.append(np.nan) ## ignore the large k values. 
    
Pk_classy = np.array(wrapper_pk)


## Create a comparison plot
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k\ [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k)\ [\mathrm{Mpc}/h]^3$')
plt.plot(k_bin_h, Pk_bin_h, color= "black", label= "CLASS Binary")
plt.plot(k_bin_h, Pk_classy/h**3, linestyle='--', color= "gray", label= "Classy Wrapper") ## Convert units 
plt.legend() 
plt.show()