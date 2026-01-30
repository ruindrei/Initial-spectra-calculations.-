## File for generating WDM transfer and pk output files for MP-Gadget sims. 
## CLASS should automatically output a transfer file? The pk file I can make myself fairly easily.
## Let's try a 10keV WDM species from k=1.05167e-05 to 1.10271 h/Mpc. 

import numpy as np
from classy import Class
import matplotlib.pyplot as plt
import pandas as pd

## Set up our k ranges we can actually go a bit larger than the given transfer file
## Ranges of k really depend on the appropriate resolution scales we want to probe. 
## Likely we will need to ask about this. This is something we need to build an intuition for. 
kk = np.geomspace(1e-4,50,num=500) 

wdm_mass_eV = 10.0 * 1000
params_WDM_10keV = { ## Set up a WDM universe for 10keV mass particle
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'omega_cdm': 0,
    'root':'class_wdm10_', ## Set the file root output
    'transfer_verbose': 1, ## Should create a proper header file for the transfer function.
    'h': 0.67556,
    'omega_b': 0.022032,
    'N_ncdm': 1,
    'omega_ncdm': 0.12038,
    'A_s': 2.215e-9,
    'n_s': 0.9619,
    'tau_reio': 0.0925,
    'm_ncdm': wdm_mass_eV, ## make sure this in eV not keV
    'k_per_decade_for_pk': 50, ## Try setting this to 50 like Jimmy did. 
    'P_k_max_h/Mpc': 200,
    'z_pk':99,
    'ncdm_fluid_approximation':3, ## For ncdm, add these additional parameters to improve precision at small scales
    'ncdm_quadrature_strategy':3,
    'ncdm_maximum_q':15,
    'ncdm_N_momentum_bins':60,
    'l_max_ncdm':60

    }



params_mass_neutrinos = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'N_ncdm':1,
    'm_ncdm':0.06,
    'N_ur':2.0328, ## Make sure you reduce this from the default N_ur = 3.044. Otherwise you wil be adding new DOF. 
    'z_pk':99,
    'k_per_decade_for_pk': 50, ## Try setting this to 50 like Jimmy did. 
    'P_k_max_h/Mpc': 200  ## Try 200 for now

}

params_mass_neutrinos_06 = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'N_ncdm':3,
    'm_ncdm': '0.02, 0.02, 0.02', ## 2 eV is probably too high
    'N_ur':0.00441, ## Make sure you reduce this from the default N_ur = 3.044. Otherwise you wil be adding new DOF. 
    'z_pk':99,
    'A_s':2.100549e-09,
    'h': 0.67810,
    'omega_b':0.02238280, 
    'omega_m':0.1424903, # Should this be omega_b + omega_cdm? Yes
    'T_cmb':2.7255,
    'k_per_decade_for_pk': 50, ## Try setting this to 50 like Jimmy did. 
    'P_k_max_1/Mpc': 500  ## Try 200 for now

}

params_mass_neutrinos_6 = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'N_ncdm':3,
    'h': 0.67810,
    'm_ncdm': '0.2, 0.2, 0.2', ## 2 eV is probably too high
    'N_ur':0.00441, ## Make sure you reduce this from the default N_ur = 3.044. Otherwise you wil be adding new DOF. 
    'z_pk':99,
    'A_s':2.100549e-09,
    'omega_b':0.02238280,
    'T_cmb':2.7255,
    'omega_m':0.1424903, # This will set the total mass density to the correct value. omega_b + omega_cdm + omega_ncdm
    'k_per_decade_for_pk': 50, ## Try setting this to 50 like Jimmy did. 
    'P_k_max_1/Mpc': 500  ## Try 200 for now

}

params_default = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'z_pk':99,
    'h':0.6781,
    'A_s':2.100549e-09,
    'omega_b':0.02238280,
    'T_cmb':2.7255,
    'omega_cdm':0.1201075, # This will set the total mass density to the correct value
    'k_per_decade_for_pk': 50, ## Try setting this to 50 like Jimmy did. 
    'P_k_max_1/Mpc': 500  ## WARNING: MUST BE SET HIGHER THAN kk(1/Mpc)*h!!!
}

## We care about redshift 99 as we want to forward model our power spectrum. 
## However, we should also generate spectra at each specific snapshot for comparison to simulation outputs. 
## We can do this later, let's just make sure we can run this for now. 


## Let's also plot the precomputed files as well:

# Read the binary CLASS outputs we just generated (these files use k in [h/Mpc] and P in [Mpc/h]^3)
Pk_mass_nu = pd.read_csv("class_mass_neutrinos_6_00_pk.dat", skiprows=3, delim_whitespace=True)
Pk_def = pd.read_csv("class_default_00_pk.dat", skiprows=3, delim_whitespace=True)

# binary files: first column is k (h/Mpc), second column is P(k) (Mpc/h)^3
k_bin = Pk_def["#"]
Pk_bin = Pk_def["1:k"]

k_bin_nu = Pk_mass_nu["#"]
Pk_bin_nu = Pk_mass_nu["1:k"]


nu_06 = Class()
nu_06.set(params_mass_neutrinos_06)
nu_06.compute() ## Does not OUTPUT a file in Python, use this code to CHECK before changing the .ini parameter and just running the C code. 

nu_6 = Class()
nu_6.set(params_mass_neutrinos_6)
nu_6.compute()

default = Class()
default.set(params_default)
default.compute()

transfer_def = default.get_transfer(z=99.0)
T_def = transfer_def['d_tot']
h_def = default.h()

## Get and plot parameters 
transfer_nu_06 = nu_06.get_transfer(z=99.0)
T_nu_06 = transfer_nu_06['d_tot']
h_nu_06 = nu_06.h()

transfer_nu_6 = nu_6.get_transfer(z=99.0)
T_nu_6 = transfer_nu_6['d_tot']
h_nu_6 = nu_6.h()

print(h_def, h_nu_06, h_nu_6) ## These should all be equal. 

kk_hmpc = np.geomspace(1.04058e-05, 309.157, num=395)
# convert k array (given in [h/Mpc]) to 1/Mpc for CLASS get_pk_all
kk = kk_hmpc / h_def

# get P(k) from the python wrapper (returns P in [Mpc]^3 for k in 1/Mpc)
Pk_nu_06 = nu_06.get_pk_all(kk, 0., nonlinear=False)
Pk_def = default.get_pk_all(kk, 0., nonlinear=False)
Pk_nu_6 = nu_6.get_pk_all(kk, 0., nonlinear=False)

# convert python-Class outputs to the same units as the binary files:
#   k_plot = k (1/Mpc) * h_model -> [h/Mpc]
#   P_plot = P_class (Mpc^3) / h_model^3 -> [Mpc/h]^3
k_nu6_h = kk * h_nu_6
P_nu6_h = Pk_nu_6 / (h_nu_6 ** 3)

k_def_h = kk * h_def
P_def_h = Pk_def / (h_def ** 3)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k\ [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P(k)\ [\mathrm{Mpc}/h]^3$')
# plot python-Class outputs (converted)
plt.plot(k_nu6_h, P_nu6_h, color='red', label='Classy 0.6 eV')
plt.plot(k_def_h, P_def_h, color='black', label='Classy CDM')
# plot the binary CLASS outputs (already in [h/Mpc] and [Mpc/h]^3)
plt.plot(k_bin, Pk_bin, color='black', linestyle='--', label='CLASS CDM (binary)')
plt.plot(k_bin_nu, Pk_bin_nu, color='red', linestyle='--', label='CLASS 0.6 eV (binary)')
plt.legend()
plt.show()
plt.savefig('pk_comparison.png', dpi=200)
plt.close()


"""
plt.xscale('log')
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P_{massive}(k)/P_{def}(k)$')
#plt.plot(kk, Pk_nu_06/Pk_def, label = "0.06 eV")
plt.plot(kk* h_def, Pk_nu_6/Pk_def, label = "0.6 eV")
plt.legend()
plt.show()
"""
"""
## Read in the transfer function binary file.
transfer_file_6 = pd.read_csv("class_mass_neutrinos_6_00_tk.dat", skiprows=10, delim_whitespace=True)
transfer_file_default = pd.read_csv("class_def_00_tk.dat", skiprows=10, delim_whitespace=True)

print(transfer_file_default.columns)



k = transfer_def['k (h/Mpc)']
T_total_matter = transfer_def['d_tot']

## Let's plot the transfer functions as well.
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k\ [h/\mathrm{Mpc}]$')
plt.ylabel(r'$|T(k)|$')
#plt.plot(transfer_file_6["#"], np.abs(transfer_file_6["12:psi"]), color='red', label='Classy 0.6 eV') 
plt.plot(transfer_file_default['1:k'], np.abs(transfer_file_default['7:d_tot']), color='red', label='CLASS CDM') 
plt.plot(k, np.abs(T_total_matter), color='black', alpha = 0.5, label='Classy CDM')
plt.legend()
plt.plot()
plt.show()
"""

## Always empy the struct to prevent memory leak
nu_06.struct_cleanup()
nu_06.empty()

default.struct_cleanup()
default.empty()

nu_6.struct_cleanup()
nu_6.empty()


## Now let's attempt to run the output simulation, with the new base CLASS power spectra. 

## Easier to just do a .ini file with the same parameters and run the binary. Use python classy wrapper for plotting only 

## This WORKS WELL BUT IT IS VERY SLOW. Let's start running the .ini file at around 3:30 pm. 
