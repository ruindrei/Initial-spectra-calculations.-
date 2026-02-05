## Attempt to fix parameters to run classy with massive neutrinos,

## default parameters. 

import numpy as np
from classy import Class
import matplotlib.pyplot as plt

z=0.0
omega_b = 0.02238280
omega_cdm = 0.1201075 ## This will change for non-default parameters. 
A_s = 2.100549e-09 
hh = 0.6781

params_default = {
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum.
    'z_pk':z,
    'z_max_pk':z,
    'h':hh,
    'A_s':2.100549e-09,
    'N_ncdm':0,
    'omega_b':omega_b, 
    'omega_cdm':omega_cdm, 
    'k_per_decade_for_pk': 50,  
    'P_k_max_1/Mpc': 300, 
    'Omega_k': 0.0

}

m_total_nu = np.linspace(0.138, 1.38, 10) # eV values of our 10 models
m_individual_nu = m_total_nu / 3.0  # eV (3 equal-mass neutrinos)
omega_ncdm = m_total_nu / 93.14  ## More accurate array of omega_ncdm for each of the 10 cases. 
omega_m = omega_b + omega_cdm ## Should be constant and the same even for non default cases. 

## These stay the same for all cases
common_params = { 
    'output': 'dTk vTk mPk', ## Want both density and velocity transfer functions as well as power spectrum. 
    'Neff':3.044, ## Should automatically scale N_ur down
    'N_ncdm':3,
    'T_ncdm':'0.71611,0.71611,0.71611',
    'A_s':2.100549e-09,
    'h':hh,
    'omega_b':0.02238280, 
    'omega_m':omega_m,
    'z_pk':z,
    'z_max_pk':z,
    'Omega_k': 0.0,
    'Omega_Lambda': 1.0 - (omega_m)/hh**2,
    'k_per_decade_for_pk': 50, 
    'P_k_max_1/Mpc': 300  
}



kk = np.geomspace(1e-5,200,num=500) ## Let's set up our k's in 1/Mpc

default=Class()
default.set(params_default)
default.compute()


Pk_default = default.get_pk_all(kk, z, nonlinear=False) ## Default P(k) in Mpc^3

Pk_list = []

for m in m_individual_nu:
    cosmo = Class()
    cosmo.set({**common_params, 'm_ncdm': f'{m},{m},{m}'})
    cosmo.compute()
    
    pk_massive = np.array([cosmo.pk_cb(k, z) for k in kk])
    Pk_list.append(pk_massive)

    cosmo.struct_cleanup()
    cosmo.empty()
    
Pk_neutrinos = np.array(Pk_list)  ## Shape (10, len(kk))
colors = ['red', 'green', 'blue', 'purple', 'cyan', 'yellow', 'black', 'orange', 'grey', 'brown'] ## colors for plotting.


## Recreate Fig 13. From Lesgourgues & Pastor 2006 review.
plt.xscale('log')
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P_{massive}(k)/P_{def}(k)$')
plt.xlim(1e-4, 1e2)
plt.ylim(0.,1.2)
plt.title("z = " + str(z))

for i in range(len(m_total_nu)):
    plt.plot(kk/hh, (Pk_neutrinos[i,:])/Pk_default, color = colors[i]) ## plot all of the neutrinos.
plt.show()


default.struct_cleanup()
default.empty()
