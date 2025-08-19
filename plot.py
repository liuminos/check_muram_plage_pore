import numpy as np
import matplotlib.pyplot as plt
import firtez_dz as frz
from astropy.io import fits
from matplotlib.colors import CenteredNorm, TwoSlopeNorm

plt.rc('text', usetex = False)

muramdir = '/dat/milic/MURaM_400G_plage_pore/'

time = np.fromfile(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_time_012000_to_044450.dat',dtype=np.float64)
nsnaps = time.size //2
time = time.reshape(2,nsnaps)

nx = 1536
ny = 1536
nz = 1056
# Continuum intensity
Ic = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_Iout_bref_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,nx,ny))
I = np.zeros_like(Ic)
# normalize
#for i in range(nsnaps):
#    mean = np.mean(Ic[i,800:1200,100:500])
#    I[i,:,:] = Ic[i,:,:] / mean

# Bz
Bz = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_tau_1.000_bx_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,nx,ny))
Bz = Bz * np.sqrt(4*np.pi)
#Bz_new = Bz.reshape(nsnaps,nx//8,8,ny//8,8).mean(axis=(2,4))
# Vz
Vz = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_tau_1.000_vx_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,nx,ny))
Vz = - Vz/1e5

for i in range(nsnaps):
    mean = np.mean(Ic[i,800:1200,100:500])
    I[i,:,:] = Ic[i,:,:] / mean
    plt.clf()
    plt.figure(figsize=(17,5))
    plt.subplot(131)
    plt.imshow(I[i,:,:].T,origin='lower',cmap='gray',vmin=0.2,vmax=2.,extent=[0,10,0,10])
    plt.colorbar()
    plt.xlabel('x [Mm]')
    plt.ylabel('y [Mm]')
    plt.title(r'$I_{500}$')
    plt.subplot(132)
    plt.imshow(Bz[i,:,:].T,origin='lower',cmap='PuOr',norm=TwoSlopeNorm(vmin=-2000,vcenter=0,vmax=2000),extent=[0,10,0,10])
    plt.colorbar()
    plt.xlabel('x [Mm]')
    plt.title(r'Bz [G] at $\tau=1$')
    plt.subplot(133)
    plt.imshow(Vz[i,:,:].T,origin='lower',cmap='bwr',vmin=-6,vmax=6,extent=[0,10,0,10])
    plt.colorbar()
    plt.xlabel('x [Mm]')
    plt.title(r'Vz [km/s] at $\tau=1$')
    plt.suptitle('t = '+ str("%.2f" %time[1,i]) + 's')
    plt.tight_layout()
    plt.savefig('/dat/xenosh/muram_plage_pore/plots/img'+str(i)+'.png',bbox_inches='tight')
'''
for i in range(nsnaps):
    plt.clf()
    plt.figure(figsize=(12,10))
    plt.imshow(Bz_new[i,:,:].T,origin='lower',cmap='PuOr',norm=TwoSlopeNorm(vmin=-500,vcenter=0,vmax=2000),extent=[0,10,0,10])
    plt.colorbar()
    plt.xlabel('x [Mm]')
    plt.ylabel('y [Mm]')
    plt.title('t = '+ str("%.2f" %time[1,i]) + 's')

    plt.tight_layout()
    plt.savefig('/dat/xenosh/muram_plage_pore/plots_asym/Bznew'+str(i)+'.png',bbox_inches='tight')
'''