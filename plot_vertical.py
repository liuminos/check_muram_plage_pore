import numpy as np
import matplotlib.pyplot as plt
import firtez_dz as frz
from astropy.io import fits
from matplotlib.colors import CenteredNorm, TwoSlopeNorm
from tqdm import tqdm

plt.rc('text', usetex = False)

muramdir = '/dat/milic/MURaM_400G_plage_pore/'

time = np.fromfile(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_time_012000_to_044450.dat',dtype=np.float64)
nsnaps = time.size //2
time = time.reshape(2,nsnaps)

nx = 1536
ny = 1536
nz = 1056

bx = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_z_0768_by_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,ny,nz))
bz = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_z_0768_bx_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,ny,nz))

vx = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_z_0768_vy_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,ny,nz))
vz = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_z_0768_vx_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,ny,nz))

tau = np.memmap(muramdir+'2D/Pore_10Mm_6x6km_res_Bz400G_z_0768_tau_012000_to_044450.dat',dtype=np.float32,mode='r',shape=(nsnaps,ny,nz))

# x 768
for t in tqdm(range(nsnaps)):
    plt.clf()
    plt.figure(figsize=(20,12))
    plt.subplot(211)
    plt.imshow(bz[t,:,700:900].T*np.sqrt(4*np.pi),origin='lower',cmap='PuOr',norm=TwoSlopeNorm(vmin=-500,vcenter=0,vmax=2000))
    plt.colorbar(aspect=40)
    stride = 10
    #xsize,ysize = bz[0,:,700:900].shape
    xmax = nx//stride+1
    ymax = 200//stride
    X,Y = np.meshgrid(stride*np.arange(0,xmax),stride*np.arange(0,ymax))
    #plt.quiver(X,Y,bx[t,::stride,700:900:stride].T,bz[t,::stride,700:900:stride].T,width=.004,color='black')
    plt.streamplot(X, Y, bx[t,::stride,700:900:stride].T,bz[t,::stride,700:900:stride].T,density=[10,1],color='black')
    plt.contour(bz[t,:,700:900].T*np.sqrt(4*np.pi),levels=[-50],colors=['brown'])
    plt.contour(tau[t,:,700:900].T,levels=[0.01,0.1,1],colors=['blue','green','red'],alpha=0.6)
    plt.xlim([800,1100])
    plt.ylim([50,170])
    plt.subplot(212)
    plt.imshow(-vz[t,:,700:900].T/1e5,origin='lower',cmap='bwr',norm=TwoSlopeNorm(vmin=-10,vcenter=0,vmax=10))
    plt.colorbar(aspect=40)
    #plt.quiver(X,Y,vx[t,::stride,700:900:stride].T,vz[t,::stride,700:900:stride].T,width=.004,color='black')
    plt.streamplot(X, Y, vx[t,::stride,700:900:stride].T,vz[t,::stride,700:900:stride].T,density=[10,1],color='black')
    plt.contour(tau[t,:,700:900].T,levels=[0.01,0.1,1],colors=['blue','green','red'],alpha=0.6)
    plt.xlim([800,1100])
    plt.ylim([50,170])
    plt.suptitle('t = '+ str("%.2f" %time[1,t]) + 's')
    plt.tight_layout()
    plt.savefig('/dat/xenosh/muram_plage_pore/vertical/bzvz'+str(t)+'.png',bbox_inches='tight')