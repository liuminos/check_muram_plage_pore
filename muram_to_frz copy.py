import numpy as np 
from astropy.io import fits
import sys
import muram as mio
import firtez_dz as frz

muramdir = '/dat/milic/muram_from_rebecca/' # can be changed

#iter = int(sys.argv[1]) #
iter = 63000


#select height range 
hh = 430 #
hl = 350 #
xh = 600
yh = 600

# Height z
z = np.arange(hh-hl) * 12 #km

# Temp
T = mio.MuramCube(muramdir, iter, 'Temp')
T = T.transpose(1,2,0)
T = T[:xh,:yh,hl:hh]

print(T.shape)
nx = T.shape[0]//3
ny = T.shape[1]//3
nz = T.shape[2]

#frz model
model3D = frz.atm_model3D(nx,ny,nz)
model3D.z = z[None,None,:]
model3D.tem = T.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
del(T)
print('T')

# Gas pressure
P = mio.MuramCube(muramdir, iter, 'Pres')
P = P.transpose(1,2,0)
P = P[:xh,:yh,hl:hh]
model3D.pg = P.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
del(P)
print('P')

# Vlos
vz = mio.MuramCube(muramdir, iter, 'vx')
vz = vz.transpose(1,2,0)
vz = vz[:xh,:yh,hl:hh]
model3D.vz = -vz.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
del(vz)
print('vz')

# Magnetic field
Bx = mio.MuramCube(muramdir, iter, 'By')
Bx = Bx.transpose(1,2,0)
Bx = Bx[:xh,:yh,hl:hh]
Bx = Bx * np.sqrt(4 * np.pi)
By = mio.MuramCube(muramdir, iter, 'Bz')
By = By.transpose(1,2,0)
By = By[:xh,:yh,hl:hh]
By = By * np.sqrt(4 * np.pi)
Bz = mio.MuramCube(muramdir, iter, 'Bx')
Bz = Bz.transpose(1,2,0)
Bz = Bz[:xh,:yh,hl:hh]
Bz = Bz * np.sqrt(4 * np.pi)

model3D.bx = Bx.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
model3D.by = By.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
model3D.bz = Bz.reshape(nx,3,ny,3,nz).mean(axis=(1,3))
del(Bx,By,Bz)
print('B')

# tau
#tau = mio.MuramCube(muramdir, iter, 'tau')
#tau = tau.transpose(1,2,0)
#tau = tau[:,:,hl:hh]

#model3D.tau = tau
#del(tau)

#write
model3D.write_model('/dat/xenosh/MiHi_Fe_I_plage/simulation/qs_nx200_ny200_nz80_dz12.bin')

print('END')
