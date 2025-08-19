import numpy as np 
from astropy.io import fits
import sys
import muram as mio
import firtez_dz as frz

muramdir = '/dat/milic/MURaM_400G_plage_pore/3D/' # can be changed

#iter = int(sys.argv[1]) #
iter = 20000


#select height range 
hh = 940 #
hl = 700 #

# Height z
z = np.arange(hh-hl) * 6 #km

# Temp
T = mio.MuramCube(muramdir, iter, 'Temp')
T = T.transpose(1,2,0)
T = T[:,:,hl:hh]

print(T.shape)
nx = T.shape[0]
ny = T.shape[1]
nz = T.shape[2]

#frz model
model3D = frz.atm_model3D(nx,ny,nz)
model3D.z = z[None,None,:]
model3D.tem = T
del(T)
print('T')

# Gas pressure
P = mio.MuramCube(muramdir, iter, 'Pres')
P = P.transpose(1,2,0)
P = P[:,:,hl:hh]
model3D.pg = P
del(P)
print('P')

# horizontal velocity
vx = mio.MuramCube(muramdir, iter, 'vy')
vx = vx.transpose(1,2,0)
vx = vx[:,:,hl:hh]
model3D.vx = vx
del(vx)
vy = mio.MuramCube(muramdir, iter, 'vz')
vy = vy.transpose(1,2,0)
vy = vy[:,:,hl:hh]
model3D.vy = vy
del(vy)

# Vlos
vz = mio.MuramCube(muramdir, iter, 'vx')
vz = vz.transpose(1,2,0)
vz = vz[:,:,hl:hh]
model3D.vz = -vz
del(vz)
print('vz')

# Magnetic field
Bx = mio.MuramCube(muramdir, iter, 'By')
Bx = Bx.transpose(1,2,0)
Bx = Bx[:,:,hl:hh]
Bx = Bx * np.sqrt(4 * np.pi)
model3D.bx = Bx
del(Bx)

By = mio.MuramCube(muramdir, iter, 'Bz')
By = By.transpose(1,2,0)
By = By[:,:,hl:hh]
By = By * np.sqrt(4 * np.pi)
model3D.by = By
del(By)

Bz = mio.MuramCube(muramdir, iter, 'Bx')
Bz = Bz.transpose(1,2,0)
Bz = Bz[:,:,hl:hh]
Bz = Bz * np.sqrt(4 * np.pi)
model3D.bz = Bz
del(Bz)

print('B')

# tau
tau = mio.MuramCube(muramdir, iter, 'tau')
tau = tau.transpose(1,2,0)
tau = tau[:,:,hl:hh]
tau = np.log10(tau)

model3D.tau = tau
del(tau)

#write
model3D.write_model('/dat/xenosh/muram_plage_pore/20000/muram_020000_nx1536_ny1536_nz240_d6.bin')

print('END')
