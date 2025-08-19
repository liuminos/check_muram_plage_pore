import numpy as np 
from astropy.io import fits
import sys
import muram as mio
import firtez_dz as frz

# bin atmos
atm = frz.read_model('/dat/xenosh/muram_plage_pore/20000/muram_020000_nx1536_ny1536_nz240_d6.bin')
# bin
b = 8
nx = atm.shape[0] // b
ny = atm.shape[1] // b
nz = atm.shape[2] // 2
cb = frz.atm_model3D(nx,ny,nz)
cb.tem[:,:,:] = atm.tem[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
cb.pg[:,:,:] = atm.pg[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
cb.vz[:,:,:] = atm.vz[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
cb.bx[:,:,:] = atm.bx[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
cb.by[:,:,:] = atm.by[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
cb.bz[:,:,:] = atm.bz[:,:,:].reshape(nx,8,ny,8,nz,2).mean(axis=(1,3,5))
print('DONE')
z = np.arange(nz) * 12
cb.z[None,None,:] = z

cb.write_model('/dat/xenosh/muram_plage_pore/20000/muram_020000_nx192_ny192_nz120_dz12.bin')
'''
# bin stokes
stokes = frz.read_profile('/dat/xenosh/muram_plage_pore/20000/stokes.bin')
# bin
nx = stokes.shape[1] // 8
ny = stokes.shape[2] // 8
nw = stokes.shape[0]
cb = frz.stk_profile3D(nw,nx,ny)
cb.stki[:,:,:] = stokes.stki[:,:,:].reshape(nw,nx,8,ny,8).mean(axis=(2,4))
cb.stkq[:,:,:] = stokes.stkq[:,:,:].reshape(nw,nx,8,ny,8).mean(axis=(2,4))
cb.stku[:,:,:] = stokes.stku[:,:,:].reshape(nw,nx,8,ny,8).mean(axis=(2,4))
cb.stkv[:,:,:] = stokes.stkv[:,:,:].reshape(nw,nx,8,ny,8).mean(axis=(2,4))
print('DONE')
cb.wave = stokes.wave
# add noise
cb.stki = cb.stki + np.random.normal(loc=0, scale=1e-3, size=cb.shape)
cb.stkq = cb.stkq + np.random.normal(loc=0, scale=1e-3, size=cb.shape)
cb.stku = cb.stku + np.random.normal(loc=0, scale=1e-3, size=cb.shape)
cb.stkv = cb.stkv + np.random.normal(loc=0, scale=1e-3, size=cb.shape)

cb.write_profile('/dat/xenosh/muram_plage_pore/20000/stokes_bin_noise0.001.bin')
'''