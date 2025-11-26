import numpy as np
import matplotlib.pyplot as plt 
import sys
from astropy.io import fits 
import firtez_dz as frz

#path = sys.argv[1]
start = 0
N = 40

NP = 10 # tau, T, p, vx, vy, vz, Bx, By, Bz, rho

cube = 0

for i in range(start,start+N):
	
	temp = frz.read_model("/dat/xenosh/muram_plage_pore/smaller_cubes/"+str(20000+i*500)+"/muram_0"+str(20000+i*500)+"_nx300_ny300_nz200_d605.bin")

	if (i == start):
		NX = temp.shape[0]
		NY = temp.shape[1]
		NZ = temp.shape[2]
		cube = np.zeros([N, NP, NX, NY, NZ])

	cube[i,0,:,:,:] = temp.tau[:,:,:]
	cube[i,1,:,:,:] = temp.tem[:,:,:]
	cube[i,2,:,:,:] = temp.pg[:,:,:]
	cube[i,3,:,:,:] = temp.vx[:,:,:]
	cube[i,4,:,:,:] = temp.vy[:,:,:]
	cube[i,5,:,:,:] = temp.vz[:,:,:]
	cube[i,6,:,:,:] = temp.bx[:,:,:]
	cube[i,7,:,:,:] = temp.by[:,:,:]
	cube[i,8,:,:,:] = temp.bz[:,:,:]
	cube[i,9,:,:,:] = temp.rho[:,:,:]

	print("info::time step", i, " done")

kek = fits.PrimaryHDU(cube)
kek.writeto('/dat/xenosh/muram_plage_pore/smaller_cubes/muram_nx300_ny300_nz200_d605_020000_to_039000_addrho.fits', overwrite=True)
	

	
