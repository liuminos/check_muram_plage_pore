import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as interpolate
from astropy.io import fits
import sys

#input_atmosphere = sys.argv[1]
#output_atmosphere = sys.argv[2]
input_atmosphere = '/dat/xenosh/muram_plage_pore/smaller_cubes/muram_nx300_ny300_nz200_d605_020000_to_039000_addrho.fits'
output_atmosphere = '/dat/xenosh/muram_plage_pore/smaller_cubes/muram_nx300_ny300_nz200_d605_020000_to_039000_tau_addrho.fits'

#NZ = int(sys.argv[3])
#index = int(sys.argv[4]) # which index we use to interpolate the atmosphere 
                         # typically it is 0 or 1 (tau or height) 
NZ = 41-5

atmos_in = fits.open(input_atmosphere)[0].data
print("info::read atmosphere with shape: ", atmos_in.shape)

dims = atmos_in.shape
NZ_old = dims[-1]
NP = dims[1]
NT = dims[0]
NX = dims[2]
NY = dims[3]

print (NZ_old, NP)

atmos_out = np.zeros([NT, NP+1, NX, NY, NZ])

z = np.arange(NZ_old) * 12.0

#start by making an independent variable, which is, of course, h
tau = np.linspace(-2.5,1,NZ)
atmos_out[None,0,None,None,:] = tau[:]


# take log of pressure
atmos_in[:,2,:,:,:] = np.log10(atmos_in[:,2,:,:,:])
#print (atmos_in[:,2,:,:,:])
atmos_in[:,9,:,:,:] = np.log10(atmos_in[:,9,:,:,:])



for t in range (0,NT):
	print(t)
	for p in  range(1,NP):
		for i in range(0,NX):
			for j in range(0,NY):

				#print(atmos_in[t,0,i,j,::-1], atmos_in[t,p,i,j,::-1])

				f = interpolate.interp1d(atmos_in[t,0,i,j,::-1], atmos_in[t,p,i,j,::-1])
				atmos_out[t,p,i,j,:] = f(tau)

	for i in range(0,NX):
		for j in range(0,NY):
			f = interpolate.interp1d(atmos_in[t,0,i,j,::-1], z[::-1])
			atmos_out[t,NP,i,j,:] = f(tau)

atmos_out[:,2,:,:,:] = 10.**atmos_out[:,2,:,:,:]
atmos_out[:,9,:,:,:] = 10.**atmos_out[:,9,:,:,:]
 
atmos_out = atmos_out[:,:,:,:,::-1]

kek = fits.PrimaryHDU(atmos_out)
kek.writeto(output_atmosphere, overwrite=True)
