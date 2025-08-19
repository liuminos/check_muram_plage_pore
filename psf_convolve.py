import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.convolution as apconv
import sys
import firtez_dz as frz

stokes = frz.read_profile('/dat/xenosh/muram_plage_pore/20000/stokes.bin')
print(stokes.shape)

#prepare the PSF
#D = 0.20 #m
D = 1.0
llambda = 500e-9
pix_size = 6e3 #m
pixel_scale_arcsec = pix_size / 725e3
diff_limit_arcsec = 1.22 * llambda / D * 206265
diff_limit_pixels = diff_limit_arcsec / pixel_scale_arcsec
PSF = apconv.AiryDisk2DKernel(diff_limit_pixels, mode='oversample')
print("info::PSF prepared")

#convolve PSF with data and bin
nl, nx, ny = stokes.shape
convolve = frz.stk_profile3D(nl, nx, ny)
cb = frz.stk_profile3D(nl, int(nx/8), int(ny/8))

for i in range(int(nl)):
    convolve.stki[i,:,:] = apconv.convolve_fft(stokes.stki[i,:,:],PSF,boundary='wrap',normalize_kernel=True)
    convolve.stkq[i,:,:] = apconv.convolve_fft(stokes.stkq[i,:,:],PSF,boundary='wrap',normalize_kernel=True)
    convolve.stku[i,:,:] = apconv.convolve_fft(stokes.stku[i,:,:],PSF,boundary='wrap',normalize_kernel=True)
    convolve.stkv[i,:,:] = apconv.convolve_fft(stokes.stkv[i,:,:],PSF,boundary='wrap',normalize_kernel=True)

    cb.stki[i,:,:] = np.sum(convolve.stki[i,:,:].reshape(int(nx/8),8,int(ny/8),8),axis=(1,3))/(8*8)
    cb.stkq[i,:,:] = np.sum(convolve.stkq[i,:,:].reshape(int(nx/8),8,int(ny/8),8),axis=(1,3))/(8*8)
    cb.stku[i,:,:] = np.sum(convolve.stku[i,:,:].reshape(int(nx/8),8,int(ny/8),8),axis=(1,3))/(8*8)
    cb.stkv[i,:,:] = np.sum(convolve.stkv[i,:,:].reshape(int(nx/8),8,int(ny/8),8),axis=(1,3))/(8*8)

cb.wave = stokes.wave

# add noise
cb.stki = cb.stki + np.random.normal(loc=0, scale=1e-2, size=cb.shape)
cb.stkq = cb.stkq + np.random.normal(loc=0, scale=1e-2, size=cb.shape)
cb.stku = cb.stku + np.random.normal(loc=0, scale=1e-2, size=cb.shape)
cb.stkv = cb.stkv + np.random.normal(loc=0, scale=1e-2, size=cb.shape)
cb.write_profile('/dat/xenosh/muram_plage_pore/20000/psf_bin_noise0.01.bin')
