import muram as mio
import numpy as np
from astropy.io import fits
from tqdm import tqdm

muramdir = '/dat/milic/MURaM_400G_plage_pore/3D/'
iter = 12000
z = 795
x1 = 850
x2 = 1150
y1 = 650
y2 = 950
slice = np.zeros((64,8,1536,1536), dtype=np.float32)
for i in tqdm(range(64)):
    test = mio.MuramSnap(muramdir,iter)
    slice[i,0,:,:] = np.copy(test.Temp[z,x1:x2,y1:y2])
    slice[i,1,:,:] = np.copy(test.vy[z,x1:x2,y1:y2])
    slice[i,2,:,:] = np.copy(test.vz[z,x1:x2,y1:y2])
    slice[i,3,:,:] = np.copy(test.vx[z,x1:x2,y1:y2])
    slice[i,4,:,:] = np.copy(test.By[z,x1:x2,y1:y2])
    slice[i,5,:,:] = np.copy(test.Bz[z,x1:x2,y1:y2])
    slice[i,6,:,:] = np.copy(test.Bx[z,x1:x2,y1:y2])
    slice[i,7,:,:] = np.copy(test.tau[z,x1:x2,y1:y2])
    iter += 500
    del(test)

hdu = fits.PrimaryHDU(slice)
hdu.writeto('/dat/xenosh/muram_plage_pore/slice/x_'+str(x1)+'_to_'+str(x2)+'_y_'+str(y1)+'_to_'+str(y2)+'_z_'+str(z)+'.fits',overwrite=True)

