import numpy as np 
from astropy.io import fits
import sys
import muram as mio
import firtez_dz as frz

muramdir = '/dat/milic/MURaM_400G_plage_pore/3D/' # can be changed

#iter = int(sys.argv[1]) #
iter = 20500

# x,y range
x1 = 850
x2 = 1150
y1 = 650
y2 = 950

# select height range 
hh = 900 #
hl = 700 #

# Height z
z = np.arange(hh-hl) * 6.5 #km

for i in range(20):
    print('Running iter', iter)
    # Temp  
    T = mio.MuramCube(muramdir, iter, 'Temp')   
    T = T.transpose(1,2,0)
    T = T[x1:x2,y1:y2,hl:hh]

    #print(T.shape)
    nx = T.shape[0]
    ny = T.shape[1]
    nz = T.shape[2]

    #frz model
    model3D = frz.atm_model3D(nx,ny,nz)
    z3d = np.zeros([nx, ny, nz])
    z3d[:,:,:] = z[None,None,:]
    model3D.set_z(z3d)
    model3D.set_tem(T)
    del(T)
    #print('T')

    # Gas pressure
    P = mio.MuramCube(muramdir, iter, 'Pres')
    P = P.transpose(1,2,0)
    P = P[x1:x2,y1:y2,hl:hh]
    model3D.set_pg(P)
    del(P)
    #print('P')

    # rho
    rho = mio.MuramCube(muramdir, iter, 'rho')
    rho = rho.transpose(1,2,0)
    rho = rho[x1:x2,y1:y2,hl:hh]
    model3D.set_rho(rho)
    del(rho)
    #print('rho')

    # horizontal velocity
    vx = mio.MuramCube(muramdir, iter, 'vy')
    vx = vx.transpose(1,2,0)
    vx = vx[x1:x2,y1:y2,hl:hh]
    model3D.set_vx(vx)
    del(vx)
    vy = mio.MuramCube(muramdir, iter, 'vz')
    vy = vy.transpose(1,2,0)
    vy = vy[x1:x2,y1:y2,hl:hh]
    model3D.set_vy(vy)
    del(vy)

    # Vlos
    vz = mio.MuramCube(muramdir, iter, 'vx')
    vz = vz.transpose(1,2,0)
    vz = vz[x1:x2,y1:y2,hl:hh]
    model3D.set_vz(-vz)
    del(vz)
    #print('vz')

    # Magnetic field
    Bx = mio.MuramCube(muramdir, iter, 'By')
    Bx = Bx.transpose(1,2,0)
    Bx = Bx[x1:x2,y1:y2,hl:hh]
    Bx = Bx * np.sqrt(4 * np.pi)
    model3D.set_bx(Bx)
    del(Bx)

    By = mio.MuramCube(muramdir, iter, 'Bz')
    By = By.transpose(1,2,0)
    By = By[x1:x2,y1:y2,hl:hh]
    By = By * np.sqrt(4 * np.pi)
    model3D.set_by(By)
    del(By)

    Bz = mio.MuramCube(muramdir, iter, 'Bx')
    Bz = Bz.transpose(1,2,0)
    Bz = Bz[x1:x2,y1:y2,hl:hh]
    Bz = Bz * np.sqrt(4 * np.pi)
    model3D.set_bz(Bz)
    del(Bz)

    #print('B')

    # tau
    tau = mio.MuramCube(muramdir, iter, 'tau')
    tau = tau.transpose(1,2,0)
    tau = tau[x1:x2,y1:y2,hl:hh]
    tau = np.log10(tau)
    model3D.set_tau(tau)
    del(tau)

    #write
    model3D.write_model('/dat/xenosh/muram_plage_pore/smaller_cubes/'+str(iter)+'/muram_0'+str(iter)+'_nx300_ny300_nz200_d605.bin')
    #
    iter += 1000

print('END')
