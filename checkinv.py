import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import matplotlib.pyplot as plt
import firtez_dz as frz
from astropy.io import fits
from matplotlib.colors import CenteredNorm, TwoSlopeNorm
plt.rc('text', usetex = False)
#read atmosphere
atm = frz.read_model('/dat/xenosh/MiHi_Fe_I_plage/Time_10s/0/c1/MHS/out_out_out_out_out_out_hsra_nx116_ny129_nz64_dz12_ME_Bh.bin')


# Create main figure
fig, ax = plt.subplots()
im=ax.imshow(atm.bz[:,:,19].T, cmap="PuOr", origin="lower",vmin=-500,vmax=500) # !!! z not determined
fig.colorbar(im,ax=ax)
#ax.set_title("")

# Create a second figure for detail plot
fig_detail, ax_detail = plt.subplots(2,2, sharex=True)

#boundary
tau_left = -4
tau_right = 1

def onclick(event):
    if event.inaxes is not None:
        ix, iy = int(event.xdata), int(event.ydata)
        #print(f"Clicked at x={ix}, y={iy}")

        # Clear and redraw the detail figure
        for a in ax_detail.flatten():
            a.clear()
        ax_detail[0,0].plot(atm.tau[ix,iy,:],atm.tem[ix,iy,:])
        ax_detail[0,0].set_xlim(tau_left, tau_right)
        ax_detail[0,0].set_title('T [K]')
        ax_detail[0,1].plot(atm.tau[ix,iy,:],atm.vz[ix,iy,:]/1e5)
        ax_detail[0,1].set_title('Vz [km/s]')
        ax_detail[1,0].plot(atm.tau[ix,iy,:],atm.bz[ix,iy,:])
        ax_detail[1,0].set_title('Bz [G]')
        ax_detail[1,1].plot(atm.tau[ix,iy,:],np.sqrt(atm.bx[ix,iy,:]**2+atm.by[ix,iy,:]**2))
        ax_detail[1,1].set_title('Bh [G]')
        fig_detail.suptitle(f"({ix}, {iy})")

        # Refresh the figure
        fig_detail.canvas.draw()

# Connect click event
cid = fig.canvas.mpl_connect("button_press_event", onclick)
plt.tight_layout()
plt.show(block=True)