import h5py
import numpy as np
from scipy.constants import epsilon_0, e, c
from rich.progress import track
import matplotlib.pyplot as plt
import scipy.io
# from cmcrameri import cm

#plt.rcParams["font.family"] = "Nimbus Roman"
plt.rcParams["font.size"] = 20
plt.rcParams["lines.linewidth"] = 2
plt.rcParams["figure.figsize"] = 9, 10

sigmaz = 9e-2
zz = np.linspace(-0.8, 0.4, 1000)

bs = scipy.io.loadmat('LHC_chm_ver.mat')
bsx = bs['Vx'][0]
bsy = bs['Vy'][0]

bsx = np.array(list(bsx) + [bsx[0]])
bsy = np.array(list(bsy) + [bsy[0]])

pinch = 'Pinch.h5'

h5p = h5py.File(pinch,'r')

xg = h5p['grid/xg'][()]
yg = h5p['grid/yg'][()]
zg = h5p['grid/zg'][()]
dh = xg[1] - xg[0]

XX, YY = np.meshgrid(xg, yg, indexing='ij')

# sigx = h5p['grid/sigx'][()]*1e3
# sigy = h5p['grid/sigy'][()]*1e3
# xc = h5p['grid/xc'][()]*1e3
# yc = h5p['grid/yc'][()]*1e3

blen = 0.09
rms_tt = blen/c*1e9

xr = np.cos(np.linspace(0,2*np.pi, 1000))
yr = np.sin(np.linspace(0,2*np.pi, 1000))

vmin=0.1
vmax=2e10

phi = h5p[f'slices/slice0/phi'][()]
rho0 = np.zeros_like(phi)
rho0[1:-1, 1:-1] = -epsilon_0*(phi[:-2,1:-1] + phi[2:,1:-1] + phi[1:-1,2:] + phi[1:-1,:-2] - 4 * phi[1:-1, 1:-1])/dh**2

fig = plt.figure(1)
for ii in track(range(len(zg))):
# for ii in [0,1]:
    fig.clear()
    fig, (ax0, ax) = plt.subplots(2,1, height_ratios=[1,2])
    ax0.plot(zz, np.exp(-zz**2/2./sigmaz**2), 'k')
    #ax0 = fig.add_subplot(211)
    #ax = fig.add_subplot(212)
    phi = h5p[f'slices/slice{ii}/phi'][()]
    rho = np.zeros_like(phi)
    rho[1:-1, 1:-1] = -epsilon_0*(phi[:-2,1:-1] + phi[2:,1:-1] + phi[1:-1,2:] + phi[1:-1,:-2] - 4 * phi[1:-1, 1:-1])/dh**2
    # rho = rho - rho0
    #rho[rho >= 0] = 10**vmin * (-e)
    rho = rho/(-e)
    #rho = np.log10(rho/(-e))
    rho[rho < vmin] = np.nan
    print(f"min: {np.min(phi)}, std: {np.std(phi)},  mean: {np.mean(phi)}, max: {np.max(phi)}")
    mpl = ax.pcolormesh(XX*1e3, YY*1e3, rho, vmin=vmin, vmax=vmax, shading="nearest", cmap='viridis')
    ax.plot(bsx*1e3, bsy*1e3, 'k-', lw=5)
    ax0.axvline(zg[ii], ls='--', c='k', lw=2)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax0.set_xlim(-0.8, 0.4)
    ax0.set_ylabel('Long. bunch profile')
    ax0.set_xlabel('z [m]')

    # cb = plt.colorbar(mpl)
    # cb.set_label("electron density, log ρ [e]")

    tt = zg[ii]/c*1e9
    # ax.set_title(f"t = {(tt):.3f} ns, {((tt)/rms_tt):.2f}σ")
    fig.tight_layout()
    fig.savefig(f"frames/frame_{ii:d}.png")

