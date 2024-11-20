import json
import scipy.io
from scipy.constants import c,e
import matplotlib.pyplot as plt
import numpy as np
# from cmcrameri import cm
from rich.progress import Progress
import glob
import os

import sys
sys.path.append("../../tools")
import filemanager as exfm
from scipy.constants import e

from PyECLOUD.buildup_simulation import BuildupSimulation
import argparse

plt.rcParams["font.size"] = 20
plt.rcParams["lines.linewidth"] = 2
plt.rcParams["figure.figsize"] = 12, 6


radius = 23
theta = np.linspace(0,2*np.pi, 1000)
bs = scipy.io.loadmat('LHC_chm_ver.mat')
bsx = bs['Vx'][0]
bsy = bs['Vy'][0]

vmin=0.1
vmax=7e8

mod_factor = 4


mp_states = []
mp_state_files = ["../MP_state_0"]
for name in mp_state_files:
    mp_system = scipy.io.loadmat(name)
    mp_system['nel_mp'] = mp_system['nel_mp']/len(mp_state_files)
    mp_states.append(mp_system)

t_min = 0# 2.5e-9 - 0.4 / c
t_max = 50e-9 #2.5e-9 + 0.8 / c
t_end_sim = t_max + 0.02 / c

total_N_particles = sum([mp_system["N_mp"][0][0] for mp_system in mp_states])
total_N_particles *= 2 #put more particles to max allowed
total_N_particles += 10000 #put more particles to max allowed

sim = BuildupSimulation(extract_sey=False, N_mp_max=total_N_particles,
                        init_unif_edens_flag=0, Dt_sc=None)


sim.cloud_list[0].MP_e.init_from_dict(mp_states[0])
sim.cloud_list[0].MP_e.nel_mp_ref *= len(mp_states)
prev_mp = sim.cloud_list[0].MP_e.N_mp
for mp_system in mp_states[1:]:
    sim.cloud_list[0].MP_e.add_from_file(mp_system)
    now_mp = sim.cloud_list[0].MP_e.N_mp
    print(now_mp, mp_system["N_mp"][0][0], now_mp - prev_mp)
    prev_mp = now_mp

N_mp = sim.cloud_list[0].MP_e.N_mp
max_nel_mp = np.max(sim.cloud_list[0].MP_e.nel_mp[:N_mp])
sim.cloud_list[0].MP_e.nel_mp_split = 3*(max_nel_mp)

print(sim.cloud_list[0].MP_e.nel_mp_split)

print("Start timestep iter")

## simulation
def time_step(sim, t_end_sim=None):
    beamtim = sim.beamtim
    if t_end_sim is not None and beamtim.tt_curr is not None:
        if beamtim.tt_curr >= t_end_sim:
            print("Reached user defined t_end_sim --> Ending simulation")
            return 1
 
    beamtim.next_time_step()
 
    if sim.flag_presence_sec_beams:
        for sec_beam in sim.sec_beams_list:
            sec_beam.next_time_step()
 
    sim.sim_time_step(force_reinterp_fields_at_substeps=True, skip_MP_cleaning=True, skip_MP_regen=True)
 
    if beamtim.flag_new_bunch_pass:
        print(
            "**** Done pass_numb = %d/%d\n"
            % (beamtim.pass_numb, beamtim.N_pass_tot)
        )
    return 0



xg = sim.spacech_ele.xg
yg = sim.spacech_ele.yg

XX, YY = np.meshgrid(xg, yg, indexing='ij')

zg = []

ii = 0

time_array = []
number_of_electrons = []

with Progress() as progress:
    task = progress.add_task("PyECLOUD tracking", total=t_end_sim)

    while not time_step(sim, t_end_sim=t_end_sim):
#        print(sim.cloud_list[0].MP_e.N_mp/total_N_particles)
        ii += 1
        tt = sim.beamtim.tt_curr
        progress.update(task, completed=tt)
        if tt > t_min and tt < t_max:
            if ii%mod_factor == 0.:
                saver = sim.cloud_list[0].pyeclsaver
                ilast = saver.i_last_save
                fig, (ax1, ax2) = plt.subplots(1,2)
                rho = sim.spacech_ele.rho / (-e)
                rho[rho < vmin] = np.nan
                ax1.pcolormesh(XX*1e3, YY*1e3, rho, vmin=vmin, vmax=vmax, shading="nearest", cmap='viridis')
                ax1.plot(radius*np.cos(theta), radius*np.sin(theta), 'k-', lw=5)
                #ax1.plot(bsx*1e3, bsy*1e3, 'k-', lw=5)
                ax2.plot(saver.t[:ilast]/1e-9, saver.Nel_timep[:ilast]/1e5, 'b-')
                ax2.plot(saver.t[:ilast]/1e-9, saver.lam_t_array[:ilast]/1e12*6, 'k-')

                ax1.set_xlabel("x [mm]")
                ax1.set_ylabel("y [mm]")
                ax1.set_aspect('equal')
                
                ax2.set_xlabel("Time [ns]")
                ax2.set_ylabel("Number of e$^-$ [10$^{5}$ / m]")
                ax2.set_xlim(0, 50)
                ax2.set_ylim(2, 6)
                fig.tight_layout()
                fig.savefig(f"frames/frame_{ilast:d}.png")
                plt.close(fig)
