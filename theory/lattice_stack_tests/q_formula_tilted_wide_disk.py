import sys; sys.path.insert(0,'theory/lattice_stack_tests')
import numpy as np
from lattice_stack_common import *
pz=13.483; E=np.hypot(1,pz); th=0.9*np.pi; a0=0.2; d=1.4e-2
n=np.array([np.sin(th),np.cos(th)])
base={'laser_a0':a0,'laser_direction_theta':'0.9 pi','circular_detector_N_R':1,'circular_detector_N_phi':1,'N_harmonics':1,'beam_lattice_spacing':'10 lambda'}
for label,extra in (('LG w0=25',{}),('plane wave',{'laser_type':'plane_wave','beam_cylinder_radius':'46.5 lambda'})):
    F={}
    for name,kv in {'0':{}, 'zp':{'beam_center_z':f'{d} lambda'}, 'zm':{'beam_center_z':f'{-d} lambda'}, 'xp':{'beam_center_x':f'{d} lambda'}, 'xm':{'beam_center_x':f'{-d} lambda'}}.items():
        r=run_variant('phase1_slope_check_tilt_point',f"{label.split()[0]}_{name}",dict(base,**extra,**kv)); F[name]=read_field(r)[1][0]
        if name=='0': k1=float(read_log_value(r,r"k_s\s+\(.*?\[\(\S+, \S+, ([-0-9.e+]+)\) k_L\]"))
    ratio=k1*(E-pz)/(E-np.cos(th)*pz)
    print(f'{label}: k1={k1:.6f} ratio={ratio:.6f} |F|={np.linalg.norm(F["0"]):.3e}')
    A=np.array([0,k1])-n; C=np.array([0,k1])-ratio*n
    for ax,(p,m),j in (('x',('xp','xm'),0),('z',('zp','zm'),1)):
        cp,e1=fit_scale(F['0'],F[p]); cm,e2=fit_scale(F['0'],F[m])
        ph=2*np.pi*C[j]*d
        dd=(np.angle(cp*np.exp(1j*ph))-np.angle(cm*np.exp(-1j*ph)))/2
        s=C[j]-dd/(2*np.pi*d)
        print(f'  {ax}: measured {s:+.6f} | dressed {C[j]:+.6f} (diff {s-C[j]:+.1e}) | k_s-k_L {A[j]:+.6f} (diff {s-A[j]:+.1e}) | eta {max(e1,e2):.1e}',flush=True)
