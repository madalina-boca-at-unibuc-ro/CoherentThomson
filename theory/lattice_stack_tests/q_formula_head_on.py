import sys; sys.path.insert(0,'theory/lattice_stack_tests')
import numpy as np
from lattice_stack_common import *
pz=13.483; b0=pz/np.hypot(1,pz)
for a0,group in ((0.1,'phase1_slope_check_a01'),(0.2,'phase1_slope_check_a02')):
    F={}
    for dz in (0.0,1.4e-2,-1.4e-2):
        run=run_variant(group,f'dz_{dz:+.1e}',{'beam_center_z':f'{dz:.10e} lambda','laser_a0':a0,'beam_lattice_spacing':'10 lambda','circular_detector_N_R':6,'circular_detector_N_phi':4})
        if dz==0: k1,q1=laser_kL_units(run)
        F[dz]=read_field(run)[1]
    print(f'a0={a0}: k1={k1:.6f}')
    for N,i in ((1,0),(2,1)):
        cp,_=fit_scale(F[0.0][i],F[1.4e-2][i]); cm,_=fit_scale(F[0.0][i],F[-1.4e-2][i])
        A=N*(k1+1); C=2*N*k1/(1+b0); D=N*k1+1
        pa=2*np.pi*A*1.4e-2
        # slope from the +/- pair, unwrapped against N(k1+1)
        dd=(np.angle(cp*np.exp(1j*pa))-np.angle(cm*np.exp(-1j*pa)))/2
        s=A-dd/(2*np.pi*1.4e-2)
        print(f'  N={N}: measured s={s:.6f} | K_s.P/K_L.P form {C:.6f} (diff {(s-C)/s:+.1e}) | N(k1+kL) {A:.6f} (diff {(s-A)/s:+.1e}) | Nk1+kL {D:.6f} (diff {(s-D)/s:+.1e})')
