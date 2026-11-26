import sys, numpy as np
sys.argv=[sys.argv[0]]; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import fig_spectrum as fs, plot_spectrum_map as psm
from scipy.signal import find_peaks
for c in ['a0_1p0_Nc8_N1','a0_1p0_Nc8_N2']:
    r=fs.load('10um_'+fs.CASES[c][1]); t=r['theta']/r['theta_star']; w=r['omega']; P=r['photons']
    dA=psm.annulus_areas(r['theta']); fl=np.nansum(P*dA[None,:],axis=0); cf=np.cumsum(fl)/fl.sum()
    q=lambda x:t[np.searchsorted(cf,x)]
    print(f"== {c}: flux in theta/theta*: 10%-50%-90% at {q(.1):.1f}, {q(.5):.1f}, {q(.9):.1f}")
    for tt in [0.5,1,2,3,4,6,8,10]:
        j=np.argmin(abs(t-tt)); s=P[:,j]/np.nanmax(P[:,j])
        pk,_=find_peaks(s,prominence=0.05)
        print(f"  theta={t[j]:4.1f}: {len(pk)} peaks at w/w1 = "+" ".join(f"{w[k]:.3f}" for k in pk))
    dN=r['dN']/r['dN'].max(); pk,pr=find_peaks(dN,prominence=0.01)
    print("  integrated: peaks", " ".join(f"{w[k]:.3f}" for k in pk), " prominences", " ".join(f"{x:.2f}" for x in pr['prominences']))
