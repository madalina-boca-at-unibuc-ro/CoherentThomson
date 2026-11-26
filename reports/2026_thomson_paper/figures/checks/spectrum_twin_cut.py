import sys, numpy as np
sys.argv=[sys.argv[0]]; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import fig_spectrum as fs, plot_spectrum_map as psm
for c in ['a0_1p0_Nc8_N1','a0_1p0_Nc8_N2','a0_0p5_Nc30_N2']:
    s=fs.CASES[c][1]; a,b=fs.load('10um_'+s),fs.load('800nm_'+s)
    ta,tb=a['theta']/a['theta_star'],b['theta']/b['theta_star']
    w=a['omega']; out=[]
    for cut in [3,6,12,1e9]:
        def spec(r,t):
            dA=psm.annulus_areas(r['theta']); m=(t<=cut)
            return np.nansum((r['photons']*dA[None,:])[:,m],axis=1)
        na,nb=spec(a,ta),spec(b,tb); na/=na.max(); nb/=nb.max()
        mean=lambda n:(n*w).sum()/n.sum()
        out.append(f"<={cut if cut<1e8 else 'all':>3}: {np.abs(na-nb).max():.1e} ({(mean(nb)-mean(na))/mean(na):+.1e})")
    print(f"{c:16s} theta_max/theta*={ta.max():.1f}  "+"  ".join(out))
