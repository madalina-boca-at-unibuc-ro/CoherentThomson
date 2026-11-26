import sys, numpy as np
sys.argv=[sys.argv[0]]
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import fig_spectrum as fs
print(f"{'case':16s} {'max|dN|':>8s} {'rms dN':>8s} {'d<w>/w':>9s} {'dwidth':>8s} {'map max':>8s} {'map rms':>8s} {'D10/D800 N ratio-1':>8s}")
for c in ['a0_0p1_N1','a0_0p1_N2','a0_1p0_Nc8_N1','a0_1p0_Nc8_N2','a0_0p5_Nc30_N1','a0_0p5_Nc30_N2']:
    s=fs.CASES[c][1]; a,b=fs.load('10um_'+s),fs.load('800nm_'+s)
    assert np.allclose(a['omega'],b['omega']), 'omega grids differ'
    ta,tb=a['theta']/a['theta_star'],b['theta']/b['theta_star']
    gr=np.nanmax(np.abs(tb[1:]/ta[1:]-1))
    na,nb=a['dN']/a['dN'].max(),b['dN']/b['dN'].max()
    w=a['omega']
    mean=lambda n:(n*w).sum()/n.sum()
    sd=lambda n:np.sqrt((n*(w-mean(n))**2).sum()/n.sum())
    mb_i=np.array([np.interp(ta,tb,row) for row in b['photons']])
    ma,mb=a['photons']/np.nanmax(a['photons']),mb_i/np.nanmax(mb_i)
    d=np.abs(ma-mb)
    print(f"{c:16s} {np.abs(na-nb).max():8.1e} {np.sqrt(((na-nb)**2).mean()):8.1e} {(mean(nb)-mean(na))/mean(na):+9.1e} {sd(nb)/sd(na)-1:+8.1e} {np.nanmax(d):8.1e} {np.sqrt(np.nanmean(d**2)):8.1e}  grid {gr:.1e}")
