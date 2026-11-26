import sys, re, numpy as np
sys.argv=[sys.argv[0]]; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
import fig_spectrum as fs, plot_spectrum_map as psm
def info(run):
    log=open(run+'run_log.txt').read(); cfg=open(run+'config.cfg').read()
    w0=float(re.search(r'LG w0\s+([0-9.e+-]+) a\.u\.',log).group(1))
    zw0=float(re.search(r'^rectangular_detector_distance\s+([0-9.e+-]+)\s+w0',cfg,re.M).group(1))
    a0=float(re.search(r'^laser_a0\s+([0-9.e+-]+)',cfg,re.M).group(1))
    r1=float(re.search(r'fundamental / laser omega\s+([0-9.e+-]+)',log).group(1))
    ne=float(re.search(r'represented charge\s+([0-9.e+-]+)',log).group(1))
    return zw0*w0, np.sqrt(r1*(1+a0**2/2)), ne
print(f"{'case':16s} {'N':>2s} {'D10/D800':>8s} {'expect':>8s} {'N_tot ratio':>11s} {'/expect-1':>9s} {'peak ratio':>10s} {'/expect-1':>9s}")
for c in ['a0_0p1_N1','a0_0p1_N2','a0_1p0_Nc8_N1','a0_1p0_Nc8_N2','a0_0p5_Nc30_N1','a0_0p5_Nc30_N2']:
    s=fs.CASES[c][1]; N=int(s[-1]); a,b=fs.load('10um_'+s),fs.load('800nm_'+s)
    (Za,Da,na),(Zb,Db,nb)=info(a['run']),info(b['run']); assert na==nb
    # dN/dx per solid angle: dPz (per screen area, per absolute omega) * Z^2 * dOmega / x ; x = omega/omega_1 (omega_1 cancels)
    da=a['dN']*Za**2; db=b['dN']*Zb**2
    x=a['omega']; dx=np.gradient(x)
    Ta,Tb=(da*dx).sum(),(db*dx).sum(); exp=(Db/Da)**(2*N)
    print(f"{c:16s} {N:2d} {Da/Db:8.4f} {exp:8.5f} {Ta/Tb:11.6f} {Ta/Tb/exp-1:+9.1e} {da.max()/db.max():10.6f} {da.max()/db.max()/exp-1:+9.1e}")
