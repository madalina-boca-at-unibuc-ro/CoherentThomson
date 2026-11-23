"""trajectory_NT convergence vs gamma (theory/future_dev/numerical_parameters_and_cost.md, section 1).

Generates configs from one_electron_on_axis.cfg (16 electrons on a 1 w0 disk, Gaussian pulse sigma = 25 cycles,
cutoff 4), runs bin/coherent_thomson_solver for 4 momenta x 6 trajectory_NT x 2 formulas, and prints the relative L2
difference of the total F (6 components, N = 1 and 2) against NT = 1600 and between the two formulas. Configs and
the run index go to ~/output_tests/nt_scan/, runs to ~/output_tests/. About 3 minutes on 16 threads.

    python3 theory/future_dev/tests/nt_scan.py run    # run + analyze
    python3 theory/future_dev/tests/nt_scan.py        # analyze the existing index (runs if there is none)
"""
import os, re, subprocess, sys, time, glob, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
S = os.path.expanduser('~/output_tests/nt_scan')
os.makedirs(S, exist_ok=True)
sys.path.insert(0, os.path.join(REPO, 'py_scripts', 'radiation'))
sys.argv = sys.argv[:1]
from plot_point_spectrum import read_radiation_data, get_point_component

BASE = os.path.join(HERE, 'one_electron_on_axis.cfg')
OUT = os.path.expanduser('~/output_tests')
PZ = [0.0, 2.0, 13.483, 50.0]
NTS = [25, 50, 100, 200, 400, 1600]
FORMULAS = ['long_distance_simplified_approx', 'long_distance_direct_approx']
A2 = 0.1**2 / 2
Z_LAMBDA = 4.0e6
W0 = 25.0


def set_keys(text, kv):
  for k, v in kv.items():
    text, n = re.subn(rf'^({re.escape(k)}[ \t]+)[^#\n]*', rf'\g<1>{v}  ', text, flags=re.M)
    assert n == 1, k
  return text


def ratio(pz):
  p0 = np.sqrt(1 + pz**2)
  return (p0 + pz) / (p0 - pz + A2 / (p0 + pz))


def run_all():
  base = open(BASE).read()
  index = {}
  for pz in PZ:
    r_max = 3 * Z_LAMBDA / (np.pi * W0 * ratio(pz))
    for f in FORMULAS:
      for nt in NTS:
        name = f'nt_pz{pz}_{f}_{nt}'
        cfg = set_keys(base, {
            'beam_particle_count': '16', 'beam_cylinder_radius': '1.0 w0', 'num_threads': '16',
            'circular_detector_N_R': '16', 'circular_detector_N_phi': '17',
            'circular_detector_R_max': f'{r_max:.6g} lambda', 'circular_detector_distance': f'{Z_LAMBDA:.6g} lambda',
            'average_pz': f'{pz} mc', 'trajectory_NT': str(nt), 'radiation_formula': f,
            'laser_flat_duration': '0.0 cycles_adim', 'laser_wing_sigma': '25.0 cycles_adim',
            'laser_wing_sigma_cutoff': '4'})
        path = os.path.join(S, name + '.cfg')
        open(path, 'w').write(cfg)
        before = set(glob.glob(OUT + '/*/'))
        t = time.time()
        subprocess.run([REPO + '/bin/coherent_thomson_solver', path], stdout=open(path + '.log', 'w'),
                       stderr=subprocess.STDOUT, check=True)
        time.sleep(1.1)  # run folders are timestamped to the second
        new = sorted(set(glob.glob(OUT + '/*/')) - before)
        assert len(new) == 1, new
        index[name] = new[0]
        print(f'{name}: {time.time() - t:.1f} s', flush=True)
  json.dump(index, open(os.path.join(S, 'nt_index.json'), 'w'), indent=1)


def fields(path):
  d = read_radiation_data(path + 'radiation_field.dat')
  out = {}
  for N in (1, 2):
    sub = d[np.isclose(d['omega'], N)]
    out[N] = np.concatenate([get_point_component(sub, 'total', mu, nu)
                             for mu, nu in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
  return out


def rel(a, b):
  return np.linalg.norm(a - b) / np.linalg.norm(b)


def analyze():
  index = json.load(open(os.path.join(S, 'nt_index.json')))
  F = {k: fields(v) for k, v in index.items()}
  for pz in PZ:
    print(f'\naverage_pz = {pz} mc  (omega_1/omega_L = {ratio(pz):.4g})')
    print('   NT | N=1: simpl vs ref  direct vs ref  simpl vs direct | N=2: simpl vs ref  direct vs ref  simpl vs direct')
    for nt in NTS:
      row = f'{nt:5d} |'
      for N in (1, 2):
        s = F[f'nt_pz{pz}_{FORMULAS[0]}_{nt}'][N]
        dr = F[f'nt_pz{pz}_{FORMULAS[1]}_{nt}'][N]
        s_ref = F[f'nt_pz{pz}_{FORMULAS[0]}_{NTS[-1]}'][N]
        d_ref = F[f'nt_pz{pz}_{FORMULAS[1]}_{NTS[-1]}'][N]
        row += f'      {rel(s, s_ref):9.2e}     {rel(dr, d_ref):9.2e}       {rel(s, dr):9.2e} |'
      print(row)


if __name__ == '__main__':
  if 'run' in sys.argv[1:] or not os.path.exists(os.path.join(S, 'nt_index.json')):
    run_all()
  analyze()
