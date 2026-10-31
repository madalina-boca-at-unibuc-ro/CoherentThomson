"""Builds the report's composite figures from the runs' PNGs (after their post-processing): figs. 1-7 from runs A-D,
fig. 9 from the bunch-length runs (Part II).

Run from this folder: python3 make_figures.py (the runs are in ~/output/progress_report_20261029/).
"""
import os

from PIL import Image, ImageDraw, ImageFont

OUT = 'figures'
RUNS_DIR = '~/output/progress_report_20261029'
RUNS = {'LG00_a0.1': 'runs_A-D/20261029_132340', 'LG22_a0.1': 'runs_A-D/20261029_140030',
        'LG00_a0.2': 'runs_A-D/20261029_143739', 'LG22_a0.2': 'runs_A-D/20261029_151457',
        # Part II: run B's setup on one transverse sample, Gaussian bunch length with q_1 sigma = 0 (ref), 0.5, 1
        'bunch_ref': 'bunch_length/20261030_143857', 'bunch_qs05': 'bunch_length/20261030_155222',
        'bunch_qs1': 'bunch_length/20261030_162911'}
# Panel kind -> PNG path inside <run>/png_folder/radiation/, with {i} the frequency index (N - 1).
KINDS = {'field_F01_radial_phase': 'processed/radial_phase/radiation_field_total_F01_omega{i}.png',
         'field_F01_raw': 'emitted/radiation_field_total_F01_omega{i}.png',
         'energy_flux': 'emitted/observable_energy_flux_omega{i}.png',
         'flux_spin': 'emitted/observable_flux_spin_omega{i}.png',
         'flux_orbital': 'emitted/observable_flux_orbital_omega{i}.png',
         'flux_total': 'emitted/observable_flux_total_omega{i}.png',
         'lg_coefficients': 'lg_decomposition/emitted/lg_coefficients_omega{i}_minus.png',
         'lg_w0_z0_map': 'lg_decomposition/emitted/lg_w0_z0_map_omega{i}_minus.png'}
CASES = [('LG00_a0.1', 'LG(0,0), a0 = 0.1'), ('LG00_a0.2', 'LG(0,0), a0 = 0.2'),
         ('LG22_a0.1', 'LG(2,2), a0 = 0.1'), ('LG22_a0.2', 'LG(2,2), a0 = 0.2')]
FONT = ImageFont.load_default(size=34)
LABEL_H = 50
MAX_WIDTH = 1600  # px: sharp on a phone or a printed page, small enough to email


def grid(panels, labels, n_cols, panel_w, out_name):
    """Pastes `panels` ((case, N, kind) triples) row by row, each with its label above it, at width `panel_w`."""
    ims = []
    for case, N, kind in panels:
        path = os.path.expanduser(f'{RUNS_DIR}/{RUNS[case]}/png_folder/radiation/' + KINDS[kind].format(i=N - 1))
        im = Image.open(path).convert('RGB')
        ims.append(im.resize((panel_w, round(im.height * panel_w / im.width)), Image.LANCZOS))
    n_rows = -(-len(ims) // n_cols)
    row_h = [max(im.height for im in ims[r * n_cols:(r + 1) * n_cols]) + LABEL_H for r in range(n_rows)]
    sheet = Image.new('RGB', (n_cols * panel_w, sum(row_h)), 'white')
    draw = ImageDraw.Draw(sheet)
    for i, (im, label) in enumerate(zip(ims, labels)):
        r, c = divmod(i, n_cols)
        y = sum(row_h[:r])
        draw.text((c * panel_w + 20, y + 8), label, fill='black', font=FONT)
        sheet.paste(im, (c * panel_w, y + LABEL_H))
    if sheet.width > MAX_WIDTH:
        sheet = sheet.resize((MAX_WIDTH, round(sheet.height * MAX_WIDTH / sheet.width)), Image.LANCZOS)
    sheet.save(f'{OUT}/{out_name}', quality=85, optimize=True, progressive=True)


def all_runs(kind, out_name, panel_w):
    panels = [(c, N, kind) for c, _ in CASES for N in (1, 2)]
    labels = [f'{name}, N = {N}' for _, name in CASES for N in (1, 2)]
    grid(panels, labels, 2, panel_w, out_name)


grid([('LG22_a0.1', 2, 'field_F01_raw'), ('LG22_a0.1', 2, 'field_F01_radial_phase')],
     ['exported field (with the e^{ikR} carrier)', 'carrier removed: F e^{-ik(R-Z)}'], 2, 1100,
     'fig1_field_carrier_removal.jpg')
all_runs('field_F01_radial_phase', 'fig2_field_all_runs.jpg', 1100)
all_runs('energy_flux', 'fig3_energy_flux_all_runs.jpg', 800)
grid([('LG22_a0.1', 1, f'flux_{k}') for k in ('spin', 'orbital', 'total')],
     ['spin flux', 'orbital flux', 'total flux'], 3, 800, 'fig4_am_flux_split.jpg')
all_runs('flux_orbital', 'fig5_orbital_flux_all_runs.jpg', 800)
grid([('LG00_a0.1', 2, 'lg_coefficients'), ('LG22_a0.1', 2, 'lg_coefficients')],
     ['LG(0,0), a0 = 0.1, N = 2: close to one mode', 'LG(2,2), a0 = 0.1, N = 2: several radial modes'], 1, 2000,
     'fig6_lg_coefficients.jpg')
all_runs('lg_w0_z0_map', 'fig7_lg_w0_z0_all_runs.jpg', 1100)
grid([(case, N, 'field_F01_radial_phase') for case in ('bunch_ref', 'bunch_qs05', 'bunch_qs1') for N in (1, 2)],
     [f'{name}, N = {N}: q_N sigma = {N * qs:g}'  # ASCII: PIL's default font has no Greek
      for name, qs in (('flat disk', 0), ('sigma = 1.1 nm', 0.5), ('sigma = 2.2 nm', 1)) for N in (1, 2)], 2, 1100,
     'fig9_bunch_length_field.jpg')
