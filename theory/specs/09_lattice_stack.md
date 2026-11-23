# Longitudinal lattice stacks and ponderomotive microbunching

Implementation brief for the CoherentThomson coding agent. Work in the existing C++ project at `~/Dropbox/work/bin/c++/codes/CoherentThomson`. Read the scientific report in `theory/scientific_report/progress_report.md`, and `theory/thomson_form_factor_diagnostic.md`. The bunching model comes from Schaap et al., *Phys. Rev. Accel. Beams* **26**, 074401 (2023), especially Eqs. (7), (9), (12), (15), and (24). The local PDF is in `~/Zotero/storage/RQ9MJ5TU/`.

Implement one phase at a time. At each checkpoint, leave a runnable configuration and the exact commands for reproducing the checks, report the measured quantities and any discrepancies, and stop so the user can run and review that phase before authorizing the next one.

## Shared conventions

- Use the report’s head-on geometry for validation. The beam and detector point along lab $+z$, and the scattering laser propagates along $-z$. Start with the Gaussian laser and the Thomson fundamental. Keep the laser, screen, transverse lattice, and total represented charge fixed across paired comparisons.
- Distinguish the scattered wavelength $\lambda_s=2\pi/k_s$, the phase-matching period $\Lambda=2\pi/|q_z|$, and the modulation period $\Lambda_{\rm mod}$. The existing initial-position diagnostic uses $\mathbf q=\mathbf k_s-\mathbf k_L$, giving $q_z=k_s+k_L$ in the baseline geometry. Read the actual solver frequency from the chosen configuration or run log. The report gives $\lambda_s$ of roughly 13.8–14.0 nm. The paper’s example uses a 272 nm modulation and targets about 13.6 nm at the twentieth density harmonic. The density harmonic index is independent of the radiation harmonic index.
- Verify the effective longitudinal phase slope in the radiation solver by translating a single sheet over several small, known offsets and fitting the complex field ratio at the screen center. Account for the Fourier sign convention, check the fitted magnitude for changes, and compare the slope with the diagnostic $q_z$. The latter is an initial-position diagnostic, so use the measured slope for radiation phase-matching if the difference matters over the stack length. Log both values and explain which one sets the modulation period.
- Start at $d=1$ nm and refine to $d=0.5$ nm at the **same physical length and total charge**. This gives roughly 14 and 28 sheets per $\Lambda$. The first sampling replica is at $2\pi/d$, well beyond the target $q_z$ at these spacings. Check convergence of the complex bunching factor, not just the plotted density. A rough bin-centering estimate $1-\operatorname{sinc}(q_zd/2)$ gives 0.8% at 1 nm and 0.2% at 0.5 nm. Treat these as estimates and measure the error for the actual density profile. Keep the transverse lattice spacing and its independent diffraction constraints as in the report.
- Use zero momentum spread for the initial validation. Preserve legacy behavior and output normalization for all existing configurations. Use the existing beam frame, unit conversion, and lattice geometry. Resolve and log the physical spacing, number of sheets, represented charge, and trajectory count. Estimate cost before full-screen runs. At 1 nm, even one 280 nm period on a 1,076-site transverse lattice costs about 300,000 trajectories. Use cheap geometry and form-factor checks first, then a detector-center run, then a modest screen if warranted.

## Phase 1. Uniform stack of parallel polar lattices

Add an opt-in stack configuration to the existing polar-lattice beam generator. Give it an explicit sheet count $M$ and spacing $d$ with units. Duplicate the same transverse lattice at beam-frame coordinates
$$
z_j=\left(j-\frac{M-1}{2}\right)d,\qquad j=0,\ldots,M-1.
$$
Apply the existing beam rotation and offset convention. For a sheet interpreted as the center of a longitudinal bin, the stack covers $H=Md$ and its first-to-last center span is $(M-1)d$. Reject conflicting settings that would smear or ambiguously position the sheets. In particular, document how the stack interacts with cylinder height and Gaussian longitudinal smearing.

Give each sheet equal represented charge at fixed total charge. Carry each site's real, nonnegative weight into the **complex field sum**, including all long-range, short-range, and boundary contributions. Keep unit weights in the legacy path. Inspect `src/core/particle/electron_factory.cpp`, `src/core/simulation/simulation.cpp`, and `src/core/radiation/radiation.cpp` for the narrowest correct weighting point. Apply weights exactly once and keep the thread-local reduction. Update `src/core/simulation/form_factor.cpp` and the associated logging to use the same weights for phase averaging and mean momentum. Distinguish simulated trajectory count from total represented charge.

The uniform stack has a normalized longitudinal amplitude factor
$$
B_M(q)=\frac{1}{M}\sum_j e^{iqz_j}=\frac{\sin(Mqd/2)}{M\sin(qd/2)}
$$
with removable singularities evaluated by the finite sum or its limit. In the continuous uniform limit it approaches $\operatorname{sinc}(qH/2)$, with $\operatorname{sinc}(x)=\sin x/x$. In the regime where translating a sheet changes mainly its phase, the stack-to-single-sheet field ratio should follow this factor up to the solver's Fourier sign convention. Power follows $|B_M|^2$.

**Checkpoint.** Show that $M=1$ reproduces the original flat lattice. At fixed charge, run stacks with increasing $H$ from much less than $\Lambda$ through the first minimum near $H=\Lambda$ and into the next lobe. Compare the calculated complex form factor and the center field with the discrete prediction. Check the artificial rephasing at $qd=2\pi$ in the diagnostic only, as a control. Include the sheet-translation slope measurement and a spacing-refinement comparison at fixed $H$. For near-zero predictions, report absolute errors instead of relative errors.

## Phase 2. Derive and validate the density model

Prepare a short derivation and a standalone numerical reference or diagnostic before enabling modulation in production runs. Use the ideal longitudinal, impulsive, cold-beam model in Schaap et al. The beat wave imprints a sinusoidal energy kick. During the drift, faster electrons overtake slower ones. In the paper's phase convention the map is
$$
\phi=\phi_0+\zeta\sin\phi_0,\qquad \zeta=L/f,\qquad
f=\frac{c\gamma^3\beta^3}{|\omega_1-\omega_2|\,|\Delta\gamma|},
$$
where $\omega_1,\omega_2$ are the **buncher** laser frequencies and $\Delta\gamma$ is the energy-modulation amplitude in units of $mc^2$.

The periodic density is the pushforward of uniform initial phase. Away from folds, its value relative to the mean is
$$
\frac{\rho(\phi)}{\bar\rho}=\sum_{\phi_0:\,\phi_0+\zeta\sin\phi_0=\phi\;({\rm mod}\;2\pi)}\frac{1}{|1+\zeta\cos\phi_0|}.
$$
Instead of evaluating the density at sheet centers, integrate the charge assigned to each longitudinal bin. One option is to map a converged deterministic quadrature of uniform initial phases into final bins, counting all branches and treating periodic boundaries consistently. Bin charges must remain nonnegative and conserve the total charge. For the first runs, use a flat finite envelope containing an integer number of periods and an explicit phase origin. If a numerical approximation is used for binning, refine its resolution independently of sheet spacing.

With the paper’s negative-exponent convention, the cold-beam benchmark is
$$
b_n=\frac{1}{2\pi}\int_0^{2\pi}e^{-in(\phi_0+\zeta\sin\phi_0)}\,d\phi_0=(-1)^nJ_n(n\zeta).
$$
A change of phase origin adds a predictable complex phase. Convert to the code's convention when comparing complex values.

**Checkpoint.** Check charge conservation and compare the direct phase integral and sheet-weight Fourier sum with the Bessel result for $n=1$ and $20$, including $\zeta=0$, $\zeta=1$, and a value above 1 where the map folds. Report the sheet-spacing error separately from the binning-integration error. Leave the production solver in its Phase 1 state at this checkpoint.

## Phase 3a. Restore coherence with the fundamental density component

Add an opt-in modulated stack using Phase 2 bin charges. Share the uniform stack geometry and radiation weighting from Phase 1. Expose the modulation period, drift ratio $\zeta$, and phase origin. Use $\Lambda_{\rm mod}=\Lambda$ and start at $\zeta=1$, for which $|b_1|=|J_1(1)|\simeq0.440$. Use several complete periods for the physics comparison. Choose the integer number of bins so the envelope boundaries and total length are identical for the uniform and modulated runs. Report the resulting actual spacing if it differs slightly from the requested spacing.

**Checkpoint.** Compare a flat reference, a uniform stack, and a $b_1$-modulated stack, each at the same total charge. Show that the uniform control has negligible longitudinal bunching at the matched wavenumber for a complete-period envelope and that the modulated stack recovers the predicted complex factor, subject to the measured translation approximation. Compare 1 nm and 0.5 nm grids. Record the corresponding power ratios and, if running a screen, the residual after fitting a single complex scale to the flat-reference field.

## Phase 3b. Restore coherence with the twentieth density component

Reuse Phase 3a with $\Lambda_{\rm mod}=20\Lambda$ and the same target scattered frequency. Use the same physical stack length and total charge for the uniform, $b_1$, and $b_{20}$ comparisons when practical. Choose that length to contain an integer number of the longer periods. Keep $d$ near 1 nm or finer, since $b_{20}$ still oscillates on the $\Lambda$ scale. For an initial cheap check, one complete modulation period is sufficient to measure the discrete coefficient. Use multiple complete periods for the finite-train radiation comparison if the cost permits.

Start with the same $\zeta$ as in Phase 3a, then scan around the first maximum of $|J_{20}(20\zeta)|$ in the cheap longitudinal diagnostic. Use the resulting value in a selected radiation run. A scan of the diagnostic versus $q_z$ around the target can check the finite-train peak and distinguish it from sheet-spacing replicas.

**Checkpoint.** Report the measured complex $b_{20}$ and its cold-beam prediction, the difference between 1 nm and 0.5 nm, and the recovered field and power relative to the flat and uniform references. Compare $b_1$ and $b_{20}$ at the same target $q_z$. State the stack length, actual plane spacing, and whether the detector calculation was center-only or used a screen.

The paper’s $b_{20}\simeq0.10$ is a proposed finite-beam result. Energy and angular spread, finite pulse overlap, diffraction of the **bunching lasers**, and correlations between density and momentum can alter that value. Treat these as a later physics extension. The phases above prescribe the longitudinal density and test the coherent radiation it produces. They require no explicit two-laser particle tracking or prediction of the incoherent background.
