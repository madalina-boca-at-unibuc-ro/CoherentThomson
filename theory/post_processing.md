# post_processing of the emitted field, stage 1: 

the overall goal of this approach is to study the structure of the emitted radiation, and to establish if it can be written as a superposition of LG modes / identify the parameters of the LG modes in the superposition.

## 1. elimination of the radial phase.

In the output folder, make a new subfolder png_folder/radiation/processed 
there, plot in the same format as in radiation/emitted/ but plotted quantity should be 

$$ F^{\mu \nu} e^{-i k_{N} (R-Z)}$$ 

were $k_N=\omega_N/c$ the nominal wavevector, $Z$ is the distance to the screen, and 

$$R=\sqrt{X^2+Y^2+Z^2}$$ 

with $X$, $Y$,  the coordinates on the screen.

