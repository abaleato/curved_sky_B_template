# curved_sky_B_template

Method for evaluating the lensing B-modes on the curved-sky to leading order in their perturbative expansion. This fast implementation drinks heavily from D.Hanson's [Quicklens](https://github.com/dhanson/quicklens) code used in the [Planck 2015 lensing analysis](https://arxiv.org/pdf/1502.01591.pdf) and itself based on the position-space approach of [Dvorkin & Smith](https://arxiv.org/pdf/0812.1566.pdf).

#### Usage
```
B_template_vlm = curved_sky_B_template.integration_functions.weights.B_template(np.ones(lmax)).eval_fullsky(wiener_filtered_phi, wiener_filtered_e_alm)

# Extract the gradient and curl modes. Discard the latter.
g_B_alm, c_B_alm = curved_sky_B_template.shts.util.vlm2alm(B_template_vlm)
```

#### Installation
The code can be run either from this directory, or installed by 
running "python setup.py install".

The code is primarily written in Python, although some low level 
spherical harmonic transform (SHT) routines are written in Fortran 
for speed. In order to use these SHT routines, the code must be 
compiled. This is done automatically when installing, however if 
running the code without installing it needs to be built with

"python setup.py build_ext --inplace"

Depending on system, you may need to specify a fortran compiler. 
For example

"python setup.py build_ext --inplace --fcompiler=gnu95"
