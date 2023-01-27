!include tutorials/tutorial_header.md

# Tutorial 1a: Simple minimization of double well free energy

This tutorial (and others) covers the basic usage of the ferroelectric phase field method implemented in FERRET.

This relatively simple problem considers solving the time-dependent LGD equation in the presence of a weak electric field along the $\hat{z}$ direction.

Consider a computational domain with a geometry $(10\times 10\times 6)$ which we define with the Mesh block.

!listing tutorial/multidomain.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. Since we will minimize an energy corresponding to a homogeneous problem with arbitrary units in the energy and system variable (polarization), the length scale will be ignored.


