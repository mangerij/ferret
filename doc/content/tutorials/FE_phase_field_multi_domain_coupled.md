!include tutorials/tutorial_header.md

# Tutorial 3: Ferroelectric thin film

This tutorial (and others) covers the basic usage of the ferroelectric (FE) phase field method implemented in FERRET. This specific example focuses on the thin film problem for $\mathrm{PbTiO}_3$ (PTO) and $\mathrm{BaTiO}_3$ (BTO) where periodicity is enforced along $x,y$ and $z$ corresponds to the film/substrate interface plane normal.

As in Tutorial 2, we use the fully-coupled problem with electrostrictive free energy density. The total free energy density for this simulation can be written as,

\begin{equation}
  \begin{aligned}
    f = f_\mathrm{bulk} + f_\mathrm{elastic} + f_{\nabla P} + f_\mathrm{electostr} - \mathbf{P}\cdot\mathbf{E}
  \end{aligned}
\end{equation}

for the bulk, elastic, gradient, electrostrictive, and electrostatic free energies respectively. We choose a computational geometry as follows,

!listing tutorial/film.i
         block=Mesh
         link=False
         language=python

where $n_x, n_y,$ and $n_z$ are chosen accordingly such that the mesh spacing $\Delta = 1.0$ nm. In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through `Materials`, `Kernels`, or other MOOSE objects. For this problem, the length scale is introduced through the units in the `Materials` objects that connect to the `Kernels`. Note that this discretization is around the upper bound for these types of calculations since the thickness of the PTO wall is about $\Delta$. Next we see that we use `SubdomainBoundingBoxGenerator` and `SideSetsBetweenSubdomainsGenerator` to split the rectilinear computational volume into two `block`s where the film sits on top of the substrate. For this example input file, we let the FE film have a thickness of 20 nm. 
