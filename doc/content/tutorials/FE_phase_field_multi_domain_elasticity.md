!include tutorials/tutorial_header.md

# Tutorial 1a: Simple minimization of double well free energy

This tutorial (and others) covers the basic usage of the ferroelectric phase field method implemented in FERRET.

This relatively simple problem considers solving the time-dependent LGD equation in the presence of a weak electric field along the $\hat{z}$ direction.

Consider a computational domain with a geometry $(10\times 10\times 6)$ which we define with the Mesh block.

!listing tutorial/problem1.i
         block=Mesh
         link=False
         language=python

In general, the geometry defined in the 'Mesh' block *never* carries units. The length scale is introduced through Materials, Kernels, or other MOOSE objects. Since we will minimize an energy corresponding to a homogeneous problem with arbitrary units in the energy and system variable (polarization), the length scale will be ignored.


Throughout FERRET, we make extensive use of the 'GlobalParams' block. For this problem,

!listing tutorial/problem1.i
         block=GlobalParams
         link=False
         language=python

Note that this is just a convenient convention for the input files to reduce length. Many of the 'GlobalParams' are required to be passed to each of the MOOSE Objects (Materials, Kernels, etc). For this problem, we also pass the Landau coefficients of a simple fourth order potential through the 'GlobalParams' block. Therefore, these are real numbers (scalars) and not Materials coefficients as we will see in later tutorial problems/examples.


The total system energy, $F$, is a fourth order Landau potential which is *uniaxial* $(P_x = P_y = 0, P_z \neq 0)$ plus a contribution from a constant electric field.

\begin{equation}\label{4thOrdLandau}
  \begin{aligned}
    F &= F_\mathrm{bulk} - P_z E_z \\
      &= \underbrace{\alpha_1 P_z^2 + \alpha_{11} P_z^4}_{F_\mathrm{bulk}} - P_z E_z
  \end{aligned}
\end{equation}

Higher order terms are zeroed in this example by setting the appropriate coefficients ($\alpha_{111}$, etc) to zero. The right-hand side of the time-dependent equation involves a variation of the energy with respect to the variables. In this case, 

\begin{equation}\label{partial4thOrdLandau}
  \begin{aligned}
    \frac{\delta F}{\delta P_z} &= \frac{\partial F}{\partial P_z} \\
                                &= 2 \alpha_1 P_z + 4 \alpha_{11} P_z^3 - E_z
  \end{aligned}
\end{equation}


