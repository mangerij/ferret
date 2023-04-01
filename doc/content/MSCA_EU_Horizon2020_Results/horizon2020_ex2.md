!config navigation breadcrumbs=False

# H2020-MSCA-IF-2020, Example: FE DWs

This page details how to obtain inhomogeneous solutions of BFO for the structural order (i.e. FE DWs). The relevant input file is `BFO_dwP1A1_100.i` located in the FERRET tutorials subdirectory. The name of the input file is already suggestive, as it identifies that one component of $\mathbf{P}$ and one component $\mathbf{A}$ will switch indicating that this is for the 1/1 DW (see [!cite](Mangeri2023) on arXiv for more description of this notation). The DW plane is oriented $(100)$ which suggests the plane normal is along $x$. The 1/1 $(100)$-oriented DW only has large gradients in the components of $\partial A_z / \partial x$ and $\partial P_z / \partial x$. As such, we can simulate this DW profile in quasi 1D with the following `Mesh` block,

!listing tutorial/BFO_dwP1A1_100.i
         block=Mesh
         link=False
         language=python

where `xMax`, `yMax`, and `zMax` define the spatial dimensions of our box. The numbers `Nx`, `Ny`, and `Nz` define the number of finite elements in the calculation. Both of these quantities are listed at the beginning of the input file. We should note that since we are now looking at inhomogeneous structure of $\mathbf{P}$ and $\mathbf{A}$, then the number of finite elements is important. We find that if the spatial discretization is at least $0.1$ nm per element, then the solution is converged to smoothness. The `Functions` block is important as it allows one to set the initial conditions of $\mathbf{P}$ and $\mathbf{A}$ which in the end define the final DW profile,

!listing tutorial/BFO_dwP1A1_100.i
         block=Functions
         link=False
         language=python

 We pick sin(x) dependence for the initial conditions where the frequency corresponds to $2\pi/L_x$ where $L_x$ is the length of the box along the long direction. These leads to exactly two DWs in the simulation box. These initial conditions are set explicitly in the `Variables` block,

 !listing tutorial/BFO_dwP1A1_100.i
          block=Variables
          link=False
          language=python

with the `FunctionIC` MOOSE `Object`. We evolve the time-dependent Landau-Ginzburg equations,

\begin{equation}\label{eqn:TDLG}
  \begin{aligned}
    \frac{\partial \mathbf{P} }{\partial t} &= -\Gamma_P \frac{\delta f}{\delta \mathbf{P}}, \\
    \frac{\partial \mathbf{A} }{\partial t} &= -\Gamma_A \frac{\delta f}{\delta \mathbf{A}},
  \end{aligned}
\end{equation}

where $f$ is the total free energy density from that of BFO. As we are not looking for real dynamics of $\mathbf{P}$ and $\mathbf{A}$, we set $\Gamma_A = \Gamma_P$ to unity. We also solve the equation for mechanical equilibrium (Einstein summation convention implied),

\begin{equation}
  \begin{aligned}
    \frac{\partial \sigma_{ij}}{\partial x_j} = \frac{\partial}{\partial x_j}\left[C_{ijkl} \left(\varepsilon_{kl} + Q_{klmn} P_m P_n + R_{klmn} A_m A_n \right)\right] = 0,
  \end{aligned}
\end{equation}

which includes the electro- and rotostrictive coupling. The set of these three equations is cast into the weak formulation sufficient for finite element analysis. The set of residual and jacobian contributions are communicated to MOOSE via the `Kernels` block,

!listing tutorial/BFO_dwP1A1_100.i
         block=Kernels
         link=False
         language=python

We direct the reader to the relevant `Syntax` page for derivations of each of these `Kernels` and how they relate to the different physics. The total free energy density, $f$, is a fourth-order expansion and contains both contributions from homogeneous phases and also inhomogeneous phases. The latter is introduced as,

\begin{equation}
  \begin{aligned}
  F_\mathrm{\nabla A} = \int\limits_\Omega d^3 \mathbf{r} \left\{ \frac{H_{11}}{2}   \left( A_{x,x}^2 + A_{y,y}^2 + A_{z,z}^2 \right) +  H_{12}  \left(A_{x,x} A_{y,y} + A_{y,y} A_{z,z} + A_{x,x} A_{z,z} \right) + \frac{H_{44}}{2} \left[\left(A_{x,y} + A_{y,x} \right)^2+ \left(A_{y,z} + A_{z,y} \right)^2 + \left(A_{x,z} + A_{z,x}\right)^2\right] \right\}
  \end{aligned}
\end{equation}
and

\begin{equation}
  \begin{aligned}
  F_\mathrm{\nabla A} = \int\limits_\Omega d^3 \mathbf{r} \left\{ \frac{H_{11}}{2}   \left( A_{x,x}^2 + A_{y,y}^2 + A_{z,z}^2 \right) +  H_{12}  \left(A_{x,x} A_{y,y} + A_{y,y} A_{z,z} + A_{x,x} A_{z,z} \right) + \frac{H_{44}}{2} \left[\left(A_{x,y} + A_{y,x} \right)^2+ \left(A_{y,z} + A_{z,y} \right)^2 + \left(A_{x,z} + A_{z,x}\right)^2\right] \right\}
  \end{aligned}
\end{equation}

The relevant materials coefficients (in base S.I. units of `nanometers`, `seconds`, `attocoulombs`, and `kilograms`) are seeded through the `Materials` block as,

!listing tutorial/BFO_dwP1A1_100.i
         block=Materials
         link=False
         language=python

Note that, in principle, our methodology allows for eighth order thermodynamic potentials but the potential relevant for BFO, as parameterized in [!cite](Fedorova2022) is of fourth-order. Therefore, there are many zeroes in this `Materials` block. We also have the `Postprocessors` block,

!listing tutorial/BFO_dwP1A1_100.i
         block=Postprocessors
         link=False
         language=python

which tracks a volume integral over the following free energy density contributions to $f$,

\begin{equation}
  \begin{aligned}
  f &= f_P + f_A + f_{AP} + f_{P\varepsilon} + f_{A\varepsilon} + f_\varepsilon + f_{\nabla P} + f_{\nabla A}.
  \end{aligned}
\end{equation}

where $f_P$ is `FbP`, $f_A$ is `FbA`, $f_{AP}$ is `FcPA`, $f_{\nabla P}$ is `FgP`, $f_{\nabla A}$ is `FgA`, $f_{P\varepsilon}$ is `FcPu`, $f_{A\varepsilon}$ is `FcAu`, and $f_\varepsilon$ is `Felu`. The `Fele` quantitity is due to the interaction energy $-\mathbf{P}\cdot\mathbf{E}$. Also within this block is the summation of all of these energies and a tracker `elapsed` for wall clock time. The boundary conditions are given by the `BCs` block,

!listing tutorial/BFO_dwP1A1_100.i
         block=BCs
         link=False
         language=python

This set of BCs along with the `GlobalStrain` system in MOOSE allows us to have periodicity in the strain tensor components $\varepsilon_{ij}$ as well as the primary order parameters $\mathbf{P}$ and $\mathbf{A}$. We refer the reader to the tutorial problem `Ferroelectric domain wall` for more information on the `GlobalStrain` system. We have the following `Executioner` options,

!listing tutorial/BFO_dwP1A1_100.i
         block=Executioner
         link=False
         language=python

where we implement the backwards- finite difference time integration `bdf2` and the preconditioned Jacobi-free Newton-Krylov (`PJFNK`) method for the solve. A possible visualization of the DW profile in the `Exodus` output is provided below

!media media/tut_FE_DW.png style=display:block;margin:auto;width:60%; caption=Left: components of $\mathbf{A}$, Middle: $\mathbf{P}$, and Right: the strain tensor $\varepsilon_{ij}$ across the DW. id=tut_FE_DW

This calculation runs in 68.3 seconds on 6 processors using the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) distribution of MOOSE. Other DWs can be obtained by switching out the relevant `Mesh` and `Function` blocks accordingly.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This project [SCALES - 897614](https://cordis.europa.eu/project/id/897614) was funded for 2021-2023 at the [Luxembourg Institute of Science and Technology](https://www.list.lu/) under principle investigator [Jorge Íñiguez](https://sites.google.com/site/jorgeiniguezresearch/). The research was carried out within the framework of the [Marie Skłodowska-Curie Action (H2020-MSCA-IF-2019)](https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/opportunities/topic-details/msca-if-2020) fellowship.

!media media/euflag.png style=display:block;margin-left:auto;margin-right:auto;width:12%;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!content pagination previous=MSCA_EU_Horizon2020_Results/horizon2020_ex1.md next=MSCA_EU_Horizon2020_Results/horizon2020_ex3.md
