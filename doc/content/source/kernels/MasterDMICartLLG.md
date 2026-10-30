# MasterDMICartLLG

!syntax description /Kernels/MasterDMICartLLG

## Overview

Residual contribution of the Dzyaloshinskii-Moriya interaction (DMI) of a ferromagnet to the LLG equation, in the same
weak form and prefactor convention as [MasterExchangeCartLLG](MasterExchangeCartLLG.md). Two energy densities are
available through `dmi_type`,

\begin{equation}
  \begin{aligned}
    f_\mathrm{DMI}^\mathrm{bulk} &= D\, \mathbf{m} \cdot \left(\nabla \times \mathbf{m}\right), \\
    f_\mathrm{DMI}^\mathrm{interfacial} &= D \left[ m_z \left(\nabla \cdot \mathbf{m}\right) - \left(\mathbf{m}\cdot\nabla\right) m_z \right],
  \end{aligned}
\end{equation}

where the interfacial (Cnv) form assumes the interface normal along $+\hat{z}$. With the torque operator
$T_\mathbf{m}\mathbf{H} = \mathbf{m}\times\mathbf{H} + \alpha\, \mathbf{m}\times(\mathbf{m}\times\mathbf{H})$ and the test
function $\mathbf{v} = \psi\, \hat{e}_k$ for component $k$, the residual is

\begin{equation}
  R_k = -\frac{\gamma_0 D}{\mu_0 (1+\alpha^2) M_s} \int\limits_\Omega d^3\mathbf{r}\, B\left(\mathbf{m}, T_\mathbf{m}^T \mathbf{v}\right),
\end{equation}

where $B$ is the symmetric bilinear form with $B(\mathbf{m},\mathbf{m}) = 2 f_\mathrm{DMI}/D$, i.e. the first variation of
$f_\mathrm{DMI}$. The gradient of $T_\mathbf{m}^T\mathbf{v}$ carries both $\nabla\psi$ and $\nabla\mathbf{m}$, and $|\mathbf{m}|^2$ is
not replaced by 1. The form is NOT integrated by parts, so the DMI surface term is retained: assembled with
[MasterExchangeCartLLG](MasterExchangeCartLLG.md) it imposes the natural boundary condition

\begin{equation}
  \begin{aligned}
    2A_e\, \partial_\mathbf{n}\mathbf{m} + D\, \mathbf{m}\times\mathbf{n} &= 0 \quad \text{(bulk)}, \\
    2A_e\, \partial_\mathbf{n}\mathbf{m} + D \left[ m_z \mathbf{n}_\parallel - (\mathbf{m}\cdot\mathbf{n}_\parallel)\hat{z} \right] &= 0 \quad \text{(interfacial)},
  \end{aligned}
\end{equation}

without a BC object. The material property `D` carries the length scale in the same way as `Ae`: for a mesh in nm
and $D$ in J/m$^2$, supply $D/10^{-9}$. The sign convention of the interfacial form is that of MuMax3 (`Dind`).
The matching energy is [MasterMagneticDMIEnergy](MasterMagneticDMIEnergy.md).

## Landau-Lifshitz-Bloch route

With `llb = true` (bulk and interfacial `dmi_type`) the kernel solves the LLB equation ([Finite-temperature micromagnetics](theory/intro_LLB.md)); `D` is $D(T)/m_e^2$ and the material properties `alpha_perp` and `alpha_par` are required. As for [MasterExchangeCartLLG](MasterExchangeCartLLG.md), the optional `grad_alpha_perp` and `grad_alpha_par` properties add the terms that integration by parts produces when the damping varies in space.

!listing test/tests/magnets/llb_alpha_ramp_dmi.i block=Kernels/MasterDMICartLLG_x

## Example Input File Syntax

The DMI class is set with `dmi_type`, here through `[GlobalParams]`. Interfacial DMI on a 1D wire
(SP6 of the DMI standard problems, Cortés-Ortuño et al., New J. Phys. 20, 113015 (2018)):

!listing test/tests/magnets/sp6_1d.i block=GlobalParams

!listing test/tests/magnets/sp6_1d.i block=Kernels/dmi_x

Bulk DMI with a Zeeman field on a cylinder (the 3D case of the same benchmark):

!listing test/tests/magnets/sp6_3d.i block=GlobalParams

The material property `D` is supplied with the other constants:

!listing test/tests/magnets/sp6_3d.i block=Materials

!syntax parameters /Kernels/MasterDMICartLLG

!syntax inputs /Kernels/MasterDMICartLLG

!syntax children /Kernels/MasterDMICartLLG
