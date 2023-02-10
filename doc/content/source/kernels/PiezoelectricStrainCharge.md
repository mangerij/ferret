# PiezoelectricStrainCharge

!alert construction title=Undocumented Class
The PiezoelectricStrainCharge has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /Kernels/PiezoelectricStrainCharge

## Overview

The Poisson equation includes contributions from the (converse piezo) strain-charge,

\begin{equation}
  \begin{aligned}
    \frac{\partial}{\partial x_r} \epsilon_b \frac{\partial \Phi_\mathrm{E} }{\partial x_r} + \frac{\partial (d_{jkl} \sigma_{kl})}{\partial x_j} = 0,
  \end{aligned}
\end{equation}

 with $\sigma_{kl} = C_{ijkl} \varepsilon_{kl}$ being the components of the elastic stress tensor and $\epsilon_b$ being a relative static permittivity.
 Multiplying by a test function and integrating by both sides,

\begin{equation}
  \begin{aligned}
   \left(\psi_h, \frac{\partial}{\partial x_r} \epsilon_b \frac{\partial \Phi_\mathrm{E} }{\partial x_r}\right) + \left(\psi_h,\frac{\partial (d_{jkl} \sigma_{kl})}{\partial x_j}\right) = 0,
  \end{aligned}
\end{equation}

If $\epsilon_b$ does not depend on space (as is usually the case), we can rewrite the above equation by integrating by parts,

\begin{equation}
  \begin{aligned}
    \underbrace{-\left(\epsilon_b \frac{\partial \psi_h }{\partial x_r}, \frac{\partial \Phi_\mathrm{E}}{\partial x_r}\right)
    + \langle \epsilon_b \frac{\partial \psi_h }{\partial x_r}, \frac{\partial \Phi_\mathrm{E}}{\partial x_r} \rangle}_\mathrm{Electrostatics} \,\,\,\,\, - \,\,\,\,\,
    \underbrace{ \left(\frac{\partial \psi_h}{\partial x_j}, d_{jkl} \sigma_{kl} \right)
    + \langle \frac{\partial \psi_h}{\partial x_j}, d_{jkl} \sigma_{kl} \rangle}_\mathrm{PiezoelectricStrainCharge} = 0,
  \end{aligned}
\end{equation}

The surface terms $\langle . \rangle$ vanish due to the properties of the test function. Note that the `Electrostatics` must be included in the `Kernels` block of the input file in order for this term to couple properly to $\Phi_\mathrm{E}$.

## Example Input File Syntax

!! Describe and include an example of how to use the PiezoelectricStrainCharge object.

!syntax parameters /Kernels/PiezoelectricStrainCharge

!syntax inputs /Kernels/PiezoelectricStrainCharge

!syntax children /Kernels/PiezoelectricStrainCharge
