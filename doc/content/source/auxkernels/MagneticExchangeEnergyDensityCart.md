# MagneticExchangeEnergyDensityCart

!alert construction title=Undocumented Class
The MagneticExchangeEnergyDensityCart has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/MagneticExchangeEnergyDensityCart

## Overview

Computes the magnetic free energy density due to the exchange interaction as,

\begin{equation}
  \begin{aligned}
  f_\mathrm{exch} = A_e M_s^2 \left\{\left(\nabla m_x\right)^2 + \left(\nabla m_y\right)^2 + \left(\nabla m_z\right)^2\right\}
  \end{aligned}
\end{equation}

where $M_s$ is the saturation magnetization density and $A_e$ is the exchange stiffness constant.

## Example Input File Syntax

!! Describe and include an example of how to use the MagneticExchangeEnergyDensityCart object.

!syntax parameters /AuxKernels/MagneticExchangeEnergyDensityCart

!syntax inputs /AuxKernels/MagneticExchangeEnergyDensityCart

!syntax children /AuxKernels/MagneticExchangeEnergyDensityCart
