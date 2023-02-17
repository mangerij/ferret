# MagneticExchangeEnergyDensityCart

!syntax description /AuxKernels/MagneticExchangeEnergyDensityCart

## Overview

Computes the magnetic free energy density due to the exchange interaction as,

\begin{equation}
  \begin{aligned}
  f_\mathrm{exch} = A_e M_s^2 \left\{\left(\frac{\partial m_x}{\partial x}\right)^2 + \left(\frac{\partial m_x}{\partial y}\right)^2 + \left(\frac{\partial m_x}{\partial z}\right)^2 + \left(\frac{\partial m_y}{\partial x}\right)^2 + \left(\frac{\partial m_y}{\partial y}\right)^2 + \left(\frac{\partial m_y}{\partial z}\right)^2 + \left(\frac{\partial m_z}{\partial x}\right)^2 + \left(\frac{\partial m_z}{\partial y}\right)^2 + \left(\frac{\partial m_z}{\partial z}\right)^2\right\}
  \end{aligned}
\end{equation}

where $M_s$ is the saturation magnetization density and $A_e$ is the exchange stiffness constant.

## Example Input File Syntax

!! Describe and include an example of how to use the MagneticExchangeEnergyDensityCart object.

!syntax parameters /AuxKernels/MagneticExchangeEnergyDensityCart

!syntax inputs /AuxKernels/MagneticExchangeEnergyDensityCart

!syntax children /AuxKernels/MagneticExchangeEnergyDensityCart
