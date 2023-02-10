# RandomConstrainedVectorFieldIC

!syntax description /ICs/RandomConstrainedVectorFieldIC

## Overview

Seeds a vector field that is constrained to the unit sphere. Typically useful as `ICs` for micromagnetic simulations where $|\mathbf{m}| = 1$ with random or pseudo-random angles. The usual spherical coordinate prescription is used with $\phi$ and $\theta$ the azimuthal and polar angles respectively.

\begin{equation}
  \begin{aligned}
    m_x &= M_s^0 \cos{\left(\phi\right)} \sin{\left(\theta\right)} \\
    m_y &= M_s^0 \sin{\left(\phi\right)} \sin{\left(\theta\right)} \\
    m_z &= M_s^0 \cos{\left(\theta\right)}.
  \end{aligned}
\end{equation}

The parameter $M_s^0$ can be used to constrain the `ICs` to something other than the unit sphere if needed.

## Example Input File Syntax

!listing tutorial/ringdown.i
         block=Variables
         link=False
         language=python

Note that $\phi$ and $\theta$ are fields that are seeded via `AuxVariables` in this example.

!syntax parameters /ICs/RandomConstrainedVectorFieldIC

!syntax inputs /ICs/RandomConstrainedVectorFieldIC

!syntax children /ICs/RandomConstrainedVectorFieldIC
