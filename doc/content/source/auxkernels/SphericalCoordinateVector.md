# SphericalCoordinateVector

!alert construction title=Undocumented Class
The SphericalCoordinateVector has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/SphericalCoordinateVector

## Overview

Computes the spherical coordinates of a vector $\mathbf{v}$ in radians with,

\begin{equation}
  \begin{aligned}
    \phi = \tan{\left(\frac{v_y}{v_x}\right)}
  \end{aligned}
\end{equation}

and

\begin{equation}
  \begin{aligned}
    \theta = \cos^{-1}{\left(\frac{v_z}{v_x^2+v_y^2+v_z^2}\right)}
  \end{aligned}
\end{equation}

where $\phi$ and $\theta$ are the azimuthal and polar angles respectively.

## Example Input File Syntax

!! Describe and include an example of how to use the SphericalCoordinateVector object.

!syntax parameters /AuxKernels/SphericalCoordinateVector

!syntax inputs /AuxKernels/SphericalCoordinateVector

!syntax children /AuxKernels/SphericalCoordinateVector
