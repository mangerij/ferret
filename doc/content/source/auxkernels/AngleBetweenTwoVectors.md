# AngleBetweenTwoVectors

!alert construction title=Undocumented Class
The AngleBetweenTwoVectors has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/AngleBetweenTwoVectors

## Overview

Computes the angle (in degrees) between two vectors $\mathbf{v}_{1}$ and $\mathbf{v}_{2}$ defined by,

\begin{equation}
  \begin{aligned}
    \theta = 180^\circ - \left(\frac{180^\circ}{\pi}\right) \cos^{-1}{\left(\frac{\mathbf{v}_1 \cdot \mathbf{v}_2}{|\mathbf{v}_{1}||\mathbf{v}_{2}|}\right)}.
  \end{aligned}
\end{equation}

## Example Input File Syntax

!! Describe and include an example of how to use the AngleBetweenTwoVectors object.

!syntax parameters /AuxKernels/AngleBetweenTwoVectors

!syntax inputs /AuxKernels/AngleBetweenTwoVectors

!syntax children /AuxKernels/AngleBetweenTwoVectors
