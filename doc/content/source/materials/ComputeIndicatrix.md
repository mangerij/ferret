# ComputeIndicatrix

!syntax description /Materials/ComputeIndicatrix

## Overview

Computes the indicatrix components

\begin{equation}
  \begin{aligned}
    \tilde{B}_{ij} = R_{i\beta}R_{j\gamma} B_{\beta\gamma} = R_{i\beta}R_{j\gamma} \frac{1}{n_i n_j}
  \end{aligned}
\end{equation}

with an option to rotate due to a rotation operator $R_{ij}$. The rotation operator $R_{ij}$ accepts Euler angles in the standard Bunge sequence ($\mathbf{ZXZ}$). Currently, this implementation only accepts diagonal (principle) values for the refractive indices $n_1 = n_a$, $n_2 = n_b$ $n_3 = n_g$. We plan to generalize this further at a later date.

## Example Input File Syntax

!! Describe and include an example of how to use the ComputeIndicatrix object.

!syntax parameters /Materials/ComputeIndicatrix

!syntax inputs /Materials/ComputeIndicatrix

!syntax children /Materials/ComputeIndicatrix
