# CoupledDirichletBC

!syntax description /BCs/CoupledDirichletBC

## Overview

Computes a `NodalBC` which computes a residual contribution for a variable $u$ of the form,

\begin{equation}
  \begin{aligned}
    \mathcal{r}_{u} = \left(u[\_\mathbf{qp}] - v[\_\mathbf{qp}]\right)
  \end{aligned}
\end{equation}

for some variable $v$.

TODO:

## Example Input File Syntax

!! Describe and include an example of how to use the CoupledDirichletBC object.

!syntax parameters /BCs/CoupledDirichletBC

!syntax inputs /BCs/CoupledDirichletBC

!syntax children /BCs/CoupledDirichletBC
