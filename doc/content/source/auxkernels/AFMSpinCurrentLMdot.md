# AFMSpinCurrentLMdot

!alert construction title=Undocumented Class
The AFMSpinCurrentLMdot has not been documented. The content listed below should be used as a starting point for
documenting the class, which includes the typical automatic documentation associated with a
MooseObject; however, what is contained is ultimately determined by what is necessary to make the
documentation clear for users.

!syntax description /AuxKernels/AFMSpinCurrentLMdot

## Overview

Computes the $k^\mathrm{th}$ component of the antiferromagnetic (AFM) spin current contribution due to

\begin{equation}
  \begin{aligned}
    j_k^\mathrm{L\dot{m}} = f \epsilon_{ijk} L_i \frac{\partial m_j}{\partial x_j}\
  \end{aligned}
\end{equation}

where $\epsilon_{ijk}$ is the Levi-Civita pseudo-tensor and $f$ is a factor that can be used to scale the postprocessed quantity. Note that $\mathbf{L}$ and $\mathbf{m}$ are normalized quantities from the FERRET output. The total spin current is $$

\begin{equation}
  \begin{aligned}
    j_k^\mathrm{total} = j_k^\mathrm{L\dot{L}}+j_k^\mathrm{L\dot{m}}+j_k^\mathrm{m\dot{L}}+j_k^\mathrm{m\dot{m}}
  \end{aligned}
\end{equation}

as described in [!cite](Cheng2014).

## Example Input File Syntax

!! Describe and include an example of how to use the AFMSpinCurrentLMdot object.

!syntax parameters /AuxKernels/AFMSpinCurrentLMdot

!syntax inputs /AuxKernels/AFMSpinCurrentLMdot

!syntax children /AuxKernels/AFMSpinCurrentLMdot
