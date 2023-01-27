# Constitutive theory of thermoelectrics

Thermoelectric materials are characterized by coupled, or interdependent conduction of heat and electricity. When heat is flowing in a thermoelectric material, an electric current arises. The converse effect also happens when a voltage is applied thus changing the temperature. The static governing equations of a thermoelectric material are given by, 

\begin{equation}
  \begin{aligned}
    \nabla \cdot \mathbf{j} = 0,
  \end{aligned}
\end{equation}

and 

\begin{equation}
  \begin{aligned}
    \nabla \cdot \mathbf{q} + \mathbf{j} \cdot \nabla \Phi_\mathrm{E} = 0,
  \end{aligned}
\end{equation}

where 

\begin{equation}
  \begin{aligned}
    j_{i} = - \gamma_{ij} \frac{\partial \Phi_\mathrm{E}}{\partial x_{j}} - \gamma_{ij} S_{jk} \frac{\partial T}{\partial x_{k}},
  \end{aligned}
\end{equation}

and 

\begin{equation}
  \begin{aligned}
    q_{i} = - \kappa _{ij} \frac{\partial T}{\partial x_{j}} - \gamma_{ij} S_{jk} T \frac{\partial V}{\partial x_{k}} - \gamma_{ij} S_{jk} S_{kl} T \frac{\partial T}{\partial x_{l}}.
  \end{aligned}
\end{equation}

The tensor components $\kappa_{ij}, S_{ij}$, and $\gamma_{ij}$ correspond to the thermal conductivity, Seebeck, and electrical conductivity respectively which require information of the material crystallographic directions. 

Within our block-restricted polycrystal approach, this allows for the thermoelectric properties to be evaluated in a computational box with a \textit{real} grain structure. More details of this method can be found in the recent work of Basuala et al Acta Materialia 228, 117743 (2022) published using FERRET. 
