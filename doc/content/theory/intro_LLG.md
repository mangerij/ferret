# Micromagnetic theory

!alert construction title=Documentation in-progress
This section requires some work before it will be finalized online. Please contact the developers if you need assistance with this aspect of the module.

The micromagnetic approach to modeling magnetic systems is similar to the phase field method in the same sense that a continuum description of an order parameter is provided. In ferromagnets, below the critical temperature, the spins order forming a nonzero and net spin magnetization $\mathbf{M}$. Additionally, the spins see a specific symmetry characterized by the presence (or absence) of magnetocrystalline anisotropy. This drives the tendency to form ordered domains but competitions from the ferromagnetic exchange or the demagnetizing fields (internal or external) can favor other spin textures such as skyrmions or vortices. To predict the magnetic ground states (or their dynamics), the Landau-Lifshitz-Gilbert equation is used. This coupled nonlinear equation reads,


\begin{equation}\label{LLG}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H}
  \end{aligned}
\end{equation}

where $\gamma$ is the electron gyromagnetic ratio and $\alpha$ is a phenomenological damping coefficient. The second term in the governing equation drives the system towards its energy minima (or system equilibrium). The vector $\mathbf{H}$ is the effective field which arises from different physics of the magnet (exchange, anisotropy, spin-torque transfer...).

In contrast to ferroelectrics where the magnitude of the electric dipole $\mathbf{P}$ can stretch or shrink, the normalized magnetization $\mathbf{m} = \mathbf{M} / M_s$ must stay on the unit sphere during the time evolution. There are two possible systems to use to enforce this constraint. One possibility is to renormalize $\mathbf{m}$ at the conclusion of every time step in the evolution of Eq. (1). Another option is to evolve the Landau-Liftshitz-Bloch equation

\begin{equation}\label{LLG-LLB}
  \begin{aligned}
    \frac{\partial \mathbf{M}}{\partial t} = -\gamma \mathbf{M}\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M} \times \mathbf{M}\times \mathbf{H} - \frac{\gamma \alpha_\mathrm{LLB}}{1+\alpha^2} \left(m^2 - 1\right) \mathbf{m}
  \end{aligned}
\end{equation}

which provides a longitudinal restoring force to $\mathbf{m}$ to map back onto the unit sphere. Both methods are essentially equivalent. We have both possibilities within the Ferret/MOOSE code and their implementation is up to the user.

If the material is an antiferromagnet, then the spins form two or more sublattices $(\mathbf{M}_1, \mathbf{M}_2, ... )$ whose interaction is mediated by the exchange. In a collinear antiferromagnet, the total net magnetization, $\mathbf{M} = \mathbf{M}_1 + \mathbf{M}_2$, is zero indicating that the spin sublattices perfectly compensating each other. To simulate the micromagnetics of an antiferromagnet, we have the two sublattice LLG equation,

\begin{equation}
  \begin{aligned}
    \frac{\partial \mathbf{M}_\eta}{\partial t} = -\gamma \mathbf{M}_\eta\times \mathbf{H} - \frac{1}{M_s} \frac{\gamma \alpha}{1+\alpha^2} \mathbf{M}_\eta \times \mathbf{M}_\eta \times \mathbf{H},
  \end{aligned}
\end{equation}

where $\eta = 1, 2$. Note as mentioned before, one can choose to use the LLB approximation and/or also let the modulus be renormalized at every time step. Both methods are equivalent (more discussion in the recent results section of the MSCA here link).
