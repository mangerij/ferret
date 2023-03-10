!include tutorials/tutorial_header.md

# Tutorial 5: Antiferromagnetic dynamics

!alert construction title=To be finished
The tutorial for antiferromagnetic dynamics has not been completed and uploaded to the website. We expect to do this in the near future.

This tutorial covers the basic implementation of antiferromagnetic (AFM) dynamics implemented in FERRET. We refer the reader to the first example in [!cite](Rezende2019) which considers the unaxial AFM material $\mathrm{MnF}_2$. The magnetic ion is $\mathrm{Mn}$ leading to two sublattices arranged collinearly. The free energy density is,

\begin{equation}
  \begin{aligned}
    f = - H_0 \left(M_{1z} + M_{2z}\right) + \frac{H_\mathrm{E}}{M_s} \mathbf{M}_1 \cdot \mathbf{M}_2 - \frac{H_\mathrm{A}}{2 M_s} \left(M_{1z}^2 + M_{2z}^2\right)
  \end{aligned}
\end{equation}

where $H_\mathrm{E}$ and $H_\mathrm{A}$ are the effective exchange and anisotropy fields. Only nearest-neighbors along the $z$ direction are considered in this example. We evolve the two-sublattice Landau-Lifshitz-Gilbert-Bloch (LLG-LLB) equations,

\begin{equation}\label{eqn:LLG_LLB}
  \begin{aligned}
    \frac{d \mathbf{m}_\eta}{dt} = -\frac{\gamma}{1+\alpha^2} \left(\mathbf{m}_\eta \times \mathbf{H}_\eta\right) - \frac{\gamma\alpha}{1+\alpha^2}\mathbf{m}_\eta\times\left(\mathbf{m}_\eta\times\mathbf{H}_\eta\right) + \frac{\gamma\tilde\alpha_\parallel}{(1+\alpha^2)} m_{\eta}^2 \left[ m_{\eta}^2 - 1 \right]\mathbf{m}_\eta,
  \end{aligned}
\end{equation}

where the constant $\gamma$ corresponds to the electron gyromagnetic factor equal to $2.8$ GHz, $\alpha$ the Gilbert damping, and $\tilde\alpha_\parallel$ the longitudinal damping constant from the LLB approximation. The effective fields are defined as $\mathbf{H}_\eta = - \mu_0^{-1} M_s^{-1} \delta f / \delta \mathbf{m}_\eta$ with $\mu_0$ the permeability of vacuum. In this example, we choose units of $H_0 = 80$ kOe, $H_\mathrm{E} = 526$ kOe, and $H_\mathrm{A} = 8.2$ kOe. One can see that the factor of $M_s$ cancels so we do not need to set its value other than $1.0$ in the input files. For this problem, we only consider the conservative dynamics given by $\alpha = 0$. This is the linear spin wave excitation limit and this problem is solved analytically in the above reference. Note that $H_0$ is the strength of an applied field which interacts with the non-equilibrium AFM spins via the Zeeman interaction. If $H_0$ is sufficiently large, a spin flop transition can occur. We choose $H_0$ below this critical value.

This tutorial problem is set up as a macrospin simulation (no spatial gradients). The computational box is simply a $2x2x2$ finite element grid with $0.01$ $\mu$m size as shown in the below Mesh block,

We seed the coefficients ($\gamma, H_0, H_\mathrm{E},...$) through the Materials block,




!content pagination previous=tutorials/magnetic_ringdown.md next=tutorials/piezoelectric.md
