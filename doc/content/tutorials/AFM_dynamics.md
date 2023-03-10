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

where the constant $\gamma$ corresponds to the electron gyromagnetic factor equal to $2.8$ GHz, $\alpha$ the Gilbert damping, and $\tilde\alpha_\parallel$ the longitudinal damping constant from the LLB approximation. The effective fields are defined as $\mathbf{H}_\eta = - \mu_0^{-1} M_s^{-1} \delta f / \delta \mathbf{m}_\eta$ with $\mu_0$ the permeability of vacuum. In this tutorial problem, we look for conservative ($\alpha = 0$) time dependent solutions of $\mathbf{m}_\eta$. This is the linear spin wave excitation limit and this problem is solved analytically in the above reference. The variables are defined as,

!listing tutorial/AFMR_MnF2_ex.i
         block=Variables
         link=False
         language=python

where we have used [`RandomConstrainedVectorFieldIC`](source/ics/RandomConstrainedVectorFieldIC.md) to set the initial condition such that $\mathbf{m}_\eta$ is on the unit sphere  $(|\mathbf{m}_\eta| = 1)$. We choose units of $H_0 = 80$ kOe, $H_\mathrm{E} = 526$ kOe, and $H_\mathrm{A} = 8.2$ kOe. One can see that the factor of $M_s$ cancels in the above equation so we do not need to set its value other than $1.0$ in the input files. Note that $H_0$ is the strength of an applied field which interacts with the non-equilibrium AFM spins via the Zeeman interaction. If $H_0$ is sufficiently large, a spin flop transition can occur. We choose $H_0$ below this critical value.

This tutorial problem (`tutorials/AFMR_MnF2_ex.i`) is set up as a macrospin simulation (no spatial gradients). The computational box is simply a $2 \times 2 \times 2$ finite element grid with $0.01$ $\mu$m size as shown in the below `Mesh` block,

!listing tutorial/AFMR_MnF2_ex.i
         block=Mesh
         link=False
         language=python

We seed the coefficients ($\gamma, H_0, H_\mathrm{E},...$) through the `Materials` block,

!listing tutorial/AFMR_MnF2_ex.i
         block=Materials
         link=False
         language=python

The value of $\tilde\alpha_\parallel$ is set in the `ParsedFunction` called `bc_func_1`. In principle, $\tilde\alpha_\parallel$ could be time-dependent or depend on spatial coordinates, but for this problem, we just set $\tilde\alpha_\parallel = 0.5$ which is sufficient to constrain $|\mathbf{m}_\eta| = 1$. To construct the `Kernels` relevant for the partial differential equation, we have the following `Kernels` block,

!listing tutorial/AFMR_MnF2_ex.i
         block=Kernels
         link=False
         language=python

There are three FERRET `Objects` here, which correspond to `LongitudinalLLB` handling the LLB damping to constrain $|\mathbf{m}_\eta| = 1$, `TimeDerivative` for $\partial \mathbf{m}_\eta / \partial t$ and [`UniaxialAFMSublattice`](source/kernels/UniaxialAFMSublattice.md) which computes the residual and jacobian contributions corresponding to $\mathbf{m}_\eta \times \mathbf{H}_\eta$. Note that the LLG-LLB equation for this tutorial problem is hardcoded with Gilbert damping $\alpha = 0$ (linear spin wave limit). We refer the reader to the hyperlinks to see the `Kernels` computations. We also compute $\mathbf{L} = \mathbf{m}_1 - \mathbf{m}_2$ and $\mathbf{m} = \mathbf{m}_1 + \mathbf{m}_2$ using the `VectorDiffOrSum` object in the `AuxKernels` block,

!listing tutorial/AFMR_MnF2_ex.i
         block=AuxKernels
         link=False
         language=python

A possible output of this problem is shown below for the components of $\mathbf{L}$ as a function of time (in $\mu$s),

!media media/tut_AFMR.png style=display:block;margin:auto;width:60%; caption=Undamped ($\alpha = 0$) evolution of the components of the Neel vector $L_k$ as a function of time of $\mathrm{MnF}_2$. id=tut_AFMR

It can be seen that the Neel vector components along the direction of the field is relatively constant identifying that the sublattices are precessing about their initial position as expected in the linear limit of the LLG equation. In the article by [!cite](Rezende2019), the resonance frequency of the spin sublattice system is calculated analytically as a function of the field strength $H_0$. We show below good agreement with the results from FERRET for $\mathbf{H}_0 \parallel \hat{\mathbf{z}}$ and $\mathbf{H}_0 \perp \hat{\mathbf{z}}$.

!media media/tut_AFMR2.png style=display:block;margin:auto;width:60%; caption=AFM resonance frequency dependence on field strength $H_0$ calculated by FFT from FERRET calculations. The dashed lines correspond to analytical solutions for different modes. See [!cite](Rezende2019) for more information. id=tut_AFMR2

!content pagination previous=tutorials/magnetic_ringdown.md next=tutorials/piezoelectric.md
