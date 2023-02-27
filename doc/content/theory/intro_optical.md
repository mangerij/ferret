# Optical properties in FERRET

!alert construction title=Documentation in-progress
This section requires some work before it will be finalized online. Please contact the developers if you need assistance with this aspect of the module.

In this section, we describe the postprocessing capability of FERRET to compute refractive indices for an arbitrary anisotropic medium. The medium can have internal structure (in the case of ferroelectrics) or can be nominally isotropic but under the influence of elastic or electric fields (i.e. elastoptics, electrooptics, or piezooptics).

The refractive index is defined by the quantity known as the indicatrix

\begin{equation}
    B_{ij} x_{i} x_{j} = 1,
\end{equation}

in which $B_{ij} = \kappa_0 \partial E_i / \partial D_j$ with $\kappa_0$ the permittivity of vacuum. In an isotropic medium, $B_{11} = B_{22} = B_{33}$ and therefore the indicatrix can be considered to be a perfect sphere in coordinate $x_1, x_2, x_3$ space. The quantities $B_{ij}$ are reciprocal to the relative dielectric permeability or $B_{ij} = 1/\kappa_{ij}$ and $B_{ij} = 1/n_{ij}^2$ where $n_{ij}$ is the refective index tensor. Usually the off-diagonal elements are very small or zero. We can see an example of the cubic medium in Fig. \ref{fig_indicatrix}.

!media media/indicatrix.png style=display:block;margin:auto;width:50%; caption=Left; Isotropic (cubic) medium. Right; Lower symmetry indicatrix characterized by a change of $n_z$ with respect to $n_x = n_y$.  id=fig_indicatrix

If the crystal symmetry is broken, then the indicatrix will also be modified thus losing the $n_x = n_y = n_z$ relationship as shown in the right panel of Fig. \ref{fig_indicatrix}. This can be due to a phase transition or an applied external field (i.e. electric or mechanical). Since there are many possible ways to change a crystal structure, we will only describe a few that make sense in the context of FERRET simulations. There are two assumptions to the phenomological theory prescribed to the topical properties of crystals and they are:

\begin{itemize}
    \item In a homogeneously deformed solid the effect of deformation is only to alter the parameters of the optical indicatrix.
    \item When strain is within elastic limits, the changes of an optical parameter (polarization constant) due to deformation can be expressed as a homogeneous linear function of the nine strain or stress components.
    \item The electromagnetic (EM) field is weak and does not induce secondary changes to the refractive index.
\end{itemize}

Typically, this last assumption is usually what is applicable to experiments involving transmission or scattering. To lowest order, the relative change in the indicatrix components due to electric and mechanical fields are given by [!cite](NyeBook),

\begin{equation}
    \Delta B_{ij} = r_{ijk} E_k + \pi_{ijkl} \sigma_{kl}
\end{equation}

where $r_{ijk}$ and $\pi_{ijkl}$ is the electrooptic and photoelastic tensors. Typical orders of magnitude and units of $r_{ijk}$ and $\pi_{ijkl}$ are $10^{-12}$ m/V and $10^{-12}$ m/N respectively. The quantity $\sigma_{kl}$ is the elastic stress tensor component. The relative change of the indicatrix relationship can also be expressed in terms of strain $\varepsilon_{rs}$,

\begin{equation}
    \Delta B_{ij} = r_{ijk} E_k + \pi_{ijkl} C_{klrs} \varepsilon_{rs}
\end{equation}

with $C_{klrs}$ is the elastic stiffness tensor. As an example, we can consider an applied uniaxial stress $\sigma_{xx}^0$ on a cubic crystal (\mathbf{E} = 0). Before the deformation,

\begin{equation}
    B_{0} \left(x^2_1 + x^2_2 + x^2_3\right) = 1
\end{equation}

with $B_{0} = 1/n_0^2$. After deformation, by using the above expressions for $\Delta B_{ij}$ and a cubic symmetry for $\pi_{ijkl}$, we find (in Voight notation),

\begin{equation}
\begin{aligned}
    \Delta n_1 &= - \frac{1}{2} n_0^3 \pi_{11} \sigma_{xx}^0 \\
    \Delta n_2 &= - \frac{1}{2} n_0^3 \pi_{13} \sigma_{xx}^0 \\
    \Delta n_3 &= - \frac{1}{2} n_0^3 \pi_{12} \sigma_{xx}^0
\end{aligned}
\end{equation}

Note that if the crystal symmetry before deformation is cubic symmetry of $\bar{4}3m, 432$, or $m3m$, then $\pi_{12} = \pi_{13}$ and the resulting indicatrix relationship is $n_x = n_y \neq n_z$ as in Fig. \ref{fig_indicatrix}. However, if the cubic symmetry is $23$ or $m3$, then $\pi_{12} \neq \pi_{13}$ resulting in $n_x \neq n_y \neq n_z$ which has much consequence on an observed birefringence.

In principle any symmetry crystal can be investigated in FERRET simulations and that these expressions also work for inhomogeneously applied fields breaking the symmetry locally leading to spatially varying refractive indices.

Another example of possible postprocessing ability in FERRET is to compute the relative changes in the refractive index due to the ferroelectric phase transition. Similar to an applied stress which breaks the symmetry leading to birefringence, the onset of the electric polarization moment arising from a structural distortion of the ionic lattice also breaks the symmetry. In the case of canonical perovskite ferroelectric (FE) $\mathrm{BaTiO}_3$, the symmetry initially goes from cubic $(m3m)$ to tetragonal $(4mm)$ as the temperature is lowered below the Curie temperature $T_C$. We consult the work of [!cite](Bernasconi1995) for a description.

The spontaneous polarization $P_S$ enters into the relative changes to the indicatrix coefficients via,

\begin{equation}
\begin{aligned}
    \Delta B_1 = \Delta B_2 &= g_{1122} P_S^2 + \left(p_{1111} + p_{1122}\right)\varepsilon_{xx} + p_{1122} \varepsilon_{zz} \\
    \Delta B_3 &= g_{1111} P_S^2 + 2 p_{1122} \varepsilon_{xx} + p_{1111} \varepsilon_{zz}
\end{aligned}
\end{equation}

where $\varepsilon_{ij}$ is the spontaneous strain arising from the phase transition. One can appreciate here that both the polarization and the elastic field adjust the refractive indices leading to a birefringence $(n_1 = n_2 \neq n_3)$. The parameters $(g_{ijkl})$ are known as the polar-optic coefficients in units of $10^{-2}$ $\mathrm{m}^4$ $\mathrm{C}^{-2}$. These materials constants have been shown to be strongly temperature dependent (as is for $P_S$) as well as dependent on the wavelength of light $\lambda_0$. Using room temperature values of $P_S$ and $\lambda_0 = 633$ nm, we have $n_3^{-2} - n_1^{-2} \approx 10^{-1}$ due to the FE phase transition. For weak EM fields, and at frequencies of the visible spectrum, the FE polarization is expected to be fixed in time which is a significant approximation to be taken into consideration. We should also mention that the values of $g_{ijkl}$ are not known for many FE materials but they should be able to be computed from first-principles methodology or measured carefully as in the case of $\mathrm{BaTiO}_3$ in the cited work.

As the birefringence is one way to observe and measure the domain topology, we provide this feature in FERRET for users. By first using the phase field method to calculate the ground state structure (in arbitrary 3D geometries), the local refractive indices can then be computed within this approach. We aim to provide an example in our tutorials listed on this website (in construction) and we plan to expand this capability in the future to handle lower symmetries and other FE materials such as those displaying magnetic order (multiferroics).
