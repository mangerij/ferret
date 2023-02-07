# Optical properties in FERRET

!alert construction title=Documentation in-progress
This section requires some work before it will be online. Please contact the developers if you need assistance with this aspect of the module.

In this section a description of the physics governing light interactions with matter are discussed.  First, an explanation of electromagnetic theory in crystals is developed based on primarily on Maxwell's equations.  This is followed by a discussion on the phenomenological description of wave propagation in crystals, mainly by connecting refractive indices to the optical indicatrix.  Lastly, a description of the photoelastic and electro-optic effects, which involve changes to the optical properties of crystals under applied external stresses and electric fields.

We start by considering Maxwell's equations, in an isotropic medium for simplicity sake, presented in tensor form,
\begin{align}
    \frac{\partial D_i}{\partial r_i} &= \rho \\
    \frac{\partial B_i}{\partial r_i} &= 0 \\
    \epsilon_{ijk} \frac{\partial E_k}{r_j} &= - \frac{\partial B_i}{dt}\\
    \epsilon_{ijk} \frac{\partial H_k}{r_j} &= j_i + \frac{\partial D_i}{\partial t}
\end{align}
where $\epsilon_{ijk}$ is the Levi-Civita symbol, $\rho$ is charge density, $j_i$ is current density, $E_i$ and $B_i$ are electric and magnetic field strength, and $D_i$ and $M_i$ are the electric displacement vector and magnetization respectively.  By assuming a medium with no free charges or currents, and replacing $D_i$ and $B_i$ with the relations $D_i = \kappa_o \kappa_r E_i$ and $B_i = \mu_o \mu_r H_i$,\footnote{It should be noted here that $\kappa_r$ and $\mu_r$ are the relative permittivity and permeability respectively.  It is important to remember that they are actually $\kappa = \kappa_r/\kappa_o$ and $\mu = \mu_r/\mu_o$, where the naught subscript is representative of the free space (vacuum) values. } we arrive at the simplified equations
\begin{align}
    \epsilon_{ijk} \frac{\partial E_k}{r_j} &= - \mu_o \mu_r \frac{\partial H_i}{dt} \\
    \epsilon_{ijk} \frac{\partial H_k}{r_j} &= \kappa_o \kappa_r \frac{\partial E_i}{\partial t}
\end{align}
where $\kappa_o$ and $\mu_o$ are the permittivity and permeability of vacuum and $\kappa_r$ and $\mu_r$ are the relative permittivity and permeability of the medium.  By combining these equations to eliminate $H_i$, and utilizing a useful vector identity, the final result is
\begin{equation}
    \frac{\partial^2 E_i}{\partial r_i^2} = \mu \kappa \frac{\partial^2 E_i}{dt^2}
\end{equation}
which has the same form as the wave equation and has a solution of the form
\begin{equation}\label{wave_prop}
    E_i (r_k, t) = E_i^0 e^{i(kr_k - \omega t)}
\end{equation}
where $k$ is the wave vector, $r_k$ is the direction of propagation, and $\omega$ is the angular frequency of the wave.  The phase velocity of the wave is given by
\begin{equation}\label{wave_velocity}
    \frac{1}{v^2} = \mu \kappa
\end{equation}
from which the speed of light in free space can be determined by setting $\kappa_r = \mu_r = 1$ to find that $c = 2.998 \times 10^8$ m/s.\footnote{The standard values of $\kappa_o = 8.854 \times 10^{-12}$ C/Vm and $\mu_o = 1.256 \times 10^{-6}$ N/A$^2$ were used to calculate this.}. From Eq.~\ref{wave_velocity} it is clear that the velocity of a wave in a medium is $v = c/ \sqrt{\kappa_r \mu_r}$ and therefore we can define that the index of refraction is
\begin{equation}
    n = \frac{c}{v}
\end{equation}
or more simply $n = \sqrt{\kappa}$ at optical frequencies (setting $\mu = 1$).  With this definition in place, as well the formula that defines the light wave propagation (Eq.~\ref{wave_prop}, we can carry on to developing the phenomenological description of optical properties in crystals.



First developed by Pockels, there are two assumptions to the phenomological theory prescribed to the topical properties of crystals and they are:
\begin{enumerate}
    \item In a homogeneously deformed solid the effect of deformation is only to alter the optical parameters of the optical indicatrix.
    \item When strain is within elastic limits, the changes of an optical parameter (polarization constant) due to deformation can be expressed as a homogeneous linear function of the nine strain components.
\end{enumerate}. The optical properties of crystals can be defined in terms of the refractive index ellipsoid, pictured above/below, and expressed as
\begin{equation}
     \frac{x^2}{n_x^2} +  \frac{y^2}{n_y^2} +  \frac{z^2}{n_z^2} = 1
\end{equation}
where x, y, and z are the principle axes, and $n_x$,$n_y$, and $n_z$ are the associated principal refractive indices.  The refractive index unfortunately does not transform like a tensor, however, the dielectric tensor $\kappa_{ij}$ does.  For an anistropic medium, the previous derivation is expanded by considering the tensorial properties of the medium; specifically the permittivity becomes $\kappa_{ij}$ with a directional dependence.  Since this second-rank tensor is symmetric (the \textit{i} and \textit{j} suffixes are interchangeable) it is possible to multiply it by the product of two coordinate vectors to create a quadric.\footnote{A quadric is a representation surface that can be used to describe any second-rank tensor, and hence any crystal property that is given by such a tensor.}. By using the relation between permittivity and index of refraction, a full tensorial description of the optical properties of a crystal can be defined.  Using all we know now, we can can define the optical indicatrix equation as
\begin{equation}
    B_{11}x^2 + B_{22}y^2 +B_{33}z^2 = 1,
\end{equation}
or more generally as
\begin{equation}\label{polarization_indicatrix}
    \sum_{i,j} B_{ij} r_i r_j = 1
\end{equation}
where $B_{ij} = 1/{n_{ij}^2 = 1 /\kappa_{ij}}$, and $r_i$ are the principle axes.  Also known as polarization constants, these parameters can be conveniently utilized to describe induced changes of the crystal optical properties.
