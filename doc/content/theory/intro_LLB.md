# Finite-temperature micromagnetics: the Landau-Lifshitz-Bloch equation

The LLG equation keeps $|\mathbf{m}| = 1$, so temperature can only enter through its parameters. Near the Curie temperature $T_c$, or when the temperature varies in space, the magnitude of the magnetization changes as well: domain walls become elliptic and then linear [!cite](Hinzke2008), and walls or skyrmions in a thermal gradient are pushed by the gradient of their own free energy [!cite](Schlickeiser2014). The Landau-Lifshitz-Bloch (LLB) equation [!cite](Garanin1997) describes these effects. It is valid below and above $T_c$. Reviews: [!cite](Atxitia2017), [!cite](Evans2012).

FERRET solves the LLB equation with the same kernels as LLG. Each torque kernel has a route switch `llb = true`, and [MasterLongitudinalLLB](MasterLongitudinalLLB.md) supplies the longitudinal field with `form = llb`. With the default settings every kernel behaves exactly as before. (The `form = penalty` term of [MasterLongitudinalLLB](MasterLongitudinalLLB.md), mentioned on the [micromagnetics page](intro_LLG.md), is a $T = 0$ norm penalty, not this equation.)

## Equation

The variable is $\mathbf{m} = \mathbf{M}/M_{s0}$, with $M_{s0}$ the saturation magnetization at $T = 0$. $|\mathbf{m}|$ is free and relaxes toward the equilibrium value $m_e(T)$. FERRET writes the LLB equation in the regular form

\begin{equation}\label{LLB}
  \frac{\partial \mathbf{m}}{\partial t} = -\gamma' \mathbf{m}\times\mathbf{H} + \gamma' \alpha_\perp \mathbf{H} - \gamma'\left(\alpha_\perp - \alpha_\parallel\right) \frac{(\mathbf{m}\cdot\mathbf{H})\,\mathbf{m}}{|\mathbf{m}|^2},
  \qquad \gamma' = \frac{g_0}{1+\lambda^2}.
\end{equation}

This is identical to Garanin's form $-\gamma'\mathbf{m}\times\mathbf{H} + \gamma'\alpha_\parallel (\mathbf{m}\cdot\mathbf{H})\mathbf{m}/|\mathbf{m}|^2 - \gamma'\alpha_\perp \mathbf{m}\times(\mathbf{m}\times\mathbf{H})/|\mathbf{m}|^2$. The only $1/|\mathbf{m}|^2$ multiplies $\alpha_\perp - \alpha_\parallel$, which vanishes at and above $T_c$. The parameter `llb_m2_eps` regularizes it at the centre of a linear wall, where $|\mathbf{m}| \to 0$.

The effective field $\mathbf{H}$ is the sum of the usual micromagnetic fields (exchange, DMI, anisotropy, Zeeman, demagnetization) and the longitudinal field

\begin{equation}
  \mathbf{H}_L = \begin{cases}
    \dfrac{1}{2\chi_\parallel}\left(1 - \dfrac{|\mathbf{m}|^2}{m_e^2}\right)\mathbf{m}, & T < T_c, \\[2ex]
    -\dfrac{1}{\chi_\parallel}\left(1 + \dfrac{3 T_c |\mathbf{m}|^2}{5 (T - T_c)}\right)\mathbf{m}, & T > T_c.
  \end{cases}
\end{equation}

$\mathbf{H}_L$ acts only through the $\alpha_\parallel$ term, so it relaxes $|\mathbf{m}|$ to $m_e$ on the time scale $\tau_\parallel = \chi_\parallel/(\gamma'\alpha_\parallel)$. This time is short: about 10 fs at 300 K for typical parameters.

### Temperature dependence goes into the material properties

| Property | Meaning |
|---|---|
| `Ms` | $M_{s0}$, a constant. Do not make it temperature dependent: $\mathbf{m}$ is already scaled by $M_{s0}$. |
| `me` | $m_e(T)$ |
| `chi_par` | $\chi_\parallel = \partial m/\partial H$ in m/A ($\mu_0$ times the usual $\tilde\chi_\parallel$ in 1/T) |
| `alpha_perp`, `alpha_par` | transverse and longitudinal damping. In mean-field theory $\alpha_\parallel = 2\lambda T/(3T_c)$, and $\alpha_\perp = \lambda(1 - T/(3T_c))$ below $T_c$; above $T_c$, $\alpha_\perp = \alpha_\parallel$. |
| `Ae`, `D`, `K1` (`K2`, `Kh`) | $A(T)/m_e^2$, $D(T)/m_e^2$, $K(T)/m_e^2$. The energies are written in $\mathbf{m}$, so these are the temperature-dependent constants divided by $m_e^2$. |
| `alpha` | $\lambda$, the coupling to the bath. In the LLB route it only sets the prefactor $g_0/(1+\lambda^2)$, that is, the time scale. The damping itself comes from `alpha_perp` and `alpha_par`. |

With a temperature variable `T`, each property can be a `PiecewiseLinearInterpolationMaterial` table in `T`, for example from a mean-field or atomistic calculation. These tables clamp outside their range, so cover every temperature the run can reach. Writing the constants this way is exact for a spatially varying $T(\mathbf{x})$. Unit-length LLG with $M_s(T(\mathbf{x}))$ is not: it adds a spurious $\nabla\cdot[(A/M_s)\nabla\mathbf{u}]$ force comparable to the entropic force itself.

In the LLB route the uniaxial anisotropy uses the transverse form $-K_1[(\mathbf{m}\cdot\mathbf{n})^2 - |\mathbf{m}|^2] - K_2[(\mathbf{m}\cdot\mathbf{n})^4 - |\mathbf{m}|^4]$. It has the same torque as the usual form but does not shift $|\mathbf{m}|$ inside the domains. The optional hard-axis term $K_h(\mathbf{m}\cdot\mathbf{n}_h)^2$ (properties `Kh`, `nhx`, `nhy`, `nhz`) makes biaxial walls with a Walker breakdown.

## Weak form and spatially varying damping

Every LLB kernel computes the residual $\mathcal{R}_c = \gamma' (\mathbf{H}, \mathbf{q})$, with the test vector

\begin{equation}
  \mathbf{q}(\mathbf{v}) = \mathbf{v}\times\mathbf{m} - \alpha_\perp \mathbf{v} + (\alpha_\perp - \alpha_\parallel)\frac{(\mathbf{v}\cdot\mathbf{m})\,\mathbf{m}}{|\mathbf{m}|^2}, \qquad \mathbf{v} = \psi\,\mathbf{e}_c .
\end{equation}

For gradient energies (exchange, DMI) the field contains $\nabla^2\mathbf{m}$ or $\nabla\mathbf{m}$, and integration by parts moves the derivative onto $\mathbf{q}$. When $\alpha_\perp$ and $\alpha_\parallel$ vary in space, $\nabla\mathbf{q}$ contains $\nabla\alpha_\perp$ and $\nabla\alpha_\parallel$. [MasterExchangeCartLLG](MasterExchangeCartLLG.md) and [MasterDMICartLLG](MasterDMICartLLG.md) include these terms when the material properties `grad_alpha_perp` and `grad_alpha_par` (`RealVectorValue`) are present. Supply them whenever the damping varies in space, in particular with $T(\mathbf{x})$. Without them the damping is treated as locally uniform. A wall then drifts at $v = -\gamma' (2A_e/\mu_0 M_{s0})\, \partial_x\alpha_\perp$, and a skyrmion settles about 7 degrees off the gradient axis. Stock MOOSE objects provide the gradients: `PiecewiseLinearInterpolationMaterial` declares the derivative `dalpha_perp/dT`, and `CoupledGradientMaterial` multiplies it by $\nabla T$:

```text
[Materials]
  [alpha_perp]                    # table in T; also declares dalpha_perp/dT
    type = PiecewiseLinearInterpolationMaterial
    property = alpha_perp
    variable = T
    x = '...'
    y = '...'
  []
  [grad_alpha_perp]
    type = CoupledGradientMaterial
    gradient_material_name = grad_alpha_perp
    coupled_variable = T
    scalar_property_factor = dalpha_perp/dT
  []
  # the same pair for alpha_par
[]
```

A prescribed $\alpha(\mathbf{x})$ can use `GenericConstantVectorMaterial` or `GenericFunctionVectorMaterial` for the gradients instead, as in the test below. The prefactor $g_0/(1+\lambda^2)$ and $M_{s0}$ are assumed uniform; a spatially varying `alpha` ($\lambda$) or `Ms` is not covered by these terms.

!listing test/tests/magnets/llb_alpha_ramp_wall.i block=Materials

## Stochastic LLB

The stochastic form sLLB-II of [!cite](Evans2012) adds a transverse thermal field $\boldsymbol{\eta}_\perp$, which enters only through the $\alpha_\perp$ damping term, and an additive longitudinal noise $\boldsymbol{\eta}_\parallel$:

- $\boldsymbol{\eta}_\perp$: a third instance of [MasterInteractionCartLLGHConst](MasterInteractionCartLLGHConst.md) with `thermal_field = true`, the material property `sigma_perp` and three noise user objects `noise_x`, `noise_y`, `noise_z`.
- $\boldsymbol{\eta}_\parallel$: the stock `LocalConservedLangevinNoise` kernel on each component, with `amplitude` $= \sigma_\parallel$ and `T_bath` set to an auxiliary variable equal to 1.
- The noise user objects are [LocalConservedNormalNoise](LocalConservedNormalNoise.md) with `conserved = false` (keep the uniform mode) and `white_noise = true` (space-time white noise, so the load vector has the fluctuation-dissipation covariance with the consistent mass matrix).

In FERRET units (nm, ns, A/m, $g_0$ in m/(A ns))

\begin{equation}
  \sigma_\perp^2 = \frac{2 k_B T (\alpha_\perp - \alpha_\parallel)(1+\lambda^2)}{g_0 \mu_0 M_{s0}\, \alpha_\perp^2 \cdot 10^{-27}}, \qquad
  \sigma_\parallel^2 = \frac{2 g_0 k_B T \alpha_\parallel}{(1+\lambda^2)\,\mu_0 M_{s0} \cdot 10^{-27}},
\end{equation}

where the factor $10^{-27}$ converts the cell volume to nm$^3$. Practical points:

- MOOSE reseeds every random object at each step with `seed + timestep`. Seeds of different noise objects must therefore differ by more than the number of time steps (the test uses $10^7$ apart). Otherwise the streams are copies of each other, shifted in time.
- Use `line_search = bt` for stochastic runs.
- An element is a macrospin only if `mag_x/y/z` are CONSTANT MONOMIAL, or if exchange locks the nodes together.

!listing test/tests/magnets/sllb_ensemble.i block=UserObjects

## Solver settings

- Do not use the renormalization user object or the tangent-plane integrators: $|\mathbf{m}|$ is a degree of freedom.
- Use `implicit-euler` or `bdf2`. `implicit-midpoint` is not L-stable, so the fast longitudinal mode rings when $\Delta t \gg \tau_\parallel$.
- Recommended: `nl_abs_tol = 1e-10`, `nl_rel_tol = 1e-8`, `nl_rel_step_tol = 1e-12`, `line_search = none` (deterministic), `nl_forced_its = 1`. On fine meshes the round-off floor of the residual exceeds `nl_abs_tol`; the step tolerance lets Newton stop there instead of cutting $\Delta t$.
- For static states, limit `dtmax`. A very large implicit step is a Newton solve for the static equations, and near a bifurcation (for example the elliptic-to-linear wall transition) it can land on the unstable branch.

## Regression tests

All are in `test/tests/magnets/`; each asserts its physics with a `Terminator`, so a wrong result stops the run.

The reduced-unit tests use $g_0 = \mu_0 = M_s = A_e = K_1 = 1$.

| Test | Setup | Asserted (tolerance) |
|---|---|---|
| `llb_tau_par` | macrospin, classical MFA, $T_c = 660$ K, $T = 620$ K | $(|\mathbf{m}| - m_e)/(|\mathbf{m}_0| - m_e) = e^{-t/\tau_\parallel}$, $\tau_\parallel = \chi_\parallel/(\gamma'\alpha_\parallel) = 7.510184\times 10^{-4}$ ns ($2\times 10^{-3}$) |
| `llb_tau_par_above` | as `llb_tau_par`, $T = 700$ K | the same relaxation, $\tau_\parallel = 1.467862\times 10^{-3}$ ns ($2\times 10^{-3}$) |
| `llb_cone` | macrospin at 300 K, $K_1 = 2\times 10^5$, $K_2 = -2\times 10^5$ J/m$^3$ (LLB inputs $K_1/m_e^2$, $K_2/m_e^4$) | the minimum of the LLB free energy (anisotropy + longitudinal): cone angle $45.005576^\circ$ and $|\mathbf{m}|/m_e = 1.0000973$ |
| `llb_helix` | bulk-DMI helix at 600 K, one period $2\pi/q$, $q = D/(2A_e)$, periodic | exchange + DMI energy density $-D^2|\mathbf{m}|^2/(4A_e)$ ($10^{-4}$); $|\mathbf{m}|$ pulled above $m_e$ to $m^2 = m_e^2\,(1 + \chi_\parallel D^2/(\mu_0 M_{s0} A_e))$, shift $6.010529\times 10^{-4}$ (2 %, the mesh error is 0.9 %) |
| `llb_demag_film` | 1-D film at 600 K along its normal $x$, vacuum on both sides, natural $\phi$ BCs, precessing in an in-plane field along $z$ | in the film, once the transient has decayed: $H_d = -\partial_x\phi = -M_{s0} m_x$ ($10^{-6}$) and the para-process $|\mathbf{m}| = m_e + \chi_\parallel H$ ($2\times 10^{-6}$); the gold CSV holds the Kittel ring-down |
| `llb_biaxial_chi` | uniform LLB state, reduced units, $m_e = 0.8$, easy $x$, hard $z$ with $K_h = 0.5$, $H_y = H_z = 0.1$ | the hard-axis term has no $|\mathbf{m}|^2$ part, so $|\mathbf{m}| = m_e$, $m_y = \mu_0 M_s H_y/(2K_1) = 0.05$, $m_z = \mu_0 M_s H_z/(2(K_1 + K_h)) = 0.0333\ldots$ ($10^{-10}$, at $t = 90$) |
| `walker_llg` | LLG 1-D wall, reduced units, easy $x$, hard $z$, $K_h = 0.5$, $\alpha = 0.1$, $H = 0.03 < H_W = \alpha K_h/(\mu_0 M_s) = 0.05$, started as the Walker profile | at $t = 5$, the exact Walker solution: $\sin 2\varphi = H/H_W$, $\varphi = 18.4349^\circ$ ($0.2^\circ$), and $v = g_0\Delta(\varphi)H/\alpha = 0.292770$ with $\Delta(\varphi) = \sqrt{A_e/(K_1 + K_h\sin^2\varphi)}$ ($2\times 10^{-3}$) |
| `walker_llg_above` | as `walker_llg`, $H = 0.1 > H_W$ | the wall plane precesses at $\dot\varphi = (H - H_W\sin 2\varphi)/(1+\alpha^2)$ at every step up to $t = 20$ (1 %); the rate is positive at every $\varphi$, which proves the breakdown |
| `llb_alpha_ramp_wall` | LLB Bloch wall, reduced units, $m_e = 0.8$, easy $z$, uniform $T$, $\alpha_\perp = 0.25 - 0.002x$, $\alpha_\parallel = 0.15 + 0.002x$ | a static LLB state does not depend on $\alpha(x)$: after $t = 10$, $|q| < 0.05\, v_\mathrm{sp} t$ at every step. Without `grad_alpha_*` the wall drifts at $v_\mathrm{sp} = -2\,\partial_x\alpha_\perp/(1+\lambda^2) = 0.003960$; the start-up ring-down leaves a displacement $|q| \sim 10^{-3}$, not a drift |
| `llb_alpha_ramp_dmi` | as `llb_alpha_ramp_wall`, Néel wall with interfacial $D = 0.3$ | the same bound, now for the exchange + DMI weak form |
| `sllb_ensemble` | sLLB-II, 200 independent macrospins (`CONSTANT MONOMIAL` $\mathbf{m}$, 5.3 nm cells) at 620 K | after $15\,\tau_\parallel$: ensemble $\langle|\mathbf{m}|\rangle$ within 0.008 of the Boltzmann value 0.31308, and its spread within 20 % of 0.02159 |

Beyond the tests, the implementation reproduces:

- the $T \to 0$ LLG limit;
- the Bloch, elliptic and linear walls of [!cite](Hinzke2008), including the transition temperature;
- the entropic wall velocities of [!cite](Schlickeiser2014);
- the Thiele motion of a skyrmion in a thermal gradient (the Hall ratio to 0.4 %, and the force equal to $-\partial F_\mathrm{sk}/\partial T\, \nabla T$ to about 1 %).
