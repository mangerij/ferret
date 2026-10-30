# Spin waves from the eigensolver

Magnon modes of a magnetic film can be obtained directly as the eigenpairs of the linearized
Landau-Lifshitz-Gilbert (LLG) equation, rather than by exciting the system and transforming a
long time trace. This tutorial follows
[`test/tests/spin_waves/sw_eigen_film.i`](https://github.com/mangerij/ferret/blob/devel/test/tests/spin_waves/sw_eigen_film.i),
which computes the magnon bands of a perpendicularly magnetized film using MOOSE's `Eigenvalue`
executioner. It uses only the standard FERRET LLG kernels - no special eigenmode objects are
required.

!media media/sw_modes_de_k1.gif style=display:block;margin-left:auto;margin-right:auto;width:95%; caption=Damon-Eshbach spin waves travelling perpendicular to the magnetization in an in-plane film. Left: the lowest band, a surface mode. Right: the first perpendicular standing spin wave (PSSW). Arrows show the precessing magnetization (tilt exaggerated), coloured by precession phase. id=fig-sw-modes

## The idea

Linearizing the LLG equation about a static ground state $\mathbf{m}_0$ gives

\begin{equation}
  \begin{aligned}
    \mathsf{M} \frac{\partial \delta\mathbf{m}}{\partial t} = -\mathsf{J}\, \delta\mathbf{m},
  \end{aligned}
\end{equation}

where $\mathsf{J}$ is the Jacobian of the LLG kernels evaluated at $\mathbf{m}_0$ and $\mathsf{M}$
collects the time-derivative terms. Writing $\delta\mathbf{m} \sim e^{-\mu t}$ turns this into the
generalized eigenproblem

\begin{equation}
  \begin{aligned}
    \mathsf{J}\, \delta\mathbf{m} = \mu\, \mathsf{M}\, \delta\mathbf{m},
  \end{aligned}
\end{equation}

and for zero damping the eigenvalues $\mu = \pm i\omega$ come in purely imaginary pairs. The
magnon frequencies are therefore the imaginary parts of the eigenvalues.

MOOSE builds $\mathsf{J}$ for free, because it is the same Jacobian the Newton solver already
assembles. The only thing the input has to supply is $\mathsf{M}$.

## Supplying the mass matrix

The transverse components $m_x$ and $m_z$ are the dynamic variables, so $\mathsf{M}$ is built from
one `CoefReaction` per component, tagged into the `eigen` vector:

!listing test/tests/spin_waves/sw_eigen_film.i
         block=Kernels
         link=False
         language=python

The `extra_vector_tags = 'eigen'` entry is what places these terms on the right-hand side of the
generalized eigenproblem instead of in the residual. `MagHStrongCart` supplies the magnetostatic
source term, so the mode frequencies include the dipolar field.

The static physics comes from `film_common.i`, which is shared with the ring-down test and uses
only stock objects - `MasterExchangeCartLLG` for exchange, `MasterInteractionCartLLGHConst` for
the Zeeman term, and `Electrostatics` for the magnetostatic potential on a graded vacuum region
surrounding the film.

## The eigensolver

!listing test/tests/spin_waves/sw_eigen_film.i
         block=Executioner
         link=False
         language=python

A shift-and-invert spectral transform (`-st_type sinvert`) targeted at zero picks out the lowest
modes, which are the ones of physical interest. `initial_eigenvalue = 1e12` is required: the
default initial guess puts SLEPc in the wrong part of the spectrum and the solve fails to
converge. `MUMPS` is used for the inner factorization, and `n_eigen_pairs = 8` returns four
$\pm i\omega$ pairs.

## Reading the frequencies out

Eigenvalues arrive as a `VectorPostprocessor`, and individual bands are pulled out with
`VectorPostprocessorComponent`. Because the eigenvalues come in conjugate pairs, band $n$ sits at
index $2n$:

!listing test/tests/spin_waves/sw_eigen_film.i
         block=Postprocessors
         link=False
         language=python

The last two objects are worth highlighting. `rel_kittel` compares the $k = 0$ frequency against
the analytic Kittel result $\omega = \gamma_0 (H_0 - M_s)$, and a `Terminator` aborts the run if
the two disagree by more than one part in $10^{10}$:

!listing test/tests/spin_waves/sw_eigen_film.i
         block=UserObjects
         link=False
         language=python

The input therefore asserts its own correctness. This is the recommended pattern for FERRET
inputs that have an analytic limit - it turns a silent numerical error into a hard failure.

## Results

For in-plane magnetization and propagation perpendicular to $\mathbf{M}$, the modes are
Damon-Eshbach waves. The lowest band is the surface mode and the higher bands are perpendicular
standing spin waves quantized across the film thickness, as shown in [fig-sw-modes].

!media media/sw_film_de_bands.png style=display:block;margin-left:auto;margin-right:auto;width:70%; caption=Damon-Eshbach band structure of the film, compared against the Kalinikos-Slavin analytic dispersion. id=fig-sw-bands

## Related inputs

Three further inputs in the same directory build on this one:

- `sw_ringdown_film.i` excites the film and extracts the same frequencies from a time trace, as an
  independent check of the eigensolve.
- `sw_eigen_vs_ringdown.i` compares the two and fails if they disagree.
- `afm_skyr_modes.i` computes modes about a *non-uniform* ground state - an antiferromagnetic
  skyrmion - rather than about a uniformly magnetized state.

!content pagination previous=tutorials/piezoelectric.md next=tutorials/flexoelectric_domain_wall.md
