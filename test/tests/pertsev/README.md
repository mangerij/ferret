# `pertsev` — misfit-strain phase diagrams of (001) single-domain perovskite films

Seven tests, one per **phase** of the misfit–temperature phase diagram, for three materials
with progressively harder physics:

| test | material | T | u_m | phase | order parameters |
|---|---|---|---|---|---|
| `bto_c`  | BaTiO3 | 298.15 K | −0.005 | **c**  (P ∥ z) | P |
| `bto_ac` | BaTiO3 | 298.15 K |  0.000 | **ac** (P₁, P₃) | P |
| `bto_aa` | BaTiO3 | 298.15 K | +0.005 | **aa** (P₁ = P₂) | P |
| `pto_c`  | PbTiO3 | 298.15 K | −0.005 | **c** | P |
| `pto_r`  | PbTiO3 | 298.15 K | +0.006 | **r** (all three) | P |
| `sto_c`  | SrTiO3 | 100 K | −0.010 | **(p₃; q₃)** | P **and AFD tilt A** |
| `sto_a`  | SrTiO3 | 100 K | +0.010 | **(p₁; q₂)** | P **and AFD tilt A** |

References:
- BaTiO3, PbTiO3 — Pertsev, Zembilgotov & Tagantsev, *PRL* **80**, 1988 (1998), Fig. 1 and
  the parameter set in their footnote [17].
- SrTiO3 — Y. L. Li *et al.*, *PRB* **73**, 184112 (2006), Table VII. SrTiO3 carries the
  antiferrodistortive oxygen-octahedron tilt as a second order parameter, coupled to P
  biquadratically, which is what makes it the useful stress case for the multi-order-parameter
  machinery.

## What these tests actually assert

Each runs **4 master steps** (`num_steps = 4`) from a seed placed in the relevant well.
That makes them **floating-point regression tests of the multiapp machinery**, fast enough
for CI (~1 s each, 7.3 s for the suite) — they are *not* converged-physics tests.

The physics validation is recorded below instead: every point has been relaxed to
convergence and checked against an independent 0-D minimisation of the same free energy.

**To reproduce the physics, delete `num_steps` from the input.** The Terminator then relaxes
the point properly, and you should get:

| test | converged Ferret | 0-D reference | wall |
|---|---|---|---|
| `bto_c`  | P = (0.000000, 0.000000, 0.290051) | (0, 0, 0.29005) | 5.3 s |
| `bto_ac` | P = (0.146966, 0.000000, 0.212726) | (0.14697, 0, 0.21273) | 11.4 s |
| `bto_aa` | P = (0.221249, 0.221249, 0.000000) | (0.22125, 0.22125, 0) | 5.3 s |
| `pto_c`  | P = (0.000000, 0.000000, 0.686073) | (0, 0, 0.68607) | 2.2 s |
| `pto_r`  | P = (0.341851, 0.341851, 0.309278) | (0.34184, 0.34184, 0.30934) | 4.3 s |
| `sto_c`  | P = (0, 0, 0.080062), A = (0, 0, 12.4314) pm | matches | 2.6 s |
| `sto_a`  | P = (0.143513, 0, 0), A = (0, 4.17618, 0) pm | matches | 1.4 s |

P in C/m²; the AFD displacement A in **picometres**. Agreement is to ~1e-5 C/m² or better —
see `benchmarks_ferret/misfit_mono_BTO_multiapp` (400-point sweep, 397/400 phases, median
|ΔP| = 3.1e-5) and `benchmarks_ferret/misfit_mono_STO` (153/153 phases, max |ΔP| = 6.3e-5).

## Scripting the full phase diagram

Every input is self-contained and parameterised at the top. Delete `num_steps`, then sweep:

```bash
for T in 200 250 300 350 400; do
  for um in -0.010 -0.005 0.000 0.005 0.010; do
    ferret-opt -i BTO_pertsev_c.i T=$T um=$um Outputs/file_base=out_T${T}_um${um}
  done
done
```

Two things to get right when you do:

1. **`T` is in KELVIN in all seven inputs.** Pertsev's BaTiO3/PbTiO3 formulas are written in
   Celsius, so those masters carry `TC = ${fparse T - 273.15}` internally; SrTiO3's Barrett
   forms are natively in Kelvin. One knob, one unit, either way.
2. **Seed each candidate well and keep the lowest energy.** A single seed lands in whichever
   basin it starts in, and several of these wells are genuinely metastable — in the BaTiO3
   sweep a random seed chose `aa` over the true `ac` ground state at 8 of 400 points. Run
   each point from a `c`, `aa`, `ac` and `r` seed (`p0x/p0y/p0z`, plus `a0x/a0y/a0z` for
   SrTiO3) and compare `Ftotal`.

## Structure

Master + two children per material, following `tests/polar-elastic`:

```
<MAT>_pertsev_<phase>.i   master: owns the clock, parameters, Terminator, CSV
<MAT>_pertsev_polar.i     child : order parameters, ActuallyExplicitEuler (lumped)
<MAT>_pertsev_mech.i      child : one linear elastic solve per master step
```

Per master step: previous `u` → polar; polar sub-cycles `dt_mech/dt_polar` explicit steps;
fresh order parameters → mech; mech does one linear solve. The film is the film ALONE — no
substrate block. The misfit enters as a constant `global_strain` added to `total_strain`,
with `u` periodic in-plane and the top face free, which is exactly Pertsev's constraint set
(ε₁₁ = ε₂₂ = u_m, σ₃₃ = σ₁₃ = σ₂₃ = 0).

## Elastic energy with more than one order parameter

`CubicParentElasticEnergy` and `CubicParentElasticAEnergy` each compute
½(ε − ε⁰):C:(ε − ε⁰) using **only their own** eigenstrain. With two order parameters they
therefore do **not** sum to the elastic energy:

```
F_elP + F_elA  =  F_elastic  +  ½ ε:C:ε  −  (Q·PP):C:(R·AA)
```

They double-count ½ε:C:ε and miss the P–A cross term. So these inputs report

```
Felastic_true = ∫ ½ σ : ε_elastic        (ElasticEnergyAux + ElementIntegralVariablePostprocessor)
```

which is correct for **any** number of order parameters: `ComputeSmallStrain` subtracts every
eigenstrain in `eigenstrain_names`, so `stress` and `elastic_strain` already carry all of
them, with no cross term left to chase. For BaTiO3 and PbTiO3 (one order parameter) it
reproduces `CubicParentElasticEnergy` to the last bit, which is the check that it is right.

Use **`Ftotal`** to rank competing states. SrTiO3 additionally carries `Ftot`
(= Fbulk + Froto + Fcouple, order parameters only, no elastic term) for continuity with
`benchmarks_ferret/misfit_mono_STO`; it is **not** the total and must not be used to choose
between seeds — at `sto_c` the elastic term is 0.0282 against an `Ftot` of 0.0047.

The same asymmetry exists on the driving-force side: `CubicParentElasticPDerivative` and
`CubicParentElasticADerivative` each subtract only their own eigenstrain, so the
strain-mediated cross term is folded into `t1111/t1122/t1212` in the SrTiO3 inputs (see
`benchmarks_ferret/misfit_mono_STO/coeffs.py`). The general fix, if these get revisited, is
the same contraction the energy uses: ∂F_el/∂OPᵢ = −σ_jk ∂ε⁰_jk/∂OPᵢ, which needs no cross
term because σ already contains every eigenstrain.

## Conventions

- Ferret units: nm, kg, s, aC. An energy density is aJ/nm³ = 1e9 J/m³, so a coefficient is
  its SI value × 1e-9. SrTiO3's AFD displacement is carried in picometres, so each power of
  it brings a further 1e-12.
- `Q44` is entered as the **tensor** component (half the published engineering value):
  0.059/2 = 0.0295 for BaTiO3, 0.0675/2 = 0.03375 for PbTiO3, because
  `ComputeCubicParentElectrostrictiveStrain` fills the tensor shear.
- Stiffnesses are the inverse of Pertsev's published compliances: BaTiO3
  C11/C12/C44 = 175.549/84.639/108.225, PbTiO3 = 174.6032/79.3651/111.1111 GPa.
- `max_parallel = 1`: the BaTiO3/PbTiO3 polar children seed with `RandomIC`, whose values
  depend on the mesh partitioning.
