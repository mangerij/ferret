
# Core-shell nanoparticle, 25 nm in diameter:
# Single-crystal round ZnO core, 15 nm in diameter.
# Single-crystal rutile TiO2 shell, 5 nm in thickness.
# Here, crystalline rTiO2 elastic parameters (both bulk and surface ones) in the shell were
# averaged out to isotropic symmetry, which effectively makes the rTiO2 shell amorphous.
# See Influence of Elastic and Surface Strains on the Optical Properties of Semiconducting Core-Shell Nanoparticles
# J. Mangeri, O. Heinonen, D. Karpeyev, and S. Nakhmanson, Phys. Rev. Applied 4, 014001 (2015)


[Mesh]
  file = core_shell_exodus.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]

  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./EgZnO]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = false
  [../]
[]


[AuxKernels]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx_elastic
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy_elastic
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz_elastic
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy_elastic
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz_elastic
    index_i = 1
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx_elastic
    index_i = 2
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx_elastic
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy_elastic
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz_elastic
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy_elastic
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz_elastic
    index_i = 1
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx_elastic
    index_i = 2
    index_j = 0
    use_displaced_mesh = false
  [../]
[]

[BCs]

# Surface elastic parameters for w-ZnO surface taken from Table III of JAP 111, 124305 (2012)
# and averaged out to represent an isotropic surface described by one diagonal elastic const C_1111,
# one shear elastic const C_1122 and one residual stress const tau_11.
#
# Original (anisotropic, w-ZnO 10-10 surface) elastic consts and residual stresses in N/m (PBE):
# C_1111 = 49.1, C_2222 = 34.9, C_1122 = 15.1, C_1212 = 13.7, tau_11 = -2.1, tau_22 = -1.4
#
# Isotropic C_1111 = [C_1111 + C_2222]/2 = 42.0, C_1122 = 15.0, C_1212 is computed automatically
# from C_1111 and C_1122, tau_11 = [tau_11 + tau_22]/2 = -1.7
# Conversion to [N/nm] multiplies everything by 10^{-9}

  [./surface_elasticity_X]
    type = SurfaceMechanicsBC
    variable = disp_x
    boundary = '1'
    S_k = '0.0 0.0 0.0'
    # Surface elastic tensor C_1111, C_1122
    Cs_ijkl = '42.0e-09 15.0e-09'
    # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]

  [./surface_elasticity_Y]
    type = SurfaceMechanicsBC
    variable = disp_y
    boundary = '1'
    S_k = '0.0 0.0 0.0'
    # Surface elastic tensor C_1111, C_1122
    Cs_ijkl = '42.0e-09 15.0e-09'
    # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]
  [./surface_elasticity_Z]
    type = SurfaceMechanicsBC
    variable = disp_z
    boundary = '1'
    S_k = '0.0 0.0 0.0'
    # Surface elastic tensor C_1111, C_1122
    Cs_ijkl = '42.0e-09 15.0e-09'
    # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]
[]


[Materials]
    # ZnO shell (full [hexagonal] crystalline symmetry)
    # Single crystal (see Table 1.6 and Ref. 86) in "Zinc Oxide: Fundamentals, Materials and Device Technology"
    # by Hadis Morkoc and Umit Ozgur (2009 WILEY-VCH Verlag GmbH & Co. KGaA, Weinheim ISBN: 978-3-527-40813-9), Chapter 1.
    # In GPa: C_ijkl = 209.7 121.1 105.1 209.7 105.1 210.9 42.47 42.47 44.29
    # In N/(nm)^2, 1 GPa = 1 * 10^{-9} N/(nm)^2:


    [./elasticity_tensor1]
      type = ComputeElasticityTensor
      block = '1'
      fill_method = symmetric9
      C_ijkl = '209.7e-09 121.1e-09 105.1e-09 209.7e-09 105.1e-09 210.9e-09 42.47e-09 42.47e-09 44.29e-09'
      euler_angle_1 = 0.0
      euler_angle_2 = 0.0
      euler_angle_3 = 0.0
    [../]
    [./strain1]
      type = ComputeSmallStrain
      block = '1'
    [../]
    [./stress1]
      type = ComputeLinearElasticStress
      block = '1'
    [../]

    # Zn core (full [hexagonal] crystalline symmetry)
    # Averaged constants for crystalline Zn, see Table 2 in H. M. Ledbetter, J. Phys. Chem. Ref. Data 6, 1181 (1977).
    # In GPa: C_ijkl = 163.0 30.6 48.1 163.0 48.1 60.3 39.4 39.4 65.9
    # In N/(nm)^2, 1 GPa = 1 * 10^{-9} N/(nm)^2:

    [./elasticity_tensor2]
       type = ComputeElasticityTensor
       block = '2'
       fill_method = symmetric9
       C_ijkl = '163.0e-09 30.6e-09 48.1e-09 163.0e-09 48.1e-09 60.3e-09 39.4e-09 39.4e-09 65.9e-09'
       euler_angle_1 = 0.0
       euler_angle_2 = 0.0
       euler_angle_3 = 0.0
    [../]
    [./strain2]
      type = ComputeSmallStrain
      block = '2'
    [../]
    [./stress2]
      type = ComputeLinearElasticStress
      block = '2'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_type -pc_type  -snes_rtol -ksp_rtol'
    petsc_options_value = '    gmres    bjacobi     1e-4       1e-4'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'       #"PJFNK, JFNK, NEWTON"
  #PJFNK IS SLOW ~1500 seconds but NEWTON introduces the rigid body modes and causes many nonlinear iterations to show up
[]

[Outputs]
  [./Exodus]
    type = Exodus
    file_base = out_Zn_ZnO_xstl000_core_iso_shell
    elemental_as_nodal = true
  [../]
[]
