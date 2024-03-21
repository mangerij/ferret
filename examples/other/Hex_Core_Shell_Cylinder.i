
# Core-shell nanoparticle, 25 nm in diameter:
# Single-crystal round ZnO core, 15 nm in diameter.
# Single-crystal rutile TiO2 shell, 5 nm in thickness.
# Here, crystalline rTiO2 elastic parameters (both bulk and surface ones) in the shell were
# averaged out to isotropic symmetry, which effectively makes the rTiO2 shell amorphous.


[Mesh]
  file = ZnO_Core_100_Shell_20_nm.e
  # uniform_refine = 1
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

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zx]
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
  [./SolidMechanics]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    use_displaced_mesh = false
  [../]
[]


[AuxKernels]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
    use_displaced_mesh = false
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx
    index_i = 2
    index_j = 0
    use_displaced_mesh = false
  [../]

### Here we use the results from Phys Rev B 88, 235210 (2013) Wagner et al
# which computes the strain-induced bandgap change of wurtzite ZnO using
# the HSE approach. Here E0, db, du, Rb and nu are material properties.
# Units will be in eV. Note we are justified in the strain_xx = strain_yy(using histo2d)
# Here we consider Gamma_7c-Gamma-7v as an example
  [./bandgap]
    type = BandGapAuxZnO
    variable = EgZnO
    relaxed_energy = 3.200
    uniaxial_strain_rate = -3.800
    biaxial_strain_rate = -0.450
    biaxial_relaxation_coeff = 0.929
    poisson_ratio = 0.31
  [../]
[]

[BCs]

# Positive pressure -- compression, negative pressure -- tension
  #active = 'back_x_bc back_y_bc back_z_bc top_z_bc'

  [./back_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 1
     value = 0
   [../]
   [./back_y_bc]
     type = DirichletBC
     variable = disp_y
     boundary = 1
     value = 0
   [../]

  [./back_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 8
    value = -15
  [../]


#Load boundary condition (-.1 um displacement along z)

#  [./top_x_bc]
#   type = DirichletBC
#    variable = disp_x
#    boundary = 2
#    value = 0
#  [../]
#
#  [./top_y_bc]
#    type = DirichletBC
#   variable = disp_y
#    boundary = 2
#    value = 0
#  [../]

  [./top_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 7
    value = 15                  #Edit this value#
 [../]

### Should I also add in the boundary conditions for my other surfaces?  According to Cubit, Surface 2
### does not exist.  Do I just go along with what it's telling me?

# Surface elastic parameters for w-ZnO surface taken from Table III of JAP 111, 124305 (2012)
# and averaged out to represent an isotropic surface described by one diagonal elastic const C_1111,
# one shear elastic const C_1122 and one residual stress const tau_11.
#
# Original (anizotropic, w-ZnO 10-10 surface) elastic consts and residual stresses in N/m (PBE):
# C_1111 = 49.1, C_2222 = 34.9, C_1122 = 15.1, C_1212 = 13.7, tau_11 = -2.1, tau_22 = -1.4
#
# Isotropic C_1111 = [C_1111 + C_2222]/2 = 42.0, C_1122 = 15.0, C_1212 is computed automatically
# from C_1111 and C_1122, tau_11 = [tau_11 + tau_22]/2 = -1.7
# Conversion to [N/nm] multiplies everything by 10^{-9}


 [./surface_elasticity_X_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '3'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]

  [./surface_elasticity_Y_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '3'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]

  [./surface_elasticity_Z_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '3'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]

 [./surface_elasticity_X_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '4'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]

  [./surface_elasticity_Y_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '4'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]

  [./surface_elasticity_Z_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '4'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]

 [./surface_elasticity_X_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '5'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]

  [./surface_elasticity_Y_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '5'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]

  [./surface_elasticity_Z_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '5'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]

  [./surface_elasticity_X_surf4]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_x
     boundary = '6'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 0
   [../]

   [./surface_elasticity_Y_surf4]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_y
     boundary = '6'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 1
   [../]

   [./surface_elasticity_Z_surf4]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_z
     boundary = '6'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 2
   [../]

  [./surface_elasticity_X_surf5]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_x
     boundary = '7'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 0
   [../]

   [./surface_elasticity_Y_surf5]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_y
     boundary = '7'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 1
   [../]

   [./surface_elasticity_Z_surf5]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_z
     boundary = '7'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 2
   [../]

  [./surface_elasticity_X_surf6]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_x
     boundary = '8'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 0
   [../]

   [./surface_elasticity_Y_surf6]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_y
     boundary = '8'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 1
   [../]

   [./surface_elasticity_Z_surf6]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_z
     boundary = '8'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 2
   [../]

   [./surface_elasticity_X_surf7]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_x
      boundary = '9'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
      taus = '-1.7e-09'
      component = 0
    [../]

    [./surface_elasticity_Y_surf7]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_y
      boundary = '9'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
      taus = '-1.7e-09'
      component = 1
    [../]

    [./surface_elasticity_Z_surf7]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_z
      boundary = '9'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
      taus = '-1.7e-09'
      component = 2
    [../]

   [./surface_elasticity_X_surf8]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_x
      boundary = '10'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
      taus = '-1.7e-09'
      component = 0
    [../]

    [./surface_elasticity_Y_surf8]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_y
      boundary = '10'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
      taus = '-1.7e-09'
      component = 1
    [../]

    [./surface_elasticity_Z_surf8]
      type = SurfaceMechanicsBC
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      variable = disp_z
      boundary = '10'
      surface_euler_angle_1 = 0.0
      surface_euler_angle_2 = 0.0
      surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
      Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
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
     petsc_options = '-snes_view -snes_converged_reason -ksp_converged_reason'
     petsc_options_iname = '-ksp_type -snes_rtol -pc_type  -pc_hypre_type '
     petsc_options_value = '   gmres    1e-9     hypre      boomeramg'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  [./Exodus]
    type = Exodus
    file_base = Out_ZnO_Core_100_Shell_20_nm
    elemental_as_nodal = true
    output_elemental_variables = 1
    output_initial = 0
  [../]
[]
