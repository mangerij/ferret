
[Mesh]
  file = exodus_disk_r8_h1.e
[]


[GlobalParams]
  # Elastic interactions are turned off for now.
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z

  displacements = 'disp_x disp_y disp_z'

  #  the negative prefactor should give compression now
  prefactor = -0.02

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

[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]
[]


[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
[]

[Postprocessors]
    [./Felastic]
      type = ElasticEnergy
      block = '1'
      execute_on = 'initial timestep_end'
    [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    120               2e-10      1e-8     1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = false
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_nanodisk_skyrm2PE
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_nanodisk_skyrm2PE
    execute_on = 'timestep_end'
  [../]
[]
