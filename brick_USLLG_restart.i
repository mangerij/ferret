
[Mesh]
  file = out_USLLG_test2_0.e
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  magnetic_x = magnetic_x
  magnetic_y = magnetic_y
  magnetic_z = magnetic_z

  alpha = 0.5 #this might need to be negative due to a typo in the original document

  Ae = 0.013
  Ms = 0.8

  g0 = 176.1

  permittivity = 1.0
  mu0 = 1256.0

[]

#theta = 0  is maybe the problem!!!

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    initial_from_file_var = polar_theta
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    initial_from_file_var = azimuth_phi
  [../]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    initial_from_file_var = potential_H_int
  [../]
  [./potential_H_ext]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
[]

[AuxVariables]
  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./Hexch_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hexch_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hexch_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./Hdemag_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hdemag_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hdemag_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]

 [./Hexch_x_c]
   type = ExchangeFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hexch_y_c]
   type =  ExchangeFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hexch_z_c]
   type =  ExchangeFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '1'
 [../]

 [./Hdemag_x_c]
   type = DemagFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hdemag_y_c]
   type =  DemagFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hdemag_z_c]
   type =  DemagFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeUSLL
    component = 0
    variable = polar_theta
    block = '1'
  [../]
  [./azimuthal_time]
    type = TimeUSLL
    component = 1
    variable = azimuth_phi
    block = '1'
  [../]

   #LLG simple

  # Exchange term

  [./d_llg_exch_th]
    type = ExchangeUSLL
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_exch_phi]
    type = ExchangeUSLL
    variable = azimuth_phi
    component = 1
  [../]

  # Anisotropy term. TURNED OFF => soft magnet


  # Magnetic interaction term

  [./d_HM_0]
    type = InteractionUSLL
    variable = polar_theta
    component = 0
  [../]
  [./d_HM_1]
    type = InteractionUSLL
    variable = azimuth_phi
    component = 1
  [../]


  # Magnetostatic Poisson equation

  [./ext_pot_lap]
    type = Electrostatics
    variable = potential_H_ext
    block = '1 2'
  [../]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '1 2'
  [../]
  [./int_bc_pot_lap]
    type = MagHStrong
    variable = potential_H_int
    block = '1'
  [../]
[]

[BCs]
   [./bc_int_pot_R]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1 2 3 4 5 6'
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
  [./Fdemag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Fdemag'
    pp_coefs = ' 1 1 ' 
    execute_on = 'initial timestep_end'
  [../]

  [./res]
    type = Residual
    residual_type = INITIAL_BEFORE_PRESET
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10      1e-8      1e-8       bjacobi'
  [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 1.0e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 1.2
    cutback_factor = 0.4
    dt = 1.0e-7
  [../]
  verbose = true
[]


[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_USLLG_test2_R0
    interval = 1
    elemental_as_nodal = true
  [../]

#0: with sqrt(1-cos^2x) // newton implc, ksp
#
#   The time to solve seems to be off. Perhaps second term in dMk/dt = - γ' Mk × δF/δMk - γ' α [Mk × (Mk × δF/δMk] needds 1/M_s thus making the damping coefficient 1/0.8 stronger...

#1: then with larger disk (*test2_0.e)
#
#   had to do a restart (*test2_R0.e)
# 
#2: next could try to flip sign on Poisson equation (this will invert the vectors at the <near> point and perhaps cause vortex nucleation).
#3:
#4:       
#5:

[]
