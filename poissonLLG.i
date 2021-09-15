[Mesh]
  file = sphere_tet_approx_size0_05.e
  dim = 3
[]

[Variables]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./mag_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Materials]
  [./constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' alpha    Ae    Ms    g0     mu0   nx ny nz   long_susc t'
    prop_values = '0.01    0.013  1.2  221010.0 1256.64   1  0  0          1.0     0'
  [../]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.0'
  [../]
[]

[Kernels]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCart
    variable = potential_H_int
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
  [../]
[]

[Executioner]
  type = Steady
  nl_abs_tol = 1e-8
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
  execute_on = timestep_end
[]
