[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 5
    xmax = 0.004
    ymax = 0.004
    zmax = 0.02
    elem_type = HEX27
  []
  [pin_a]
    type = ExtraNodesetGenerator
    input = gen
    new_boundary = pin_a
    coord = '0 0 0'
  []
  [pin_b]
    type = ExtraNodesetGenerator
    input = pin_a
    new_boundary = pin_b
    coord = '0.004 0 0'
  []
[]

[GlobalParams]
  potential_E_int = phi
  disp_x = u_x
  disp_y = u_y
  disp_z = u_z
  displacements = 'u_x u_y u_z'
[]

[Variables]
  [u_x]
    order = SECOND
  []
  [u_y]
    order = SECOND
  []
  [u_z]
    order = SECOND
  []
  [phi]
    order = SECOND
  []
[]

[Kernels]
  [SolidMechanics]
  []
  [piezocouple_0]
    type = ConversePiezoelectricStrain
    variable = u_x
    component = 0
  []
  [piezocouple_1]
    type = ConversePiezoelectricStrain
    variable = u_y
    component = 1
  []
  [piezocouple_2]
    type = ConversePiezoelectricStrain
    variable = u_z
    component = 2
  []
  [FE_E_int]
    type = Electrostatics
    variable = phi
  []
  [strain_charge]
    type = PiezoelectricStrainCharge
    variable = phi
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '127.205e9 80.212e9 84.670e9 127.205e9 84.670e9 117.436e9 22.988e9 22.988e9 23.474e9'
  []
  [piezostrictive_tensor]
    type = ComputePiezostrictiveTensor
    fill_method = general
    e_ijk = '0 0 17.0345 0 0 0 17.0345 0 0  0 0 0 0 0 17.0345 0 17.0345 0  -6.62281 0 0 0 -6.62281 0 0 0 23.2403'
  []
  [permittivity]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.26934e-08'
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[BCs]
  [slide_z]
    type = DirichletBC
    variable = u_z
    boundary = back
    value = 0
  []
  [pin_ax]
    type = DirichletBC
    variable = u_x
    boundary = pin_a
    value = 0
  []
  [pin_ay]
    type = DirichletBC
    variable = u_y
    boundary = pin_a
    value = 0
  []
  [pin_by]
    type = DirichletBC
    variable = u_y
    boundary = pin_b
    value = 0
  []
  [ground]
    type = DirichletBC
    variable = phi
    boundary = back
    value = 0
  []
  [press]
    type = Pressure
    variable = u_z
    boundary = front
    factor = 625000
    use_displaced_mesh = false
  []
[]

[Postprocessors]
  [phi_front]
    type = SideAverageValue
    variable = phi
    boundary = front
  []
  [stroke]
    type = SideAverageValue
    variable = u_z
    boundary = front
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  automatic_scaling = true
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
  file_base = open_circuit_sensor_out
[]
