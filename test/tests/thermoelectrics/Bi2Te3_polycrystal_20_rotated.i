[Mesh]
  [./basic_mesh]
    type = FileMeshGenerator
    file = 20g.e
  []

  [./add_sidesets]
    type = SideSetsFromNormalsGenerator
    input = basic_mesh
    normals = '1  0  0
              -1  0  0
               0  1  0
               0 -1  0
               0  0  1
               0  0 -1'
    fixed_normal = true
    new_boundary = 'right left front back top bottom'
    normal_tol = 0.5
  [../]
[]

[Variables]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
    initial_condition = 273
  [../]
  []

[Kernels]
  [./residualV_x]
    type = TensorDivCurrentV
    component = 0
    variable = potential_E_int
    T = T
    potential_E_int = potential_E_int
  [../]

  [./residualT_x]
    type = TensorHeatFlowElectricT
    component = 0
    variable = T
    T = T
    potential_E_int = potential_E_int
  [../]
[]

[AuxVariables]
  [./j_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./q_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  []

[AuxKernels]
  [./Electric_flux_x]
    type = ElectricFluxTensor
    variable = j_x
    T = T
    potential_E_int = potential_E_int
    component = 0
  [../]
  [./Electric_flux_y]
    type = ElectricFluxTensor
    variable = j_y
    T = T
    potential_E_int = potential_E_int
    component = 1
  [../]
  [./Electric_flux_z]
    type = ElectricFluxTensor
    variable = j_z
    T = T
    potential_E_int = potential_E_int
    component = 2
  [../]

  [./Heat_flux_x]
    type = HeatFluxTensor
    variable = q_x
    T = T
    potential_E_int = potential_E_int
    component = 0
  [../]
  [./heat_flux_y]
    type = HeatFluxTensor
    variable = q_y
    T = T
    potential_E_int = potential_E_int
    component = 1
  [../]
  [./Heat_flux_z]
    type = HeatFluxTensor
    variable = q_z
    T = T
    potential_E_int = potential_E_int
    component = 2
  [../]
[]


[Materials]
  [./thermal_conductivity_tensor1]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'#Bi2Te3-mp-568390
    euler_angle_1 = 0.0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block  = '1'
  [../]
  [./thermal_conductivity_tensor2]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 0.0
    euler_angle_2 = 15
    euler_angle_3 = 30
    block  = '2'
  [../]
  [./thermal_conductivity_tensor3]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 0.0
    euler_angle_2 = -15
    euler_angle_3 = -30
    block  = '3'
  [../]
  [./thermal_conductivity_tensor4]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 0
    euler_angle_2 = 45
    euler_angle_3 = 60
    block  = '4'
  [../]
  [./thermal_conductivity_tensor5]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'#Bi2Te3-mp-568390
    euler_angle_1 = 0
    euler_angle_2 = -45
    euler_angle_3 = -60
    block  = '5'
  [../]
  [./thermal_conductivity_tensor6]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 15
    euler_angle_2 = 30
    euler_angle_3 = 60
    block  = '6'
  [../]
  [./thermal_conductivity_tensor7]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = -15
    euler_angle_2 = -30
    euler_angle_3 = -60
    block  = '7'
  [../]
  [./thermal_conductivity_tensor8]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 20
    euler_angle_2 = 40
    euler_angle_3 = 80
    block  = '8'
  [../]
  [./thermal_conductivity_tensor9]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'#Bi2Te3-mp-568390
    euler_angle_1 = -20
    euler_angle_2 = -40
    euler_angle_3 = -80
    block  = '9'
  [../]
  [./thermal_conductivity_tensor10]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 30
    euler_angle_2 = 45
    euler_angle_3 = 60
    block  = '10'
  [../]
  [./thermal_conductivity_tensor11]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = -30
    euler_angle_2 = -45
    euler_angle_3 = -60
    block  = '11'
  [../]
  [./thermal_conductivity_tensor12]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 45
    euler_angle_2 = 60
    euler_angle_3 = 90
    block  = '12'
  [../]
  [./thermal_conductivity_tensor13]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'#Bi2Te3-mp-568390
    euler_angle_1 = -45
    euler_angle_2 = -60
    euler_angle_3 = -90
    block  = '13'
  [../]
  [./thermal_conductivity_tensor14]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 50
    euler_angle_2 = 75
    euler_angle_3 = 100
    block  = '14'
  [../]
  [./thermal_conductivity_tensor15]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = -50
    euler_angle_2 = -75
    euler_angle_3 = -100
    block  = '15'
  [../]
  [./thermal_conductivity_tensor16]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 10
    euler_angle_2 = 20
    euler_angle_3 = 30
    block  = '16'
  [../]
  [./thermal_conductivity_tensor17]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'#Bi2Te3-mp-568390
    euler_angle_1 = -10
    euler_angle_2 = -20
    euler_angle_3 = -30
    block  = '17'
  [../]
  [./thermal_conductivity_tensor18]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 30
    euler_angle_2 = 60
    euler_angle_3 = 90
    block  = '18'
  [../]
  [./thermal_conductivity_tensor19]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = -30
    euler_angle_2 = -60
    euler_angle_3 = -90
    block  = '19'
  [../]
  [./thermal_conductivity_tensor20]
    type = ComputeThermalConductivityTensor
    k_ij = '2.58 2.4e-17 0.0077 2.4e-17 2.523 1.8e-17 0.0077 1.8e-17 2.524'
    euler_angle_1 = 90
    euler_angle_2 = 90
    euler_angle_3 = 0
    block  = '20'
  [../]
  [./electrical_conductivity_tensor1]
    type = ComputeElectricalConductivityTensor
    g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
    euler_angle_1 = 0.0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block  = '1'
  [../]
   [./electrical_conductivity_tensor2]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = 15
     euler_angle_3 = 30
     block  = '2'
   [../]
   [./electrical_conductivity_tensor3]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = -15
     euler_angle_3 = -30
     block  = '3'
   [../]
   [./electrical_conductivity_tensor4]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = 45
     euler_angle_3 = 60
     block  = '4'
   [../]
   [./electrical_conductivity_tensor5]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = -45
     euler_angle_3 = -60
     block  = '5'
   [../]
   [./electrical_conductivity_tensor6]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 15
     euler_angle_2 = 30
     euler_angle_3 = 60
     block  = '6'
   [../]
   [./electrical_conductivity_tensor7]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -15
     euler_angle_2 = -30
     euler_angle_3 = -60
     block  = '7'
   [../]
   [./electrical_conductivity_tensor8]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 20
     euler_angle_2 = 40
     euler_angle_3 = 80
     block  = '8'
   [../]
   [./electrical_conductivity_tensor9]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -20
     euler_angle_2 = -40
     euler_angle_3 = -80
     block  = '9'
   [../]
   [./electrical_conductivity_tensor10]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 30
     euler_angle_2 = 45
     euler_angle_3 = 60
     block  = '10'
   [../]
   [./electrical_conductivity_tensor11]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -30
     euler_angle_2 = -45
     euler_angle_3 = -60
     block  = '11'
   [../]
   [./electrical_conductivity_tensor12]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 45
     euler_angle_2 = 60
     euler_angle_3 = 90
     block  = '12'
   [../]
   [./electrical_conductivity_tensor13]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -45
     euler_angle_2 = -60
     euler_angle_3 = -90
     block  = '13'
   [../]
   [./electrical_conductivity_tensor14]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 50
     euler_angle_2 = 75
     euler_angle_3 = 100
     block  = '14'
   [../]
   [./electrical_conductivity_tensor15]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -50
     euler_angle_2 = -75
     euler_angle_3 = -100
     block  = '15'
   [../]
   [./electrical_conductivity_tensor16]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 10
     euler_angle_2 = 20
     euler_angle_3 = 30
     block  = '16'
   [../]
   [./electrical_conductivity_tensor17]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -10
     euler_angle_2 = -20
     euler_angle_3 = -30
     block  = '17'
   [../]
   [./electrical_conductivity_tensor18]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 30
     euler_angle_2 = 60
     euler_angle_3 = 90
     block  = '18'
   [../]
   [./electrical_conductivity_tensor19]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = -30
     euler_angle_2 = -60
     euler_angle_3 = -90
     block  = '19'
   [../]
   [./electrical_conductivity_tensor20]
     type = ComputeElectricalConductivityTensor
     g_ij = '2.62e4 123e-15 432 123e-15 2.29e4 61.03e-15 432 61.03e-15 2.3e4'#mp-568390
     euler_angle_1 = 90
     euler_angle_2 = 90
     euler_angle_3 = 0
     block  = '20'
   [../]



   [./seebeck_tensor1]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = 0
     euler_angle_3 = 0
     block  = '1'
   [../]
   [./seebeck_tensor2]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = 15
     euler_angle_3 = 30
     block  = '2'
   [../]
   [./seebeck_tensor3]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = -15
     euler_angle_3 = -30
     block  = '3'
   [../]
   [./seebeck_tensor4]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = 45
     euler_angle_3 = 60
     block  = '4'
   [../]
   [./seebeck_tensor5]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 0.0
     euler_angle_2 = -45
     euler_angle_3 = -60
     block  = '5'
   [../]
   [./seebeck_tensor6]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 15
     euler_angle_2 = 30
     euler_angle_3 = 60
     block  = '6'
   [../]
   [./seebeck_tensor7]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -15
     euler_angle_2 = -30
     euler_angle_3 = -60
     block  = '7'
   [../]
   [./seebeck_tensor8]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 20
     euler_angle_2 = 40
     euler_angle_3 = 80
     block  = '8'
   [../]
   [./seebeck_tensor9]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -20
     euler_angle_2 = -40
     euler_angle_3 = -80
     block  = '9'
   [../]
   [./seebeck_tensor10]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 30
     euler_angle_2 = 45
     euler_angle_3 = 60
     block  = '10'
   [../]
   [./seebeck_tensor11]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -30
     euler_angle_2 = -45
     euler_angle_3 = -60
     block  = '11'
   [../]
   [./seebeck_tensor12]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 45
     euler_angle_2 = 60
     euler_angle_3 = 90
     block  = '12'
   [../]
   [./seebeck_tensor13]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -45
     euler_angle_2 = -60
     euler_angle_3 = -90
     block  = '13'
   [../]
   [./seebeck_tensor14]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 50
     euler_angle_2 = 75
     euler_angle_3 = 100
     block  = '14'
   [../]
   [./seebeck_tensor15]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -50
     euler_angle_2 = -75
     euler_angle_3 = -100
     block  = '15'
   [../]
   [./seebeck_tensor16]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 10
     euler_angle_2 = 20
     euler_angle_3 = 30
     block  = '16'
   [../]
   [./seebeck_tensor17]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -10
     euler_angle_2 = -20
     euler_angle_3 = -30
     block  = '17'
   [../]
   [./seebeck_tensor18]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 30
     euler_angle_2 = 60
     euler_angle_3 = 90
     block  = '18'
   [../]
   [./seebeck_tensor19]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = -30
     euler_angle_2 = -60
     euler_angle_3 = -90
     block  = '19'
   [../]
   [./seebeck_tensor20]
     type = ComputeSeebeckTensor
     a_ij = '59e-6 1.03e-21 5.8e-6 9.4e-22 14.4e-6 1.17e-21 5.8e-6 1.18e-21 15.18e-6'#mp-568390
     euler_angle_1 = 90
     euler_angle_2 = 90
     euler_angle_3 = 0
     block  = '20'
   [../]
 []


[BCs]
  [./sideset_1T]
    type = DirichletBC
    variable = T
    boundary = 'left'
    value = 273
  [../]
  [./sideset_2T]
    type = DirichletBC
    variable = T
    boundary = 'right'
    value = 273
  [../]

  [./sideset_1V]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'right'
    value = 0
  [../]
  [./sideset_2V]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'left'
    value = 0.05
  [../]
[]

[Postprocessors]
  [./potential_x]
    type = PointValue
    point = '0 0 0'
    variable = potential_E_int
  [../]

  [./T]
    type = PointValue
    point = '0 0 0'
    variable = T
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160               1e-10      1e-8      1e-6          bjacobi       allreduce'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_Bi2Te3_test
    elemental_as_nodal = true
  [../]
[]
