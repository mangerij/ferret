[Mesh]
  file = poissonstripe_coarse.e
[]
[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./potential]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=3.8e5  # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha12=7.5e8
    alpha111=2.6e8
    alpha112=6.1e8
    alpha123=-3.7e9
    block='interior'
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=3.8e5  # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha12=7.5e8
    alpha111=2.6e8
    alpha112=6.1e8
    alpha123=-3.7e9
    block='interior'
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=3.8e5  # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha12=7.5e8
    alpha111=2.6e8
    alpha112=6.1e8
    alpha123=-3.7e9
    block='interior'
  [../]
  [./walled_x]
     type=WallEnergyDerivative
     variable=polar_x
     component=0
     polar_x=polar_x
     polar_y=polar_y
     polar_z=polar_z
     G110=1.73e-10
     G11/G110=0.6
     G12/G110=0.0
     G44/G110=0.3
     G44P/G110=0.3
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component=1
     polar_x=polar_x
     polar_y=polar_y
     polar_z=polar_z
     G110=1.73e-10
     G11/G110=0.6
     G12/G110=0.0
     G44/G110=0.3
     G44P/G110=0.3
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component=2
     polar_x=polar_x
     polar_y=polar_y
     polar_z=polar_z
     G110=1.73e-10
     G11/G110=0.6
     G12/G110=0.0
     G44/G110=0.3
     G44P/G110=0.3
  [../]
  [./polar_electric_E]
     type=PolarElectricE
     variable=potential
     polar_x = polar_x
     polar_y = polar_y
     polar_z = polar_z
     permittivity=8.85e-12
     block='interior'
  [../]
  [./diffusion_E]
     type=ElectricStatics
     variable=potential
     permittivity=8.85e-12
     block='exterior interior'
  [../]
  [./polar_electric_px]
     type=PolarElectricP
     variable=polar_x
     component=0
     potential=potential
  [../]
  [./polar_electric_py]
     type=PolarElectricP
     variable=polar_y
     component=1
     potential=potential
  [../]
  [./polar_electric_px]
     type=PolarElectricP
     variable=polar_z
     component=2
     potential=potential
  [../]
[]

[ICs]
  [./polar_x]
    type=RandomIC
    variable=polar_x
  [../]
  [./polar_y]
    type=RandomIC
    variable=polar_y
  [../]
  [./polar_z]
    type=RandomIC
    variable=polar_z
  [../]
[]

[BCs]
  [./potential_upz]
    type = DirichletBC
    variable = potential
    boundary = 'upz'
    value = 0
  [../]
  [./potential_downz]
    type = DirichletBC
    variable = potential
    boundary = 'downz'
    value = 0
  [../]
  [./Periodic]
    [./potential_x]
       variable = potential
       primary = 'downx'
       secondary = 'upx'
       translation = '1 0 0'
    [../]
    [./polar_x_x]
       variable = polar_x
       primary = 'downx'
       secondary = 'upx'
       translation = '1 0 0'
    [../]
    [./polar_y_x]
       variable = polar_y
       primary = 'downx'
       secondary = 'upx'
       translation = '1 0 0'
    [../]
    [./polar_z_x]
       variable = polar_z
       primary = 'downx'
       secondary = 'upx'
       translation = '1 0 0'
    [../]
    [./potential_y]
       variable = potential
       primary = 'downy'
       secondary ='upy'
       translation = '0 1 0'
    [../]
    [./polar_x_y]
       variable = polar_x
       primary = 'downy'
       secondary ='upy'
       translation = '0 1 0'
    [../]
    [./polar_y_y]
       variable = polar_y
       primary = 'downy'
       secondary ='upy'
       translation = '0 1 0'
    [../]
    [./polar_z_y]
       variable = polar_z
       primary = 'downy'
       secondary ='upy'
       translation = '0 1 0'
    [../]
  [../]
[]

[Postprocessors]
  [./BulkEnergy]
   type=BulkEnergy
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   alpha1=3.8e5  # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha12=7.5e8
   alpha111=2.6e8
   alpha112=6.1e8
   alpha123=-3.7e9
  [../]
[]

[Executioner]
  type = Steady
  nl_max_its=1000
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-snes_monitor -snes_converged_reason'
[]

[Output]
  file_base = out
  #interval = 1
  exodus = true
  perf_log = true
[]
