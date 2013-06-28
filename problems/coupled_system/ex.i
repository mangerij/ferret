[Mesh]
  file = poissonstripe_coarse.e
  uniform_refine=1
[]
# [Mesh]
#   #file = mug.e
#   type=GeneratedMesh
#   dim=3
#   nx=16
#   ny=16
#   nz=16
#   xmin=0.0
#   xmax=1.0
#   ymin=0.0
#   ymax=1.0
#   zmin=0.0
#   zmax=1.0
#   #uniform_refine=1
# []
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
  # [./polar_x]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   #block='interior'
  # [../]
  # [./polar_y]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   #block='interior'
  # [../]
  # [./polar_z]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   #block='interior'
  # [../]
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha111=2.6e8
    alpha12=7.5e8
    alpha112=6.1e8
    alpha123=-3.7e9
    #alpha12=0.0
    #alpha112=0.0
    #alpha123=0.0
    #block='interior'
   # alpha112=0.0
   # alpha123=0.0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha111=2.6e8
    alpha12=7.5e8
    alpha112=6.1e8
    alpha123=-3.7e9
    #alpha12=0.0
    #alpha112=0.0
    #alpha123=0.0
    #block='interior'
   # alpha112=0.0
   # alpha123=0.0

  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
    alpha11=-7.3e7
    alpha111=2.6e8
    alpha12=7.5e8
    alpha112=6.1e8
    alpha123=-3.7e9
    #alpha12=0.0
    #alpha112=0.0
    #alpha123=0.0
    #block='interior'
   # alpha112=0.0
   # alpha123=0.0
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
    #type=RandomIC
    #type=ConstantIC
    type=PerturbedIC
    variable=polar_x
    mean=0.0
    factor=0.1
  [../]
  [./polar_y]
    #type=ConstantIC
    #type=RandomIC
    type=PerturbedIC
    variable=polar_y
    mean=0.0
    factor=0.1
  [../]
  [./polar_z]
    #type=ConstantIC
    #type=RandomIC
    type=PerturbedIC
    variable=polar_z
    mean=1.0
    factor=0.1
  [../]
[]

[BCs]
  [./potential_upz]
    type = DirichletBC
    variable = potential
    boundary = 'upz'
    value = 1.0
  [../]
  [./potential_downz]
    type = DirichletBC
    variable = potential
    boundary = 'downz'
    value = -1.0
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
[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     #off_diag_row='var_name'
     #off_diag_column='var_name'
     full=true   #to use every off diagonal block
     #petsc_options='snes_mf_operator'
     #petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'
     #petsc_options_value = 'lu       1e-6                 ds'
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
  petsc_options='-snes_monitor -snes_view -snes_converged_reason'
  #petsc_options_iname='-snes_linearsearch_type -snes_rtol'
  #petsc_options_value='basic                   1e-16'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]

[Output]
  file_base = out
  #interval = 1
  exodus = true
  perf_log = true
[]
