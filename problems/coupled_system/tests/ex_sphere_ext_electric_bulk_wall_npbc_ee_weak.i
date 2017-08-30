[Mesh]
  file = concentric_sphere.e
  uniform_refine=0
[]
[Variables]
#active='polar_x polar_y polar_z'
#active='potential_int potential_ext'
  [./polar_x]
    #scaling=1e-3
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./polar_y]
    #scaling=1e-3
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./polar_z]
    #scaling=1e-3
    order = FIRST
    family = LAGRANGE
    block='interior'
  [../]
  [./potential_int]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
[./auxv_es_energy_density_e] #es for electrostatic
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./auxv_es_energy_density]
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./auxv_bulk_energy_density]
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./auxv_wall_energy_density]
      order=CONSTANT
      family=MONOMIAL
  [../]
[]

[GlobalParams]
   len_scale=1e-7
   #len_scale=1.0
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=1.73e4
   #G110=0.0
   G11_G110=0.6
   G12_G110=0.0
   G44_G110=0.3
   G44P_G110=0.3
   permittivity=8.85e-12
   #permittivity=1.0
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
   potential_int=potential_int
   potential_ext=potential_ext
[]

[AuxKernels]
  #active='diff'
  [./auxk_electrostatic_energy_density_e]
    type =ElectrostaticEnergyDensityE
    variable =auxv_es_energy_density_e
    potential=potential_int
  [../]
  [./auxk_electrostatic_energy_density]
    type =ElectrostaticEnergyDensity
    variable =auxv_es_energy_density
  [../]
  [./auxk_bulk_energy_density]
    type =BulkEnergyDensity
    variable =auxv_bulk_energy_density
  [../]
  [./auxk_wall_energy_density]
     type=WallEnergyDensity
     variable=auxv_wall_energy_density
  [../]
[]

[Kernels]
  active='diffusion_E diffusion_E_Ext  polar_electric_E polar_electric_px polar_electric_py polar_electric_pz polar_x_time polar_y_time polar_z_time bed_x bed_y bed_z walled_x walled_y walled_z'
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
    implicit=false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
    implicit=false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
    implicit=false
  [../]
  [./walled_x]
     type=WallEnergyDerivative
     variable=polar_x
     component=0
     implicit=false
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component=1
     implicit=false
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component=2
     implicit=false
  [../]
  [./polar_electric_E]
     type=PolarElectricE
     variable=potential_int
     block='interior'
     implicit=false
  [../]
  [./diffusion_E]
     type=ElectricStatics
     #type=Diffusion
     variable=potential_int
     block='exterior interior'
  [../]
  [./diffusion_E_Ext]
     type=ElectricStatics
     #type=Diffusion
     variable=potential_ext
     block='exterior interior'
  [../]
  [./polar_electric_px]
     type=PolarElectricP
     variable=polar_x
     component=0
     implicit=false
  [../]
  [./polar_electric_py]
     type=PolarElectricP
     variable=polar_y
     component=1
     implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricP
     variable=polar_z
     component=2
     implicit=false
  [../]
  [./polar_x_time]
     type=TimeDerivative
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative
     variable=polar_z
  [../]
[]

[ICs]
  active='polar_x_constic polar_y_constic polar_z_constic'
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     value=0.6
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value=0.3
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value=0.3
  [../]
[]

[BCs]
 # active ='Periodic'
  [./potential_int]
    type = DirichletBC
    variable = potential_int
    boundary = 'outter'
    value = 0.0
    #implicit=false
  [../]

   [./potential_ext_upz]
    type = DirichletBC
    variable = potential_ext
    boundary = 'outter'
    value = 1.0
    #implicit=false
  [../]
  # [./Periodic]
  #   #active='polar_x_x polar_y_x polar_z_x polar_x_y polar_y_y polar_z_y'
  #   #active='potential_int_x potential_ext_x potential_int_y potential_ext_y'
  #   [./potential_int_x]
  #      variable = potential_int
  #      primary = 'downx'
  #      secondary = 'upx'
  #      translation = '1 0 0'
  #      #implicit=false
  #   [../]
  #   [./potential_ext_x]
  #      variable = potential_ext
  #      primary = 'downx'
  #      secondary = 'upx'
  #      translation = '1 0 0'
  #      #implicit=false
  #   [../]
  #   [./polar_x_x]
  #      variable = polar_x
  #      primary = 'downx'
  #      secondary = 'upx'
  #      translation = '1 0 0'
  #      #implicit=false
  #   [../]
  #   [./polar_y_x]
  #      variable = polar_y
  #      primary = 'downx'
  #      secondary = 'upx'
  #      translation = '1 0 0'
  #      #implicit=false
  #   [../]
  #   [./polar_z_x]
  #      variable = polar_z
  #      primary = 'downx'
  #      secondary = 'upx'
  #      translation = '1 0 0'
  #      #implicit=false
  #   [../]
  #   [./potential_int_y]
  #      variable = potential_int
  #      primary = 'downy'
  #      secondary ='upy'
  #      translation = '0 1 0'
  #      #implicit=false
  #   [../]
  #   [./potential_ext_y]
  #      variable = potential_ext
  #      primary = 'downy'
  #      secondary ='upy'
  #      translation = '0 1 0'
  #      #implicit=false
  #   [../]
  #   [./polar_x_y]
  #      variable = polar_x
  #      primary = 'downy'
  #      secondary ='upy'
  #      translation = '0 1 0'
  #      #implicit=false
  #   [../]
  #   [./polar_y_y]
  #      variable = polar_y
  #      primary = 'downy'
  #      secondary ='upy'
  #      translation = '0 1 0'
  #      #implicit=false
  #   [../]
  #   [./polar_z_y]
  #      variable = polar_z
  #      primary = 'downy'
  #      secondary ='upy'
  #      translation = '0 1 0'
  #      #implicit=false
  #   [../]
  # [../]
[]
[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     #off_diag_row='var_name'
     #off_diag_column='var_name'
     full=true   #to use every off diagonal block
     pc_side=left
     #petsc_options='snes_mf_operator'
     #petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'
     #petsc_options_value = 'lu       1e-6                 ds'
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1
  nl_max_its=100
  num_steps=800
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-5        gmres       1e-8      jacobi       basic                1e-50'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]
[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
   [../]
   [./wall_energy]
    type=WallEnergy
   [../]
   [./electric_energy]
    type=ElectrostaticEnergy
    [../]
    [./total_energy]
    type=TotalEnergy
    bulk_energy=bulk_energy
    wall_energy=wall_energy
    electric_energy=electric_energy
    [../]
[]

[Output]
  #file_base = out
  output_initial=1
  #interval = 1
  exodus = true
  perf_log = true
[]
