#######################
#In this case, we only consider the electrostatic energy:
# we can compute the exactly analytic solution and we can show the minimizer is px=py=0 component and pz=-1.1e-4 inside the material domain. And we indeed obtain this solution by 800 step wich step size 1e-3.
# However, the interesting fact is that, the electrostatic energy is almost unchanged after 250, however, at 250 steps, px and py is still quite large. This indicates the energy surface actually is quite very flat on x and y-direction and it demonstrate a need of adaptive time step. It's also a sign of ill-conditioning. 
[Mesh]
  file=slab.e
  #uniform_refine=1
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
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
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
[]

[GlobalParams]
   len_scale=1e-7
   energy_scale=1e12
   polar_elecctric_scale=1e-2
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   #G110=1.73e4
   G110=0.0
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   permittivity=8.85e-12
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
[]

[Kernels]
  active='diffusion_E diffusion_E_Ext  polar_electric_E polar_electric_px polar_electric_py polar_electric_pz polar_x_time polar_y_time polar_z_time'
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
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component=2
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
  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = 'upz_outter'
    value = 0.0
    #implicit=false
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = 'downz_outter'
    value = 0.0
    #implicit=false
  [../]
   [./potential_ext_upz]
    type = DirichletBC
    variable = potential_ext
    boundary = 'upz_outter'
    value = 1.0
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = DirichletBC
    variable = potential_ext
    boundary = 'downz_outter'
    value = 0.0
    #implicit=false
  [../]
[]
[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     full=true   
     pc_side=left
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e-3
  nl_max_its=100
  num_steps=800
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-5        gmres       1e-8      jacobi       basic                1e-50'
[]

[Postprocessors]
   [./electric_energy]
    type=ElectrostaticEnergy
    [../]
[]

[Output]
  #file_base = out
  output_initial=1
  #interval = 1
  exodus = true
  perf_log = true
[]