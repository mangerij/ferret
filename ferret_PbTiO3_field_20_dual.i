[Mesh]
  file = slab_dual_100.e
  uniform_refine=0
[]

[Variables]
#active='polar_x polar_y polar_z'
#active='potential_int potential_ext'
  [./polar_x]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2 3'
  [../]
  [./polar_y]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2 3'
  [../]
  [./polar_z]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2 3'
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
  [./fieldEz]
     order=CONSTANT
     family=MONOMIAL
  [../]
[]

[GlobalParams]
   len_scale=1e-9
   #len_scale=1.0
   alpha1=-1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=1e-10
#   #G110=0.0
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
   time_scale = 1.0e-31

# this ostensibly defines dt = 1e-14. Note if time_scale=1 then dt=1e17 to see domain relaxation
# Olle suggests that the dynamics should occur on the 10 ps time scale. This parameter allows the time scale to be
# user-input and most likely depends on material and problem in question, maybe even temperature.
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
  [./Ez_fieldaux]
    type = Ez_fieldAux
    variable = fieldEz
    potential_ext = potential_ext
    potential_int = potential_int
  [../]
#  [./Ex_fieldaux]
#    type = Ex_fieldAux
#    variable = fieldEx
#    potential_ext = potential_ext
#    potential_int = potential_int
#  [../]
#  [./Ey_fieldaux]
#    type = Ey_fieldAux
#    variable = fieldEy
#    potential_ext = potential_ext
#    potential_int = potential_int
#  [../]
[]

[Kernels]
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
     type=WallEnergyDerivative_scaled
     variable=polar_x
     component=0
  [../]
  [./walled_y]
     type=WallEnergyDerivative_scaled
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative_scaled
     variable=polar_z
     component=2
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'
     implicit=false
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     permittivity=8.85e-12
     block='2 3'
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     permittivity=8.85e-12
     block='1'
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
     permittivity=8.85e-12
     block='1'
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
     permittivity=8.85e-12
     block='2 3'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     component=0
     implicit=false
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     component=1
     implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     component=2
     implicit=false
  [../]
  [./polar_x_time]
     type=TimeDerivative_scaled
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative_scaled
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative_scaled
     variable=polar_z
  [../]
[]

[ICs]
  [./polar_x_constic2]
     type=ConstantIC
     variable=polar_x
     value = 0.6
     block = '2'
  [../]
  [./polar_y_constic2]
     type=ConstantIC
     variable=polar_y
     value = 0.6
     block = '2'
  [../]
  [./polar_z_constic2]
     type=ConstantIC
     variable=polar_z
     value = 0.1
     block = '2'
  [../]
  [./polar_x_constic3]
     type=ConstantIC
     variable=polar_x
     value = -0.6
     block = '3'
  [../]
  [./polar_y_constic3]
     type=ConstantIC
     variable=polar_y
     value = -0.6
     block = '3'
  [../]
  [./polar_z_constic3]
     type=ConstantIC
     variable=polar_z
     value = -0.1
     block = '3'
  [../]
[]

[BCs]

  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0
    #implicit=false
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0.0
    #implicit=false
  [../]
# Applied field: for zero field use NeumannBC on the external potential = 0. A
# Note that \nabla^2 \Phi_{ext} = 0 is satisfied if \Phi_{ext} = 0, ie Dirichlet and Neumann BC classes are equivalent
   [./potential_ext_upz]
    type = NeumannBC
    variable = potential_ext
    boundary = '1'
    #value = 5.0e-11 #1e-12 gives a field of ~1e8 V/m ~ 1000 kV/cm for size 1, what about size 100?
                    # seems that this doesn't matter here? dt needs to change for the size however...
                    # lowering this value, to about 1e-16 is the limit here. 5e-11 gives 5.6e9 V/m for 100 size
    value = 0.0
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    #value = -5.0e-11
    value = 0.0
    #implicit=false
  [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme = explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e-16    #adjustable within a few orders of magnitude.
  num_steps=2500
  dtmin=1e-19
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value='1e-6        gmres       1e-8      hypre       basic                1e-50                     0.25' #may need thresold to 0.5 for 3D problems?
[]

[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
    block = '2 3'
   [../]
   [./wall_energy]
    type=WallEnergy

   [../]
   [./electric_energy]
    type=ElectrostaticEnergy
    block = '2 3'
    [../]
    [./total_energy]
    type=TotalEnergy
    bulk_energy=bulk_energy
    wall_energy=wall_energy
    electric_energy=electric_energy
     block = '2 3'
    [../]
    [./dE(P)dPx]
     type=TotalEnergyGradient
     component = 0
     execute_on = timestep_end
     block = '2 3'
    [../]
    [./dE(P)dPy]
     type=TotalEnergyGradient
     component = 1
     execute_on = timestep_end
     block = '2 3'
    [../]
    [./dE(P)dPz]
     type=TotalEnergyGradient
     component = 2
     execute_on = timestep_end
     block = '2 3'
    [../]
    [./R(i)]
     type=Residual
     execute_on = timestep_end
     block = '1 2 3'
    [../]
[]


[Outputs]
  file_base = out_PbTiO3_bulk_len-9_icConst_n700_field_10wall_dual_100
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    interval = 1
  [../]
  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
