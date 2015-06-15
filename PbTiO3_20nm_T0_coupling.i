[Mesh]
  file = slab_exodus_coarse_20.e
  uniform_refine=0
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
  [../]

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

[AuxVariables]
 #[./auxv_es_energy_density_e] #es for electrostatic
 #    order=CONSTANT
 #    family=MONOMIAL
 # [../]
 # [./auxv_es_energy_density]
 #    order=CONSTANT
 #    family=MONOMIAL
 # [../]
 # [./auxv_bulk_energy_density]
 #    order=CONSTANT
 #    family=MONOMIAL
 # [../]

  #[./Ez]
  #   order=CONSTANT
  #   family=MONOMIAL
  #[../]
  #[./Ex]
  #   order=CONSTANT
  #   family=MONOMIAL
  #[../]
  #[./Ey]
  #   order=CONSTANT
  #   family=MONOMIAL
  #[../]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zx]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[GlobalParams]
   len_scale=1e-9
   #len_scale=1.0
   alpha1=-1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=1.0e-10
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
   polar_scale = 1.0
   disp_x = disp_x
   disp_y = disp_y
   disp_z = disp_z

# this ostensibly defines dt = 1e-14. Note if time_scale=1 then dt=1e17 to see domain relaxation
# Olle suggests that the dynamics should occur on the 10 ps time scale. This parameter allows the time scale to be
# user-input and most likely depends on material and problem in question, maybe even temperature.
[]

[AuxKernels]
  #[./auxk_electrostatic_energy_density_e]
  #  type = ElectrostaticEnergyDensityE
  #  variable =auxv_es_energy_density_e
  #  potential=potential_int
  #[../]
  #[./auxk_electrostatic_energy_density]
  #  type = ElectrostaticEnergyDensity
  #  variable =auxv_es_energy_density
  #[../]
  #[./auxk_bulk_energy_density]
  #  type = BulkEnergyDensity
  #  variable =auxv_bulk_energy_density
  #[../]
  #[./Ez_fieldaux]
  #  type = Ez_fieldAux
  #  variable = Ez
  #[../]
  #[./Ex_fieldaux]
  #  type = Ex_fieldAux
  #  variable = Ex
  #[../]
  #[./Ey_fieldaux]
  #  type = Ey_fieldAux
  #  variable = Ey
  #[../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx
    index_i = 2
    index_j = 0
  [../]
[]




[Kernels]
#Bulk energy density
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

#Wall energy penalty
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
#Ferroelectric-strain coupling
  [./ferroelectriccouplingu_x]
     type = FerroelectricCouplingU
     variable=polar_x
     component=0
     block = '2'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingU
     variable=polar_y
     component=1
     block = '2'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingU
     variable=polar_z
     component=2
     block = '2'
  [../]

  [./ferroelectriccouplingp_x]
     type = FerroelectricCouplingP
     variable=polar_x
     component=0
  [../]
  [./ferroelectriccouplingp_y]
     type = FerroelectricCouplingP
     variable=polar_y
     component=1
  [../]
  [./ferroelectriccouplingp_z]
     type = FerroelectricCouplingP
     variable=polar_z
     component=2
  [../]


#Tensor mechanics--Hooke's law
  [./TensorMechanicsScaled]
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
  [../]

#Electrostatics
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'
     implicit=false
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     block='1'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     block='2'
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
     block='1'
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
     block='2'
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

[Materials]
#  [./slab_elastic]
#    type=LinearElasticMaterial
#    block = '2'
#    disp_x = disp_x
#    disp_y = disp_y
#    disp_z = disp_z
#in GPA. from N. Pandech et al. Ceramics International,
# just multiply by 1e9 to convert to N/m^2
# C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
#    euler_angle_1 = 0.0
#    euler_angle_2 = 0.0
#    euler_angle_3 = 0.0
#  [../]

  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
#in GPA. from N. Pandech et al. Ceramics International,
# just multiply by 1e9 to convert to N/m^2
# C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
#in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
# Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.034 0.034 0.034'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./vacuum]
    type=LinearElasticMaterial
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
#vacuum elasticity
    C_ijkl = '0 0 0 0 0 0 0 0 0'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
[]


[ICs]
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     value = 0.6
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value = 0.6
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value = 0.1
  [../]



  [./disp_x_constic_vacuum]
     type=ConstantIC
     variable=disp_x
     block = '1'
     value = 0
  [../]
  [./disp_y_constic_vacuum]
     type=ConstantIC
     variable=disp_y
     block = '1'
     value = 0
  [../]
  [./disp_z_constic_vacuum]
     type=ConstantIC
     variable=disp_z
     block = '1'
     value = 0
  [../]
[]

[BCs]

  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0.0
  [../]



  [./disp_z_slab]
    type = DirichletBC
    variable = disp_z
    boundary = '12'
    value = -1.0
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
  [../]
  [./potential_ext_downz]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    #value = -5.0e-11
    value = 0.0
  [../]
[]

[Executioner]
  type=Transient
  solve_type=newton
  scheme = explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e-25  #adjustable within a few orders of magnitude. It seems that we need small dt steps but large num_steps.
  dtmin=1e-30
#NOTE: First time step calculates the depolarization field due to the unphysical initial condition. Energy may increase, which is allowed.
  num_steps=3500
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-snes_rtol -ksp_type -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1e-6        gmres       1e-8      ilu       basic                1e-50       '
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
    permittivity=8.85e-12
   [../]
   [./total_energy]
    type=TotalEnergy
    bulk_energy=bulk_energy
    wall_energy=wall_energy
    electric_energy=electric_energy
   [../]
   [./gradx]
     type=TotalEnergyGradient
     component = 0
   [../]
   [./grady]
     type=TotalEnergyGradient
     component = 1
   [../]
   [./gradz]
     type=TotalEnergyGradient
     component = 2
   [../]
   [./GradE-L2]
     type=TotalEnergyGradientL2
     gradx = gradx
     grady = grady
     gradz = gradz
   [../]
   [./R(i)]
     type=Residual
     block = '1 2'
   [../]
[]


[Outputs]

print_linear_residuals = true
print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_PbTiO3_T0_20nm_icConst111_0field_1wall_coupling
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]

  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
