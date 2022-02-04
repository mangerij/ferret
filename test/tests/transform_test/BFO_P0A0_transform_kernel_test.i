Nx = 3
Ny = 3
Nz = 3

xMax = 1.0
yMax = 1.0
zMax = 1.0

    ############################################
    ##
    ##   This input file aims to solve the BFO
    ##   problem where one of the components of 
    ##   the order parameters is aligned along 
    ##   [111] || global z
    ##
    ##   As such, all of the residuals and
    ##   jacobians need to be transformed (rotated)
    ##   This will be done in general at a later date
    ##   However, since this the Aux system is used
    ##   extensively to do this, indices will appear
    ##   as p0,p1,p2,a0,a1,a2,u0,u1,u2 which indicate
    ##   derivatives w.r.t the components of various
    ##   order parameters in the transformed coords.
    ##   
    ##   Indices of microforces and jacobians
    ##   will also have bp, br, rp, els, ros, sels, sros
    ##   which denote the bulk polar, bulk roto, rotopolar, 
    ##   electrostrictive, rotostrictive, 
    ##   stress- electrostrictive, and stress- rotostrictive 
    ##   terms.
    ##
    ##   i.e, Jbp_p0p0 corresoponds to the on-diagonal
    ##   jacobian associated with the bulk energy for the
    ##   polarization.
    ## 
    ##   any instance of 1_x,1_y,1_z will denote the current 
    ##   global cartesian coordinate system which is
    ##   orthonormal and defined by the S matrix.
    ##
    ##   any instance of o_x,o_y,o_z denote the original
    ##   cartesian coordinate system.
    ##
    ##   i.e, P1 = S Po => Po = S^{-1} P1
    ##
    ##   Finally, some indices are q0,q1,q2,q3,
    ##   which follow q = <p0,p1,p2,a0,a1,a2,u0,u1,u2>
    ##   ordering.
    ##
    ##
    ############################################

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${Nx}
    ny = ${Ny}
    nz = ${Nz}
    xmin = 0.0
    xmax = ${xMax}
    ymin = 0.0
    ymax = ${yMax}
    zmin = 0.0
    zmax = ${zMax}
    elem_type = HEX8
  []
  [./cnode]
    input = gen

    ############################################
    ##
    ##   additional boundary sideset (one node) 
    ##   to zero one of the elastic displacement vectors 
    ##   vectors and eliminates rigid body translations 
    ##   from the degrees of freedom
    ##
    ##   NOTE: This must conform with the about
    ##         [Mesh] block settings
    ##
    ############################################

    type = ExtraNodesetGenerator
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0


  displacements = 'u1_x u1_y u1_z'
[]


[Functions]
  [./constPp]
    type = ParsedFunction
    value = 0.54
  [../]
  [./constAp]
    type = ParsedFunction
    value = 7.37
  [../]
[]


[Variables]
  [./u1_x]
  [../]
  [./u1_y]
  [../]
  [./u1_z]
  [../]
 # [./global_strain]
 #   order = SIXTH
 #   family = SCALAR
 # [../]

  [./P1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
    [../]
  [../]
  [./P1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
    [../]
  [../]
  [./P1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constPp
    [../]
  [../]
  [./A1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
    [../]
  [../]
  [./A1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
    [../]
  [../]
  [./A1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constAp
    [../]
  [../]
[]

[AuxVariables]

    ############################################
    ##
    ## Transformed coordinates
    ##  
    ##   These follow Po = Inv[S] P1
    ##   In the AuxKernels that calculate this
    ##
    ##   We will flag this as 'inverse = true'
    ##
    ############################################

  [./Po_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Po_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Po_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Ao_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Ao_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Ao_z]
    order = FIRST
    family = LAGRANGE
  [../]


    ############################################
    ##
    ## Microforces
    ##  
    ##  We have bulk energy P, bulk energy A 
    ##  and a RP coupling.
    ##  This means we should have 3+3+6 forces. 
    ##  
    ############################################

  [./fbp_p0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fbp_p1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fbp_p2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fba_a0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fba_a1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fba_a2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./frp_p0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./frp_p1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./frp_p2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./frp_a0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./frp_a1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./frp_a2]
    order = CONSTANT
    family = MONOMIAL
  [../]



  [./Jbp_p0p0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jbp_p1p1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jbp_p2p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]

  [./Jbp_p0p1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jbp_p1p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jbp_p0p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]

  [./Jba_a0a0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jba_a1a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jba_a2a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jba_a0a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jba_a1a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jba_a0a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]


  [./Jrp_p0p0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p1p1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p2p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p0p1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p1p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p0p2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]

  [./Jrp_a0a0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_a1a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_a2a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_a0a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_a1a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_a0a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]

  [./Jrp_p0a0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p0a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p0a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p1a0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p1a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p1a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p2a0]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p2a1]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./Jrp_p2a2]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
[]

[AuxKernels]
  [./p1]
    type = Transformed111Order
    variable = Po_x
    inverse = true
    component = 0
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
  [../]
  [./p2]
    type = Transformed111Order
    variable = Po_y
    inverse = true
    component = 1
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
  [../]
  [./p3]
    type = Transformed111Order
    variable = Po_z
    inverse = true
    component = 2
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
  [../]

  [./a1]
    type = Transformed111Order
    variable = Ao_x
    inverse = true
    component = 0
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z
  [../]
  [./a2]
    type = Transformed111Order
    variable = Ao_y
    inverse = true
    component = 1
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z
  [../]
  [./a3]
    type = Transformed111Order
    variable = Ao_z
    inverse = true
    component = 2
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z
  [../]

  ##################################################
  ##
  ##
  ##   Microforces
  ##
  ##
  ##################################################

  [./fbpp1]
    type = MicroforceBulkEnergy
    variable = fbp_p0
    component = 0
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./fbpp2]
    type = MicroforceBulkEnergy
    variable = fbp_p1
    component = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./fbpp3]
    type = MicroforceBulkEnergy
    variable = fbp_p2
    component = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]

  [./fbaa1]
    type = MicroforceRotoBulkEnergy
    variable = fba_a0
    component = 0
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./fbaa2]
    type = MicroforceRotoBulkEnergy
    variable = fba_a1
    component = 1
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./fbaa3]
    type = MicroforceRotoBulkEnergy
    variable = fba_a2
    component = 2
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]


  [./frpp0]
    type = MicroforceRotopolarCoupledPolarEnergy
    variable = frp_p0
    component = 0
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./frpp1]
    type = MicroforceRotopolarCoupledPolarEnergy
    variable = frp_p1
    component = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./frpp2]
    type = MicroforceRotopolarCoupledPolarEnergy
    variable = frp_p2
    component = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
    outputs = 'none'
  [../]

  [./frpa0]
    type = MicroforceRotopolarCoupledDistortEnergy
    variable = frp_a0
    component = 0
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./frpa1]
    type = MicroforceRotopolarCoupledDistortEnergy
    variable = frp_a1
    component = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./frpa2]
    type = MicroforceRotopolarCoupledDistortEnergy
    variable = frp_a2
    component = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]

  ##################################################
  ##
  ##
  ##   Jacobians
  ##
  ##
  ##################################################


  [./Jbpp0p0]
    type = JacobiansBulkEnergy
    variable = Jbp_p0p0
    index_i = 0
    index_j = 0
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbpp1p1]
    type = JacobiansBulkEnergy
    variable = Jbp_p1p1
    index_i = 1
    index_j = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbpp2p2]
    type = JacobiansBulkEnergy
    variable = Jbp_p2p2
    index_i = 2
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbpp0p1]
    type = JacobiansBulkEnergy
    variable = Jbp_p0p1
    index_i = 0
    index_j = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbpp1p2]
    type = JacobiansBulkEnergy
    variable = Jbp_p1p2
    index_i = 1
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbpp0p2]
    type = JacobiansBulkEnergy
    variable = Jbp_p0p2
    index_i = 0
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]


  [./Jbaa0a0]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a0a0
    index_i = 0
    index_j = 0
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbaa1a1]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a1a1
    index_i = 1
    index_j = 1
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbaa2a2]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a2a2
    index_i = 2
    index_j = 2
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbaa0a1]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a0a1
    index_i = 0
    index_j = 1
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbaa1a2]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a1a2
    index_i = 1
    index_j = 2
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jbaa0a2]
    type = JacobiansRotoBulkEnergy
    variable = Jba_a0a2
    index_i = 0
    index_j = 2
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]




  [./Jrpp0p0]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0p0
    index_i = 0
    index_j = 0
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp1p1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p1p1
    index_i = 1
    index_j = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp2p2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p2p2
    index_i = 2
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa0a0]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a0a0
    index_i = 3
    index_j = 3
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa1a1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a1a1
    index_i = 4
    index_j = 4
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa2a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a2a2
    index_i = 5
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]

  [./Jrpp0p1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0p1
    index_i = 0
    index_j = 1
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp1p2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p1p2
    index_i = 1
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp0p2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0p2
    index_i = 0
    index_j = 2
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp0a0]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0a0
    index_i = 0
    index_j = 3
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp0a1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0a1
    index_i = 0
    index_j = 4
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp0a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p0a2
    index_i = 0
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]

  [./Jrpp1a0]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p1a0
    index_i = 1
    index_j = 3
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp1a1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p1a1
    index_i = 1
    index_j = 4
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp1a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p1a2
    index_i = 1
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp2a0]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p2a0
    index_i = 2
    index_j = 3
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp2a1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p2a1
    index_i = 2
    index_j = 4
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpp2a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_p2a2
    index_i = 2
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa0a1]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a0a1
    index_i = 3
    index_j = 4
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa0a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a0a2
    index_i = 3
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]
  [./Jrpa1a2]
    type = JacobiansRotopolarCoupledEnergy
    variable = Jrp_a1a2
    index_i = 4
    index_j = 5
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
    execute_on = 'INITIAL LINEAR NONLINEAR'
  [../]

[]

#[ScalarKernels]
#  [./global_strain]
#    type = GlobalStrain
#    variable = global_strain
#    global_strain_uo = global_strain_uo
#  [../]
#[]

[Kernels]

  [./TensorMechanics]
  [../]


  ### Operators for the polar field: ###

  [./bed_x]
    type = Transformed111KernelOp3
    variable = P1_x
    component = 0
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z

    f_q0 = fbp_p0
    f_q1 = fbp_p1
    f_q2 = fbp_p2

    J_q0q0 = Jbp_p0p0
    J_q1q1 = Jbp_p1p1
    J_q2q2 = Jbp_p2p2
    J_q0q1 = Jbp_p0p1
    J_q1q2 = Jbp_p1p2
    J_q0q2 = Jbp_p0p2
  [../]
  [./bed_y]
    type = Transformed111KernelOp3
    variable = P1_y
    component = 1
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z

    f_q0 = fbp_p0
    f_q1 = fbp_p1
    f_q2 = fbp_p2

    J_q0q0 = Jbp_p0p0
    J_q1q1 = Jbp_p1p1
    J_q2q2 = Jbp_p2p2
    J_q0q1 = Jbp_p0p1
    J_q1q2 = Jbp_p1p2
    J_q0q2 = Jbp_p0p2
  [../]
  [./bed_z]
    type = Transformed111KernelOp3
    variable = P1_z
    component = 2
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z

    f_q0 = fbp_p0
    f_q1 = fbp_p1
    f_q2 = fbp_p2

    J_q0q0 = Jbp_p0p0
    J_q1q1 = Jbp_p1p1
    J_q2q2 = Jbp_p2p2
    J_q0q1 = Jbp_p0p1
    J_q1q2 = Jbp_p1p2
    J_q0q2 = Jbp_p0p2
  [../]


  # Operators for the AFD field

  [./rbed_x]
    type = Transformed111KernelOp3
    variable = A1_x
    component = 0
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z

    f_q0 = fba_a0
    f_q1 = fba_a1
    f_q2 = fba_a2

    J_q0q0 = Jba_a0a0
    J_q1q1 = Jba_a1a1
    J_q2q2 = Jba_a2a2
    J_q0q1 = Jba_a0a1
    J_q1q2 = Jba_a1a2
    J_q0q2 = Jba_a0a2
  [../]
  [./rbed_y]
    type = Transformed111KernelOp3
    variable = A1_y
    component = 1
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z

    f_q0 = fba_a0
    f_q1 = fba_a1
    f_q2 = fba_a2

    J_q0q0 = Jba_a0a0
    J_q1q1 = Jba_a1a1
    J_q2q2 = Jba_a2a2
    J_q0q1 = Jba_a0a1
    J_q1q2 = Jba_a1a2
    J_q0q2 = Jba_a0a2
  [../]
  [./rbed_z]
    type = Transformed111KernelOp3
    variable = A1_z
    component = 2
    order_param_x = A1_x
    order_param_y = A1_y
    order_param_z = A1_z

    f_q0 = fba_a0
    f_q1 = fba_a1
    f_q2 = fba_a2

    J_q0q0 = Jba_a0a0
    J_q1q1 = Jba_a1a1
    J_q2q2 = Jba_a2a2
    J_q0q1 = Jba_a0a1
    J_q1q2 = Jba_a1a2
    J_q0q2 = Jba_a0a2
  [../]

  [./rpp_x]
    type = Transformed111KernelOp6
    variable = P1_x
    component = 0
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]

  [./rpp_y]
    type = Transformed111KernelOp6
    variable = P1_y
    component = 1
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]

  [./rpp_z]
    type = Transformed111KernelOp6
    variable = P1_z
    component = 2
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]

 [./rpa_x]
    type = Transformed111KernelOp6
    variable = A1_x
    component = 3
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]

 [./rpa_y]
    type = Transformed111KernelOp6
    variable = A1_y
    component = 4
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]


 [./rpa_z]
    type = Transformed111KernelOp6
    variable = A1_z
    component = 5
    order_param_x = P1_x
    order_param_y = P1_y
    order_param_z = P1_z
    order_param2_x = A1_x
    order_param2_y = A1_y
    order_param2_z = A1_z

    f_q0 = frp_p0
    f_q1 = frp_p1
    f_q2 = frp_p2
    f_q3 = frp_a0
    f_q4 = frp_a1
    f_q5 = frp_a2

    J_q0q0 = Jrp_p0p0
    J_q1q1 = Jrp_p1p1
    J_q2q2 = Jrp_p2p2
    J_q3q3 = Jrp_a0a0
    J_q4q4 = Jrp_a1a1
    J_q5q5 = Jrp_a2a2
    J_q0q1 = Jrp_p0p1
    J_q1q2 = Jrp_p1p2
    J_q0q2 = Jrp_p0p2
    J_q0q3 = Jrp_p0a0
    J_q0q4 = Jrp_p0a1
    J_q0q5 = Jrp_p0a2
    J_q1q3 = Jrp_p1a0
    J_q1q4 = Jrp_p1a1
    J_q1q5 = Jrp_p1a2
    J_q2q3 = Jrp_p2a0
    J_q2q4 = Jrp_p1a1
    J_q2q5 = Jrp_p1a2
    J_q3q4 = Jrp_a0a1
    J_q3q5 = Jrp_a0a2
    J_q4q5 = Jrp_a1a2

  [../]


  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = P1_x
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable = P1_y
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_z_time]
    type = TimeDerivativeScaled
    variable = P1_z
    time_scale = 1.0
    block = '0'
  [../]

  [./a_x_time]
    type = TimeDerivativeScaled
    variable = A1_x
    time_scale = 0.01
    block = '0'
  [../]
  [./a_y_time]
    type = TimeDerivativeScaled
    variable = A1_y
    time_scale = 0.01
    block = '0'
  [../]
  [./a_z_time]
    type = TimeDerivativeScaled
    variable = A1_z
    time_scale = 0.01
    block = '0'
  [../]
[]


[Materials]
  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-2.81296 1.72351 2.24147 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./Landau_A]
    type = GenericConstantMaterial
    prop_names = 'beta1 beta11 beta12 beta111 beta112 beta123 beta1111 beta1112 beta1122 beta1123'
    prop_values = '-0.0137763 0.0000349266 0.0000498846 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./P_A_couple]
    type = GenericConstantMaterial
    prop_names = 't1111 t1122 t1212 t42111111 t24111111 t42111122 t24112222 t42112233 t24112233 t42112211 t24111122 t42111212   t42123312 t24121112 t24121233 t6211111111 t2611111111 t6211111122 t2611222222 t4411111111 t4411112222'
    prop_values = '0.012516 0.0180504 -0.036155 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '295.179 117.567 117.567 295.179 117.567 295.179 74.0701 74.0701 74.0701'

    ########################################################################################
    ##
    ## The below Euler rotation below should rotate the [001]-oriented elasticity tensor to 
    ## the [111]-orientation as it is equivalent to our S matrix. 
    ## Note that the MOOSE convention is slightly different than Mathematica's so we use 
    ## three different angles.
    ## Therefore, all of our strains will be calculated using ComputeSmallStrain in primed 
    ## coordinates [i.e. e_{||,||}, e_{1,||}, ...]
    ##
    ##   ...but then, it is important that our strains talk to our primed variables correctly. 
    ##
    ########################################################################################

    euler_angle_1 = 135.0
    euler_angle_2 = -54.735610317245346   
    euler_angle_3 = -90.0

  [../]

  [./strain]
    type = ComputeSmallStrain
    #global_strain = global_strain
  [../]

  #[./global_strain]
  #  type = ComputeGlobalStrain
  #  scalar_global_strain = global_strain
  #  global_strain_uo = global_strain_uo
  #[../]

  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./FbP1]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    polar_x = P1_x
    polar_y = P1_y
    polar_z = P1_z
  [../]

  [./FbPo]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
  [../]
  [./FbA1]
    type = RotoBulkEnergyEighth
    execute_on = 'timestep_end'
    antiphase_A_x = A1_x
    antiphase_A_y = A1_y
    antiphase_A_z = A1_z
  [../]
  [./FbAo]
    type = RotoBulkEnergyEighth
    execute_on = 'timestep_end'
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
  [../]

  [./FcPA1]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'timestep_end'
    polar_x = P1_x
    polar_y = P1_y
    polar_z = P1_z
    antiphase_A_x = A1_x
    antiphase_A_y = A1_y
    antiphase_A_z = A1_z
  [../]
  [./FcPAo]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'timestep_end'
    polar_x = Po_x
    polar_y = Po_y
    polar_z = Po_z
    antiphase_A_x = Ao_x
    antiphase_A_y = Ao_y
    antiphase_A_z = Ao_z
  [../]


  [./Felu]
    type = ElasticEnergy
    execute_on = 'timestep_end'
  [../]

  [./Ftoto]
    type = LinearCombinationPostprocessor
    pp_names = 'FbPo FbAo FcPAo'
    pp_coefs = ' 1 1 1'
    execute_on = 'timestep_end'
  
    ##########################################
    #
    # NOTE: Ferret output is in attojoules
    #
    ##########################################
  [../]
  [./Ftot1]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP1 FbA1 FcPA1'
    pp_coefs = ' 1 1 1'
    execute_on = 'timestep_end'
  
    ##########################################
    #
    # NOTE: Ferret output is in attojoules
    #
    ##########################################
  [../]
  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot1
    execute_on = 'timestep_end'
    dt = dt
  [../]
[]


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y z'
      variable = 'u1_x u1_y u1_z P1_x P1_y P1_z A1_x A1_y A1_z'
    [../]
  [../]


  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u1_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u1_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 100
    variable = u1_z
    value = 0
  [../]
[]

[UserObjects]
  #[./global_strain_uo]
  #  type = GlobalBFOMaterialRVEUserObject
  #  execute_on = 'Initial Linear Nonlinear'
  #  polar_x = P1_x
  #  polar_y = P1_y
  #  polar_z = P1_z
  #  antiphase_A_x = A1_x
  #  antiphase_A_y = A1_y
  #  antiphase_A_z = A1_z
  #[../]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 5.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -build_twosided'
    petsc_options_value = '    121            1e-10          1e-10       1e-6     bjacobi    allreduce'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.08
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 10.0

  num_steps = 10
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_P0A0_transformed_kernel_test
    elemental_as_nodal = true
  [../]
[]
