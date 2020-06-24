import numpy as np
import math
import random
import csv

# grains = int(input('How many grains are there?: '))
grains = 100
q = 0
# inputs_needed = int(input('How many input files are needed?:')) + 1
inputs_needed = 1

while q < inputs_needed:
    list_of_orientations = []

    for line in open("/Users/lukaszkuna/projects/ferret/trigonal_grains.txt"):
        csv_row = line.split()
        list_of_orientations.append(csv_row)


    blocks = np.arange(1,grains + 1)

    file = open("Alumina_{}-100_z.i".format(grains),"w")

    file.write("""

    [Mesh]
      file = 100grains.e
    []

    [MeshModifiers]
      [./add_sidesets]
        type = SideSetsFromNormals
        normals = '1  0  0
                  -1  0  0
                   0  1  0
                   0 -1  0
                   0  0  1
                   0  0 -1'
        fixed_normal = true
        new_boundary = 'right left front back top bottom'
        variance = 0.5
      [../]
    []

    [GlobalParams]
      potential_E_int = potential_int
      displacements = 'u_x u_y u_z'
      n_a = 1.76530
      n_b = 1.76530
      n_g = 1.75730

    []

    [Variables]
      [./u_x]
      [../]
      [./u_y]
      [../]
      [./u_z]
      [../]
      [./potential_int]
        order = FIRST
        family = LAGRANGE
      [../]
    []


    [AuxVariables]
      [./disp_x]
      [../]
      [./disp_y]
      [../]
      [./disp_z]
      [../]
      [./stress_xx]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./stress_yy]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./stress_xy]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./stress_xz]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./stress_zz]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./stress_yz]
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
      [./strain_xy]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./strain_xz]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./strain_zz]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./strain_yz]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dn_1]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dn_2]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dn_3]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dne_1]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dne_2]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./dne_3]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./n_1]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./n_2]
        order = CONSTANT
        family = MONOMIAL
      [../]

      [./n_3]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [./n_avg]
        order = CONSTANT
        family = MONOMIAL
      [../]
    []


    [AuxKernels]
      [./matl_e11]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 0
        index_j = 0
        variable = strain_xx
        execute_on = 'timestep_end'
      [../]
      [./matl_e12]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 0
        index_j = 1
        variable = strain_xy
        execute_on = 'timestep_end'
      [../]
      [./matl_e13]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 0
        index_j = 2
        variable = strain_xz
        execute_on = 'timestep_end'
      [../]
      [./matl_e22]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 1
        index_j = 1
        variable = strain_yy
        execute_on = 'timestep_end'
      [../]
      [./matl_e23]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 1
        index_j = 2
        variable = strain_yz
        execute_on = 'timestep_end'
      [../]
      [./matl_e33]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        index_i = 2
        index_j = 2
        variable = strain_zz
        execute_on = 'timestep_end'
      [../]
      [./matl_s11]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 0
        index_j = 0
        variable = stress_xx
        execute_on = 'timestep_end'
      [../]
      [./matl_s12]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 0
        index_j = 1
        variable = stress_xy
        execute_on = 'timestep_end'
      [../]
      [./matl_s13]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 0
        index_j = 2
        variable = stress_xz
        execute_on = 'timestep_end'
      [../]
     [./matl_s22]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 1
        index_j = 1
        variable = stress_yy
        execute_on = 'timestep_end'
      [../]
      [./matl_s23]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 1
        index_j = 2
        variable = stress_yz
        execute_on = 'timestep_end'
      [../]
      [./matl_s33]
        type = RankTwoAux
        rank_two_tensor = stress
        index_i = 2
        index_j = 2
        variable = stress_zz
        execute_on = 'timestep_end'
      [../]

      [./dn_s1]
        type = ChangeInRefractiveIndex
        index_i = 0
        index_j = 0
        index_k = 0
        index_l = 0
        variable = dn_1
        execute_on = 'timestep_end'
      [../]

      [./dn_s2]
        type = ChangeInRefractiveIndex
        index_i = 1
        index_j = 1
        index_k = 1
        index_l = 1
        variable = dn_2
        execute_on = 'timestep_end'
      [../]

      [./dn_s3]
        type = ChangeInRefractiveIndex
        index_i = 2
        index_j = 2
        index_k = 2
        index_l = 2
        variable = dn_3
        execute_on = 'timestep_end'
      [../]

      [./dne_s1]
        type = ChangeInRefractiveIndexElectro
        index_i = 0
        index_j = 0
        index_k = 0
        index_l = 0
        variable = dne_1
        execute_on = 'timestep_end'
      [../]

      [./dne_s2]
        type = ChangeInRefractiveIndexElectro
        index_i = 1
        index_j = 1
        index_k = 1
        index_l = 1
        variable = dne_2
        execute_on = 'timestep_end'
      [../]

      [./dne_s3]
        type = ChangeInRefractiveIndexElectro
        index_i = 2
        index_j = 2
        index_k = 2
        index_l = 2
        variable = dne_3
        execute_on = 'timestep_end'
      [../]

      [./n_1_c]
        type = RefractiveIndex
        variable = n_1
        elasto = true
        electro = false
        index_j = 0
        index_k = 0
        var1 = dn_1
        execute_on = 'timestep_end'
      [../]

      [./n_2_c]
        type = RefractiveIndex
        variable = n_2
        elasto = true
        electro = false
        index_j = 1
        index_k = 1
        var1 = dn_2
        execute_on = 'timestep_end'
      [../]

      [./n_3_c]
        type = RefractiveIndex
        variable = n_3
        elasto = true
        electro = false
        index_j = 2
        index_k = 2
        var1 = dn_3
        execute_on = 'timestep_end'
      [../]

    []

    [Materials]
          [./strain_1]
            type = ComputeSmallStrain
          [../]
          [./stress_1]
            type = ComputeLinearElasticStress
          [../]

    """.format(grains))

    for x in blocks:
        angle1 = list_of_orientations[x][0]
        print(angle1)
        angle2 = list_of_orientations[x][1]
        print(angle2)
        angle3 = list_of_orientations[x][2]
        print(angle3)
        file.write("""
        [./elasticity_tensor_{}]
          type = ComputeElasticityTensor
          fill_method = symmetric21
          #w-ZnO    C1111 C1122 C1133 C2222 C2233  C3333 C2323 C1313 C1212
          C_ijkl = '452 150 107 20 0 0 452 107 -20 0 0 454 0 0 0 132 0 0 132 20 151'
          block = '{}'
          euler_angle_1 = {}
          euler_angle_2 = {}
          euler_angle_3 = {}
        [../]

        [./Piezo_block_{}]
            type = ComputePiezostrictiveTensor
            fill_method = general
            block = '{}'
            e_ijk = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
            euler_angle_1 = {}
            euler_angle_2 = {}
            euler_angle_3 = {}
        [../]

        [./elastooptic_tensor_{}]
          type = ComputeElastoopticTensor
          fill_method = symmetric9
          #        C1111 C1122 C1133 C2222 C2233 C3333 C2323 C1313 C1212
          P_mnkl = '-0.23 -0.03 0.00 -0.23 0.00 -0.20 -0.10 -0.10 -0.10'
          block = '{}'
          euler_angle_1 = {}
          euler_angle_2 = {}
          euler_angle_3 = {}
        [../]


        [./electrooptic_tensor_{}]
          type = ComputeElectroopticTensor
          fill_method = general #For ZnO (32)
          # Fill Method (see Nye page 113)
          # r111, r112, r113, r121, r122, r123, r131, r132, r133,
          # r211, r212, r213, r221, r222, r223, r231, r232, r233,
          # r311, r312, r313, r321, r322, r323, r331, r332, r333,
          r_ijk ='0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
          block = '{}'
          euler_angle_1 = {}
          euler_angle_2 = {}
          euler_angle_3 = {}
        [../]

        [./beta_tensor_{}]
          type = ComputeIndicatrix
          block = '{}'
          euler_angle_1 = {}
          euler_angle_2 = {}
          euler_angle_3 = {}
        [../]

        [./delta_beta_tensor_{}]
          type = ComputeDeltaIndicatrix
          block = '{}'
        [../]

        [./delta_beta_tensor_e_{}]
          type = ComputeDeltaIndicatrixElectro
          block = '{}'
          potential_int = potential_int
        [../]
            """.format(x,x,angle1,angle2,angle3,
                       x,x,angle1,angle2,angle3,
                       x,x,angle1,angle2,angle3,
                       x,x,angle1,angle2,angle3,
                       x,x,angle1,angle2,angle3,
                       x,x,
                       x,x))


    file.write("""
    []


    [Kernels]
      #Elastic problem
      [./TensorMechanics]
      #This is an action block
      [../]

      [./ElectroStats1]
        type = Electrostatics
        variable = potential_int
        permittivity = 0.101
      [../]
    []


    [BCs]
      # [./top_pot]
      #   type = DirichletBC
      #   variable = u_z
      #   boundary = 'top'
      #   value = -4.0    ## Negative Compression Positive Tension
      #   #value = 0.0
      # [../]
      #
      # [./bot_pot]
      #   type = DirichletBC
      #   variable = u_z
      #   boundary = 'bottom'
      #   value = 4.0
      #   #value = 0.0
      # [../]

    [./left_pot]
      type = DirichletBC
      variable = u_x
      boundary = 'left'
      value = 0.5
    [../]

    [./right_pot]
      type = DirichletBC
      variable = u_x
      boundary = 'right'
      value = -0.5 ## Negative Compression Positive Tension
      #value = 0
    [../]

    # [./front_pot]
    #   type = DirichletBC
    #   variable = u_y
    #   boundary = 'front'
    #   value = 0
    # [../]
    #
    # [./back_pot]
    #   type = DirichletBC
    #   variable = u_y
    #   boundary = 'back'
    #   value = 0.0
    # [../]
    []

    [Postprocessors]

    """)

    for y in blocks:
        file.write("""
    [./n_x_{}]
      type = ElementAverageValue
      variable = n_1 #Post#
      execute_on = 'timestep_end'
      block = '{}'
    [../]

        """.format(y,y))
    for c in blocks:
        file.write("""
    [./n_y_{}]
      type = ElementAverageValue
      variable = n_2 #Post#
      execute_on = 'timestep_end'
      block = '{}'
    [../]
        """.format(c,c))

    for w in blocks:
        file.write("""
    [./n_z_{}]
      type = ElementAverageValue
      variable = n_3 #Post#
      execute_on = 'timestep_end'
      block = '{}'
    [../]
        """.format(w,w))




    file.write("""

        [./n1]
          type = ElementAverageValue
          variable = n_1
          execute_on = 'timestep_end'
        [../]
        [./n2]
          type = ElementAverageValue
          variable = n_2
          execute_on = 'timestep_end'
        [../]
        [./n3]
          type = ElementAverageValue
          variable = n_3
          execute_on = 'timestep_end'
        [../]

        [./n_avg_0]
          type = ElementAverageValue
          variable = n_avg
          execute_on = 'timestep_end'
        [../]

        [./dn33_elasto]
          type = ElementAverageValue
          variable = dn_3
          execute_on = 'timestep_end'
        [../]
        [./dn22_elasto]
          type = ElementAverageValue
          variable = dn_2
          execute_on = 'timestep_end'
        [../]
        [./dn11_elasto]
          type = ElementAverageValue
          variable = dn_1
          execute_on = 'timestep_end'
        [../]
        [./strain_xx]
           type = ElementAverageValue
           variable = strain_xx
           execute_on = 'timestep_end'
        [../]
        [./strain_yy]
           type = ElementAverageValue
           variable = strain_yy
           execute_on = 'timestep_end'
        [../]
        [./strain_zz]
           type = ElementAverageValue
           variable = stress_zz
           execute_on = 'timestep_end'
        [../]
        [./strain_xy]
           type = ElementAverageValue
           variable = strain_xy
           execute_on = 'timestep_end'
        [../]
        [./strain_xz]
           type = ElementAverageValue
           variable = strain_xz
           execute_on = 'timestep_end'
        [../]
        [./strain_yz]
           type = ElementAverageValue
           variable = strain_yz
           execute_on = 'timestep_end'
        [../]
        []

      [Problem]
        null_space_dimension = 6
      []



      [Preconditioning]
        [./smp]
          type = SMP
          full = true
          petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -pc_hypre_type'
          petsc_options_value = '    250              1e-10      1e-8      1e-6      hypre       boomeramg '
        [../]
      []


      [Executioner]
        type = Steady
        solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
      []


      [Outputs]
        print_linear_residuals = false
        print_perf_log = true
        [./out]
          type = Exodus
          execute_on = 'timestep_end'
          file_base = out_alumina_x_{}
          elemental_as_nodal = true
        [../]
        [./CSV]
          type = CSV
          file_base = out_alumina_x_{}
        [../]
      []

    """.format(grains,grains))


    file.close()
    q = q + 1
