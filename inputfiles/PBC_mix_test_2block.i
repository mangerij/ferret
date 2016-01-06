#debug file to find how to implement PBC at end caps of cylinder
# translation vector doesnt work for regular quasi-2D pbc but auto_direction does
# for irregular auto_direction complains about being irregular so try translation
# should check auto_direction for quasi-2D for a non-generated mesh (perhaps a cheap superlattice) 
#-- doesn't work.. auto_direction NEEDs generated mesh

#distribution serial does not help -- perhaps it has to do with the LEVEL 1 test.. ie it seems that is what 
# libmesh is complaining about

# check element type...nope... error still happens for quads or tets

#Error:
#Assertion `s_neigh != libMesh::invalid_uint' failed.
#s_neigh = 4294967295
#libMesh::invalid_uint = 4294967295

#This means that we need s_neigh to not equal libMesh::invalid_uint but it does so we exit the program.. 
#what does this mean? What are these two objects?

#one note is that perhaps we need to use the transform functions as in the trapezoid example.. this is more akin to what we are doing -> but perhaps 3D isn't accessible there?

# transform functions dont work either. 


[Mesh]
  file = qexodus_twoblock.e
  distribution = SERIAL
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./dot]
    type = TimeDerivative
    variable = u
  [../]
[]


[BCs]
  [./Periodic]
    [./zz]
      variable = u
      auto_direction = 'z'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.5
  num_steps = 10
  solve_type = NEWTON
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
