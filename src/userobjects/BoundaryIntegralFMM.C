/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

   The FMM boundary condition is contributed by X. Jiang <xikaij@imech.ac.cn> 
**/


#include "BoundaryIntegralFMM.h"
#include <iostream>
#include <fstream>

// MOOSE headers
#include "SystemBase.h"
#include "FerretConfig.h"

// libMesh headers
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/node.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/boundary_volume_solution_transfer.h"
#include "libmesh/mesh_tools.h"

// PETSc headers
#include "petscsys.h"

// ScalFMM headers
#ifdef FERRET_HAVE_SCALFMM
#include "Utils/FPoint.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Core/FFmmAlgorithmTsm.hpp"
#endif

// Global consts
const double PI_4 = 4.*M_PI;

using namespace libMesh;

registerMooseObject("FerretApp", BoundaryIntegralFMM);

template<>
InputParameters validParams<BoundaryIntegralFMM>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.set<std::string>("built_by_action") = "add_user_object";

//  params.addRequiredCoupledVar("coupled_var", "Potential from Poisson solver");

  params.addRequiredParam<Real>("cx","x-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("cy","y-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("cz","z-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("boxWidth","Width of the FMM cubic box");
  params.addRequiredParam<unsigned int>("TreeHeight","Octree height or level");

  return params;
}

BoundaryIntegralFMM::BoundaryIntegralFMM(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _cx(getParam<Real>("cx")),
    _cy(getParam<Real>("cy")),
    _cz(getParam<Real>("cz")),
    _boxWidth(getParam<Real>("boxWidth")),
    _TreeHeight(getParam<unsigned int>("TreeHeight"))
{
}

void
BoundaryIntegralFMM::initialize()
{
}

void
BoundaryIntegralFMM::execute()
{
#ifdef FERRET_HAVE_SCALFMM
  // Get a constant reference to the mesh object
  MeshBase & mesh = _subproblem.mesh().getMesh();

  // Element dimension of volume mesh
  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int dim_boundary = dim -1;

  // Extract boundary mesh
  BoundaryMesh boundary_mesh(mesh.comm() , dim_boundary);
  mesh.get_boundary_info().sync(boundary_mesh);
  boundary_mesh.print_info();

  // Re-partition boundary mesh, before adding EquationSystems
  MetisPartitioner bs_partitioner;
  bs_partitioner.partition(boundary_mesh);

  // Equation system on boundary mesh
  EquationSystems boundary_system (boundary_mesh);
  ExplicitSystem & boundary_potential = boundary_system.add_system<ExplicitSystem> ("BoundaryPotential");
  unsigned int bs_phi1 = boundary_system.get_system("BoundaryPotential").add_variable("phi1", FIRST);
  unsigned int bs_phi2 = boundary_system.get_system("BoundaryPotential").add_variable("phi2", FIRST);
  boundary_system.init();

  // Construct SolutionTransfer object. The communicator won't be
  // used directly, it just needs to be passed in to satisfy the
  // ParallelObject interface.
  BoundaryVolumeSolutionTransfer transfer(mesh.comm());

  // Get a reference to the EquationSystem
  EquationSystems & bs = _subproblem.es();
  bs.print_info();

  // Keep track of which variables are coupled so we know what we depend on
//  const std::vector<MooseVariable *> & coupled_vars = getCoupledMooseVars();

  // System on the volume mesh, 'System' is a libMesh class
  System & sys = bs.get_system("nl0");
  System & sys_aux = bs.get_system("aux0");

  // Transfer nodal values from VolumeMesh to BoundaryMesh.
//  transfer.transfer(sys.variable(0), // Note: there is only one variable in volume mesh, so variable number is 0
  transfer.transfer(sys_aux.variable(0), // Note: there is only one variable in volume mesh, so variable number is 0
                    boundary_potential.variable(bs_phi1));

  boundary_system.get_system("BoundaryPotential").solution->close();
  boundary_system.get_system("BoundaryPotential").update();

  std::vector<Number> global_potential(boundary_potential.solution->size());

  // A reference to the DofMap object for this system.
  const DofMap & dof_map_bs = boundary_potential.get_dof_map();

  // Define dof_indices holder for phi1
  std::vector<dof_id_type> dof_indices_phi1;

  // Get system number and variable numbers
  const unsigned short int        system_number = boundary_potential.number();
  const unsigned short int variable_number_phi1 = boundary_potential.variable_number("phi1");
  const unsigned short int variable_number_phi2 = boundary_potential.variable_number("phi2");

  // Get a constant reference to variable phi2, get their number of components
  const Variable&                 variable_phi1 = boundary_potential.variable(variable_number_phi1);
  const unsigned short int   variable_comp_phi1 = variable_phi1.n_components();

  const Variable &                variable_phi2 = boundary_potential.variable(variable_number_phi2);
  const unsigned short int   variable_comp_phi2 = variable_phi2.n_components();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType bs_type = dof_map_bs.variable_type(variable_number_phi1);

  // Build a Finite Element object of the specified type. Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as an UniquePtr<FEBase>. This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> bs_face (FEBase::build(dim_boundary, bs_type));

  // Quadraure rule for surface integration with dimensionality
  // one less than the dimensionality of the element.
  QGauss qface(dim_boundary, FIFTH);

  // Tell FE object to use quadrature rule.
  bs_face->attach_quadrature_rule (&qface);

  // ScalFMM
  typedef double FReal;
  // In order to be used in a template, a constant value must be initialized
  const unsigned int ORDER = 5;
  FPoint<FReal> centerOfBox( _cx, _cy, _cz );
  const unsigned int SubTreeHeight = 1;

  // Particle, Leaf, Cell, Octree
  typedef FP2PParticleContainerIndexed<FReal>                 ContainerClass;
  typedef FTypedLeaf<FReal,ContainerClass>                    LeafClass;
  typedef FTypedChebCell<FReal,ORDER>                         CellClass;
  typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>   OctreeClass;

  // Kernel D(1/r)/Dx2,D(1/r)/Dy2,D(1/r)/Dz2
  typedef FInterpMatrixKernelRx2<FReal> MatrixKernelClass1;
  typedef FInterpMatrixKernelRy2<FReal> MatrixKernelClass2;
  typedef FInterpMatrixKernelRz2<FReal> MatrixKernelClass3;

  // Three kernel classes
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass1,ORDER> KernelClass1;
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass2,ORDER> KernelClass2;
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass3,ORDER> KernelClass3;

  // Three kernels
  const   MatrixKernelClass1 MatrixKernel1;
  const   MatrixKernelClass2 MatrixKernel2;
  const   MatrixKernelClass3 MatrixKernel3;

  // Octrees 
  OctreeClass tree1(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);
  OctreeClass tree2(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);
  OctreeClass tree3(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);

  FPoint<FReal> particlePosition;
  FSize indexPart = 0;

  // Node iterator for global mesh, because here we don't partition target points.
  MeshBase::const_node_iterator           nd = boundary_mesh.nodes_begin();
  const MeshBase::const_node_iterator end_nd = boundary_mesh.nodes_end();

  // Loop over all nodes - targets
  for ( ; nd != end_nd ; ++nd){   
      const Node* node_bs = *nd;

      // Target point coords
      const Real xt = (*node_bs)(0);
      const Real yt = (*node_bs)(1);
      const Real zt = (*node_bs)(2);
      particlePosition.setPosition( xt, yt, zt );

      // Insert into trees             particleType    index     physicalValue  pot forces 
      tree1.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);
      tree2.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);
      tree3.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);

      indexPart += 1;
  }

  const FSize nbTargets = indexPart;

  // Iterator el will iterate from the first to the last element on
  // this processor. Because here we partition the source points.
  const MeshBase::const_element_iterator end_el = boundary_mesh.local_elements_end();
  MeshBase::const_element_iterator           el = boundary_mesh.local_elements_begin();

  // Loop over all sources at quadrature points in every elements.
  // ++el requires an unnecessary temporary object.
  for ( ; el != end_el ; ++el){
      // Store a pointer to the element
      const Elem* elem_bs = *el;
  
      // The Jacobian * Quadrature Weight at the quadrature points on the face.
      const std::vector<Real>& JxW_face = bs_face->get_JxW();

      // The shape function at quadrature points
      const std::vector<std::vector<Real> >& phi = bs_face->get_phi();
  
      // The XYZ locations (in physical space) of the quadrature points on the face.
      const std::vector<Point >& qface_point = bs_face->get_xyz();
  
      // Tangent direction of xi and eta, cross product to get normal
      const std::vector<RealGradient >& qface_dxyzdxi  = bs_face->get_dxyzdxi();
      const std::vector<RealGradient >& qface_dxyzdeta = bs_face->get_dxyzdeta();
 
      // Compute the shape function values on the element face.
      bs_face->reinit(elem_bs);

      // Calculate the normal vector, we only need normal vector at one quadrature point because surface element is 2D
      // Cross product
      Real nx = qface_dxyzdxi[0](1)*qface_dxyzdeta[0](2)-qface_dxyzdxi[0](2)*qface_dxyzdeta[0](1);
      Real ny = qface_dxyzdxi[0](2)*qface_dxyzdeta[0](0)-qface_dxyzdxi[0](0)*qface_dxyzdeta[0](2);
      Real nz = qface_dxyzdxi[0](0)*qface_dxyzdeta[0](1)-qface_dxyzdxi[0](1)*qface_dxyzdeta[0](0);

      // Normalize the normal vector
      Real nunit = sqrt (nx*nx + ny*ny +nz*nz);
      nx = nx / nunit;
      ny = ny / nunit;
      nz = nz / nunit;

      // Global dof_indices for this element and variable phi1
      dof_map_bs.dof_indices(elem_bs, dof_indices_phi1, variable_number_phi1);
      // Number of dof indices for phi1 on this element
      // used to loop through all nodes to calculate phi1 on quadrature point
      const unsigned int n_phi1_dofs = dof_indices_phi1.size();

      // Loop over the face quadrature points for integration.
      for (unsigned int qp=0; qp<qface.n_points(); qp++){
          // Location of source points
          const Real x_qp = qface_point[qp](0);
          const Real y_qp = qface_point[qp](1);
          const Real z_qp = qface_point[qp](2);
 
          // Value of phi1 at quadrature point
          Real phi1_qp = 0.0;
          for (unsigned int l=0; l < n_phi1_dofs; l++){
              phi1_qp += phi[l][qp] * (*boundary_potential.current_local_solution)(dof_indices_phi1[l]);
          }

          Real phys = JxW_face[qp]*phi1_qp;

          particlePosition.setPosition( x_qp , y_qp , z_qp );
          // Insert into trees             particleType    index      physicalValue  pot forces
          tree1.insert(particlePosition, FParticleType(0), indexPart,       phys*nx, 0., 0., 0., 0.);
          tree2.insert(particlePosition, FParticleType(0), indexPart,       phys*ny, 0., 0., 0., 0.);
          tree3.insert(particlePosition, FParticleType(0), indexPart,       phys*nz, 0., 0., 0., 0.);

          indexPart += 1;
      }
  // End of looping over elements
  }

  const FSize nbTotal = indexPart;

  //std::cout << "Particle insertion into octree done." << std::endl
  //          << nbTargets << " target particles, " << nbTotal-nbTargets << " source particles." << std::endl;

  // Apply kernels, here performs the compression and set M2L operators
  KernelClass1 kernel1(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel1);
  KernelClass2 kernel2(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel2);
  KernelClass3 kernel3(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel3);

  // FMM algorithm combine octree and kernel
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass1,LeafClass> FmmClass1;
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass2,LeafClass> FmmClass2;
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass3,LeafClass> FmmClass3;

  FmmClass1 algorithm1(&tree1, &kernel1);
  FmmClass2 algorithm2(&tree2, &kernel2);
  FmmClass3 algorithm3(&tree3, &kernel3);

  algorithm1.execute();
  algorithm2.execute();
  algorithm3.execute();

  // Store particle information
  struct TestParticle{
     FReal potential;
  };
  TestParticle* const particles = new TestParticle[nbTotal];

  // Get surface integral at each target point
  tree1.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] = potentials[idxPart];
      }
  });

  tree2.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] += potentials[idxPart];
      }
  });

  tree3.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] += potentials[idxPart];
      }
  });

  // Get contribution from all source points among all processors.
  boundary_mesh.comm().sum(global_potential);

  // Get BoundingBox of the boundary mesh
  BoundingBox bounding_box = MeshTools::create_bounding_box(boundary_mesh);
  Point p_min = bounding_box.min();
  Point p_max = bounding_box.max();

  // Reset nd and particle index, assign values to nodes at this processor
  indexPart = 0;
  nd = boundary_mesh.local_nodes_begin();
  const MeshBase::const_node_iterator end_nd_local = boundary_mesh.local_nodes_end();

  for ( ; nd != end_nd_local ; ++nd){
      const Node* node_bs = *nd;

      // Node's coordinates.
      const Real xt = (*node_bs)(0);
      const Real yt = (*node_bs)(1);
      const Real zt = (*node_bs)(2);

      // Dof_index for each node
      const dof_id_type node_dof_index_phi1 = node_bs->dof_number(system_number,
                                                                  variable_number_phi1,
                                                                  variable_comp_phi1-1);

      const dof_id_type node_dof_index_phi2 = node_bs->dof_number(system_number,
                                                                  variable_number_phi2,
                                                                  variable_comp_phi2-1);

      // Coefficient depends on solid angle.
      Real omega;

      // Determine point is on surface, edge or corner of the cuboid.
      if ( (xt == p_min(0))||(xt == p_max(0)) ){
         if ( (yt == p_min(1))||(yt == p_max(1)) ){
            if ( (p_min(2)<zt)&&(zt<p_max(2)) ) {omega=1./4.;} // Edge.
            else {omega=1./8.;} // Corner.
         } else {
            if ( (p_min(2)<zt)&&(zt<p_max(2)) ) {omega=1./2.;} // Surface.
            else {omega=1./4.;} // Edges.
         }
      } else {
         if ( ((yt == p_min(1))||(yt == p_max(1))) && ((zt == p_min(2))||(zt == p_max(2))) ){omega=1./4.;} // Edge.
         else {omega=1./2.;} // Surface.
      }

      // Integral value, only works for smooth surface,
      // minus SIGN get (R_target-R_source) in ScalFMM
      Real bi_value = -global_potential[node_bs->id()]/PI_4
                      -(1.-omega)*(*boundary_potential.solution)(node_dof_index_phi1);
      //std::cout << node_bs->id() << ", BI = "<< global_potential[node_bs->id()] << ". bi_value = " << bi_value
      //          << ". Original solution is " << (*boundary_potential.solution)(node_dof_index_phi1) << std::endl;

      // Boundary integral value to solution vector
      boundary_potential.solution->set(node_dof_index_phi2, bi_value);
 
      indexPart += 1;
  }

  // End of BI-FMM
  boundary_potential.solution->close();
  boundary_potential.update();

  // Transfer nodal values from BoundaryMesh to VolumeMesh
  transfer.transfer(boundary_potential.variable(bs_phi2),
                    sys_aux.variable(0));

  sys_aux.solution->close();
  sys_aux.update();
#endif
}

void
BoundaryIntegralFMM::finalize()
{
}
