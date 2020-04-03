// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE =====
// @FUSE_MPI
// @FUSE_BLAS
// ================

#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>


#include "ScalFmmConfig.h"
#include "Containers/FOctree.hpp"
#include "Utils/FMpi.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "BalanceTree/FLeafBalance.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#ifdef SCALFMM_USE_EZTRACE
#include "eztrace.h"
#endif
/// \file
//!
//! \brief This program runs the MPI FMM with Chebyshev interpolation of 1/r kernel
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the FMM Algorithm Proc with Chebyshev Interpolation for the 1/r kernel


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    ///////// PARAMETERS HANDLING //////////////////////////////////////
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev Interpolation kernel using MPI  (1/r kernel). "
                         "Usully run using : mpirun -np nb_proc_needed ./ChebyshevInterpolationAlgorithm [params].",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
                         FParameterDefinitions::NbThreads);

    const std::string defaultFile("../Data/test20k.fma");
    const std::string filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

    //
    std::cout << "Parameters"<< std::endl
              << "      Octree Depth      " << TreeHeight    << std::endl
              << "      SubOctree depth   " << SubTreeHeight << std::endl
              << "      Input file  name: " << filename      << std::endl
              << "      Thread count :    " << NbThreads     << std::endl
              << std::endl;


    ///////// VAR INIT /////////////////////////////////////////////////

    // Initialize values for MPI
#ifdef SCALFMM_USE_EZTRACE
   eztrace_start();
#endif
    FMpi app(argc,argv);
#ifdef SCALFMM_USE_EZTRACE
   eztrace_pause();
#endif  
    //
    // Initialize timer
    FTic time;

    typedef double FReal;

    // Creation of the particle loader
    FMpiFmaGenericLoader<FReal> loader(filename,app.global());
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;


    // Begin spherical kernel
    // Accuracy
    const unsigned int ORDER = 7;
    // Typedefs
    using ContainerClass = FP2PParticleContainerIndexed<FReal>;
    using LeafClass      = FSimpleLeaf<FReal, ContainerClass>;
    using CellClass      = FChebCell<FReal,ORDER>;
    using OctreeClass    = FOctree<FReal,CellClass,ContainerClass,LeafClass>;

    using MatrixKernelClass = FInterpMatrixKernelR<FReal>;
    const MatrixKernelClass MatrixKernel;

    using KernelClass    = FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER>;
    using FmmClassProc   = FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;

    // Initialize empty oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
        if(app.global().processId() == 0){
            std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                      << " particles ..." << std::endl;
            std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
        }
        time.tic();

        /* Mock particle structure to balance the tree over the processes. */
        struct TestParticle{
            FSize index;             // Index of the particle in the original file.
            FPoint<FReal> position;  // Spatial position of the particle.
            FReal physicalValue;     // Physical value of the particle.
            /* Returns the particle position. */
            const FPoint<FReal>& getPosition(){
                return position;
            }
        };

        // Temporary array of particles read by this process.
        TestParticle* particles = new TestParticle[loader.getMyNumberOfParticles()];
        memset(particles, 0, (sizeof(TestParticle) * loader.getMyNumberOfParticles()));

        // Index (in file) of the first particle that will be read by this process.
        FSize idxStart = loader.getStart();
        std::cout << "Proc:" << app.global().processId() << " start-index: " << idxStart << std::endl;

        // Read particles from parts.
        for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
            // Store the index (in the original file) the particle.
            particles[idxPart].index = idxPart + idxStart;
            // Read particle from file
            loader.fillParticle(&particles[idxPart].position,
                                &particles[idxPart].physicalValue);
        }

        // Final vector of particles
        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;

        // Redistribute particules between processes
        FMpiTreeBuilder< FReal, TestParticle >::
            DistributeArrayToContainer(app.global(),
                                       particles,
                                       loader.getMyNumberOfParticles(),
                                       tree.getBoxCenter(),
                                       tree.getBoxWidth(),
                                       tree.getHeight(),
                                       &finalParticles,
                                       &balancer);

        // Free temporary array memory.
        delete[] particles;

        // Insert final particles into tree.
        for(FSize idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].position,
                        finalParticles[idx].index,
                        finalParticles[idx].physicalValue);
        }

        time.tac();
        double timeUsed = time.elapsed();
        double minTime,maxTime;
        std::cout << "Proc:" << app.global().processId()
                  << " "     << finalParticles.getSize()
                  << "particles have been inserted in the tree. (@Reading and Inserting Particles = "
                  << time.elapsed() << " s)."
                  << std::endl;

        MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
        MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
        if(app.global().processId() == 0){
            std::cout << "readinsert-time-min:" << minTime
                      << " readinsert-time-max:" << maxTime
                      << std::endl;
        }
    } // -----------------------------------------------------

    { // -----------------------------------------------------
        std::cout << "\nChebyshev Interpolation  FMM Proc (P="<< ORDER << ") ... " << std::endl;

        time.tic();

        // Kernels to use (pointer because of the limited size of the stack)
        KernelClass *kernels = new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        // MPI FMM algorithm
        FmmClassProc algorithm(app.global(),&tree, kernels);
        // FMM exectution
        algorithm.execute();

        time.tac();
        double timeUsed = time.elapsed();
        double minTime,maxTime;
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s)." << std::endl;
        MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
        MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
        if(app.global().processId() == 0){
            std::cout << "exec-time-min:" << minTime
                      << " exec-time-max:" << maxTime
                      << std::endl;
        }

        // Free kernels from memory
        delete kernels;
    }
    // -----------------------------------------------------
    //
    // Some output
    //
    //
    { // -----------------------------------------------------
        FSize N1=0, N2= loader.getNumberOfParticles()/2, N3= (loader.getNumberOfParticles()-1); ;
        FReal energy =0.0 ;
        //
        //   Loop over all leaves
        //
        std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
        std::cout << std::scientific;
        std::cout.precision(10) ;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const posX = leaf->getTargets()->getPositions()[0];
            const FReal*const posY = leaf->getTargets()->getPositions()[1];
            const FReal*const posZ = leaf->getTargets()->getPositions()[2];

            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
            const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

            const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const FSize indexPartOrig = indexes[idxPart];
                if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  ) {
                    std::cout << "Proc "<< app.global().processId() << " Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
                                 << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
                                    << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
                }
                energy += potentials[idxPart]*physicalValues[idxPart] ;
            }
        });
        FReal gloEnergy = app.global().reduceSum(energy);
        if(0 == app.global().processId()){
            std::cout <<std::endl << "Proc "<< app.global().processId() << " Energy: "<< gloEnergy <<std::endl;
            std::cout <<std::endl <<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;
        }
    }
    // -----------------------------------------------------

    return 0;
}
