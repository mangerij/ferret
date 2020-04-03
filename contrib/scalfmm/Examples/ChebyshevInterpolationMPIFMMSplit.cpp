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
#include "Files/FMpiSplitFmaLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "BalanceTree/FLeafBalance.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"


/// \file ChebyshevInterpolationMPIFMM
//!
//! \brief This program runs the MPI FMM with Chebyshev interpolation of 1/r kernel
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the FMM Algorithm Proc with Chebyshev Interpolation for the 1/r kernel


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    ///////// PARAMETERS HANDLING //////////////////////////////////////
    FHelpDescribeAndExit(
        argc, argv,
        "Driver for Chebyshev Interpolation kernel using MPI  (1/r kernel).\n\n"
        "usage: mpirun -np nb_proc_needed ./ChebyshevInterpolationAlgorithm [params]\n",
        FParameterDefinitions::InputFile,
        FParameterDefinitions::OctreeHeight,
        FParameterDefinitions::OctreeSubHeight,
        FParameterDefinitions::NbThreads);

    const std::string defaultFile("../Data/test20k.main.fma");
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

    using FReal = double;
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


    // Initialize values for MPI
    FMpi app(argc,argv);
    //
    // Initialize timer
    FTic time;


    // Creation of the particle loader
    FMpiSplitFmaLoader<FReal> loader(filename,app.global().processId());
    if(!loader.isOpen())
        throw std::runtime_error("Particle file couldn't be opened!") ;

    // Initialize empty oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
        if(app.global().processId() == 0){
            std::cout << "Creating & Inserting "
                      << loader.getNumberOfParticles()
                      << " particles..."
                      << std::endl;
            std::cout << "\tHeight : "
                      << TreeHeight
                      << "\tsub-height : "
                      << SubTreeHeight
                      << std::endl;
        }
        time.tic();


        FPoint<FReal> position;  // Spatial position of the particle.
        FReal physicalValue;     // Physical value of the particle.

        // Read particles from parts.
        for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
            // Read particle from file
            loader.fillParticle(&position,
                                &physicalValue);

            tree.insert(position, idxPart, physicalValue);
        }


        time.tac();
        double timeUsed = time.elapsed();
        double minTime,maxTime;
        std::cout << "Proc:" << app.global().processId()
                  << " "     << loader.getMyNumberOfParticles()
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
                if ((indexPartOrig == N1)
                    || (indexPartOrig == N2)
                    || (indexPartOrig == N3)  ) {
                    std::cout << "Proc "<< app.global().processId()
                              << " Index "<< indexPartOrig
                              <<"  potential  " << potentials[idxPart]
                              << " Pos " << posX[idxPart]
                              << " "     << posY[idxPart]
                              << " "     << posZ[idxPart]
                              << "   Forces: " << forcesX[idxPart]
                              << " "           << forcesY[idxPart]
                              << " "           << forcesZ[idxPart]
                              << std::endl;
                }
                energy += potentials[idxPart]*physicalValues[idxPart] ;
            }
        });
        FReal gloEnergy = app.global().reduceSum(energy);
        std::cout << std::endl
                  << "Proc " << app.global().processId()
                  << " Energy: " << gloEnergy
                  << std::endl;
        std::cout << std::endl
                  << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "
                  << std::endl
                  << std::endl;
    }
    // -----------------------------------------------------

    return 0;
}
