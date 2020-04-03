// ===================================================================================
// Copyright ScalFmm 2013 INRIA
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
/// @author Pierre Blanchard
// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// @FUSE_FFT
// ================
// Keep in private GIT


#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Kernels/Uniform/FUnifCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"
#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the Uniform kernel
 */

// Simply create particles and try the kernel
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Test Uniform kernel with MPI and compare it with the direct computation.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::NbThreads,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile);

    typedef double FReal;
    const unsigned int ORDER = 7;
    const FReal epsilon = FReal(1e-7);

    typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;

    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    typedef FUnifCell<FReal,ORDER> CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

    typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

    FMpi app(argc,argv);

    const char* const filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

    std::cout << ">> This executable has to be used to test Proc Uniform Algorithm. \n";


#ifdef _OPENMP
    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
    std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
    
    std::cout << "Opening : " <<filename << "\n" << std::endl;
    // init timer
    FTic time;

    // Create Matrix Kernel
    const MatrixKernelClass MatrixKernel;

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    // open particle file
    FMpiFmaGenericLoader<FReal> loader(filename,app.global());
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

    OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());
    time.tic();
    TestParticle* particles = new TestParticle[loader.getMyNumberOfParticles()];
    memset(particles,0,(unsigned int) (sizeof(TestParticle)* loader.getMyNumberOfParticles()));
    for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
    }
    FVector<TestParticle> finalParticles;
    FLeafBalance balancer;

    FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                loader.getMyNumberOfParticles(),
                                                                tree.getBoxCenter(),
                                                                tree.getBoxWidth(),tree.getHeight(),
                                                                &finalParticles, &balancer);
    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getMyNumberOfParticles() << " particles ..." << std::endl;
        std::cout << "For a total of " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
        time.tic();

        for(FSize idxPart = 0 ; idxPart < finalParticles.getSize() ; ++idxPart){
            // put in tree
            tree.insert(finalParticles[idxPart].position, idxPart, finalParticles[idxPart].physicalValue);
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
        std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ",EPS="<< epsilon <<") ... " << std::endl;
        time.tic();
        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        FmmClass algorithm(app.global(),&tree, &kernels);
        algorithm.execute();
        time.tac();
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    { // -----------------------------------------------------
        std::cout << "\nError computation ... " << std::endl;
        FReal potential;
        FReal fx, fy, fz;
        { // Check that each particle has been summed with all other
            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    potential += potentials[idxPart];
                    fx += forcesX[idxPart];
                    fy += forcesY[idxPart];
                    fz += forcesZ[idxPart];
                }
            });
        }

        // Print for information
        std::cout << "Potential " << potential << std::endl;
        std::cout << "Fx " << fx << std::endl;
        std::cout << "Fy " << fy << std::endl;
        std::cout << "Fz " << fz << std::endl;

    } // end Unif kernel


    return 0;
}
