// Keep in private GIT

// @FUSE_STARPU
// @FUSE_OPENCL

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"

#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"

#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUOpenClWrapper.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#include "../../Src/GroupTree/TestKernel/FTestOpenCLCode.hpp"
#include "../../Src/GroupTree/TestKernel/FTestCellPOD.hpp"

#include "../../Src/GroupTree/OpenCl/FOpenCLDeviceWrapper.hpp"

int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","0",1);
    setenv("STARPU_NOPENCL","1",1);
    setenv("STARPU_OPENCL_ONLY_ON_CPUS","1",1);

    setenv("STARPU_DISABLE_ASYNCHRONOUS_OPENCL_COPY","1",1);
    setenv("STARPU_OPENCL_PIPELINE","0",1);
    setenv("STARPU_OPENCL_ON_CPUS","1",1);

    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles."
                         "Usually run with STARPU_NCPU=0 STARPU_NOPENCL=1 STARPU_OPENCL_ONLY_ON_CPUS=1 ./Tests/Release/testBlockedWithOpenCLAlgorithm",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize);

    typedef double FReal;

    typedef FTestCellPODCore  GroupCellSymbClass;
    typedef FTestCellPODData  GroupCellUpClass;
    typedef FTestCellPODData  GroupCellDownClass;
    typedef FTestCellPOD      GroupCellClass;

    typedef FGroupTestParticleContainer<FReal>                                     GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
    typedef FStarPUAllCpuOpenCLCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper
    #ifdef SCALFMM_ENABLE_CUDA_KERNEL
        , FStarPUCudaWrapper<KernelClass, FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<>>
    #endif
        , FStarPUOpenClWrapper<GroupKernelClass, FOpenCLDeviceWrapper<GroupKernelClass, FTestOpenCLCode<FReal> > >
         > GroupAlgorithm;

    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

//#define LOAD_FILE
#ifndef LOAD_FILE
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), 0);
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader<FReal> loader(filename);
#endif
    FAssertLF(loader.isOpen());

    // Usual octree
    OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());

    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
#ifndef LOAD_FILE
        loader.fillParticle(&particlePosition);
#else
        FReal ph;
        loader.fillParticle(&particlePosition, &ph);
#endif
        allParticles.push(particlePosition);
        tree.insert(particlePosition);
    }

    // Put the data into the tree
    //GroupOctreeClass groupedTree(NbLevels, groupSize, &tree);
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();

    // Check tree structure at leaf level
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass gcell, GroupContainerClass* gleaf){
        const ContainerClass* src = tree.getLeafSrc(gcell.getMortonIndex());
        if(src == nullptr){
            std::cout << "[PartEmpty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
        }
        else {
            if(src->getNbParticles() != gleaf->getNbParticles()){
                std::cout << "[Part] Nb particles is different at index " << gcell.getMortonIndex() << " is " << gleaf->getNbParticles() << " should be " << src->getNbParticles() << "\n";
            }
        }
    });

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);
    groupalgo.execute(); // FFmmP2M TODO

    // Usual algorithm
    KernelClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    // Validate the result
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        if(cell.getDataUp() != nbPartsInLeaf){
            std::cout << "[P2M] Error a Cell has " << cell.getDataUp() << " (it should be " << nbPartsInLeaf << ")\n";
        }
    });
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != loader.getNumberOfParticles()-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (loader.getNumberOfParticles()-1) << ") at index " << cell.getMortonIndex() << "\n";
            }
        }
    });
    // Compare the results
    groupedTree.forEachCellWithLevel([&](GroupCellClass gcell, const int level){
        const CellClass* cell = tree.getCell(gcell.getMortonIndex(), level);
        if(cell == nullptr){
            std::cout << "[Empty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
        }
        else {
            if(gcell.getDataUp() != cell->getDataUp()){
                std::cout << "[Up] Up is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataUp() << " should be " << cell->getDataUp() << "\n";
            }
            if(gcell.getDataDown() != cell->getDataDown()){
                std::cout << "[Down] Down is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataDown() << " should be " << cell->getDataDown() << "\n";
            }
        }
    });

    return 0;
}
