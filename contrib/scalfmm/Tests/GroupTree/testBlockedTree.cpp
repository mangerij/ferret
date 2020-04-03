
// Keep in private GIT

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "../../Src/GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#endif
#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/Rotation/FRotationCellPOD.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv,
                         "Test the blocked tree.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::InputFile, LocalOptionBlocSize);

    typedef double FReal;
    static const int P = 3;
    typedef FRotationCell<FReal,P>               CellClass;
    typedef FP2PParticleContainer<FReal>          ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    typedef FRotationCellPODCore     GroupCellSymbClass;
    typedef FRotationCellPODPole<FReal,P>  GroupCellUpClass;
    typedef FRotationCellPODLocal<FReal,P> GroupCellDownClass;
    typedef FRotationCellPOD<FReal,P>      GroupCellClass;

    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;


    FTic counter;
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.bin.fma.double");

    FFmaGenericLoader<FReal> loader(filename);
    FAssertLF(loader.isOpen());

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    FP2PParticleContainer<FReal> allParticles;

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
        allParticles.push(particlePosition, physicalValue);
    }

    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.tacAndElapsed() << "s)." << std::endl;

    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    counter.tic();
    GroupOctreeClass groupedTree2(NbLevels, groupSize, &tree);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    counter.tic();
    GroupOctreeClass groupedTree3(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    groupedTree2.printInfoBlocks();
    groupedTree3.printInfoBlocks();


#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>>   GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>  GroupKernelClass;
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass,
            GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>  GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

    GroupKernelClass kernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    GroupAlgorithm algo(&groupedTree2,&kernel);

    algo.execute();

    return 0;
}
