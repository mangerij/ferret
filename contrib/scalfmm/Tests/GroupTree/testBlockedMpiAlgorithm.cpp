
// Keep in private GIT
// @FUSE_MPI
// @FUSE_STARPU

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "../../Src/GroupTree/TestKernel/FTestCellPOD.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Containers/FCoordinateComputer.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"




int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize);
    typedef double FReal;
    // Initialize the types
    typedef FTestCellPODCore  GroupCellSymbClass;
    typedef FTestCellPODData  GroupCellUpClass;
    typedef FTestCellPODData  GroupCellDownClass;
    typedef FTestCellPOD      GroupCellClass;

    typedef FGroupTestParticleContainer<FReal>                                     GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper> GroupAlgorithm;


    FMpi mpiComm(argc, argv);
    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
    const FSize totalNbParticles = (NbParticles*mpiComm.global().processCount());

    // Load the particles
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

    // Fill the particles
    struct TestParticle{
        FPoint<FReal> position;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    std::unique_ptr<TestParticle[]> particles(new TestParticle[loader.getNumberOfParticles()]);
    memset(particles.get(), 0, sizeof(TestParticle) * loader.getNumberOfParticles());
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position);
    }
    // Sort in parallel
    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder<FReal, TestParticle >::DistributeArrayToContainer(mpiComm.global(),
                                                                particles.get(),
                                                                loader.getNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),
                                                                NbLevels,
                                                                &myParticles,
                                                                &balancer);

    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        allParticles.push(myParticles[idxPart].position);
    }

    // Each proc need to know the righest morton index
    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                NbLevels,
                myParticles[myParticles.getSize()-1].position );
    const MortonIndex myLeftLimite = host.getMortonIndex(NbLevels-1);
    MortonIndex leftLimite = -1;
    if(mpiComm.global().processId() != 0){
        FMpi::Assert(MPI_Recv(&leftLimite, sizeof(leftLimite), MPI_BYTE,
                              mpiComm.global().processId()-1, 0,
                              mpiComm.global().getComm(), MPI_STATUS_IGNORE), __LINE__);
    }
    if(mpiComm.global().processId() != mpiComm.global().processCount()-1){
        FMpi::Assert(MPI_Send(const_cast<MortonIndex*>(&myLeftLimite), sizeof(myLeftLimite), MPI_BYTE,
                              mpiComm.global().processId()+1, 0,
                              mpiComm.global().getComm()), __LINE__);
    }
    FLOG(std::cout << "My last index is " << leftLimite << "\n");
    FLOG(std::cout << "My left limite is " << myLeftLimite << "\n");


    // Put the data into the tree
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 &allParticles, true, leftLimite);
    groupedTree.printInfoBlocks();

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel);
    groupalgo.execute();

    std::cout << "Wait Others... " << std::endl;
    mpiComm.global().barrier();

    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != totalNbParticles-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (totalNbParticles-1) << ") at index " << cell.getMortonIndex() << "\n";
            }
        }
    });



    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    {
        // Usual octree
        OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxProc = 0 ; idxProc < mpiComm.global().processCount() ; ++idxProc){
            FRandomLoader<FReal> loaderAll(NbParticles, 1.0, FPoint<FReal>(0,0,0), idxProc);
            for(FSize idxPart = 0 ; idxPart < loaderAll.getNumberOfParticles() ; ++idxPart){
                FPoint<FReal> pos;
                loaderAll.fillParticle(&pos);
                tree.insert(pos);
            }
        }
        // Usual algorithm
        KernelClass kernels;            // FTestKernels FBasicKernels
        FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
        algo.execute();

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
    }

    return 0;
}

