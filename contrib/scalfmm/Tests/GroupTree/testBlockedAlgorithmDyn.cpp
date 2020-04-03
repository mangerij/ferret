
// Keep in private GIT


#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTreeDyn.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "../../Src/GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#endif
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/GroupTree/TestKernel/FTestCellPOD.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"


template <class FReal>
class FTestAttachedLeafDyn : public FGroupAttachedLeafDyn<FReal> {
    typedef FGroupAttachedLeafDyn<FReal> Parent;
public:
    using Parent::Parent;

    void init(const MortonIndex inIndex, const UnknownDescriptor<FReal> /*inParticles*/[],
                  const FSize inNbParticles, const size_t inSymbSize, const size_t inDownSize) override {
        memset(Parent::symbPart, 0, inSymbSize);
        memset(Parent::downPart, 0, inDownSize);
        *(FSize*)Parent::symbPart = inNbParticles;
        *(MortonIndex*)(Parent::symbPart+sizeof(FSize))= inIndex;
    }

    static void GetSizeFunc(const MortonIndex /*inIndex*/, const UnknownDescriptor<FReal> /*inParticles*/[],
                            const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
        *inSymbSize = sizeof(FSize) + sizeof(MortonIndex);
        *inDownSize = (inNbParticles*sizeof(long long int));
    }

    template <class ParticleClassContainer>
    static void GetSizeContainerFunc(const MortonIndex /*inIndex*/, const void* container,
                            size_t* inSymbSize, size_t* inDownSize){
        *inSymbSize = sizeof(FSize) + sizeof(MortonIndex);
        *inDownSize = (((const ParticleClassContainer*)container)->getNbParticles()*sizeof(long long int));
    }

    template<class ParticleClassContainer>
    void copyFromContainer(const MortonIndex inMindex, const ParticleClassContainer* particles){
        FAssertLF(Parent::isAttachedToSomething());
        *(FSize*)Parent::symbPart = particles->getNbParticles();
        *(MortonIndex*)(Parent::symbPart+sizeof(FSize))= inMindex;
        FAssertLF(getNbParticles() == particles->getNbParticles());
        memcpy(Parent::downPart, particles->getDataDown(), particles->getNbParticles()*sizeof(long long int));
    }

    FSize getNbParticles() const{
        return *(FSize*)Parent::symbPart;
    }

    MortonIndex getMortonIndex() const{
        return *(MortonIndex*)(Parent::symbPart+sizeof(FSize));
    }

    long long int* getDataDown(){
        return (long long int*) Parent::downPart;
    }

    const long long int* getDataDown()const{
        return (const long long int*) Parent::downPart;
    }
};


class FTestCellPODDyn {
protected:
    FTestCellPODCore* symb;
    long long int* up;
    long long int* down;

public:
    FTestCellPODDyn(unsigned char* inSymb, unsigned char* inUp,
                 unsigned char* inDown)
        : symb(reinterpret_cast<FTestCellPODCore*>(inSymb)), up(reinterpret_cast<long long int*>(inUp)),
          down(reinterpret_cast<long long int*>(inDown)){
    }
    FTestCellPODDyn()
        : symb(nullptr), up(nullptr), down(nullptr){
    }

    void release(){
        // nothing
    }

    /** To get the morton index */
    MortonIndex getMortonIndex() const {
        return symb->mortonIndex;
    }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) {
        symb->mortonIndex = inMortonIndex;
    }

    /** To get the position */
    FTreeCoordinate getCoordinate() const {
        return FTreeCoordinate(symb->coordinates[0],
                symb->coordinates[1], symb->coordinates[2]);
    }

    /** To set the position */
    void setCoordinate(const FTreeCoordinate& inCoordinate) {
        symb->coordinates[0] = inCoordinate.getX();
        symb->coordinates[1] = inCoordinate.getY();
        symb->coordinates[2] = inCoordinate.getZ();
    }

    /** To set the position from 3 FReals */
    void setCoordinate(const int inX, const int inY, const int inZ) {
        symb->coordinates[0] = inX;
        symb->coordinates[1] = inY;
        symb->coordinates[2] = inZ;
    }

    /** When doing the upward pass */
    long long int getDataUp() const {
        return (*up);
    }
    /** When doing the upward pass */
    void setDataUp(const long long int inData){
        (*up) = inData;
    }
    /** When doing the downard pass */
    long long int getDataDown() const {
        return (*down);
    }
    /** When doing the downard pass */
    void setDataDown(const long long int inData){
        (*down) = inData;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        (*down) = 0;
        (*up)   = 0;
    }

};



int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","1",1);
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles, LocalOptionBlocSize);

    typedef double FReal;

    // Initialize the types
    typedef FTestCellPODCore  GroupCellSymbClass;
    typedef FTestCellPODData  GroupCellUpClass;
    typedef FTestCellPODData  GroupCellDownClass;
    typedef FTestCellPODDyn   GroupCellClass;

    typedef FTestAttachedLeafDyn<FReal> GroupContainerClass;

    typedef FGroupTreeDyn< FReal, GroupCellClass, GroupContainerClass>  GroupOctreeClass;
#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FTestKernels< GroupCellClass, GroupContainerClass >  GroupKernelClass;
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass,
            unsigned char, unsigned char, unsigned char, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FTestKernels< GroupCellClass, GroupContainerClass >  GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

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

    std::unique_ptr<UnknownDescriptor<FReal>[]> allParticles(new UnknownDescriptor<FReal>[loader.getNumberOfParticles()]);
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
#ifndef LOAD_FILE
        loader.fillParticle(&particlePosition);
#else
        FReal ph;
        loader.fillParticle(&particlePosition, &ph);
#endif
        allParticles[idxPart].pos = (particlePosition);
        tree.insert(particlePosition);
    }

    std::unique_ptr<size_t[]> cellSymbSizePerLevel(new size_t[NbLevels]);
    std::unique_ptr<size_t[]> cellUpSizePerLevel(new size_t[NbLevels]);
    std::unique_ptr<size_t[]> cellDownSizePerLevel(new size_t[NbLevels]);
     for(int idx = 0 ; idx < NbLevels ; ++idx){
         cellSymbSizePerLevel[idx] = sizeof(GroupCellSymbClass);
         cellUpSizePerLevel[idx] = sizeof(GroupCellUpClass);
         cellDownSizePerLevel[idx] = sizeof(GroupCellDownClass);
     }

    // Put the data into the tree
//    GroupOctreeClass groupedTree(NbLevels, groupSize, &tree,
//              cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                 [](const MortonIndex inIndex, const void* inParticles,
//                                    size_t* inSymbSize, size_t* inDownSize) {
//                                        GroupContainerClass::GetSizeContainerFunc<ContainerClass>(
//                                                    inIndex, inParticles, inSymbSize, inDownSize);
//                                 },
//     [](const MortonIndex /*mindex*/,
//                        unsigned char* symbBuff, const size_t /*symbSize*/,
//                        unsigned char* upBuff, const size_t /*upSize*/,
//                        unsigned char* downBuff, const size_t /*downSize*/,
     //                         const int /*inLevel*/){
//         GroupCellClass cell(symbBuff, upBuff, downBuff);
//     });
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
                                 allParticles.get(), loader.getNumberOfParticles(),
                                 [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
                                    const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
        GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
    },
    [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
                  const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
    unsigned char* downBuffer, const size_t inDownSize){
        GroupContainerClass leaf(symbBuffer, downBuffer);
        leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
    },
    [](const MortonIndex /*mindex*/,
                       unsigned char* symbBuff, const size_t /*symbSize*/,
                       unsigned char* upBuff, const size_t /*upSize*/,
                       unsigned char* downBuff, const size_t /*downSize*/,
                       const int /*inLevel*/){
        GroupCellClass cell(symbBuff, upBuff, downBuff);
    });
//    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
//                                 cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                  allParticles.get(), loader.getNumberOfParticles(),
//                                  [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                                     const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
//                                         GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//                                  },
//    [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                  const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
//    unsigned char* downBuffer, const size_t inDownSize){
//        GroupContainerClass leaf(symbBuffer, downBuffer);
//        leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//    },
//    [](const MortonIndex /*mindex*/,
//                       unsigned char* symbBuff, const size_t /*symbSize*/,
//                       unsigned char* upBuff, const size_t /*upSize*/,
//                       unsigned char* downBuff, const size_t /*downSize*/,
    //                         const int /*inLevel*/){
//        GroupCellClass cell(symbBuff, upBuff, downBuff);
//    }
//                                 false, true);
//    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
//                                cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                 allParticles.get(), loader.getNumberOfParticles(),
//                                 [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                                    const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
//                                        GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//                                 },
//    [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                  const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
//    unsigned char* downBuffer, const size_t inDownSize){
//        GroupContainerClass leaf(symbBuffer, downBuffer);
//        leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//    },
//    [](const MortonIndex /*mindex*/,
//                       unsigned char* symbBuff, const size_t /*symbSize*/,
//                       unsigned char* upBuff, const size_t /*upSize*/,
//                       unsigned char* downBuff, const size_t /*downSize*/,
//                         const int /*inLevel*/){
//        GroupCellClass cell(symbBuff, upBuff, downBuff);
//    }
//                                false, true, 0.2);
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
    groupalgo.execute();

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

