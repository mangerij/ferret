// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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
#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H

#include "../Utils/FMpi.hpp"
#include "../Utils/FQuickSortMpi.hpp"
#include "../Utils/FBitonicSort.hpp"
#include "../Utils/FTic.hpp"

#include "../Utils/FMemUtils.hpp"

#include "../Containers/FVector.hpp"

#include "../BalanceTree/FLeafBalance.hpp"
#include "../BalanceTree/FEqualize.hpp"

#include "../Containers/FCoordinateComputer.hpp"

/**
 * This class manage the loading of particles for the mpi version.
 * It work in several steps.
 * First it load the data from a file or an array and sort them amon the MPI processes.
 * Then, it carrefully manage if a leaf is shared by multiple processes.
 * Finally it balances the data using an external interval builder.
 *
 */
template<class FReal, class ParticleClass>
class FMpiTreeBuilder{
private:
    /** To keep the leaves information after the sort */
    struct LeafInfo {
        MortonIndex mindex;
        FSize nbParts;
        FSize startingPoint;
    };

public:
    /** What sorting algorithm to use */
    enum SortingType{
        QuickSort,
        BitonicSort,
    };


    /**
     * A particle may not have a MortonIndex Method (set/get morton index)
     * But in this algorithm they are sorted based on their morton indexes.
     * So an IndexedParticle is storing a real particle + its index.
     */
    struct IndexedParticle{
    public:
        MortonIndex index;
        ParticleClass particle;

        operator MortonIndex() const {
            return this->index;
        }
    };

    //////////////////////////////////////////////////////////////////////////
    // Methods to sort the particles
    //////////////////////////////////////////////////////////////////////////


    /** Get an array of particles sorted from their morton indexes */
    template <class LoaderClass>
    static void GetSortedParticlesFromLoader( const FMpi::FComm& communicator, LoaderClass& loader, const SortingType sortingType,
                                              const int TreeHeight, IndexedParticle**const outputSortedParticles, FSize* const outputNbParticlesSorted){
        // Allocate the particles array
        IndexedParticle*const originalParticlesUnsorted = new IndexedParticle[loader.getNumberOfParticles()];
        FMemUtils::memset(originalParticlesUnsorted, 0, sizeof(IndexedParticle) * loader.getNumberOfParticles());

        FPoint<FReal> boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
        FTreeCoordinate host;
        const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / FReal(1 << (TreeHeight - 1) );

        // Fill the array and compute the morton index
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(originalParticlesUnsorted[idxPart].particle);
            host.setX( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getX() - boxCorner.getX(), loader.getBoxWidth(), boxWidthAtLeafLevel,
                                           TreeHeight ));
            host.setY( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getY() - boxCorner.getY(), loader.getBoxWidth(), boxWidthAtLeafLevel,
                                           TreeHeight ));
            host.setZ( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), loader.getBoxWidth(), boxWidthAtLeafLevel,
                                           TreeHeight ));

            originalParticlesUnsorted[idxPart].index = host.getMortonIndex(TreeHeight - 1);
        }

        // Sort particles
        if(sortingType == QuickSort){
            FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(originalParticlesUnsorted, loader.getNumberOfParticles(), *outputSortedParticles, *outputNbParticlesSorted,communicator);
            delete [] (originalParticlesUnsorted);
        }
        else {
            FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( originalParticlesUnsorted, loader.getNumberOfParticles(), communicator );
            *outputSortedParticles = originalParticlesUnsorted;
            *outputNbParticlesSorted = loader.getNumberOfParticles();
        }
    }

    /** Get an array of particles sorted from their morton indexes */
    static void GetSortedParticlesFromArray( const FMpi::FComm& communicator, const ParticleClass inOriginalParticles[], const FSize originalNbParticles, const SortingType sortingType,
                                             const FPoint<FReal>& centerOfBox, const FReal boxWidth,
                                             const int TreeHeight, IndexedParticle**const outputSortedParticles, FSize* const outputNbParticlesSorted){
        // Allocate the particles array
        IndexedParticle*const originalParticlesUnsorted = new IndexedParticle[originalNbParticles];
        FMemUtils::memset(originalParticlesUnsorted, 0, sizeof(IndexedParticle) * originalNbParticles);

        FPoint<FReal> boxCorner(centerOfBox - (boxWidth/2));
        FTreeCoordinate host;
        const FReal boxWidthAtLeafLevel = boxWidth / FReal(1 << (TreeHeight - 1) );

        FLOG(FTic counterTime);

        // Fill the array and compute the morton index
        for(FSize idxPart = 0 ; idxPart < originalNbParticles ; ++idxPart){
            originalParticlesUnsorted[idxPart].particle = inOriginalParticles[idxPart];
            host.setX( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel,
                                           TreeHeight ));
            host.setY( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel,
                                           TreeHeight ));
            host.setZ( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesUnsorted[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel,
                                           TreeHeight ));

            originalParticlesUnsorted[idxPart].index = host.getMortonIndex(TreeHeight - 1);
        }

        FLOG( FLog::Controller << "\tPrepare particles ("  << counterTime.tacAndElapsed() << "s)\n" );

        // Sort particles
        if(sortingType == QuickSort){
            FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(originalParticlesUnsorted, originalNbParticles, outputSortedParticles, outputNbParticlesSorted,communicator);
            delete [] (originalParticlesUnsorted);
        }
        else {
            FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( originalParticlesUnsorted, originalNbParticles, communicator );
            *outputSortedParticles = originalParticlesUnsorted;
            *outputNbParticlesSorted = originalNbParticles;
        }
    }


    //////////////////////////////////////////////////////////////////////////
    // To merge the leaves
    //////////////////////////////////////////////////////////////////////////

    static void MergeSplitedLeaves(const FMpi::FComm& communicator, IndexedParticle* workingArray, FSize* workingSize,
                                   FSize ** leavesOffsetInParticles, ParticleClass** particlesArrayInLeafOrder, FSize* const leavesSize){
        const int myRank = communicator.processId();
        const int nbProcs = communicator.processCount();

        FVector<LeafInfo> leavesInfo;
        { // Get the information of the leaves
            leavesInfo.clear();
            if((*workingSize)){
                leavesInfo.push({workingArray[0].index, 1, 0});
                for(FSize idxPart = 1 ; idxPart < (*workingSize) ; ++idxPart){
                    if(leavesInfo.data()[leavesInfo.getSize()-1].mindex == workingArray[idxPart].index){
                        leavesInfo.data()[leavesInfo.getSize()-1].nbParts += 1;
                    }
                    else{
                        leavesInfo.push({workingArray[idxPart].index, 1, idxPart});
                    }
                }
            }
        }

        if(nbProcs != 1){
            // Some leaf might be divived on several processes, we should move them to the first process
            const MortonIndex noDataFlag = std::numeric_limits<MortonIndex>::max();

            LeafInfo borderLeavesState[2] = { {noDataFlag, 0, 0}, {noDataFlag, 0, 0} };
            if( (*workingSize) != 0 ){
                borderLeavesState[0] = leavesInfo[0];
                borderLeavesState[1] = leavesInfo[leavesInfo.getSize()-1];
            }

            std::unique_ptr<LeafInfo[]> allProcFirstLeafStates(new LeafInfo[nbProcs*2]);
            FMpi::MpiAssert(MPI_Allgather(&borderLeavesState, sizeof(LeafInfo)*2, MPI_BYTE,
                                          allProcFirstLeafStates.get(), sizeof(LeafInfo)*2, MPI_BYTE, communicator.getComm()),__LINE__);

            FVector<MPI_Request> requests;

            // Find what to send/recv from who
            bool hasSentFirstLeaf = false;
            if( (*workingSize) != 0 ){
                // Find the owner of the leaf
                int idProcToSendTo = myRank;
                while(0 < idProcToSendTo &&
                      (allProcFirstLeafStates[(idProcToSendTo-1)*2 + 1].mindex == borderLeavesState[0].mindex
                       || allProcFirstLeafStates[(idProcToSendTo-1)*2 + 1].mindex == noDataFlag)){
                    idProcToSendTo -= 1;
                }
                // We found someone
                if(idProcToSendTo != myRank && allProcFirstLeafStates[(idProcToSendTo)*2 + 1].mindex == borderLeavesState[0].mindex){
                    // Post and send message for the first leaf
                    requests.push((MPI_Request)0);
                    FAssertLF(borderLeavesState[0].nbParts < std::numeric_limits<int>::max());
                    FMpi::MpiAssert(MPI_Isend(&workingArray[0], int(borderLeavesState[0].nbParts), MPI_BYTE, idProcToSendTo,
                            FMpi::TagExchangeIndexs, communicator.getComm(), &requests[0]),__LINE__);
                    hasSentFirstLeaf = true;
                }
            }

            bool hasExtendLastLeaf = false;
            std::vector<IndexedParticle> receivedParticles;

            {
                // Count all the particle of our first leaf on other procs
                FSize totalNbParticlesToRecv = 0;
                int idProcToRecvFrom = myRank;
                while(idProcToRecvFrom+1 < nbProcs &&
                      (borderLeavesState[1].mindex == allProcFirstLeafStates[(idProcToRecvFrom+1)*2].mindex
                       || allProcFirstLeafStates[(idProcToRecvFrom+1)*2].mindex == noDataFlag)){
                    idProcToRecvFrom += 1;
                    totalNbParticlesToRecv += allProcFirstLeafStates[(idProcToRecvFrom)*2].nbParts;
                }
                // If there are some
                if(totalNbParticlesToRecv){
                    // Alloc a received buffer
                    receivedParticles.resize(totalNbParticlesToRecv);
                    // Post the recv
                    FSize postPositionRecv = 0;
                    for(int postRecvIdx = (myRank+1); postRecvIdx <= idProcToRecvFrom ; ++postRecvIdx){
                        // If there are some on this proc
                        if(allProcFirstLeafStates[(postRecvIdx)*2].mindex != noDataFlag){
                            requests.push((MPI_Request)0);
                            FAssertLF(allProcFirstLeafStates[(postRecvIdx)*2].nbParts < std::numeric_limits<int>::max());
                            FMpi::MpiAssert(MPI_Irecv(&receivedParticles[postPositionRecv], int(allProcFirstLeafStates[(postRecvIdx)*2].nbParts), MPI_BYTE, postRecvIdx,
                                            FMpi::TagExchangeIndexs, communicator.getComm(), &requests[0]),__LINE__);
                            // Inc the write position
                            postPositionRecv += allProcFirstLeafStates[(postRecvIdx)*2].nbParts;
                        }
                    }
                    hasExtendLastLeaf = true;
                }
            }

            // Finalize communication
            FMpi::MpiAssert(MPI_Waitall(int(requests.getSize()), requests.data(), MPI_STATUSES_IGNORE),__LINE__);

            // IF we sent we need to remove the first leaf
            if(hasSentFirstLeaf){
                const FSize offsetParticles = borderLeavesState[0].nbParts;
                // Move all the particles
                for(FSize idxPart = offsetParticles ; idxPart < (*workingSize) ; ++idxPart){
                    workingArray[idxPart - offsetParticles] = workingArray[idxPart];
                }
                // Move all the leaf
                for(int idxLeaf = 1 ; idxLeaf < leavesInfo.getSize() ; ++idxLeaf){
                    leavesInfo[idxLeaf].startingPoint -= offsetParticles;
                    leavesInfo[idxLeaf - 1] = leavesInfo[idxLeaf];
                }
                (*workingSize) -= offsetParticles;
            }

            // If we received we need to merge both arrays
            if(hasExtendLastLeaf){
                // Allocate array
                const FSize finalParticlesNumber = (*workingSize) + receivedParticles.size();
                IndexedParticle* particlesWithExtension = new IndexedParticle[finalParticlesNumber];
                // Copy old data
                memcpy(particlesWithExtension, workingArray, (*workingSize)*sizeof(IndexedParticle));
                // Copy received data
                memcpy(particlesWithExtension + (*workingSize), receivedParticles.data(), receivedParticles.size()*sizeof(IndexedParticle));
                // Move ptr
                delete[] workingArray;
                workingArray   = particlesWithExtension;
                (*workingSize) = finalParticlesNumber;
            }
        }
        {//Filling the Array with leaves and parts //// COULD BE MOVED IN AN OTHER FUCTION

            (*leavesSize)    = 0; //init ptr
            (*particlesArrayInLeafOrder)   = nullptr; //init ptr
            (*leavesOffsetInParticles) = nullptr; //init ptr

            if((*workingSize)){
                //Copy all the particles
                (*particlesArrayInLeafOrder) = new ParticleClass[(*workingSize)];
                for(FSize idxPart = 0 ; idxPart < (*workingSize) ; ++idxPart){
                    memcpy(&(*particlesArrayInLeafOrder)[idxPart],&workingArray[idxPart].particle,sizeof(ParticleClass));
                }
                // Assign the number of leaf
                (*leavesSize) = leavesInfo.getSize();
                // Store the offset position for each leaf
                (*leavesOffsetInParticles) = new FSize[leavesInfo.getSize() + 1];
                for(int idxLeaf = 0 ; idxLeaf < leavesInfo.getSize() ; ++idxLeaf){
                    (*leavesOffsetInParticles)[idxLeaf] = leavesInfo[idxLeaf].startingPoint;
                }
                (*leavesOffsetInParticles)[leavesInfo.getSize()] = (*workingSize);
            }
        }
    }


    //////////////////////////////////////////////////////////////////////////
    // To equalize (same number of leaves among the procs)
    //////////////////////////////////////////////////////////////////////////

    /** Put the interval into a tree */
    template <class ContainerClass>
    static void EqualizeAndFillContainer(const FMpi::FComm& communicator,  ContainerClass* particlesSaver,
                                         const FSize leavesOffsetInParticles[], const ParticleClass particlesArrayInLeafOrder[],
                                         const FSize currentNbLeaves,
                                         const FSize currentNbParts, FAbstractBalanceAlgorithm * balancer){
        const FSize MAX_BYTE_PER_MPI_MESS = 2000000000;
        const FSize MAX_PARTICLES_PER_MPI_MESS = FMath::Max(FSize(1), FSize(MAX_BYTE_PER_MPI_MESS/sizeof(ParticleClass)));
        const int myRank = communicator.processId();
        const int nbProcs = communicator.processCount();

        if(nbProcs == 1){
            //Just copy each part into the Particle Saver
            for(FSize idxPart =0 ; idxPart < currentNbParts ; ++idxPart){
                particlesSaver->push(particlesArrayInLeafOrder[idxPart]);
            }
        }
        else{
            // We need to know the number of leaves per procs
            std::unique_ptr<FSize[]> numberOfLeavesPerProc(new FSize[nbProcs]);
            FMpi::MpiAssert(MPI_Allgather(const_cast<FSize*>(&currentNbLeaves), 1, MPI_LONG_LONG_INT, numberOfLeavesPerProc.get(),
                                          1, MPI_LONG_LONG_INT, communicator.getComm()), __LINE__);

            //prefix sum
            std::unique_ptr<FSize[]> diffNumberOfLeavesPerProc(new FSize[nbProcs+1]);
            diffNumberOfLeavesPerProc[0] = 0;
            for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc ){
                diffNumberOfLeavesPerProc[idxProc+1] = diffNumberOfLeavesPerProc[idxProc] + numberOfLeavesPerProc[idxProc];
            }

            const FSize totalNumberOfLeavesInSimulation  = diffNumberOfLeavesPerProc[nbProcs];
            // Compute the objective interval
            std::vector< std::pair<size_t,size_t> > allObjectives;
            allObjectives.resize(nbProcs);
            for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
                allObjectives[idxProc].first  = balancer->getLeft(totalNumberOfLeavesInSimulation,nbProcs,idxProc);
                allObjectives[idxProc].second = balancer->getRight(totalNumberOfLeavesInSimulation,nbProcs,idxProc);
            }

            // Ask for the pack to send
            std::pair<size_t, size_t> myCurrentInter = {diffNumberOfLeavesPerProc[myRank], diffNumberOfLeavesPerProc[myRank+1]};
            const std::vector<FEqualize::Package> packsToSend = FEqualize::GetPackToSend(myCurrentInter, allObjectives);
            std::unique_ptr<FSize[]> nbPartsPerPackToSend(new FSize[packsToSend.size()]);
            // Store the requests
            std::vector<MPI_Request> requestsParts;
            std::vector<MPI_Request> requestsNbParts;
            // Send every thing except for me or if size == 0
            for(unsigned int idxPack = 0; idxPack< packsToSend.size() ; ++idxPack){
                const FEqualize::Package& pack = packsToSend[idxPack];
                if(pack.idProc != myRank && 0 < (pack.elementTo-pack.elementFrom)){
                    // If not to me and if there is something to send
                    nbPartsPerPackToSend[idxPack] = leavesOffsetInParticles[pack.elementTo]-leavesOffsetInParticles[pack.elementFrom];
                    // Send the size of the data
                    requestsNbParts.emplace_back();
                    FMpi::MpiAssert(MPI_Isend(&nbPartsPerPackToSend[idxPack],1,MPI_LONG_LONG_INT,pack.idProc,
                                              FMpi::TagExchangeIndexs+1, communicator.getComm(), &requestsNbParts.back()),__LINE__);
                    // Send the data
                    for(FSize idxMess = 0 ; idxMess < nbPartsPerPackToSend[idxPack]; idxMess += MAX_PARTICLES_PER_MPI_MESS){
                        const int nbElementsInMessage = int(FMath::Min(nbPartsPerPackToSend[idxPack]-idxMess, MAX_PARTICLES_PER_MPI_MESS));
                        requestsParts.emplace_back();
                        FMpi::MpiAssert(MPI_Isend(const_cast<ParticleClass*>(&particlesArrayInLeafOrder[leavesOffsetInParticles[pack.elementFrom]+idxMess]),
                                int(sizeof(ParticleClass)*nbElementsInMessage),
                                MPI_BYTE, pack.idProc, int(FMpi::TagExchangeIndexs + 2 + idxMess), communicator.getComm(), &requestsParts.back()), __LINE__);
                    }
                }
                else {
                    // Nothing to send
                    nbPartsPerPackToSend[idxPack] = 0;
                }
            }
            // Compute the current intervals
            std::vector< std::pair<size_t,size_t> > allCurrentIntervals;
            allCurrentIntervals.resize(nbProcs);
            for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
                allCurrentIntervals[idxProc].first  = diffNumberOfLeavesPerProc[idxProc];
                allCurrentIntervals[idxProc].second = diffNumberOfLeavesPerProc[idxProc+1];
            }
            // Ask the packs to receive to fill my objective
            std::pair<size_t, size_t> myObjective = allObjectives[myRank];
            const std::vector<FEqualize::Package> packsToRecv = FEqualize::GetPackToRecv(myObjective, allCurrentIntervals);

            // Count the number of parts to receive
            std::unique_ptr<FSize[]> nbPartsPerPackToRecv(new FSize[packsToRecv.size()]);
            for(unsigned int idxPack = 0; idxPack < packsToRecv.size(); ++idxPack){
                const FEqualize::Package& pack = packsToRecv[idxPack];
                if(pack.idProc != myRank && 0 < (pack.elementTo-pack.elementFrom)){
                    // We need to know how much particles to receive
                    requestsNbParts.emplace_back();
                    FMpi::MpiAssert(MPI_Irecv(&nbPartsPerPackToRecv[idxPack], 1, MPI_LONG_LONG_INT, pack.idProc,
                                              FMpi::TagExchangeIndexs+1, communicator.getComm(), &requestsNbParts.back()), __LINE__);
                }
                else{
                    if(pack.idProc == myRank){
                        // Take my own data
                        const FSize sourcePosition = FMath::Max(myObjective.first, myCurrentInter.first) - myCurrentInter.first;
                        const FSize nbLeavesToCopy = pack.elementTo-pack.elementFrom;
                        nbPartsPerPackToRecv[idxPack] = leavesOffsetInParticles[sourcePosition+nbLeavesToCopy] - leavesOffsetInParticles[sourcePosition];
                    }
                    else{
                        // Nothing to receive from this so avoid communication
                        nbPartsPerPackToRecv[idxPack] = 0;
                    }
                }
            }

            FMpi::MpiAssert(MPI_Waitall(int(requestsNbParts.size()), requestsNbParts.data(), MPI_STATUSES_IGNORE), __LINE__);

            // Count the number of leaf to receive
            FSize totalPartsToReceive = 0;
            for(unsigned int idxPack = 0; idxPack < packsToRecv.size(); ++idxPack){
                totalPartsToReceive += nbPartsPerPackToRecv[idxPack];
            }

            std::vector<ParticleClass> particlesRecvBuffer;
            // Post all the receive and copy mine
            if(totalPartsToReceive){
                particlesRecvBuffer.resize(totalPartsToReceive);
                FSize offsetToRecv = 0;
                for(unsigned int idxPack = 0; idxPack < packsToRecv.size(); ++idxPack){
                    const FEqualize::Package& pack = packsToRecv[idxPack];
                    if(pack.idProc != myRank && 0 < (pack.elementTo-pack.elementFrom)){
                        for(FSize idxMess = 0 ; idxMess < nbPartsPerPackToRecv[idxPack]; idxMess += MAX_PARTICLES_PER_MPI_MESS){
                            const int nbElementsInMessage = int(FMath::Min(nbPartsPerPackToRecv[idxPack]-idxMess, MAX_PARTICLES_PER_MPI_MESS));
                            requestsParts.emplace_back();
                            FMpi::MpiAssert( MPI_Irecv(&particlesRecvBuffer[offsetToRecv+idxMess],
                                             int(sizeof(ParticleClass)*nbElementsInMessage), MPI_BYTE, pack.idProc,
                                             int(FMpi::TagExchangeIndexs + 2 + idxMess), communicator.getComm(), &requestsParts.back()), __LINE__);
                        }
                    }
                    else if(pack.idProc == myRank){
                        // Copy my particles
                        const FSize sourcePosition = FMath::Max(myObjective.first, myCurrentInter.first) - myCurrentInter.first;
                        memcpy(&particlesRecvBuffer[offsetToRecv], &particlesArrayInLeafOrder[leavesOffsetInParticles[sourcePosition]],
                                nbPartsPerPackToRecv[idxPack]*sizeof(ParticleClass));
                    }
                    offsetToRecv += nbPartsPerPackToRecv[idxPack];
                }
            }

            // Finalize communication
            FMpi::MpiAssert(MPI_Waitall(int(requestsParts.size()), requestsParts.data(), MPI_STATUSES_IGNORE), __LINE__);

            // Insert in the particle saver
            for(FSize idPartsToStore = 0 ; idPartsToStore < int(particlesRecvBuffer.size()) ; ++idPartsToStore){
                particlesSaver->push(particlesRecvBuffer[idPartsToStore]);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // The builder function
    //////////////////////////////////////////////////////////////////////////

    template <class ContainerClass>
    static void DistributeArrayToContainer(const FMpi::FComm& communicator, const ParticleClass originalParticlesArray[], const FSize originalNbParticles,
                                           const FPoint<FReal>& boxCenter, const FReal boxWidth, const int treeHeight,
                                           ContainerClass* particleSaver, FAbstractBalanceAlgorithm* balancer, const SortingType sortingType = QuickSort){

        FLOG( FLog::Controller << "Enter DistributeArrayToContainer\n" );
        FLOG( FTic timer );

        IndexedParticle* sortedParticlesArray = nullptr;
        FSize nbParticlesInArray = 0;
        // From ParticleClass get array of IndexedParticle sorted
        GetSortedParticlesFromArray(communicator, originalParticlesArray, originalNbParticles, sortingType, boxCenter, boxWidth, treeHeight,
                                    &sortedParticlesArray, &nbParticlesInArray);
        FLOG( FLog::Controller << "\t GetSortedParticlesFromArray is over (" << timer.tacAndElapsed() << "s)\n" );
        FLOG( timer.tic() );

        ParticleClass* particlesArrayInLeafOrder = nullptr;
        FSize * leavesOffsetInParticles = nullptr;
        FSize nbLeaves = 0;
        // Merge the leaves
        MergeSplitedLeaves(communicator, sortedParticlesArray, &nbParticlesInArray, &leavesOffsetInParticles, &particlesArrayInLeafOrder, &nbLeaves);
        delete[] sortedParticlesArray;

        FLOG( FLog::Controller << "\t MergeSplitedLeaves is over (" << timer.tacAndElapsed() << "s)\n" );
        FLOG( timer.tic() );

        // Equalize and balance
        EqualizeAndFillContainer(communicator, particleSaver, leavesOffsetInParticles, particlesArrayInLeafOrder, nbLeaves,
                                 nbParticlesInArray, balancer);
        delete[] particlesArrayInLeafOrder;
        delete[] leavesOffsetInParticles;

        FLOG( FLog::Controller << "\t EqualizeAndFillContainer is over (" << timer.tacAndElapsed() << "s)\n" );

        FLOG( FLog::Controller << "\t DistributeArrayToContainer is over (" << timer.cumulated() << "s)\n" );

#ifdef SCALFMM_USE_LOG
        /** To produce stats after the Equalize phase  */
        {
            const FSize finalNbParticles = particleSaver->getSize();

            if(communicator.processId() != 0){
                FMpi::MpiAssert(MPI_Gather(const_cast<FSize*>(&finalNbParticles),1,FMpi::GetType(finalNbParticles),nullptr,
                                           1,FMpi::GetType(finalNbParticles),0,communicator.getComm()), __LINE__);
            }
            else{
                const int nbProcs = communicator.processCount();
                std::unique_ptr<FSize[]> nbPartsPerProc(new FSize[nbProcs]);

                FMpi::MpiAssert(MPI_Gather(const_cast<FSize*>(&finalNbParticles),1,FMpi::GetType(finalNbParticles),nbPartsPerProc.get(),
                                           1,FMpi::GetType(finalNbParticles),0,communicator.getComm()), __LINE__);

                FReal averageNbParticles = 0;
                FSize minNbParticles = finalNbParticles;
                FSize maxNbParticles = finalNbParticles;

                for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
                    maxNbParticles = FMath::Max(maxNbParticles, nbPartsPerProc[idxProc]);
                    minNbParticles = FMath::Min(minNbParticles, nbPartsPerProc[idxProc]);
                    averageNbParticles += FReal(nbPartsPerProc[idxProc]);
                }
                averageNbParticles /= float(nbProcs);

                printf("End of Equalize Phase : \n \t Min number of parts : %lld \n \t Max number of parts : %lld \n \t Average number of parts : %e \n",
                       minNbParticles,maxNbParticles,averageNbParticles);
            }
        }
#endif
    }


};

#endif // FMPITREEBUILDER_H
