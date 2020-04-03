#ifndef FUSERKERNELDISTRENGINE_HPP
#define FUSERKERNELDISTRENGINE_HPP

#include "vector"
#include <algorithm>
/**
 * @file FUserKernelDistrEngine.hpp
 * @brief This file add MPI features to FUserKernelEngine, by reimplement or implement new methods.
 */


#include "Utils/FMpi.hpp"
#include "FUserKernelEngine.hpp"
#include "BalanceTree/FLeafBalance.hpp"
#include "Files/FMpiTreeBuilder.hpp"
#include "Containers/FVector.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp"

class CoreCellDist : public CoreCell, public FAbstractSendable{
    int level;

public:
    CoreCellDist() : level(-1){

    }

    /**
     * To save the current cell into a buffer in order to send it.
     */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FBasicCell::save<BufferWriterClass>();
        buffer << level;
        if(user_cell_descriptor.user_get_size){
            FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
            char * temp = new char[sizeToSave];
            user_cell_descriptor.user_copy_cell(CoreCell::getContainer(),sizeToSave,(void *) temp);
            buffer.write(temp,sizeToSave);
            delete[] temp;
        }
    }

    /**
     * To "create" a cell from an array in order to receive cells.
     */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FBasicCell::restore<BufferReaderClass>();
        buffer >> level;
        if(user_cell_descriptor.user_restore_cell){
            FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
            char * temp = new char[sizeToSave];
            buffer.fillArray(temp,sizeToSave);
            CoreCell::setContainer(user_cell_descriptor.user_restore_cell(level,temp));
            delete[] temp;
        }
    }

    /**
     * @brief Part where we reimplement {de,}serialize{Up,Down}
     */
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const {
        FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
        char * temp = new char[sizeToSave];
        user_cell_descriptor.user_copy_cell(CoreCell::getContainer(),sizeToSave,(void *) temp);
        buffer.write(temp,sizeToSave);
        delete[] temp;
    }
    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const {
        FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
        char * temp = new char[sizeToSave];
        user_cell_descriptor.user_copy_cell(CoreCell::getContainer(),sizeToSave,(void *) temp);
        buffer.write(temp,sizeToSave);
        delete[] temp;
    }
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer) const {
        FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
        char * temp = new char[sizeToSave];
        buffer.fillArray(temp,sizeToSave);
        CoreCell::setContainer(user_cell_descriptor.user_restore_cell(level,temp));

        delete[] temp;
    }
    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer) const {
        FSize sizeToSave = user_cell_descriptor.user_get_size(level,CoreCell::getContainer(),getMortonIndex());
        char * temp = new char[sizeToSave];
        buffer.fillArray(temp,sizeToSave);
        CoreCell::setContainer(user_cell_descriptor.user_restore_cell(level,temp));
        delete[] temp;
    }

    /**
     * @brief Function to get the size of a cell : (size of user datas
     * + size of internal data)
     */
    FSize getSavedSize() const {
        FSize toReturn = user_cell_descriptor.user_get_size(level,
                                                            CoreCell::getContainer(),
                                                            getMortonIndex())
            + FBasicCell::getSavedSize() //Size of storage needed for Basic Cell
            + sizeof(int)                //Size of storage needed for this class
            + sizeof(nullptr);           //Size of storage needed for the parent class
        return toReturn;
    }
    FSize getSavedSizeUp() const{
        FSize toReturn = user_cell_descriptor.user_get_size(level,
                                                            CoreCell::getContainer(),
                                                            getMortonIndex());
        return toReturn;
    }
    FSize getSavedSizeDown() const{
        FSize toReturn = user_cell_descriptor.user_get_size(level,
                                                            CoreCell::getContainer(),
                                                            getMortonIndex());
        return toReturn;
    }
};


template<class FReal, class LeafClass>
class FUserKernelDistrEngine: public FUserKernelEngine<FReal,LeafClass>{

private:
    using Parent = FUserKernelEngine<FReal,LeafClass>;

    //Same as in Parent class.
    using ContainerClass = FUserLeafContainer<FReal>;
    using OctreeClass = FOctree<FReal,CoreCellDist,ContainerClass,LeafClass>;
    using CoreKernelClass =  CoreKernel<CoreCellDist,ContainerClass>;

    /**
     * @brief Class for convert global indices to local indices.
     */
    class Indirection{
    private:
        //Keep track of number of local points
        FSize nbPoints;
        //Keep track of global indexes
        std::vector<FSize> indirection;
        //Number of part per proc before partitionning
        std::vector<FSize> oriPartPerProc;
        //Number of part per proc after Partitionning
        std::vector<FSize> partPerProc;
        //Intervals of part owned by each proc before partitionning
        std::vector< std::pair<FSize,FSize> > oriIntervalPerProc;
        std::vector< std::vector<FSize> > * mapComm;

        //In order to send again following the first partition, we
        //store everything to do the Alltoall
        std::vector<FSize> orderToSortToSend;
        //This will hold the indices of original parts in order to be
        //sent directly
        std::vector<FSize> initialArray;
        std::vector<FSize> currentArray;
        std::vector<FSize> currentArraySorted;
        std::vector<int> howManyToRecv;
        std::vector<int> howManyToSend;
        std::vector<int> displRecv;
        std::vector<int> displSend;

    public:
        Indirection() : indirection(0), oriPartPerProc(0), partPerProc(0), oriIntervalPerProc(0), mapComm(nullptr), initialArray(0){
        }

        ~Indirection(){
            nbPoints=0;
            if(mapComm){
                delete mapComm;
            }
        }

        /**
         *@brief Return the global (original) index
         */
        FSize getGlobalInd(FSize localInd) const {
            return indirection[localInd];
        }

        /**
         *@brief Return the local (to this process) index
         */
        FSize getLocalInd(FSize globalInd) const {
            for(auto ite = indirection.cbegin(); ite != indirection.cend() ; ++ite){
                if(*ite == globalInd){
                    return (ite - indirection.cbegin());
                }
            }
            //will cause a segmentation fault
            return -1;
        }
        /**
         * Get the current number of particles
         */
        FSize getLocalSize() const {
            return ((FSize) indirection.size());
        }

        /**
         * Add a global index
         */
        void addPoint(FSize globalInd){
            indirection.push_back(globalInd);
            ++nbPoints;
        }

        //Add multiple global indices
        void addArrayPoints(FSize inNbPoint, FSize * globalIndices){
            for(int ite = 0; ite < inNbPoint ; ++ite){
                indirection().push_back(globalIndices[ite]);
            }
            nbPoints += inNbPoint;
        }

        void setNbProc(int inNbProc){
            this->oriPartPerProc.resize(inNbProc);
            this->partPerProc.resize(inNbProc);
            this->oriIntervalPerProc.resize(inNbProc);
            this->mapComm = new std::vector<std::vector<FSize> > (inNbProc,std::vector<FSize>(inNbProc));
        }

        void setOriIntervals(const FMpi::FComm * comm){
            FSize acc = 0;
            for(int i=0; i<comm->processCount() ; ++i){
                oriIntervalPerProc[i] = std::make_pair(acc,acc+oriPartPerProc[i]-1);
                //std::cout << "Rank "<< comm->processId() << "(acc= "<< acc <<",acc+oriPartPerproc["<<i<<"]-1= " << acc+oriPartPerProc[i]-1 << ")\n";
                acc+=oriPartPerProc[i];
            }
        }

        void setNbPartPerProc(FSize nbPart, const FMpi::FComm * comm){
            if(partPerProc.size() < comm->processCount()){
                setNbProc(comm->processCount());
            }
            MPI_Allgather(&nbPart, 1, MPI_LONG_LONG, partPerProc.data() , 1, MPI_LONG_LONG, comm->getComm());
        }
        void setOriNbPartPerProc(FSize nbPart, const FMpi::FComm * comm){
            if(oriPartPerProc.size() < comm->processCount()){
                setNbProc(comm->processCount());
            }
            MPI_Allgather(&nbPart, 1, MPI_LONG_LONG, oriPartPerProc.data(), 1, MPI_LONG_LONG, comm->getComm());
        }

        template<class ParticleSaver>
        void setMap(const FMpi::FComm * comm,ParticleSaver& particles){
            FSize * temp_map = new FSize[comm->processCount()*comm->processCount()];
            memset(temp_map,0,sizeof(FSize)*comm->processCount()*comm->processCount());
            //Each process write from temp[my_rank*nb_proc] to
            //temp[(my_rank+1)*nb_proc - 1]
            for(int idPart=0 ; idPart<getCount(comm->processId()) ; ++idPart){
                temp_map[comm->processId()*comm->processCount() + particles[idPart].orOwner]++;
            }
            //Gather the Map
            MPI_Allgather(&temp_map[comm->processId()*comm->processCount()], comm->processCount(),
                          MPI_LONG_LONG, temp_map, comm->processCount(), MPI_LONG_LONG, comm->getComm());
            for(int idProc=0; idProc<comm->processCount() ; ++idProc){
                for(int idProcSource=0; idProcSource<comm->processCount() ; ++idProcSource){
                    (*mapComm)[idProc][idProcSource] = temp_map[idProc*comm->processCount() + idProcSource];
                }
            }
            //For outputing the map
            // for(int i=0; i<comm->processCount(); i++){
            //     for(int j=0; j<comm->processCount(); ++j){
            //         std::cout <<comm->processId()<<" LALALA " << temp_map[i*comm->processCount()+j] <<" ";
            //     }
            //     std::cout << std::endl;
            // }
            delete [] temp_map;

            buildIndirectionFromMap<ParticleSaver>(comm,particles);
        }

        template<class ParticleSaver>
        void buildIndirectionFromMap(const FMpi::FComm * comm, ParticleSaver& particles){
            //MPI_Allgatherv in order to get the indices of the parts
            //sent
            int rank = comm->processId();

            //initialise original array (to receive original part)
            initialArray.resize(oriIntervalPerProc[rank].second - oriIntervalPerProc[rank].first + 1);
            //initialise current array (to send to others)
            currentArraySorted.resize(partPerProc[rank]);
            currentArray.resize(partPerProc[rank]);

            for(int idPart=0 ;idPart < partPerProc[rank] ; ++idPart){
                currentArraySorted[idPart] = particles[idPart].orIndex;
                currentArray[idPart] = particles[idPart].orIndex;
            }

            //Sort the current array by Indice (so, there will be
            //different tile each one corresponding to one owner)
            std::sort(currentArraySorted.begin(),currentArraySorted.end());

            //Compute the array of number of indices to recv
            for(int i=0; i< comm->processCount() ; ++i){
                howManyToRecv.push_back(static_cast<int>((*mapComm)[i][rank]));
            }

            //Compute the array of number of indices to send
            howManyToSend.resize(comm->processCount());
            std::transform((*mapComm)[rank].begin(),(*mapComm)[rank].end(), howManyToSend.begin(),
                           [](const FSize & A) -> int{
                               return static_cast<int>(A);
                           });

            //Compute the array of displacement for receiving
            int acc=0;
            displRecv.resize(comm->processCount());
            std::transform(howManyToRecv.begin(),howManyToRecv.end(),displRecv.begin(),
                           [&](const FSize & A) -> int {
                               auto temp = acc;
                               acc += static_cast<int>(A);
                               return temp;
                           });
            //Same for sending
            acc=0;
            displSend.resize(comm->processCount());
            std::transform(howManyToSend.begin(),howManyToSend.end(),displSend.begin(),
                           [&](const FSize & A){
                               auto temp = acc;
                               acc += static_cast<int>(A);
                               return temp;
                           });
            //Process i will send to j its j_th block of indices
            //Then each proc will know which part he used to own
            MPI_Alltoallv(currentArraySorted.data(),howManyToSend.data(),displSend.data(),MPI_LONG_LONG,
                          initialArray.data(),howManyToRecv.data(),displRecv.data(),MPI_LONG_LONG,comm->getComm());

        }

        std::pair<FSize,FSize>& getoriIntervals(int inRank) {
            return oriIntervalPerProc[inRank];
        }

        FSize getOriCount(int inRank){
            return oriPartPerProc[inRank];
        }
        FSize getCount(int inRank){
            return partPerProc[inRank];
        }
        template<typename T = int>
        FSize getInitialArrayValue(T inInd){
            return initialArray[inInd];
        }
        template<typename T = int>
        FSize getCurrentArrayValue(T inInd){
            return currentArray[inInd];
        }
        template<typename T = int>
        FSize getCurrentArraySortedValue(T inInd){
            return currentArraySorted[inInd];
        }
        std::vector<int>& getDisplSendVector(){
            return displSend;
        }
        std::vector<int>& getDisplRecvVector(){
            return displRecv;
        }
        std::vector<int>& getHowManyToSend(){
            return howManyToSend;
        }
        std::vector<int>& getHowManyToRecv(){
            return howManyToRecv;
        }
    };

    struct PartToSort{
        FPoint<FReal> position;
        MortonIndex index;
        FSize orIndex;
        int orOwner;
        FPoint<FReal> getPosition(){
            return position;
        }

        //Copy operator for memoryLocation
        PartToSort& operator=(PartToSort const & inPart){
            this->position = inPart.position;
            this->index = inPart.index;
            this->orIndex = inPart.orIndex;
            this->orOwner = inPart.orOwner;
            return *this;
        }
    };

    void checkTree() const{
        if(octreeDist == nullptr){
            FAssertLF(0,"Tree need to be built first \n");
        }
    }

    void checkAlready() const {
        if(! alreadyPartionned){
            FAssertLF(0,"Particles need to be partitionned first \n");
        }
    };

    void checkPartType(PartType type) const {
        if(type != BOTH){
            FAssertLF(0,"No source/Target Algorithm available in distributed version \n");
        }
    };


    //Communicator
    const FMpi::FComm * comm;
    //Check if Equalize phase is already done
    bool alreadyPartionned;
    //wrapper of Array Indirection
    Indirection Ind;
    OctreeClass * octreeDist;
    CoreKernelClass* kernel;

public:

    FUserKernelDistrEngine(scalfmm_kernel_type KernelType, scalfmm_algorithm algo, const MPI_Comm inComm)
        : Ind(), octreeDist(nullptr), kernel(nullptr){
        FScalFMMEngine<FReal>::kernelType = KernelType;
        FScalFMMEngine<FReal>::Algorithm = algo;
        this->comm = new FMpi::FComm(inComm);
        Ind.setNbProc(comm->processCount());
    }

    ~FUserKernelDistrEngine(){
        delete comm;
        comm = nullptr;
    }

    //Qu'est-ce qu'il faut que je surcharge ?
    void tree_insert_particles( int NbPositions, double * X, double * Y, double * Z, PartType type){
        checkPartType(type);
        checkAlready();
        for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
            octreeDist->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),idPart);
        }
        FScalFMMEngine<FReal>::nbPart += NbPositions;
        this->init_cell();
    }

    void build_tree(int TreeHeight,double BoxWidth,double* BoxCenter,Scalfmm_Cell_Descriptor user_cell_descriptor){
        CoreCell::Init(user_cell_descriptor);
        Parent::treeHeight = TreeHeight;
        Parent::boxCenter = FPoint<FReal>(BoxCenter[0],BoxCenter[1],BoxCenter[2]);
        Parent::boxWidth = BoxWidth;
        Parent::boxCorner.setX(BoxCenter[0] - BoxWidth/2.0);
        Parent::boxCorner.setY(BoxCenter[1] - BoxWidth/2.0);
        Parent::boxCorner.setZ(BoxCenter[2] - BoxWidth/2.0);
        Parent::boxWidthAtLeafLevel = BoxWidth/(2<<TreeHeight);
        printf("Tree Height : %d \n",TreeHeight);
        octreeDist = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint<FReal>(BoxCenter));
    }

    void user_kernel_config( Scalfmm_Kernel_Descriptor userKernel, void * userDatas){
        if(!kernel){
            kernel = new CoreKernelClass(userKernel,userDatas);
        }
    }
    /**
     * The attributes arg will be partitionned, too.
     */
    void create_local_partition(int nbPoints, double * particleXYZ, double ** localArrayFilled, FSize ** indexesFilled, FSize * outputNbPoint){
        checkTree();
        // This is the array that will be sorted
        PartToSort * arrayToBeSorted = new PartToSort[nbPoints];
        int myRank = comm->processId();

        //Store in indirection how many parts per proc
        Ind.setOriNbPartPerProc(nbPoints,comm);
        Ind.setOriIntervals(comm);

        //        std::cout << myRank << " 0:: " << Ind.getoriIntervals(0).first  << " " << Ind.getoriIntervals(0).second  << " 1:: " << Ind.getoriIntervals(1).first  << " " << Ind.getoriIntervals(2).second << std::endl;

        for(int id = 0 ; id<nbPoints ; ++id){
            FTreeCoordinate host;
            arrayToBeSorted[id].position = FPoint<FReal>(particleXYZ[id+0],particleXYZ[id+1],particleXYZ[id+2]);
            arrayToBeSorted[id].orIndex = id + Ind.getoriIntervals(myRank).first;
            arrayToBeSorted[id].orOwner = myRank;
            //Evaluate Morton Index
            host.setX(FCoordinateComputer::GetTreeCoordinate<FReal>(particleXYZ[id*3+0] - this->getBoxCorner().getX(),this->getBoxWidth(),this->getBoxWidthAtLeafLevel(),
                                                                    this->getTreeHeight()));
            host.setX(FCoordinateComputer::GetTreeCoordinate<FReal>(particleXYZ[id*3+1] - this->getBoxCorner().getX(),this->getBoxWidth(),this->getBoxWidthAtLeafLevel(),
                                                                    this->getTreeHeight()));
            host.setX(FCoordinateComputer::GetTreeCoordinate<FReal>(particleXYZ[id*3+2] - this->getBoxCorner().getX(),this->getBoxWidth(),this->getBoxWidthAtLeafLevel(),
                                                                    this->getTreeHeight()));
            arrayToBeSorted[id].index = host.getMortonIndex(this->getTreeHeight());
        }
        // This is the array that will contains our particles.
        FVector<struct PartToSort> particles;
        //Balancer
        FLeafBalance balancer;

        //Call to partition sorting algorithm
        FMpiTreeBuilder<FReal,PartToSort>::DistributeArrayToContainer(*comm,arrayToBeSorted,nbPoints,
                                                                      this->getBoxCenter(),this->getBoxWidth(),this->getTreeHeight(),
                                                                      &particles,&balancer);
        printf("Finish ! output nb : %lld, rank : %d\n",particles.getSize(),comm->processId());

        delete[] arrayToBeSorted;
        printf("%d First Part Container deleted\n",comm->processId());

        //Now I have inside "particles" what part I shoud deal with.
        FSize newNumberOfPoint = particles.getSize();
        //Then, copy back inside container
        *localArrayFilled = new double[3*newNumberOfPoint];
        *indexesFilled = new FSize[newNumberOfPoint];

        Ind.setNbPartPerProc(newNumberOfPoint,comm);

        for(int id = 0 ; id<newNumberOfPoint ; ++id){
            //Store inside user's struct (C-like)
            (*localArrayFilled)[3*id+0] = particles[id].getPosition().getX();
            (*localArrayFilled)[3*id+1] = particles[id].getPosition().getY();
            (*localArrayFilled)[3*id+2] = particles[id].getPosition().getZ();
            (*indexesFilled)[id] = particles[id].orIndex;

            //Keep what we've done through Indirection's class
            Ind.addPoint(particles[id].orIndex);
        }
        *outputNbPoint = newNumberOfPoint;

        //Before erasing particles, give it to Indirection in order to
        //build the map and prepare the generic_partition method.
        Ind.setMap(comm,particles);
        //Delete everything
        particles.clear();
        //Set the bool to true
        this->alreadyPartionned = true;

    }

    void generic_partition(FSize nbThings, size_t stride, void * arrayOfThing, void ** newArray){
        checkAlready();
        int myRank = comm->processId();
        //assign args to a vector
        std::vector<char > sendBuffer ;
        //change type of array
        char * arrayToBePartitionned = reinterpret_cast<char *>(arrayOfThing);
        sendBuffer.resize(nbThings*stride);
        long unsigned int cursor = 0;
        //I put the array in an order where I can send everything
        for(int i=0; i<nbThings ; ++i){
            //The place I want to store the value is in the initial Array, minus the minimum indice.
            FSize realIndice = Ind.getInitialArrayValue(i) - Ind.getoriIntervals(myRank).first;
            memcpy(&sendBuffer.data()[cursor],&arrayToBePartitionned[realIndice*stride],stride);
            cursor += stride;
        }

        //Then i use again displ and howmany arrays. (need to transform that)
        std::vector<int> howManyToRecvByte;
        std::vector<int> howManyToSendByte;
        std::vector<int> displRecvByte;
        std::vector<int> displSendByte;

        howManyToRecvByte.resize(comm->processCount());
        howManyToSendByte.resize(comm->processCount());
        displRecvByte.resize(comm->processCount());
        displSendByte.resize(comm->processCount());
        //every one is multiple by stride
        std::transform(Ind.getHowManyToRecv().begin(),Ind.getHowManyToRecv().end(),howManyToRecvByte.begin(),
                       [&](const FSize & A){
                           return A*stride;
                       });
        std::transform(Ind.getHowManyToSend().begin(),Ind.getHowManyToSend().end(),howManyToSendByte.begin(),
                       [&](const FSize & A){
                           return A*stride;
                       });
        std::transform(Ind.getDisplSendVector().begin(),Ind.getDisplSendVector().end(),displSendByte.begin(),
                        [&](const FSize & A){
                           return A*stride;
                       });
        std::transform(Ind.getDisplRecvVector().begin(),Ind.getDisplRecvVector().end(),displRecvByte.begin(),
                        [&](const FSize & A){
                           return A*stride;
                       });
        std::vector<char > recvBuffer;
        recvBuffer.resize(Ind.getCount(myRank)*stride);

        MPI_Alltoallv(sendBuffer.data(),howManyToRecvByte.data(),displRecvByte.data(),MPI_BYTE,
                      recvBuffer.data(),howManyToSendByte.data(),displSendByte.data(),MPI_BYTE,comm->getComm());
        //Then, i need to find the order.
        std::vector<char> output;
        output.resize(recvBuffer.size());

        //A double loop will do ...

        //Loop over the recvBuffer (following its order)
        for(int idRead = 0; idRead < Ind.getCount(myRank); ++idRead){
            //get what shall be there
            FSize currentIndex = Ind.getCurrentArraySortedValue(idRead);
            //Then, i need to find the index where current index
            //is sorted inside currentArray;
            FSize idSearch = 0;
            for(idSearch=0 ; idSearch < Ind.getCount(myRank) && (Ind.getCurrentArrayValue(idSearch) != currentIndex) ; ++idSearch);
            memcpy(&output[idSearch*stride],&recvBuffer[idRead*stride],stride);
        }

        //C Part
        *newArray = malloc(output.size());
        memcpy(*newArray,output.data(),output.size());
    }

    void tree_insert_particles_xyz( int NbPositions, double * XYZ, PartType type){
        checkPartType(type);
        checkAlready();
        for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
            octreeDist->insert(FPoint<FReal>(&XYZ[3*idPart]),idPart);
        }
        FScalFMMEngine<FReal>::nbPart += NbPositions;
        this->init_cell();
    }

    void init_cell(){
        void * generic_ptr = nullptr;
        if(kernel){
            generic_ptr = kernel->getUserKernelDatas();
        }
        else{
            std::cout <<"Warning, no user kernel data set, need to call kernel config first"<< std::endl;
        }
        double boxwidth = octreeDist->getBoxWidth();
        //apply user function on each cell
        octreeDist->forEachCellWithLevel([&](CoreCellDist * currCell,const int currLevel){
                if(!(currCell->getContainer())){
                    FTreeCoordinate currCoord = currCell->getCoordinate();
                    int arrayCoord[3] = {currCoord.getX(),currCoord.getY(),currCoord.getZ()};
                    MortonIndex    currMorton = currCoord.getMortonIndex(currLevel);
                    double position[3];
                    FPoint<FReal> boxC = this->getBoxCorner();
                    position[0] = boxC.getX() + currCoord.getX()*boxwidth/double(1<<currLevel);
                    position[1] = boxC.getY() + currCoord.getY()*boxwidth/double(1<<currLevel);
                    position[2] = boxC.getZ() + currCoord.getZ()*boxwidth/double(1<<currLevel);
                    currCell->setContainer(CoreCell::GetInit()(currLevel,currMorton,arrayCoord,position,generic_ptr));
                }
            });
    }

    void apply_on_each_leaf(Callback_apply_on_leaf function){
        if(octreeDist){
            FUserKernelEngine<FReal,LeafClass>::template generic_apply_on_each_leaf<ContainerClass,CoreCellDist>(octreeDist,kernel->getUserKernelDatas(),function);
        }else{
            std::cout << "Need to Build the tree and insert the parts First\n" << std::endl;
        }
    }


    void execute_fmm(){
        FAssertLF(octreeDist,
                  "No Tree set, please use scalfmm_user_kernel_config before calling the execute routine ... Exiting \n");
        //Only one config shall work , so let's use it
        switch(FScalFMMEngine<FReal>::Algorithm){
        case 5:
            {
                typedef FFmmAlgorithmThreadProc<OctreeClass,CoreCellDist,ContainerClass,CoreKernelClass,LeafClass> AlgoProcClass;
                AlgoProcClass * algoProc = new AlgoProcClass(*comm,octreeDist,kernel);
                FScalFMMEngine<FReal>::algoTimer = algoProc;
                algoProc->execute();
                break;
            }
        default:
            break;
        }
    }

    void free_cell(Callback_free_cell user_cell_deallocator){
        octreeDist->forEachCell([&](CoreCellDist * currCell){
                if(currCell->getContainer()){
                    user_cell_deallocator(currCell->getContainer());
                    currCell->setContainer(nullptr);
                }
            });
    }


    void intern_dealloc_handle(Callback_free_cell userDeallocator){
        free_cell(userDeallocator);
    }

};

#endif //FUSERKERNELDISTRENGINE_HPP
