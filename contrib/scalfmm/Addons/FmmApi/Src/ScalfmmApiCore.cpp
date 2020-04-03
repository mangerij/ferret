// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"

#ifdef SCALFMM_USE_MPI
#include "../../Src/Utils/FMpi.hpp"
#endif

#include "FmmApi.h"

template <class ContainerClass>
class CoreCell : public FBasicCell {
    char* multipole;
    char* local;
    int level;
    FReal positions[3];
    mutable ContainerClass* container;

public:
    CoreCell() : multipole(nullptr), local(nullptr), level(0), container(nullptr) {
    }
    void createArrays(const int multipoleSize, const int localSize){
        multipole = new char[multipoleSize];
        local     = new char[localSize];

        memset(multipole, 0, multipoleSize);
        memset(local, 0, localSize);
    }
    ~CoreCell(){
        delete[] multipole;
        delete[] local;
    }

    const void* getMultipole() const{
        return multipole;
    }

    const void* getLocal() const{
        return local;
    }

    void* getMultipole(){
        return multipole;
    }

    void* getLocal(){
        return local;
    }

    void setLevel(const int inLevel){
        level = inLevel;
    }

    int getLevel() const {
        return level;
    }

    void setPosition(const FReal inPositions[3]){
        positions[0] = inPositions[0];
        positions[1] = inPositions[1];
        positions[2] = inPositions[2];
    }

    const FReal* getPosition() const {
        return positions;
    }

    void setContainer(ContainerClass* inContainer) const {
        container = inContainer;
    }

    ContainerClass* getContainer() const {
        return container;
    }
};



template< class CellClass, class ContainerClass>
class CoreKernel : public FAbstractKernels<CellClass,ContainerClass> {
    void* fmmCore;

public:
    CoreKernel(void* inFmmCore): fmmCore(inFmmCore){
    }

    /** Default destructor */
    virtual ~CoreKernel(){
    }

    /** Do nothing */
    virtual void P2M(CellClass* const cell, const ContainerClass* const container) {
        cell->setContainer(const_cast<ContainerClass*>(container));
        FmmKernel_P2M(fmmCore,(void*) cell);
    }

    /** Do nothing */
    virtual void M2M(CellClass* const FRestrict cell, const CellClass*const FRestrict *const FRestrict children, const int ) {
        for(int idx = 0 ; idx < 8 ; ++idx){
            if( children[idx] ){
                FmmKernel_M2M(fmmCore, (void*)cell, (void*)children[idx]);
            }
        }
    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict cell, const CellClass* interactions[], const int , const int ) {
        for(int idx = 0 ; idx < 343 ; ++idx){
            if( interactions[idx] ){
                FmmKernel_M2L(fmmCore, (void*)cell, (void*)interactions[idx]);
            }
        }
    }

    /** Do nothing */
    virtual void L2L(const CellClass* const FRestrict cell, CellClass* FRestrict *const FRestrict children, const int ) {
        for(int idx = 0 ; idx < 8 ; ++idx){
            if( children[idx] ){
                FmmKernel_L2L(fmmCore, (void*)cell, (void*)children[idx]);
            }
        }
    }

    /** Do nothing */
    virtual void L2P(const CellClass* const cell, ContainerClass* const container){
        cell->setContainer((ContainerClass*)container);
        FmmKernel_L2P(fmmCore, (void*)cell);
    }


    /** Do nothing */
    virtual void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
                     ContainerClass* const neighbors[27], const int ){
        CellClass* cell = (CellClass*)targets->getParentCell();
        FmmKernel_P2P_inner(fmmCore, cell);

        for(int idx = 0 ; idx < 27 ; ++idx){
            if( neighbors[idx] ){
                FmmKernel_P2P(fmmCore, (CellClass*)neighbors[idx]->getParentCell(), cell);
            }
        }
    }

    /** Do nothing */
    virtual void P2PRemote(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [27], const int ){
        printf("Error remote not implemented!!!!\n");
    }

};


class CoreVector {
    mutable void* parentCell;
    void* fields;
    int sizeOfField;
    void* potentials;
    int sizeOfPential;
    FVector<FReal> positions;
    FVector<int> indexes;
public:
    CoreVector() : parentCell(nullptr), fields(nullptr), sizeOfField(0),
        potentials(nullptr), sizeOfPential(0){
    }

    ~CoreVector(){
        delete[] (char*)fields;
        delete[] (char*)potentials;
    }

    void setParentCell(void* inCell) const{
        parentCell = inCell;
    }

    void* getParentCell() const{
        return parentCell;
    }

    void allocateFields(const int memSizePerParticles){
        if(fields) delete[](char*)fields;
        fields = new char[getNbParticles() * memSizePerParticles];
        sizeOfField = memSizePerParticles;
        memset(fields, 0, getNbParticles() * memSizePerParticles);
    }

    void* getFields() const{
        return fields;
    }

    int getSizeOfField() const{
        return sizeOfField;
    }

    void allocatePotentials(const int memSizePerParticles){
        if(potentials) delete[](char*)potentials;
        potentials = new char[getNbParticles() * memSizePerParticles];
        sizeOfPential = memSizePerParticles;
        memset(potentials, 0, getNbParticles() * memSizePerParticles);
    }

    void* getPotentials() const{
        return potentials;
    }

    int getSizeOfPotential() const{
        return sizeOfPential;
    }

    const int* getIndexes() const {
        return indexes.data();
    }

    const FReal* getPositions() const {
        return positions.data();
    }

    int getNbParticles() {
        return indexes.getSize();
    }

    void push(const FPoint<FReal>& partPosition, const int partIndex){
        positions.push(partPosition.getX());
        positions.push(partPosition.getY());
        positions.push(partPosition.getZ());
        indexes.push(partIndex);
    }
};


typedef CoreVector          CoreContainerClass;

typedef CoreCell<CoreContainerClass>      CoreCellClass;
typedef FSimpleLeaf<FReal, CoreContainerClass >                        LeafClass;
typedef FOctree<CoreCellClass, CoreContainerClass , LeafClass >     OctreeClass;
typedef CoreKernel<CoreCellClass, CoreContainerClass>         CoreKernelClass;

typedef FFmmAlgorithm<OctreeClass, CoreCellClass, CoreContainerClass, CoreKernelClass, LeafClass >     FmmClass;
typedef FFmmAlgorithmThread<OctreeClass, CoreCellClass, CoreContainerClass, CoreKernelClass, LeafClass >     FmmClassThread;


struct ScalFmmCoreHandle {
    struct ScalFmmCoreConfig {
        //paramètres en lecture/écriture :
        int treeHeight;     // hombre de niveaux de l'arbre (int)
        FReal boxWidth;    // taille de la boîte racine (FReal)
        FReal boxCenter[3]; // position du centre de la boîte racine (FReal[3])
#ifdef SCALFMM_USE_MPI
        MPI_Comm mpiCom;    // communicateur MPI (MPI_Comm)
#endif
        int nbThreads;      // nombre de threads (int)
        int rhsNumber;      // nombre de seconds membres (int)
    };

    ScalFmmCoreConfig config;
    OctreeClass* octree;
    void *kernelHandle;
};

int FmmCore_init(void **fmmCore) {
    ScalFmmCoreHandle* corehandle = new ScalFmmCoreHandle;
    memset(corehandle, 0, sizeof(corehandle));

    corehandle->config.nbThreads = omp_get_max_threads();

    *fmmCore = corehandle;

    return FMMAPI_NO_ERROR;
} /*alloue et initialise le FmmCore*/


int FmmCore_free(void *fmmCore) {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    if(corehandle->octree) delete corehandle->octree;
    delete corehandle;
    return FMMAPI_NO_ERROR;
} /*libère le FmmCore*/

int FmmCore_isParameterUsed(void */*fmmCore*/, int *name, int *flag){
    switch( *name ){
    case FMMCORE_ROOT_BOX_WIDTH :
    case FMMCORE_ROOT_BOX_CENTER :
    case FMMCORE_TREE_HEIGHT :
#ifdef SCALFMM_USE_MPI
    case FMMCORE_MPI_COMMUNICATOR:
#endif
    case FMMCORE_THREADS_NUMBER:
    case FMMCORE_THREAD_ID:
    case FMMCORE_RHS_NUMBER:
    case  FMMCORE_HANDLES_P2P:
        *flag = FMMAPI_SUPPORTED_PARAMETER;
        break;
    case FMMCORE_LEAF_BOX_WIDTH :
    case FMMCORE_POINTS_PER_LEAF :
        *flag = FMMAPI_UNSUPPORTED_PARAMETER;
        break;
    default:
        *flag = FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmCore_setParameter(void *fmmCore, int *name, void*value){
    int flag;

    FmmCore_isParameterUsed(fmmCore, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    switch( *name ){
    case FMMCORE_TREE_HEIGHT :
        corehandle->config.treeHeight = *(int*)value ;
        break;
    case FMMCORE_ROOT_BOX_WIDTH :
        corehandle->config.boxWidth = *(FReal*)value;
        break;
    case FMMCORE_ROOT_BOX_CENTER :
        memcpy(corehandle->config.boxCenter, value, sizeof(FReal)*3);
        break;
#ifdef SCALFMM_USE_MPI
    case FMMCORE_MPI_COMMUNICATOR:
        corehandle->config.mpiCom = *(MPI_Comm*)value;
        break;
#endif
    case FMMCORE_THREADS_NUMBER:
        corehandle->config.nbThreads = *(int*)value;
        break;
    case FMMCORE_RHS_NUMBER:
        corehandle->config.rhsNumber = *(int*)value;
        break;
    case FMMCORE_HANDLES_P2P:
    case FMMCORE_THREAD_ID:
        return FMMAPI_UNSUPPORTED_PARAMETER;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmCore_setParameter(void *fmmCore, int name, void*value){
    return FmmCore_setParameter(fmmCore, &name, value);
}


int FmmCore_getParameter(void *fmmCore, int *name, void*value){
    int flag;

    FmmCore_isParameterUsed(fmmCore, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    switch( *name ){
    case FMMCORE_TREE_HEIGHT :
        *(int*)value = corehandle->config.treeHeight;
        break;
    case FMMCORE_ROOT_BOX_WIDTH :
        *(FReal*)value = corehandle->config.boxWidth;
        break;
    case FMMCORE_ROOT_BOX_CENTER :
        memcpy(value,corehandle->config.boxCenter, sizeof(FReal)*3);
        break;
#ifdef SCALFMM_USE_MPI
    case FMMCORE_MPI_COMMUNICATOR:
        *(MPI_Comm*)value = corehandle->config.mpiCom;
        break;
#endif
    case FMMCORE_THREADS_NUMBER:
        *(int*)value = corehandle->config.nbThreads;
        break;
    case FMMCORE_THREAD_ID:
        *(int*)value = omp_get_thread_num();
        break;
    case FMMCORE_RHS_NUMBER:
        *(int*)value = corehandle->config.rhsNumber;
        break;
    case  FMMCORE_HANDLES_P2P:
        *(int*)value = true;
        break;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmCore_getParameter(void *fmmCore, int name, void*value){
    return FmmCore_getParameter(fmmCore, &name, value);
}


int FmmCore_getRadius(void*fmmCore, void *boxId, FReal *radius) {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    *radius = corehandle->octree->getBoxWidth() / FReal(1<<boxhandle->getLevel());
    return FMMAPI_NO_ERROR;
} /*Renvoie le rayon de la boîte*/

int FmmCore_getCentre(void*/*fmmCore*/, void *boxId, FReal **centre) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;

    (*centre)[0] = boxhandle->getPosition()[0];
    (*centre)[1] = boxhandle->getPosition()[1];
    (*centre)[2] = boxhandle->getPosition()[2];

    return FMMAPI_NO_ERROR;
} /*Renvoie dans le FReal[3] centre les coordonnées du centre*/

int FmmCore_getLevel(void*/*fmmCore*/, void *boxId, int *level) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    *level = boxhandle->getLevel();
    return FMMAPI_NO_ERROR;
} /*Renvoie dans level le niveau de la boîte boxId dans son arbre*/

int FmmCore_getMultipoleArray(void* /*fmmCore*/, void *boxId, void **F) {
    CoreCellClass* cell = (CoreCellClass*)boxId;
    *F = cell->getMultipole();
    return FMMAPI_NO_ERROR;
} /*Renvoie dans F l'adresse où stocker l'expansion multipôle associée à la boîte boxId*/

int FmmCore_getLocalArray(void* /*fmmCore*/, void *boxId, void **F) {
    CoreCellClass* cell = (CoreCellClass*)boxId;
    *F = cell->getLocal();
    return FMMAPI_NO_ERROR;
} /*Renvoie dans F l'adresse où stocker l'expansion locale associée à la boîte boxId*/

int FmmCore_getCoord(void*/*fmmCore*/, void *boxId, int *coord) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    coord[0] = boxhandle->getCoordinate().getX();
    coord[1] = boxhandle->getCoordinate().getY();
    coord[2] = boxhandle->getCoordinate().getZ();
    return FMMAPI_NO_ERROR;
} /*Renvoie dans coord la position dans l'arbre*/


/* Données potentiel/champ */
int FmmCore_getSource(void* /*fmmCore*/, void *boxId, FReal** position, void** potential, int *number) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    CoreContainerClass* sources = boxhandle->getContainer();

    *number = sources->getNbParticles();
    *position = const_cast<FReal*>(sources->getPositions());
    *potential = sources->getPotentials();

    return FMMAPI_NO_ERROR;
} /* Appelé par P2P et P2M pour obtenir le nombre, la position et le potentiel des sources.
Les différents tableaux sont (éventuellement) alloués par le FmmCore. */

int FmmCore_releaseSource(void*fmmCore, void *boxId, void* potential, FReal* position) {
    return FMMAPI_NO_ERROR;
} /* si le core veut libérer ces tableaux potentiel et position.*/

int FmmCore_getTargetPoints(void* /*fmmCore*/, void *boxId, FReal** position, int *number) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    CoreContainerClass* targets = boxhandle->getContainer();

    *number = targets->getNbParticles();
    *position = const_cast<FReal*>(targets->getPositions());

    return FMMAPI_NO_ERROR;
} /* : Appelé par P2P et L2P pour obtenir le nombre et la position des points cibles.*/

int FmmCore_releaseTargetPoints(void*fmmCore, void *boxId, FReal* position) {
    return FMMAPI_NO_ERROR;
} /* si le core veut libérer ce tableau "position".*/

int FmmCore_getTargetField(void* /*fmmCore*/, void *boxId, FReal* F) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    CoreContainerClass* targets = boxhandle->getContainer();

    memcpy(F, targets->getFields(), targets->getNbParticles() * targets->getSizeOfField());

    return FMMAPI_NO_ERROR;
} /* obtient dans un tableau F alloué/libéré par L2P/P2P les valeurs des champs aux points cible
(pour le cas où P2P et L2P doivent sommer leurs résultats).*/

int FmmCore_setTargetField(void* /*fmmCore*/, void *boxId, FReal* F) {
    CoreCellClass* boxhandle = (CoreCellClass*)boxId;
    CoreContainerClass* targets = boxhandle->getContainer();

    memcpy(targets->getFields(), F, targets->getNbParticles() * targets->getSizeOfField());

    return FMMAPI_NO_ERROR;
} /* transmets au FmmCore dans F les valeurs des champs aux points cibles mis à jour.*/



/* Entrée/sortie principale */
int FmmCore_setKernelData(void *fmmCore, void *fmmKernel)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    corehandle->kernelHandle = fmmKernel;
    return FMMAPI_NO_ERROR;
} /* stocke l'identifiant du FmmKernel dans le FmmCore. Cela permet par la suite aux
opérateurs FMM d'accéder au FmmKernel, et donc d'accéder aux données spécifiques au kernel
(p.ex. fréquence dans le cas Helmholtz, …)*/

int FmmCore_getKernelData(void*fmmCore, void **fmmKernel)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    *fmmKernel = corehandle->kernelHandle;
    return FMMAPI_NO_ERROR;
} /* récupère l'identifiant du FmmKernel. */


int FmmCore_setPositions(void *fmmCore, int *nb, FReal *position)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;

    corehandle->octree = new OctreeClass(corehandle->config.treeHeight, FMath::Min(3,corehandle->config.treeHeight-1),
                                         corehandle->config.boxWidth, FPoint<FReal>(corehandle->config.boxCenter));

    if( corehandle->config.nbThreads != 0){
        omp_set_num_threads(corehandle->config.nbThreads);
    }

    for(FSize idxPart = 0 ; idxPart < (*nb) ; ++idxPart){
        const FReal* pos = &position[idxPart * 3];
        corehandle->octree->insert(FPoint<FReal>(pos[0], pos[1], pos[2]),idxPart);
    }

    return FMMAPI_NO_ERROR;
} /* transmet au FmmCore les potentiels associés aux points sources.
Le tableau potential est alloué et libéré par la routine appelant, le FmmCore doit donc en faire une copie.*/


int FmmCore_setPotentials(void *fmmCore, void *potentials)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    OctreeClass* octree = corehandle->octree;

    void* fmmKernel;
    FmmCore_getKernelData(fmmCore, &fmmKernel);

    int sizeOfPotentials;
    FmmKernel_getParameter(fmmKernel, FMMKERNEL_POTENTIAL_DATA_SIZE, &sizeOfPotentials);

    typename OctreeClass::Iterator octreeIterator(octree);
    octreeIterator.gotoBottomLeft();
    do{
        CoreContainerClass* particles = octreeIterator.getCurrentLeaf()->getSrc();
        particles->allocatePotentials(sizeOfPotentials);

        for(int idx = 0 ; idx < particles->getNbParticles() ; ++idx){
            memcpy(((char*)particles->getPotentials())+idx*sizeOfPotentials,
                   ((char*)potentials)+particles->getIndexes()[idx]*sizeOfPotentials,
                   sizeOfPotentials);
        }

    } while( octreeIterator.moveRight() );

    return FMMAPI_NO_ERROR;
} /* transmet au FmmCore les potentiels associés aux points sources.
Le tableau potential est alloué et libéré par la routine appelant, le FmmCore doit donc en faire une copie.*/

int FmmCore_doComputation(void *fmmCore)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;

    { // Ceck if there is number of NbPart summed at level 1
        FPoint<FReal> corner = corehandle->octree->getBoxCenter();
        corner -= (corehandle->octree->getBoxWidth()/2);

        FReal boxWidth = corehandle->octree->getBoxWidth()/( 1<< (corehandle->config.treeHeight-1) );

        typename OctreeClass::Iterator octreeIterator(corehandle->octree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = corehandle->config.treeHeight - 1 ; idxLevel >= 1 ; --idxLevel ){
            int multipoleSize;
            int localSize;

            FmmKernel_getMultipoleArraySize(fmmCore, &multipoleSize);
            FmmKernel_getLocalArraySize(fmmCore, &localSize);

            do{
                octreeIterator.getCurrentCell()->createArrays(multipoleSize, localSize);
                octreeIterator.getCurrentCell()->setLevel(idxLevel);
                const FTreeCoordinate coord = octreeIterator.getCurrentCell()->getCoordinate();
                FPoint<FReal> position( coord.getX()*boxWidth + boxWidth/2.0 + corner.getX(),
                                 coord.getY()*boxWidth + boxWidth/2.0 + corner.getY(),
                                 coord.getZ()*boxWidth + boxWidth/2.0 + corner.getZ());
                octreeIterator.getCurrentCell()->setPosition(position.getDataValue());
            } while( octreeIterator.moveRight() );

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            boxWidth *= FReal(2.0);
        }
    }
    { // Ceck if there is number of NbPart summed at level 1
        void* fmmKernel;
        FmmCore_getKernelData(fmmCore, &fmmKernel);

        int sizeOfFields;
        FmmKernel_getParameter(fmmKernel, FMMKERNEL_FIELD_DATA_SIZE, &sizeOfFields);

        typename OctreeClass::Iterator octreeIterator(corehandle->octree);
        octreeIterator.gotoBottomLeft();
        do{
            octreeIterator.getCurrentLeaf()->getTargets()->setParentCell(octreeIterator.getCurrentCell());
            octreeIterator.getCurrentLeaf()->getTargets()->allocateFields(sizeOfFields);
        } while( octreeIterator.moveRight() );
    }

    if( corehandle->config.nbThreads <= 1){
        CoreKernelClass kernels(fmmCore );
        FmmClass algo(corehandle->octree,&kernels);
        algo.execute();
    }
    else{
        CoreKernelClass kernels(fmmCore );
        FmmClassThread algo(corehandle->octree,&kernels);
        algo.execute();
    }

    return FMMAPI_NO_ERROR;
} /* réalise le produit multipôle. */

/* !!! Warning use *filed and not **field */
int FmmCore_getField(void *fmmCore, void *fields)  {
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;

    typename OctreeClass::Iterator octreeIterator(corehandle->octree);
    octreeIterator.gotoBottomLeft();
    do{
        CoreContainerClass* particles = octreeIterator.getCurrentLeaf()->getTargets();
        for(int idx = 0 ; idx < particles->getNbParticles() ; ++idx){
            memcpy(((char*)fields)+particles->getIndexes()[idx]*particles->getSizeOfField(),
                   ((char*)particles->getFields())+idx*particles->getSizeOfField(),
                   particles->getSizeOfField());
        }
    } while(octreeIterator.moveRight());

    return FMMAPI_NO_ERROR;
} /* récupère après le produit multipôle la valeur des champs en chaque point.
Le tableau field doit être alloué et libéré par la routine appelante. */


