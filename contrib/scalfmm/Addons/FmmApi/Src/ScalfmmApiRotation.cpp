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
#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FAssert.hpp"

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "FmmApi.h"

////////////////////// Opérateurs FMM Kernel : //////////////////////////

class KernelCell : public FBasicCell {
    FComplex<FReal>* multipole;
    FComplex<FReal>* local;
public:
    KernelCell() : multipole(nullptr), local(nullptr){
    }
    void attachArrays(FComplex<FReal> inMultipole[], FComplex<FReal> inLocal[]){
        multipole = inMultipole;
        local = inLocal;
    }

    const FComplex<FReal>* getMultipole() const{
        return multipole;
    }

    const FComplex<FReal>* getLocal() const{
        return local;
    }

    FComplex<FReal>* getMultipole(){
        return multipole;
    }

    FComplex<FReal>* getLocal(){
        return local;
    }
};


static const int P = 5;

typedef KernelCell               KernelCellClass;
typedef FP2PParticleContainer<FReal>          ContainerClass;
typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
typedef FRotationKernel< KernelCellClass, ContainerClass , P>   KernelClass;

struct ScalFmmKernelHandle {
    KernelClass** kernel;
    int potentialDataSize;
    int fieldDataSize;
    int nbthread;
};

/******* Allocation : ******/
int FmmKernel_init(void *fmmCore, void **fmmKernel){
    ScalFmmKernelHandle* kernelhandle = new ScalFmmKernelHandle;
    memset(kernelhandle, 0, sizeof(ScalFmmKernelHandle));

    int NbLevels;
    FmmCore_getParameter(fmmCore, FMMCORE_TREE_HEIGHT, &NbLevels);
    FReal boxWidth;
    FmmCore_getParameter(fmmCore, FMMCORE_ROOT_BOX_WIDTH, &boxWidth);
    FReal centerOfBox[3];
    FmmCore_getParameter(fmmCore, FMMCORE_ROOT_BOX_CENTER, centerOfBox);

    FmmCore_getParameter(fmmCore, FMMCORE_THREADS_NUMBER, &kernelhandle->nbthread);

    KernelClass original( NbLevels, boxWidth, FPoint<FReal>(centerOfBox) );
    kernelhandle->kernel = new KernelClass*[kernelhandle->nbthread];
    for(int idxThread = 0 ; idxThread < kernelhandle->nbthread ; ++idxThread){
        kernelhandle->kernel[idxThread] = new KernelClass(original);
    }

    kernelhandle->potentialDataSize = 1;
    kernelhandle->fieldDataSize = 4;

    *fmmKernel = kernelhandle;

    return FMMAPI_NO_ERROR;
}/* : alloue et initialise le FmmKernel */
int FmmKernel_free(void *fmmKernel){
    ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;

    for(int idxThread = 0 ; idxThread < kernelhandle->nbthread ; ++idxThread){
        delete kernelhandle->kernel[idxThread];
    }

    delete[] kernelhandle->kernel;
    delete kernelhandle;

    return FMMAPI_NO_ERROR;
} /* libére le FmmKernel */




/******* Configuration : ***/

int FmmKernel_isParameterUsed(void * /*fmm*/, int *name, int *flag){
    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE:
    case FMMKERNEL_FIELD_DATA_SIZE:
    case  FMMKERNEL_HANDLES_P2P:
        *flag = FMMAPI_SUPPORTED_PARAMETER;
        break;
    case FMMKERNEL_ACCURACY :
        *flag = FMMAPI_UNSUPPORTED_PARAMETER;
        break;
    default:
        *flag = FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_setParameter(void *fmmKernel, int *name, void*value){
    /*ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;*/
    int flag;

    FmmKernel_isParameterUsed(fmmKernel, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE :
    case FMMKERNEL_FIELD_DATA_SIZE :
        return FMMAPI_SUPPORTED_PARAMETER;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_setParameter(void *fmmKernel, int name, void*value){
    return FmmKernel_setParameter( fmmKernel, &name, value);
}

int FmmKernel_getParameter(void *fmmKernel, int *name, void*value){
    ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;
    int flag;

    FmmKernel_isParameterUsed(fmmKernel, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE :
        *(int*)value = kernelhandle->potentialDataSize*sizeof(FReal);
        break;
    case FMMKERNEL_FIELD_DATA_SIZE :
        *(int*)value = kernelhandle->fieldDataSize*sizeof(FReal);
        break;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_getParameter(void *fmmKernel, int name, void*value){
    return FmmKernel_getParameter(fmmKernel, &name, value);
}


/****** Données FMM : *****/
int FmmKernel_getMultipoleArraySize(void */*fmmCore*/, int *size) {
    *size = ((P+2)*(P+1))/2 * sizeof(FComplex<FReal>);
    return FMMAPI_NO_ERROR;
} /* Renvoie dans size la taille (en octets) de l'expansion multipôle associée à la boîte boxId */

int FmmKernel_getLocalArraySize(void */*fmmCore*/, int *size){
    *size = ((P+2)*(P+1))/2 * sizeof(FComplex<FReal>);
    return FMMAPI_NO_ERROR;
} /* Renvoie dans size la taille (en octets) de l'expansion locale associée à la boîte boxId*/


/******* Opérateurs FMM : **/
int FmmKernel_P2M(void *fmmCore, void* boxId){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplex<FReal>* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxId, (void**)&multipole);

    KernelCellClass cell;
    cell.attachArrays(multipole, nullptr);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxId, coord);
    cell.setCoordinate(coord[0], coord[1], coord[2]);

    FReal* positions;
    FReal* physicalValues;
    int number;
    FmmCore_getSource(fmmCore, boxId, &positions, (void**)&physicalValues, &number);

    FP2PParticleContainer<FReal> sources;
    for(FSize idxPart = 0 ; idxPart < number ; ++idxPart){
        sources.push(FPoint<FReal>(positions[idxPart*3],positions[idxPart*3+1],positions[idxPart*3+2]),physicalValues[idxPart]);
    }

    kernelhandle->kernel[threadId]->P2M(&cell, &sources);

    FmmCore_releaseSource(fmmCore, boxId, physicalValues, positions);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_L2P(void *fmmCore, void* boxId){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplex<FReal>* local;
    FmmCore_getLocalArray(fmmCore, boxId, (void**)&local);

    KernelCellClass cell;
    cell.attachArrays(nullptr,local);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxId, coord);
    cell.setCoordinate(coord[0], coord[1], coord[2]);

    FReal* physicalValues;
    FReal* positions;
    int number;
    FmmCore_getSource(fmmCore, boxId, &positions, (void**)&physicalValues, &number);

    FReal* fields = new FReal[number*kernelhandle->fieldDataSize];
    FmmCore_getTargetField(fmmCore, boxId, fields);

    FP2PParticleContainer<FReal> targets;
    for(FSize idxPart = 0 ; idxPart < number ; ++idxPart){
        targets.push(FPoint<FReal>(&positions[idxPart*3]),physicalValues[idxPart],
                fields[idxPart*kernelhandle->fieldDataSize],
                fields[idxPart*kernelhandle->fieldDataSize+1],
                fields[idxPart*kernelhandle->fieldDataSize+2],
                fields[idxPart*kernelhandle->fieldDataSize+3]);
    }

    kernelhandle->kernel[threadId]->L2P(&cell, &targets);

    const FReal*const potentials = targets.getPotentials();
    const FReal*const forcesX = targets.getForcesX();
    const FReal*const forcesY = targets.getForcesY();
    const FReal*const forcesZ = targets.getForcesZ();

    for(FSize idxPart = 0 ; idxPart < number ; ++idxPart){
        fields[idxPart*kernelhandle->fieldDataSize] += potentials[idxPart];
        fields[idxPart*kernelhandle->fieldDataSize+1] += forcesX[idxPart];
        fields[idxPart*kernelhandle->fieldDataSize+2] += forcesY[idxPart];
        fields[idxPart*kernelhandle->fieldDataSize+3] += forcesZ[idxPart];
    }

    FmmCore_releaseTargetPoints(fmmCore, boxId, positions);
    FmmCore_setTargetField(fmmCore, boxId, fields);
    delete[] fields;

    return FMMAPI_NO_ERROR;
}

int FmmKernel_M2M(void *fmmCore, void *boxIdFather, void *boxIdSon){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplex<FReal>* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxIdFather, (void**)&multipole);

    KernelCellClass cellFather;
    cellFather.attachArrays(multipole, nullptr);
    int coordFather[3];
    FmmCore_getCoord(fmmCore, boxIdFather, coordFather);
    cellFather.setCoordinate(coordFather[0], coordFather[1], coordFather[2]);

    FmmCore_getMultipoleArray(fmmCore, boxIdSon, (void**)&multipole);

    KernelCellClass cellSon;
    cellSon.attachArrays(multipole, nullptr);
    int coordChild[3];
    FmmCore_getCoord(fmmCore, boxIdSon, coordChild);
    cellSon.setCoordinate(coordChild[0], coordChild[1], coordChild[2]);

    int level;
    FmmCore_getLevel(fmmCore,boxIdFather, &level);

    const KernelCellClass* children[8];
    memset(children, 0, sizeof(KernelCellClass*)*8);
    const int mindex = ((coordChild[0]&1) * 2 + (coordChild[1]&1)) * 2 + (coordChild[2]&1);
    children[mindex] = &cellSon;

    kernelhandle->kernel[threadId]->M2M(&cellFather, children, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_L2L(void *fmmCore, void *boxIdFather, void *boxIdSon){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplex<FReal>* local;
    FmmCore_getLocalArray(fmmCore, boxIdFather, (void**)&local);

    KernelCellClass cellFather;
    cellFather.attachArrays(nullptr, local);
    int coordFather[3];
    FmmCore_getCoord(fmmCore, boxIdFather, coordFather);
    cellFather.setCoordinate(coordFather[0], coordFather[1], coordFather[2]);

    FmmCore_getLocalArray(fmmCore, boxIdSon, (void**)&local);

    KernelCellClass cellSon;
    cellSon.attachArrays(nullptr, local);
    int coordChild[3];
    FmmCore_getCoord(fmmCore, boxIdSon, coordChild);
    cellSon.setCoordinate(coordChild[0], coordChild[1], coordChild[2]);

    int level;
    FmmCore_getLevel(fmmCore,boxIdFather, &level);

    KernelCellClass* children[8];
    memset(children, 0, sizeof(KernelCellClass*)*8);
    const int mindex = ((coordChild[0]&1) * 2 + (coordChild[1]&1)) * 2 + (coordChild[2]&1);
    children[mindex] = &cellSon;

    kernelhandle->kernel[threadId]->L2L(&cellFather, children, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_M2L(void *fmmCore, void *boxIdSrc, void *boxIdDest){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplex<FReal>* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxIdSrc, (void**)&multipole);
    KernelCellClass cellSrc;
    cellSrc.attachArrays(multipole,nullptr);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxIdSrc, coord);
    cellSrc.setCoordinate(coord[0], coord[1], coord[2]);

    FComplex<FReal>* local;
    FmmCore_getLocalArray(fmmCore, boxIdDest, (void**)&local);
    KernelCellClass cellDst;
    cellDst.attachArrays(nullptr, local);
    FmmCore_getCoord(fmmCore, boxIdDest, coord);
    cellDst.setCoordinate(coord[0], coord[1], coord[2]);

    int level;
    FmmCore_getLevel(fmmCore, boxIdDest, &level);

    const int xdiff = cellSrc.getCoordinate().getX() - cellDst.getCoordinate().getX();
    const int ydiff = cellSrc.getCoordinate().getY() - cellDst.getCoordinate().getY();
    const int zdiff = cellSrc.getCoordinate().getZ() - cellDst.getCoordinate().getZ();
    const int index = (((xdiff+3) * 7) + (ydiff+3)) * 7 + zdiff + 3;

    const KernelCellClass* inter[343];
    memset(inter, 0, sizeof(KernelCellClass*)*343);
    inter[index] = &cellSrc;

    kernelhandle->kernel[threadId]->M2L(&cellDst, inter, 1, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_P2P_inner(void *fmmCore, void *boxIdSrcDest){
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FReal* positionsTargets;
    FReal* potentialsTargets;
    int numberTargets;
    FmmCore_getSource(fmmCore, boxIdSrcDest, &positionsTargets, (void**)&potentialsTargets, &numberTargets);

    FReal* fieldsTargets = new FReal[numberTargets*kernelhandle->fieldDataSize];
    FmmCore_getTargetField(fmmCore, boxIdSrcDest, fieldsTargets);

    int coordTargets[3];
    FmmCore_getCoord(fmmCore, boxIdSrcDest, coordTargets);
    FTreeCoordinate treecoord(coordTargets[0], coordTargets[1], coordTargets[2]);

    FP2PParticleContainer<FReal> targets;
    for(FSize idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        targets.push(FPoint<FReal>(positionsTargets[idxPart*3],positionsTargets[idxPart*3+1],positionsTargets[idxPart*3+2]),
                potentialsTargets[idxPart],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+1],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+2],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+3]);
    }


    FP2PParticleContainer<FReal>* ptrs[27];
    memset(ptrs, 0, sizeof(FP2PParticleContainer<FReal>*) * 27);
    kernelhandle->kernel[threadId]->P2P(treecoord, &targets, &targets, ptrs, 0);

    const FReal*const potentials = targets.getPotentials();
    const FReal*const forcesX = targets.getForcesX();
    const FReal*const forcesY = targets.getForcesY();
    const FReal*const forcesZ = targets.getForcesZ();

    for(FSize idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        fieldsTargets[idxPart*kernelhandle->fieldDataSize] = potentials[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+1] = forcesX[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+2] = forcesY[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+3] = forcesZ[idxPart];
    }

    FmmCore_releaseTargetPoints(fmmCore, boxIdSrcDest, positionsTargets);
    FmmCore_setTargetField(fmmCore, boxIdSrcDest, fieldsTargets);
    delete[] fieldsTargets;

    return FMMAPI_NO_ERROR;
}

int FmmKernel_P2P(void *fmmCore, void *boxIdSrc, void *boxIdDest){ //return FMMAPI_NO_ERROR;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(fmmCore, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FReal* positionsTargets;
    FReal* potentialsTargets;
    int numberTargets;
    FmmCore_getSource(fmmCore, boxIdDest, &positionsTargets, (void**)&potentialsTargets, &numberTargets);

    FReal* fieldsTargets = new FReal[numberTargets*kernelhandle->fieldDataSize];
    FmmCore_getTargetField(fmmCore, boxIdDest, fieldsTargets);

    int coordTargets[3];
    FmmCore_getCoord(fmmCore, boxIdDest, coordTargets);
    FTreeCoordinate treecoord(coordTargets[0], coordTargets[1], coordTargets[2]);

    FP2PParticleContainer<FReal> targets;
    for(FSize idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        targets.push(FPoint<FReal>(positionsTargets[idxPart*3],positionsTargets[idxPart*3+1],positionsTargets[idxPart*3+2]),
                potentialsTargets[idxPart],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+1],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+2],
                fieldsTargets[idxPart*kernelhandle->fieldDataSize+3]);
    }

    FReal* positionsSrc;
    FReal* potentialsSrc;
    int numberSources;
    FmmCore_getSource(fmmCore, boxIdSrc, &positionsSrc, (void**)&potentialsSrc, &numberSources);

    FP2PParticleContainer<FReal> sources;
    for(FSize idxPart = 0 ; idxPart < numberSources ; ++idxPart){
        sources.push(FPoint<FReal>(positionsSrc[idxPart*3],positionsSrc[idxPart*3+1],positionsSrc[idxPart*3+2]),potentialsSrc[idxPart]);
    }

    FP2PParticleContainer<FReal>* ptrs[27];
    memset(ptrs, 0, sizeof(FP2PParticleContainer<FReal>*) * 27);
    ptrs[0] = &sources;
    kernelhandle->kernel[threadId]->P2PRemote(treecoord, &targets, &targets, ptrs, 1);

    const FReal*const potentials = targets.getPotentials();
    const FReal*const forcesX = targets.getForcesX();
    const FReal*const forcesY = targets.getForcesY();
    const FReal*const forcesZ = targets.getForcesZ();

    for(FSize idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        fieldsTargets[idxPart*kernelhandle->fieldDataSize] = potentials[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+1] = forcesX[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+2] = forcesY[idxPart];
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+3] = forcesZ[idxPart];
    }

    FmmCore_releaseSource(fmmCore, boxIdDest, potentialsSrc, positionsSrc);
    FmmCore_releaseTargetPoints(fmmCore, boxIdDest, positionsTargets);
    FmmCore_setTargetField(fmmCore, boxIdDest, fieldsTargets);
    delete[] fieldsTargets;

    return FMMAPI_NO_ERROR;
} /* pas mutuel, i.e. on fait seulement dans 1 sens. */


