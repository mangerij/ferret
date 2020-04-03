#ifndef FCUDAP2P_HPP
#define FCUDAP2P_HPP

#include "../Cuda/FCudaGlobal.hpp"
#include "../Cuda/FCudaGroupAttachedLeaf.hpp"
#include "../Cuda/FCudaEmptyCellSymb.hpp"
#include "../Cuda/FCudaCompositeCell.hpp"

#define DirectMacro(targetX, targetY, targetZ, targetPhys, \
    forceX, forceY, forceZ, potential,\
    sourcesX, sourcesY, sourcesZ, sourcesPhys)\
{\
    FReal dx = sourcesX - targetX;\
    FReal dy = sourcesY - targetY;\
    FReal dz = sourcesZ - targetZ;\
    \
    FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);\
    FReal inv_distance = sqrt(inv_square_distance);\
    \
    inv_square_distance *= inv_distance;\
    inv_square_distance *= targetPhys * sourcesPhys;\
    \
    dx *= inv_square_distance;\
    dy *= inv_square_distance;\
    dz *= inv_square_distance;\
    \
    forceX += dx;\
    forceY += dy;\
    forceZ += dz;\
    sourcesPhys += inv_distance * sourcesPhys;\
    }

#define Min(x,y) ((x)<(y)?(x):(y))
#define Max(x,y) ((x)>(y)?(x):(y))

/**
 * This class defines what should be a Cuda kernel.
 */
template <class FReal>
class FCudaP2P {
protected:
public:
    static double DSqrt(const double val){
        return sqrt(val);
    }

    static float FSqrt(const float val){
        return sqrtf(val);
    }

    typedef FCudaGroupAttachedLeaf<FReal,4,4,FReal> ContainerClass;
    typedef FCudaCompositeCell<FCudaEmptyCellSymb,int,int> CellClass;

    static const int SHARE_SIZE = 128;

    __device__ void P2M(CellClass /*pole*/, const ContainerClass* const /*particles*/) {
    }

    __device__ void M2M(CellClass  /*pole*/, const CellClass  /*child*/[8], const int /*level*/) {
    }

    __device__ void M2L(CellClass  /*pole*/, const CellClass /*distantNeighbors*/[343],
    const int /*size*/, const int /*level*/) {
    }

    __device__ void L2L(const CellClass  /*local*/, CellClass  /*child*/[8], const int /*level*/) {
    }

    __device__ void L2P(const CellClass  /*local*/, ContainerClass*const /*particles*/){
    }

    __device__ void P2P(const int3& pos,
                        ContainerClass* const  targets, const ContainerClass* const  sources,
                        ContainerClass* const directNeighborsParticles[27], const int counter){
        // Compute with other
        P2PRemote(pos, targets, sources, directNeighborsParticles, counter);
        // Compute inside
        const int nbLoops = (targets->getNbParticles()+blockDim.x-1)/blockDim.x;

        for(int idxLoop = 0 ; idxLoop < nbLoops; ++idxLoop){
            const int idxPart = (idxLoop*blockDim.x+threadIdx.x);
            const bool threadCompute = (idxPart < targets->getNbParticles());

            FReal targetX, targetY, targetZ, targetPhys;
            FReal forceX = 0, forceY = 0, forceZ = 0, potential = 0;

            if(threadCompute){
                targetX = targets->getPositions()[0][idxPart];
                targetY = targets->getPositions()[1][idxPart];
                targetZ = targets->getPositions()[2][idxPart];
                targetPhys = targets->getAttribute(0)[idxPart];
            }

            for(int idxCopy = 0 ; idxCopy < targets->getNbParticles() ; idxCopy += SHARE_SIZE){
                __shared__ FReal sourcesX[SHARE_SIZE];
                __shared__ FReal sourcesY[SHARE_SIZE];
                __shared__ FReal sourcesZ[SHARE_SIZE];
                __shared__ FReal sourcesPhys[SHARE_SIZE];

                const int nbCopies = Min(SHARE_SIZE, targets->getNbParticles()-idxCopy);
                if(threadIdx.x < nbCopies){
                    sourcesX[threadIdx.x] = targets->getPositions()[0][idxPart];
                    sourcesY[threadIdx.x] = targets->getPositions()[1][idxPart];
                    sourcesZ[threadIdx.x] = targets->getPositions()[2][idxPart];
                    sourcesPhys[threadIdx.x] = targets->getAttribute(0)[idxPart];
                }

                __syncthreads();

                if(threadCompute){
                    const int leftCopies = Min(idxPart, nbCopies);
                    // Left Part
                    for(int otherIndex = 0; otherIndex < leftCopies - 3; otherIndex += 4) { // unrolling x4
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+1], sourcesY[otherIndex+1], sourcesZ[otherIndex+1], sourcesPhys[otherIndex+1]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+2], sourcesY[otherIndex+2], sourcesZ[otherIndex+2], sourcesPhys[otherIndex+2]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+3], sourcesY[otherIndex+3], sourcesZ[otherIndex+3], sourcesPhys[otherIndex+3]);
                    }

                    for(int otherIndex = (leftCopies/4) * 4; otherIndex < nbCopies; ++otherIndex) { // if nk%4 is not zero
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                    }
                    // Right Part
                    for(int otherIndex = leftCopies+1; otherIndex < nbCopies - 3; otherIndex += 4) { // unrolling x4
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+1], sourcesY[otherIndex+1], sourcesZ[otherIndex+1], sourcesPhys[otherIndex+1]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+2], sourcesY[otherIndex+2], sourcesZ[otherIndex+2], sourcesPhys[otherIndex+2]);
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex+3], sourcesY[otherIndex+3], sourcesZ[otherIndex+3], sourcesPhys[otherIndex+3]);
                    }

                    for(int otherIndex = Max(leftCopies+1, (nbCopies/4) * 4); otherIndex < nbCopies; ++otherIndex) { // if nk%4 is not zero
                        DirectMacro(targetX, targetY, targetZ, targetPhys,
                                    forceX, forceY, forceZ, potential,
                                    sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                    }
                }

                __syncthreads();
            }

            if( threadCompute ){
                targets->getAttribute(1)[idxPart] += forceX;
                targets->getAttribute(2)[idxPart] += forceY;
                targets->getAttribute(3)[idxPart] += forceZ;
                targets->getAttribute(4)[idxPart] += potential;
            }

        }
    }

    __device__ void P2PRemote(const int3& ,
                              ContainerClass* const  targets, const ContainerClass* const  /*sources*/,
                              ContainerClass* const directNeighborsParticles[27], const int ){
        for(int idxNeigh = 0 ; idxNeigh < 27 ; ++idxNeigh){
            if(directNeighborsParticles[idxNeigh]){
                const int nbLoops = (targets->getNbParticles()+blockDim.x-1)/blockDim.x;

                for(int idxLoop = 0 ; idxLoop < nbLoops; ++idxLoop){
                    const int idxPart = (idxLoop*blockDim.x+threadIdx.x);
                    const bool threadCompute = (idxPart < targets->getNbParticles());

                    FReal targetX, targetY, targetZ, targetPhys;
                    FReal forceX = 0, forceY = 0, forceZ = 0, potential = 0;

                    if(threadCompute){
                        targetX = targets->getPositions()[0][idxPart];
                        targetY = targets->getPositions()[1][idxPart];
                        targetZ = targets->getPositions()[2][idxPart];
                        targetPhys = targets->getAttribute(0)[idxPart];
                    }

                    for(int idxCopy = 0 ; idxCopy < directNeighborsParticles[idxNeigh]->getNbParticles() ; idxCopy += SHARE_SIZE){
                        __shared__ FReal sourcesX[SHARE_SIZE];
                        __shared__ FReal sourcesY[SHARE_SIZE];
                        __shared__ FReal sourcesZ[SHARE_SIZE];
                        __shared__ FReal sourcesPhys[SHARE_SIZE];

                        const int nbCopies = Min(SHARE_SIZE, directNeighborsParticles[idxNeigh]->getNbParticles()-idxCopy);
                        if(threadIdx.x < nbCopies){
                            sourcesX[threadIdx.x] = directNeighborsParticles[idxNeigh]->getPositions()[0][idxPart];
                            sourcesY[threadIdx.x] = directNeighborsParticles[idxNeigh]->getPositions()[1][idxPart];
                            sourcesZ[threadIdx.x] = directNeighborsParticles[idxNeigh]->getPositions()[2][idxPart];
                            sourcesPhys[threadIdx.x] = directNeighborsParticles[idxNeigh]->getAttribute(0)[idxPart];
                        }

                        __syncthreads();

                        if(threadCompute){
                            for(int otherIndex = 0; otherIndex < nbCopies - 3; otherIndex += 4) { // unrolling x4
                                DirectMacro(targetX, targetY, targetZ, targetPhys,
                                            forceX, forceY, forceZ, potential,
                                            sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                                DirectMacro(targetX, targetY, targetZ, targetPhys,
                                            forceX, forceY, forceZ, potential,
                                            sourcesX[otherIndex+1], sourcesY[otherIndex+1], sourcesZ[otherIndex+1], sourcesPhys[otherIndex+1]);
                                DirectMacro(targetX, targetY, targetZ, targetPhys,
                                            forceX, forceY, forceZ, potential,
                                            sourcesX[otherIndex+2], sourcesY[otherIndex+2], sourcesZ[otherIndex+2], sourcesPhys[otherIndex+2]);
                                DirectMacro(targetX, targetY, targetZ, targetPhys,
                                            forceX, forceY, forceZ, potential,
                                            sourcesX[otherIndex+3], sourcesY[otherIndex+3], sourcesZ[otherIndex+3], sourcesPhys[otherIndex+3]);
                            }

                            for(int otherIndex = (nbCopies/4) * 4; otherIndex < nbCopies; ++otherIndex) { // if nk%4 is not zero
                                DirectMacro(targetX, targetY, targetZ, targetPhys,
                                            forceX, forceY, forceZ, potential,
                                            sourcesX[otherIndex], sourcesY[otherIndex], sourcesZ[otherIndex], sourcesPhys[otherIndex]);
                            }
                        }

                        __syncthreads();
                    }

                    if( threadCompute ){
                        targets->getAttribute(1)[idxPart] += forceX;
                        targets->getAttribute(2)[idxPart] += forceY;
                        targets->getAttribute(3)[idxPart] += forceZ;
                        targets->getAttribute(4)[idxPart] += potential;
                    }

                }
            }
        }
    }

    __host__ static FCudaP2P* InitKernelKernel(void*){
        return nullptr;
    }

    __host__ static void ReleaseKernel(FCudaP2P* /*todealloc*/){
        // nothing to do
    }

    __host__ static dim3 GetGridSize(const int intervalSize){
        return intervalSize;
    }

    __host__ static dim3 GetBlocksSize(){
        return SHARE_SIZE;
    }
};

#endif // FCUDAP2P_HPP

