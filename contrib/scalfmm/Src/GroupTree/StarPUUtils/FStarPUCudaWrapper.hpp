#ifndef FSTARPUCUDAWRAPPER_HPP
#define FSTARPUCUDAWRAPPER_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FAssert.hpp"
#include "../../Utils/FAssert.hpp"

#include "../Core/FOutOfBlockInteraction.hpp"

#ifdef SCALFMM_USE_MPI
#include "../../Utils/FMpi.hpp"
#endif

#include <vector>
#include <memory>

#include <omp.h>

#include <starpu.h>

#ifdef STARPU_USE_MPI
#include <starpu_mpi.h>
#endif

#include "../Cuda/FCudaDeviceWrapper.hpp"

#include "FStarPUUtils.hpp"

template <class KernelClass, class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CudaCellGroupClass, class CudaParticleGroupClass, class CudaParticleContainerClass,
          class CudaKernelClass>
class FStarPUCudaWrapper {
protected:
    typedef FStarPUCudaWrapper<KernelClass, SymboleCellClass, PoleCellClass, LocalCellClass,
        CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass> ThisClass;

    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        int otherBlockId;
        std::vector<OutOfBlockInteraction> interactions;
    };

    const int treeHeight;
    CudaKernelClass* kernels[STARPU_MAXCUDADEVS];        //< The kernels

public:
    FStarPUCudaWrapper(const int inTreeHeight): treeHeight(inTreeHeight){
        memset(kernels, 0, sizeof(CudaKernelClass*)*STARPU_MAXCUDADEVS);
    }

    void initKernel(const int workerId, KernelClass* originalKernel){
        FAssertLF(kernels[workerId] == nullptr);
        kernels[workerId] = FCuda__BuildCudaKernel<CudaKernelClass>(originalKernel);
    }

    void releaseKernel(const int workerId){
        FCuda__ReleaseCudaKernel(kernels[workerId]);
        kernels[workerId] = nullptr;
    }

    ~FStarPUCudaWrapper(){
        for(int idxKernel = 0 ; idxKernel < STARPU_MAXCUDADEVS ; ++idxKernel ){
            FAssertLF(kernels[idxKernel] == nullptr);
        }
    }

    static void bottomPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__bottomPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                    kernel, starpu_cuda_get_local_stream(),
                    FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Upward Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void upwardPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel, &intervalSize);

        FCudaParams<unsigned char*,9> subCellGroupsPtr;
        memset(&subCellGroupsPtr, 0, sizeof(subCellGroupsPtr));
        FCudaParams<std::size_t,9> subCellGroupsSize;
        memset(&subCellGroupsPtr, 0, sizeof(subCellGroupsSize));
        FCudaParams<unsigned char*,9> subCellGroupsUpPtr;
        memset(&subCellGroupsUpPtr, 0, sizeof(subCellGroupsUpPtr));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroupsPtr.values[idxSubGroup] = ((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+2]));
            subCellGroupsSize.values[idxSubGroup] = STARPU_VARIABLE_GET_ELEMSIZE(buffers[(idxSubGroup*2)+2]);
            subCellGroupsUpPtr.values[idxSubGroup] = (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+3]);
        }

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__upwardPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                subCellGroupsPtr,subCellGroupsSize,subCellGroupsUpPtr,
                nbSubCellGroups, idxLevel, kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass Mpi
    /////////////////////////////////////////////////////////////////////////////////////
#ifdef STARPU_USE_MPI
    static void transferInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInoutPassCallbackMpi< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                idxLevel, outsideInteractions->data(), outsideInteractions->size(), kernel,
                starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void transferInPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                idxLevel, kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }

    static void transferInoutPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInoutPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[3]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[4]),
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[5]),
                    idxLevel, outsideInteractions->data(), int(outsideInteractions->size()), kernel,
                    starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Downard Pass
    /////////////////////////////////////////////////////////////////////////////////////
    static void downardPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel, &intervalSize);

        FCudaParams<unsigned char*,9> subCellGroupsPtr;
        memset(&subCellGroupsPtr, 0, sizeof(subCellGroupsPtr));
        FCudaParams<std::size_t,9> subCellGroupsSize;
        memset(&subCellGroupsPtr, 0, sizeof(subCellGroupsSize));
        FCudaParams<unsigned char*,9> subCellGroupsDownPtr;
        memset(&subCellGroupsDownPtr, 0, sizeof(subCellGroupsDownPtr));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroupsPtr.values[idxSubGroup] = ((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+2]));
            subCellGroupsSize.values[idxSubGroup] = (STARPU_VARIABLE_GET_ELEMSIZE(buffers[(idxSubGroup*2)+2]));
            subCellGroupsDownPtr.values[idxSubGroup] = ((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+3]));
        }

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__downardPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                subCellGroupsPtr,subCellGroupsSize,subCellGroupsDownPtr,
                nbSubCellGroups, idxLevel, kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass MPI
    /////////////////////////////////////////////////////////////////////////////////////

#ifdef STARPU_USE_MPI
    static void directInoutPassCallbackMpi(void *buffers[], void *cl_arg){

        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInoutPassCallbackMpi< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                outsideInteractions->data(), outsideInteractions->size(),
                worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight ,kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void directInPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize);
        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                    (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight, kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }

    static void directInoutPassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInoutPassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                outsideInteractions->data(), int(outsideInteractions->size()), worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight,
                kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }


    /////////////////////////////////////////////////////////////////////////////////////
    /// Merge Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void mergePassCallback(void *buffers[], void *cl_arg){
        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__mergePassCallback< SymboleCellClass, PoleCellClass, LocalCellClass,
                CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>(
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                kernel, starpu_cuda_get_local_stream(),
                FCuda__GetGridSize(kernel,intervalSize),FCuda__GetBlockSize(kernel));
    }
};


#endif // FSTARPUCUDAWRAPPER_HPP

