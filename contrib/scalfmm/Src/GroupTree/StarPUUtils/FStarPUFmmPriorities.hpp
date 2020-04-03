#ifndef FSTARPUFMMPRIORITIES_HPP
#define FSTARPUFMMPRIORITIES_HPP

#include "../../Utils/FGlobal.hpp"
#include "FStarPUUtils.hpp"

#include "FStarPUKernelCapacities.hpp"

/**
 * @brief The FStarPUFmmPriorities class
 * This class should have an static method to be called by hetero getPrio.
 */
#ifdef STARPU_SUPPORT_SCHEDULER

#include "FStarPUHeteoprio.hpp"

class FStarPUFmmPriorities{
    int insertionPositionP2M;
    int insertionPositionM2M;

    int insertionPositionP2MSend;
    int insertionPositionM2MSend;

    int insertionPositionM2L;
    int insertionPositionM2LExtern;
    int insertionPositionM2LLastLevel;
    int insertionPositionL2L;
    int insertionPositionL2P;
    int insertionPositionP2P;
    int insertionPositionP2PExtern;

    int treeHeight;

    FStarPUKernelCapacities* capacities;

    static FStarPUFmmPriorities controller;


    FStarPUFmmPriorities(){
    }

public:
    static FStarPUFmmPriorities& Controller(){
        return controller;
    }

    static void InitSchedulerCallback(unsigned sched_ctx_id, void* heteroprio){
        Controller().initSchedulerCallback(sched_ctx_id, (struct _starpu_heteroprio_center_policy_heteroprio*)heteroprio);
    }

    void init(struct starpu_conf* conf, const int inTreeHeight,
              FStarPUKernelCapacities* inCapacities){
        capacities = inCapacities;

        conf->sched_policy = &_starpu_sched_heteroprio_policy;
        starpu_heteroprio_set_callback(&InitSchedulerCallback);

        treeHeight  = inTreeHeight;

        if(inTreeHeight > 2){
            int incPrio = 0;

            FLOG( FLog::Controller << "Buckets:\n" );

            insertionPositionP2MSend = incPrio++;
            FLOG( FLog::Controller << "\t P2M Send "  << insertionPositionP2MSend << "\n" );
            insertionPositionP2M     = incPrio++;
            FLOG( FLog::Controller << "\t P2M "  << insertionPositionP2M << "\n" );

            insertionPositionM2MSend = incPrio++;
            FLOG( FLog::Controller << "\t M2M Send "  << insertionPositionM2MSend << "\n" );
            insertionPositionM2M     = incPrio++;
            FLOG( FLog::Controller << "\t M2M "  << insertionPositionM2M << "\n" );

            insertionPositionM2L     = incPrio++;
            FLOG( FLog::Controller << "\t M2L "  << insertionPositionM2L << "\n" );
            insertionPositionM2LExtern = incPrio++;
            FLOG( FLog::Controller << "\t M2L Outer "  << insertionPositionM2LExtern << "\n" );

            insertionPositionL2L     = incPrio++;
            FLOG( FLog::Controller << "\t L2L "  << insertionPositionL2L << "\n" );

            incPrio += (treeHeight-3) - 1;   // M2L is done treeHeight-2 times
            incPrio += (treeHeight-3) - 1;   // M2L is done treeHeight-2 times
            incPrio += (treeHeight-3) - 1;   // L2L is done treeHeight-3 times

            insertionPositionP2P       = incPrio++;
            FLOG( FLog::Controller << "\t P2P "  << insertionPositionP2P << "\n" );

            insertionPositionM2LLastLevel = incPrio++;
            FLOG( FLog::Controller << "\t M2L last "  << insertionPositionM2LLastLevel << "\n" );

            insertionPositionL2P     = incPrio++;
            FLOG( FLog::Controller << "\t L2P "  << insertionPositionL2P << "\n" );

            insertionPositionP2PExtern = incPrio++;
            FLOG( FLog::Controller << "\t P2P Outer "  << insertionPositionP2PExtern << "\n" );

            assert(incPrio == 8 + (treeHeight-3) + (treeHeight-3) + (treeHeight-3));
        }
        else{
            int incPrio = 0;

            insertionPositionP2MSend = -1;
            insertionPositionP2M     = -1;

            insertionPositionM2MSend = -1;
            insertionPositionM2M     = -1;

            insertionPositionM2L     = -1;
            insertionPositionM2LExtern = -1;
            insertionPositionM2LLastLevel = -1;

            insertionPositionL2L     = -1;

            insertionPositionP2P     = incPrio++;
            insertionPositionP2PExtern = insertionPositionP2P;

            insertionPositionL2P     = -1;
            assert(incPrio == 1);
        }
    }

    void initSchedulerCallback(unsigned /*sched_ctx_id*/,
                               struct _starpu_heteroprio_center_policy_heteroprio *heteroprio){
        const bool workOnlyOnLeaves = (treeHeight <= 2);
#ifdef STARPU_USE_CPU
        // CPU follows the real prio
        {
            int cpuCountPrio = 0;

            if( !workOnlyOnLeaves && capacities->supportP2M(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio P2M Send "  << cpuCountPrio << " bucket " << insertionPositionP2MSend << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionP2MSend;
                heteroprio->buckets[insertionPositionP2MSend].valide_archs |= STARPU_CPU;

                FLOG( FLog::Controller << "\t CPU prio P2M "  << cpuCountPrio << " bucket " << insertionPositionP2M << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionP2M;
                heteroprio->buckets[insertionPositionP2M].valide_archs |= STARPU_CPU;
            }
            if(!workOnlyOnLeaves && capacities->supportM2M(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio M2M Send "  << cpuCountPrio << " bucket " << insertionPositionM2MSend << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionM2MSend;
                heteroprio->buckets[insertionPositionM2MSend].valide_archs |= STARPU_CPU;

                FLOG( FLog::Controller << "\t CPU prio M2M "  << cpuCountPrio << " bucket " << insertionPositionM2M << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionM2M;
                heteroprio->buckets[insertionPositionM2M].valide_archs |= STARPU_CPU;
            }
            for(int idxLevel = 2 ; idxLevel < treeHeight-1 ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_CPU_IDX)){
                    const int prioM2LAtLevel = getInsertionPosM2L(idxLevel);
                    FLOG( FLog::Controller << "\t CPU prio M2L "  << cpuCountPrio << " bucket " << prioM2LAtLevel << "\n" );
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_CPU;

                    const int prioM2LAtLevelExtern = getInsertionPosM2LExtern(idxLevel);
                    FLOG( FLog::Controller << "\t CPU prio M2L extern "  << cpuCountPrio << " bucket " << prioM2LAtLevelExtern << "\n" );
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2LAtLevelExtern;
                    heteroprio->buckets[prioM2LAtLevelExtern].valide_archs |= STARPU_CPU;
                }
                if(capacities->supportL2L(FSTARPU_CPU_IDX)){
                    const int prioL2LAtLevel = getInsertionPosL2L(idxLevel);
                    FLOG( FLog::Controller << "\t CPU prio L2L "  << cpuCountPrio << " bucket " << prioL2LAtLevel << "\n" );
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_CPU;
                }
            }

            if( capacities->supportP2P(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio P2P "  << cpuCountPrio << " bucket " << insertionPositionP2P << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionP2P;
                heteroprio->buckets[insertionPositionP2P].valide_archs |= STARPU_CPU;
            }
#ifndef STARPU_USE_REDUX
            if( capacities->supportP2PExtern(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio P2P Extern "  << cpuCountPrio << " bucket " << insertionPositionP2PExtern << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionP2PExtern;
                heteroprio->buckets[insertionPositionP2PExtern].valide_archs |= STARPU_CPU;
            }
#endif
            if(capacities->supportM2L(FSTARPU_CPU_IDX)){
                const int prioM2LAtLevel = getInsertionPosM2L(treeHeight-1);
                FLOG( FLog::Controller << "\t CPU prio M2L "  << cpuCountPrio << " bucket " << prioM2LAtLevel << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2LAtLevel;
                heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_CPU;
            }
            if( !workOnlyOnLeaves && capacities->supportL2P(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio L2P "  << cpuCountPrio << " bucket " << insertionPositionL2P << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionL2P;
                heteroprio->buckets[insertionPositionL2P].valide_archs |= STARPU_CPU;
            }
#ifdef STARPU_USE_REDUX
            if( capacities->supportP2PExtern(FSTARPU_CPU_IDX)){
                FLOG( FLog::Controller << "\t CPU prio P2P Extern "  << cpuCountPrio << " bucket " << insertionPositionP2PExtern << "\n" );
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = insertionPositionP2PExtern;
                heteroprio->buckets[insertionPositionP2PExtern].valide_archs |= STARPU_CPU;
            }
#endif
            heteroprio->nb_prio_per_arch_index[FSTARPU_CPU_IDX] = unsigned(cpuCountPrio);
            FLOG( FLog::Controller << "\t CPU Priorities: "  << cpuCountPrio << "\n" );
        }
#endif
#ifdef STARPU_USE_OPENCL
        {
            int openclCountPrio = 0;

            //insertionPositionP2P       = insertionPositionL2L + (treeHeight-3)*2+1 +1;
            //insertionPositionP2PExtern = insertionPositionP2P;
            //insertionPositionP2PMpi    = insertionPositionP2P;
            if(capacities->supportP2P(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionP2P;
                heteroprio->buckets[insertionPositionP2P].factor_base_arch_index = FSTARPU_OPENCL_IDX;
                heteroprio->buckets[insertionPositionP2P].valide_archs |= STARPU_OPENCL;
#ifdef STARPU_USE_CPU
                heteroprio->buckets[insertionPositionP2P].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
            }

            // insertionPositionM2L       = insertionPositionM2M+1;
            // insertionPositionM2LExtern = insertionPositionM2L;
            // insertionPositionM2LMpi    = insertionPositionM2L;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_OPENCL_IDX)){
                    const int prioM2LAtLevel = getInsertionPosM2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].factor_base_arch_index = FSTARPU_OPENCL_IDX;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_OPENCL;
#ifdef STARPU_USE_CPU
                    heteroprio->buckets[prioM2LAtLevel].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
                }
            }

            //insertionPositionP2MSend = 0;
            //insertionPositionP2M     = insertionPositionP2MSend+1;
            if( !workOnlyOnLeaves && capacities->supportP2M(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionP2MSend;
                heteroprio->buckets[insertionPositionP2MSend].valide_archs |= STARPU_OPENCL;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionP2M;
                heteroprio->buckets[insertionPositionP2M].valide_archs |= STARPU_OPENCL;
            }

            //insertionPositionM2MSend = insertionPositionP2M+1;
            //insertionPositionM2M     = insertionPositionM2MSend+1;
            if( !workOnlyOnLeaves && capacities->supportM2M(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionM2MSend;
                heteroprio->buckets[insertionPositionM2MSend].valide_archs |= STARPU_OPENCL;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionM2M;
                heteroprio->buckets[insertionPositionM2M].valide_archs |= STARPU_OPENCL;
            }

            // insertionPositionL2L     = insertionPositionM2L+1;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(idxLevel != treeHeight-1 && capacities->supportL2L(FSTARPU_OPENCL_IDX)){
                    const int prioL2LAtLevel = getInsertionPosL2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_OPENCL;
                }
            }

            //insertionPositionL2P     = insertionPositionP2PMpi+1;
            if( !workOnlyOnLeaves && capacities->supportL2P(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = insertionPositionL2P;
                heteroprio->buckets[insertionPositionL2P].valide_archs |= STARPU_OPENCL;
            }

            heteroprio->nb_prio_per_arch_index[FSTARPU_OPENCL_IDX] = unsigned(openclCountPrio);
        }
#endif
#ifdef STARPU_USE_CUDA
        {
            int openclCountPrio = 0;

            //insertionPositionP2P       = insertionPositionL2L + (treeHeight-3)*2+1 +1;
            //insertionPositionP2PExtern = insertionPositionP2P;
            //insertionPositionP2PMpi    = insertionPositionP2P;
            if(capacities->supportP2P(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionP2P;
                heteroprio->buckets[insertionPositionP2P].valide_archs |= STARPU_CUDA;
                heteroprio->buckets[insertionPositionP2P].factor_base_arch_index = FSTARPU_CUDA_IDX;
#ifdef STARPU_USE_CPU
                heteroprio->buckets[insertionPositionP2P].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
            }

            // insertionPositionM2L       = insertionPositionM2M+1;
            // insertionPositionM2LExtern = insertionPositionM2L;
            // insertionPositionM2LMpi    = insertionPositionM2L;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_CUDA_IDX)){
                    const int prioM2LAtLevel = getInsertionPosM2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_CUDA;
                    heteroprio->buckets[prioM2LAtLevel].factor_base_arch_index = FSTARPU_CUDA_IDX;
#ifdef STARPU_USE_CPU
                    heteroprio->buckets[prioM2LAtLevel].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
                }
            }

            //insertionPositionP2MSend = 0;
            //insertionPositionP2M     = insertionPositionP2MSend+1;
            if( !workOnlyOnLeaves && capacities->supportP2M(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionP2MSend;
                heteroprio->buckets[insertionPositionP2MSend].valide_archs |= STARPU_CUDA;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionP2M;
                heteroprio->buckets[insertionPositionP2M].valide_archs |= STARPU_CUDA;
            }

            //insertionPositionM2MSend = insertionPositionP2M+1;
            //insertionPositionM2M     = insertionPositionM2MSend+1;
            if( !workOnlyOnLeaves && capacities->supportM2M(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionM2MSend;
                heteroprio->buckets[insertionPositionM2MSend].valide_archs |= STARPU_CUDA;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionM2M;
                heteroprio->buckets[insertionPositionM2M].valide_archs |= STARPU_CUDA;
            }

            // insertionPositionL2L     = insertionPositionM2L+1;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(idxLevel != treeHeight-1 && capacities->supportL2L(FSTARPU_CUDA_IDX)){
                    const int prioL2LAtLevel = getInsertionPosL2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_CUDA;
                }
            }

            //insertionPositionL2P     = insertionPositionP2PMpi+1;
            if( !workOnlyOnLeaves && capacities->supportL2P(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = insertionPositionL2P;
                heteroprio->buckets[insertionPositionL2P].valide_archs |= STARPU_CUDA;
            }

            heteroprio->nb_prio_per_arch_index[FSTARPU_CUDA_IDX] = int(openclCountPrio);
        }
#endif
    }

    int getInsertionPosP2M() const {
        return insertionPositionP2M;
    }
    int getInsertionPosM2M(const int /*inLevel*/) const {
        return insertionPositionM2M;
    }
    int getInsertionPosP2M(bool willBeSend) const {
        return willBeSend?insertionPositionP2MSend:insertionPositionP2M;
    }
    int getInsertionPosM2M(const int /*inLevel*/, bool willBeSend) const {
        return willBeSend?insertionPositionM2MSend:insertionPositionM2M;
    }
    int getInsertionPosM2L(const int inLevel) const {
        return (inLevel==treeHeight-1? insertionPositionM2LLastLevel : insertionPositionM2L + (inLevel - 2)*3);
    }
    int getInsertionPosM2LExtern(const int inLevel) const {
        return (inLevel==treeHeight-1? insertionPositionM2LLastLevel : insertionPositionM2LExtern + (inLevel - 2)*3);
    }
    int getInsertionPosL2L(const int inLevel) const {
        return insertionPositionL2L + (inLevel - 2)*3;
    }
    int getInsertionPosL2P() const {
        return insertionPositionL2P;
    }
    int getInsertionPosP2P() const {
        return insertionPositionP2P;
    }
    int getInsertionPosP2PExtern() const {
        return insertionPositionP2PExtern;
    }
};

FStarPUFmmPriorities FStarPUFmmPriorities::controller;

#elif defined(SCALFMM_STARPU_USE_PRIO)// STARPU_SUPPORT_SCHEDULER

#include "FOmpPriorities.hpp"

class FStarPUFmmPriorities {
    static FStarPUFmmPriorities controller;
    FOmpPriorities ompPrio;

    FStarPUFmmPriorities() : ompPrio(0){
    }

public:
    static FStarPUFmmPriorities& Controller(){
        return controller;
    }


    void init(struct starpu_conf* /*conf*/, const int inTreeHeight,
              FStarPUKernelCapacities* /*inCapacities*/){
        ompPrio = FOmpPriorities(inTreeHeight);
    }

    int getInsertionPosP2M() const {
        return ompPrio.getInsertionPosP2M();
    }
    int getInsertionPosM2M(const int inLevel) const {
        return ompPrio.getInsertionPosM2M(inLevel);
    }
    int getInsertionPosP2M(bool /*willBeSend*/) const {
        return ompPrio.getInsertionPosP2M();
    }
    int getInsertionPosM2M(const int inLevel, bool /*willBeSend*/) const {
        return ompPrio.getInsertionPosM2M(inLevel);
    }
    int getInsertionPosM2L(const int inLevel) const {
        return ompPrio.getInsertionPosM2L(inLevel);
    }
    int getInsertionPosM2LExtern(const int inLevel) const {
        return ompPrio.getInsertionPosM2LExtern(inLevel);
    }
    int getInsertionPosL2L(const int inLevel) const {
        return ompPrio.getInsertionPosL2L(inLevel);
    }
    int getInsertionPosL2P() const {
        return ompPrio.getInsertionPosL2P();
    }
    int getInsertionPosP2P() const {
        return ompPrio.getInsertionPosP2P();
    }
    int getInsertionPosP2PExtern() const {
        return ompPrio.getInsertionPosP2PExtern();
    }
};

FStarPUFmmPriorities FStarPUFmmPriorities::controller;

#endif // SCALFMM_STARPU_USE_PRIO - STARPU_SUPPORT_SCHEDULER



#endif // FSTARPUFMMPRIORITIES_HPP
