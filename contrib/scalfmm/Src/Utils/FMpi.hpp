// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#ifndef FMPI_HPP
#define FMPI_HPP


#include <cstdio>
#include <stdexcept>

#include "FGlobal.hpp"
#ifndef SCALFMM_USE_MPI
#error The MPI header is included while SCALFMM_USE_MPI is turned OFF
#endif


#include "FNoCopyable.hpp"
#include "FMath.hpp"

//Need that for converting datas
#include "FComplex.hpp"


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


#include <mpi.h>


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief MPI context management
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 *
 * This class manages the initialization and destruction of the MPI context. It
 * also gives access to the main communicator via its `const FComm& global()`
 * member.
 *
 * @warning Do not create more that one instance of this class for the whole
 * program lifetime. MPI does not support multiple initializations.
 *
 * Please read the license.
 */

class FMpi {
public:
    ////////////////////////////////////////////////////////
    // MPI Flag
    ////////////////////////////////////////////////////////
    enum FMpiTag {
        // FMpiTreeBuilder
        TagExchangeIndexs = 100,
        TagSplittedLeaf = 200,
        TagExchangeNbLeafs = 300,
        TagSandSettling = 400,

        // FQuickSort
        TagQuickSort = 500,

        // FMM
        TagFmmM2M = 1000,
        TagFmmM2MSize = 1500,
        TagFmmL2L = 2000,
        TagFmmL2LSize = 2500,
        TagFmmP2P = 3000,

        // Bitonic,
        TagBitonicMin = 4000,
        TagBitonicMax = 5000,
        TagBitonicMinMess = 6000,
        TagBitonicMaxMess = 7000,

        // Last defined tag
        TagLast = 8000,
    };

    ////////////////////////////////////////////////////////
    // FComm to factorize MPI_Comm work
    ////////////////////////////////////////////////////////
    
    /**
     * \brief MPI comunicator abstraction
     *
     * This class is used to gather the usual methods related to identifying an
     * MPI communicator.
     */
    class FComm : public FNoCopyable {
        int rank;   ///< rank related to the comm
        int nbProc; ///< nb proc in this group

        MPI_Comm communicator;  ///< current mpi communicator
        MPI_Group group;        ///< current mpi group


        /// Updates current process rank and process count from mpi
        void reset(){
            FMpi::Assert( MPI_Comm_rank(communicator,&rank),  __LINE__ );
            FMpi::Assert( MPI_Comm_size(communicator,&nbProc),  __LINE__ );
        }

    public:
        /// Constructor : duplicates the given communicator
        explicit FComm(MPI_Comm inCommunicator ) {
            FMpi::Assert( MPI_Comm_dup(inCommunicator, &communicator),  __LINE__ , "comm dup");
            FMpi::Assert( MPI_Comm_group(communicator, &group),  __LINE__ , "comm group");

            reset();
        }

        /// Frees communicator and group
        virtual ~FComm(){
            FMpi::Assert( MPI_Comm_free(&communicator),  __LINE__ );
            FMpi::Assert( MPI_Group_free(&group),  __LINE__ );
        }

        /// Gets the underlying MPI communicator
        MPI_Comm getComm() const {
            return communicator;
        }

        /// Gets the current rank
        int processId() const {
            return rank;
        }

        /// The current number of procs in the group */
        int processCount() const {
            return nbProc;
        }

        ////////////////////////////////////////////////////////////
        // Split/Chunk functions
        ////////////////////////////////////////////////////////////

        /** Get a left index related to a size */
        template< class T >
        T getLeft(const T inSize)  const {
            const double step = (double(inSize) / double(processCount()));
            return T(FMath::Ceil(step * double(processId())));
        }

        /** Get a right index related to a size */
        template< class T >
        T getRight(const T inSize)  const {
            const double step = (double(inSize) / double(processCount()));
            const T res = T(FMath::Ceil(step * double(processId()+1)));
            if(res > inSize) return inSize;
            else return res;
        }

        /** Get a right index related to a size and another id */
        template< class T >
        T getOtherRight(const T inSize, const int other)  const {
            const double step = (double(inSize) / double(processCount()));
            const T res = T(FMath::Ceil(step * double(other+1)));
            if(res > inSize) return inSize;
            else return res;
        }

        /** Get a left index related to a size and another id */
        template< class T >
        T getOtherLeft(const T inSize, const int other) const {
            const double step = (double(inSize) / double(processCount()));
            return T(FMath::Ceil(step * double(other)));
        }

        /** Get a proc id from and index */
        template< class T >
        int getProc(const int position, const T inSize) const {
            const double step = (double(inSize) / processCount());
            return int(position/step);
        }

        ////////////////////////////////////////////////////////////
        // Mpi interface functions
        ////////////////////////////////////////////////////////////


        /** Reduce a value for proc == 0 */
        template< class T >
        T reduceSum(T data) const {
            T result(0);
            FMpi::Assert( MPI_Reduce( &data, &result, 1, FMpi::GetType(data), MPI_SUM, 0, communicator ), __LINE__);
            return result;
        }

        /** Reduce a value for all procs */
        template< class T >
        T allReduceSum(T data) const {
            T result(0);
            FMpi::Assert( MPI_Allreduce( &data, &result, 1, FMpi::GetType(data), MPI_SUM, communicator ), __LINE__);
            return result;
        }

        /** Reduce an average */
        template< class T >
        T reduceAverageAll(T data) const {
            T result[processCount()];
            FMpi::Assert( MPI_Allgather( &data, 1, FMpi::GetType(data), result, 1, FMpi::GetType(data), getComm()),  __LINE__ );

            T average = 0;
            for(int idxProc = 0 ; idxProc < processCount() ;++idxProc){
                average += result[idxProc] / processCount();
            }
            return average;
        }

        /** Change the group size */
        void groupReduce(const int from , const int to){
            int * procsIdArray = new int [to - from + 1];
            for(int idxProc = from ;idxProc <= to ; ++idxProc){
                procsIdArray[idxProc - from] = idxProc;
            }

            MPI_Group previousGroup = group;
            FMpi::Assert( MPI_Group_incl(previousGroup, to - from + 1 , procsIdArray, &group),  __LINE__ );

            MPI_Comm previousComm = communicator;
            FMpi::Assert( MPI_Comm_create(previousComm, group, &communicator),  __LINE__ );

            MPI_Comm_free(&previousComm);
            MPI_Group_free(&previousGroup);

            reset();
            delete[]  procsIdArray ;
        }

        void barrier() const {
            FMpi::Assert(MPI_Barrier(getComm()), __LINE__);
        }
    };

    ////////////////////////////////////////////////////////
    // FMpi methods
    ////////////////////////////////////////////////////////

    /// Constructor
    /**
     * We use `MPI_Init_thread` because of an openmpi error:
     *
     *     [fourmi062:15896] [[13237,0],1]-[[13237,1],1] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)
     *     [fourmi056:04597] [[13237,0],3]-[[13237,1],3] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)
     *     [fourmi053:08571] [[13237,0],5]-[[13237,1],5] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)
     *    
     * Error for process 1:
     *
     *     [[13237,1],1][btl_openib_component.c:3227:handle_wc] from fourmi062 to: fourmi056 error polling LP CQ with status LOCAL LENGTH ERROR status number 1 for wr_id 7134664 opcode 0  vendor error 105 qp_idx 3
     *                ^
     *
     * All processes raise the same error, the second 1 (see caret) is replaced by the rank.
     */
    FMpi() : communicator(nullptr) {
        if( instanceCount > 0) {
            throw std::logic_error("FMpi should not be instanciated more than once.");
        } else {
            instanceCount++;
        }

        int provided = 0;
        FMpi::Assert( MPI_Init_thread(nullptr,nullptr, MPI_THREAD_SERIALIZED, &provided), __LINE__);
        communicator = new FComm(MPI_COMM_WORLD);
    }

    /// Constructor
    FMpi(int inArgc, char **  inArgv ) : communicator(nullptr) {
        if( instanceCount > 0) {
            throw std::logic_error("FMpi should not be instanciatedmore than once.");
        } else {
            instanceCount++;
        }

        int provided = 0;
        FMpi::Assert( MPI_Init_thread(&inArgc,&inArgv, MPI_THREAD_SERIALIZED, &provided), __LINE__);
        communicator = new FComm(MPI_COMM_WORLD);
    }

    /// Delete the communicator and call `MPI_Finalize`
    ~FMpi(){
        delete communicator;
        MPI_Finalize();
    }

    /// Get the global communicator
    const FComm& global() {
        return (*communicator);
    }

    ////////////////////////////////////////////////////////////
    // Mpi Types meta function
    ////////////////////////////////////////////////////////////

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const long long&){
        return MPI_LONG_LONG;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const long int&){
        return MPI_LONG;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const double&){
        return MPI_DOUBLE;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const float&){
        return MPI_FLOAT;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const int&){
        return MPI_INT;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const char&){
        return MPI_CHAR;
    }

    /// Get the MPI datatype corresponding to a variable.
    static MPI_Datatype GetType(const unsigned char&){
        return MPI_UNSIGNED_CHAR;
    }

    /// Get the MPI datatype corresponding to a variable.
    template <class FReal>
    static MPI_Datatype GetType(const FComplex<FReal>& a){
        MPI_Datatype FMpiComplexe;
        MPI_Type_contiguous(2, GetType(a.getReal()) , &FMpiComplexe);
        return FMpiComplexe;
    }

    ////////////////////////////////////////////////////////////
    // Mpi interface functions
    ////////////////////////////////////////////////////////////

    /// Generic mpi assert function
    static void Assert(const int test, const unsigned line, const char* const message = nullptr){
        if(test != MPI_SUCCESS){
            printf("[ERROR-QS] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

    /// Compute a left index from data
    template <class T>
    static T GetLeft(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return T(ceil(step * inIdProc));
    }

    /// Compute a right index from data
    template <class T>
    static T GetRight(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const T res = T(ceil(step * (inIdProc+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    /// Compute a proc id from index & data */
    template <class T>
    static int GetProc(const T position, const T inSize, const int inNbProc) {
        const double step = double(inSize) / double(inNbProc);
        return int(double(position)/step);
    }

    /// assert if mpi error */
    static void MpiAssert(const int test, const unsigned line, const char* const message = nullptr){
        if(test != MPI_SUCCESS){
            printf("[ERROR] Test failed at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

private:
    /// The original communicator
    FComm* communicator;
    
    /// Counter to avoid several instanciations
    static int instanceCount;
};


#endif //FMPI_HPP


