
// Keep in private GIT
#ifndef FGROUPOFCELLS_HPP
#define FGROUPOFCELLS_HPP

#include "../../Utils/FAssert.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

#include <list>
#include <functional>

/**
* @brief The FGroupOfCells class manages the cells in block allocation.
*/
template <class CompositeCellClass, class SymboleCellClass, class PoleCellClass, class LocalCellClass>
class FGroupOfCells {
    /** One header is allocated at the beginning of each block */
    struct alignas(FStarPUDefaultAlign::StructAlign) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfCellsInBlock;
    };

protected:
    //< The size of the memoryBuffer
    size_t allocatedMemoryInByte;
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    MortonIndex*    cellIndexes;
    //< Pointer to the cells inside the block memory
    SymboleCellClass*      blockCells;

    //< The multipole data
    PoleCellClass* cellMultipoles;
    //< The local data
    LocalCellClass* cellLocals;

    //< To kown if the object has to delete the memory
    bool deleteBuffer;

public:
    typedef CompositeCellClass CompleteCellClass;

    FGroupOfCells()
        : allocatedMemoryInByte(0), memoryBuffer(nullptr),
          blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(false){
    }

    void reset(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
               unsigned char* inCellMultipoles, unsigned char* inCellLocals){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
                (&blockCells[idxCellPtr])->~SymboleCellClass();
                (&cellMultipoles[idxCellPtr])->~PoleCellClass();
                (&cellLocals[idxCellPtr])->~LocalCellClass();
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
        }
        // Move the pointers to the correct position
        allocatedMemoryInByte = (inAllocatedMemoryInByte);
        memoryBuffer = (inBuffer);
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(inBuffer);
        inBuffer += (blockHeader->numberOfCellsInBlock*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymboleCellClass*>(inBuffer);
        inBuffer += (sizeof(SymboleCellClass)*blockHeader->numberOfCellsInBlock);
        FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);

        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
        deleteBuffer = (false);
    }

    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    FGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                  unsigned char* inCellMultipoles, unsigned char* inCellLocals)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(false){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(inBuffer);
        inBuffer += (blockHeader->numberOfCellsInBlock*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymboleCellClass*>(inBuffer);
        inBuffer += (sizeof(SymboleCellClass)*blockHeader->numberOfCellsInBlock);
        FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);

        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
    }

    /**
 * @brief FGroupOfCells
 * @param inStartingIndex first cell morton index
 * @param inEndingIndex last cell morton index + 1
 * @param inNumberOfCells total number of cells in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfCells(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells)
        : allocatedMemoryInByte(0), memoryBuffer(nullptr), blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(true){
        FAssertLF((inEndingIndex-inStartingIndex) >= MortonIndex(inNumberOfCells));
        // Total number of bytes in the block
        const size_t memoryToAlloc = sizeof(BlockHeader) + (inNumberOfCells*sizeof(MortonIndex))
                + (inNumberOfCells*sizeof(SymboleCellClass));

        // Allocate
        FAssertLF(0 <= int(memoryToAlloc) && int(memoryToAlloc) < std::numeric_limits<int>::max());
        allocatedMemoryInByte = memoryToAlloc;
        memoryBuffer = (unsigned char*)FAlignedMemory::AllocateBytes<32>(memoryToAlloc);
        FAssertLF(memoryBuffer);
        memset(memoryBuffer, 0, memoryToAlloc);

        // Move the pointers to the correct position
        unsigned char* ptrBuff = memoryBuffer;
        blockHeader         = reinterpret_cast<BlockHeader*>(ptrBuff);
        ptrBuff += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(ptrBuff);
        ptrBuff += (inNumberOfCells*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymboleCellClass*>(ptrBuff);
        ptrBuff += (sizeof(SymboleCellClass)*inNumberOfCells);
        FAssertLF(size_t(ptrBuff-memoryBuffer) == allocatedMemoryInByte);

        // Init header
        blockHeader->startingIndex = inStartingIndex;
        blockHeader->endingIndex   = inEndingIndex;
        blockHeader->numberOfCellsInBlock  = inNumberOfCells;

        cellMultipoles = (PoleCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(PoleCellClass));
        cellLocals     = (LocalCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(LocalCellClass));
        for(int idxCell = 0 ; idxCell < inNumberOfCells ; ++idxCell){
            new (&cellMultipoles[idxCell]) PoleCellClass();
            new (&cellLocals[idxCell]) LocalCellClass();
            cellIndexes[idxCell] = -1;
        }
    }

    /** Call the destructor of cells and dealloc block memory */
    ~FGroupOfCells(){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
                (&blockCells[idxCellPtr])->~SymboleCellClass();
                (&cellMultipoles[idxCellPtr])->~PoleCellClass();
                (&cellLocals[idxCellPtr])->~LocalCellClass();
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
        }
    }

    /** Give access to the buffer to send the data */
    const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    size_t getBufferSizeInByte() const {
        return allocatedMemoryInByte;
    }

    /** Give access to the buffer to send the data */
    const PoleCellClass* getRawMultipoleBuffer() const{
        return cellMultipoles;
    }

    /** Give access to the buffer to send the data */
    PoleCellClass* getRawMultipoleBuffer() {
        return cellMultipoles;
    }

    /** The the size of the allocated buffer */
    size_t getMultipoleBufferSizeInByte() const {
        return sizeof(PoleCellClass)*blockHeader->numberOfCellsInBlock;
    }

    /** Give access to the buffer to send the data */
    LocalCellClass* getRawLocalBuffer(){
        return cellLocals;
    }

    /** Give access to the buffer to send the data */
    const LocalCellClass* getRawLocalBuffer() const{
        return cellLocals;
    }

    /** The the size of the allocated buffer */
    size_t getLocalBufferSizeInByte() const {
        return sizeof(LocalCellClass)*blockHeader->numberOfCellsInBlock;
    }

    /** To know if the object will delete the memory block */
    bool getDeleteMemory() const{
        return deleteBuffer;
    }

    /** The index of the fist cell (set from the constructor) */
    MortonIndex getStartingIndex() const {
        return blockHeader->startingIndex;
    }

    /** The index of the last cell + 1 (set from the constructor) */
    MortonIndex getEndingIndex() const {
        return blockHeader->endingIndex;
    }

    /** The number of cell (set from the constructor) */
    int getNumberOfCellsInBlock() const {
        return blockHeader->numberOfCellsInBlock;
    }

    /** The size of the interval endingIndex-startingIndex (set from the constructor) */
    MortonIndex getSizeOfInterval() const {
        return MortonIndex(blockHeader->endingIndex-blockHeader->startingIndex);
    }

    /** Return true if inIndex should be located in the current block */
    bool isInside(const MortonIndex inIndex) const{
        return blockHeader->startingIndex <= inIndex && inIndex < blockHeader->endingIndex;
    }

    /** Return the idx in array of the cell */
    MortonIndex getCellMortonIndex(const int cellPos) const{
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellIndexes[cellPos];
    }

    /** Check if a cell exist (by binary search) and return it index */
    int getCellIndex(const MortonIndex cellIdx) const{
        int idxLeft = 0;
        int idxRight = blockHeader->numberOfCellsInBlock-1;
        while(idxLeft <= idxRight){
            const int idxMiddle = (idxLeft+idxRight)/2;
            if(cellIndexes[idxMiddle] == cellIdx){
                return idxMiddle;
            }
            if(cellIdx < cellIndexes[idxMiddle]){
                idxRight = idxMiddle-1;
            }
            else{
                idxLeft = idxMiddle+1;
            }
        }
        return -1;
    }

    /** Check if a cell exist (by binary search) and return it index */
    int getFistChildIdx(const MortonIndex parentIdx) const{
        int idxLeft = 0;
        int idxRight = blockHeader->numberOfCellsInBlock-1;
        while(idxLeft <= idxRight){
            int idxMiddle = (idxLeft+idxRight)/2;
            if((cellIndexes[idxMiddle]>>3) == parentIdx){
                while(0 < idxMiddle && (cellIndexes[idxMiddle-1]>>3) == parentIdx){
                    idxMiddle -= 1;
                }
                return idxMiddle;
            }
            if(parentIdx < (cellIndexes[idxMiddle]>>3)){
                idxRight = idxMiddle-1;
            }
            else{
                idxLeft = idxMiddle+1;
            }
        }
        return -1;
    }

    /** Return true if inIndex is located in the current block and is not empty */
    bool exists(const MortonIndex inIndex) const {
        return isInside(inIndex) && (getCellIndex(inIndex) != -1);
    }

    /** Return the address of the cell if it exists (or NULL) */
    CompositeCellClass getCompleteCell(const int cellPos){
        FAssertLF(cellMultipoles && cellLocals);
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], &cellLocals[cellPos]);
    }

    /** Return the address of the cell if it exists (or NULL) */
    CompositeCellClass getUpCell(const int cellPos){
        FAssertLF(cellMultipoles);
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], nullptr);
    }

    /** Return the address of the cell if it exists (or NULL) */
    CompositeCellClass getDownCell(const int cellPos){
        FAssertLF(cellLocals);
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return CompositeCellClass(&blockCells[cellPos], nullptr, &cellLocals[cellPos]);
    }

    /** Allocate a new cell by calling its constructor */
    template<typename... CellConstructorParams>
    void newCell(const MortonIndex inIndex, const int id, CellConstructorParams... args){
        FAssertLF(isInside(inIndex));
        FAssertLF(!exists(inIndex));
        FAssertLF(id < blockHeader->numberOfCellsInBlock);
        new((void*)&blockCells[id]) SymboleCellClass(args...);
        cellIndexes[id] = inIndex;
    }

    /** Iterate on each allocated cells */
    template<typename... FunctionParams>
    void forEachCell(std::function<void(CompositeCellClass, FunctionParams...)> function, FunctionParams... args){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
            function(CompositeCellClass(&blockCells[idxCellPtr], &cellMultipoles[idxCellPtr], &cellLocals[idxCellPtr]), args...);
        }
    }

    void forEachCell(std::function<void(CompositeCellClass)> function){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
            function(CompositeCellClass(&blockCells[idxCellPtr], &cellMultipoles[idxCellPtr], &cellLocals[idxCellPtr]));
        }
    }
};


#endif // FGROUPOFCELLS_HPP
