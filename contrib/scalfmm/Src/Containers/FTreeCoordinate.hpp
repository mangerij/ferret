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
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP

#include <string>

#include "../Utils/FGlobal.hpp"
#include "../Utils/FMath.hpp"

#include "../Components/FAbstractSerializable.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTreeCoordinate
* Please read the license
*
* This class represents tree coordinate. It is used to save
* the position in "box unit" (not system/space unit!).
* It is directly related to morton index, as interleaves
* bits from this coordinate make the morton index
*/
class FTreeCoordinate : public FAbstractSerializable {
private:
    int data[3];	//< all box-th position

public:	
    /** Default constructor (position = {0,0,0})*/
    FTreeCoordinate() {
        data[0] = data[1] = data[2] = 0;
    }

    /** Default constructor (position = {0,0,0})*/
    explicit FTreeCoordinate(const MortonIndex mindex, const int indexLevel) {
        setPositionFromMorton(mindex, indexLevel);
    }

    /**
        * Default constructor
        * @param inX the x
        * @param inY the y
        * @param inZ the z
        */
    explicit FTreeCoordinate(const int inX,const int inY,const int inZ) {
        data[0] = inX;
        data[1] = inY;
        data[2] = inZ;
    }

    explicit FTreeCoordinate(const int inPosition[3]) {
        data[0] = inPosition[0];
        data[1] = inPosition[1];
        data[2] = inPosition[2];
    }

    /**
    * Copy constructor
    * @param other the source class to copy
    */
    FTreeCoordinate(const FTreeCoordinate& other) {
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
    }

    /**
        * Copy constructor
        * @param other the source class to copy
        */
    FTreeCoordinate(const FTreeCoordinate& other, const int inOffset) {
        data[0] = other.data[0] + inOffset;
        data[1] = other.data[1] + inOffset;
        data[2] = other.data[2] + inOffset;
    }

    /**
    * Copy constructor
    * @param other the source class to copy
    * @return this a reference to the current object
    */
    FTreeCoordinate& operator=(const FTreeCoordinate& other){
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
        return *this;
    }

    /**
    * Position setter
        * @param inX the new x
        * @param inY the new y
        * @param inZ the new z
    */
    void setPosition(const int inX,const int inY,const int inZ){
        data[0] = inX;
        data[1] = inY;
        data[2] = inZ;
    }

    /**
    * X Getter
        * @return data[0]
    */
    int getX() const{
        return data[0];
    }

    /**
    * Y Getter
        * @return data[1]
    */
    int getY() const{
        return data[1];
    }

    /**
    * Z Getter
        * @return data[2]
    */
    int getZ() const{
        return data[2];
    }

    /**
    * X Setter, simply change x position
    * @param the new x
    */
    void setX(const int inX){
        data[0] = inX;
    }

    /**
    * Y Setter, simply change y position
    * @param the new y
    */
    void setY(const int inY){
        data[1] = inY;
    }

    /**
    * Z Setter, simply change z position
    * @param the new z
    */
    void setZ(const int inZ){
        data[2] = inZ;
    }

    /**
    * To get the morton index of the current position
    * @complexity inLevel
    * @param inLevel the level of the component
    * @return morton index
    */
    MortonIndex getMortonIndex(const int inLevel) const{
        MortonIndex index = 0x0LL;
        MortonIndex mask = 0x1LL;
        // the ordre is xyz.xyz...
        MortonIndex mx = data[0] << 2;
        MortonIndex my = data[1] << 1;
        MortonIndex mz = data[2];

        for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
            index |= (mz & mask);
            mask <<= 1;
            index |= (my & mask);
            mask <<= 1;
            index |= (mx & mask);
            mask <<= 1;

            mz <<= 2;
            my <<= 2;
            mx <<= 2;
        }

        return index;
    }

    /** This function set the position of the current object using a morton index
          * @param inIndex the morton index to compute position
          * @param the level of the morton index
          */
    void setPositionFromMorton(MortonIndex inIndex, const int inLevel){
        MortonIndex mask = 0x1LL;

        data[0] = 0;
        data[1] = 0;
        data[2] = 0;

        for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
            data[2] |= int(inIndex & mask);
            inIndex >>= 1;
            data[1] |= int(inIndex & mask);
            inIndex >>= 1;
            data[0] |= int(inIndex & mask);

            mask <<= 1;
        }

    }

    /** Test equal operator
          * @param other the coordinate to compare
          * @return true if other & current object have same position
          */
    bool operator==(const FTreeCoordinate& other) const {
        return data[0] == other.data[0] && data[1] == other.data[1] && data[2] == other.data[2];
    }

    /** Test equal operator
          * @param other the coordinate to compare
          * @return true if other & current object have same position
          */
    bool equals(const int inX, const int inY, const int inZ) const {
        return data[0] == inX && data[1] == inY && data[2] == inZ;
    }

    /** To test difference
      *
      */
    bool operator!=(const FTreeCoordinate& other) const{
        return data[0] != other.data[0] || data[1] != other.data[1] || data[2] != other.data[2];
    }

    /**
     * Operator stream FTreeCoordinate to std::ostream
     * This can be used to simply write out a tree coordinate
     * @param[in,out] output where to write the coordinate
     * @param[in] inCoordinate the coordinate to write out
     * @return the output for multiple << operators
     */
    template <class StreamClass>
    friend StreamClass& operator<<(StreamClass& output, const FTreeCoordinate& inCoordinate){
        output << "(" <<  inCoordinate.getX() << ", " << inCoordinate.getY() << ", " << inCoordinate.getZ() <<")";
        return output;  // for multiple << operators.
    }

    /** Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        buffer << data[0] << data[1] << data[2];
    }
    /** Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        buffer >> data[0] >> data[1] >> data[2];
    }

    /** To know the size when we save it */
    FSize getSavedSize() const {
        return FSize(sizeof(data[0]) + sizeof(data[1]) + sizeof(data[2]));
    }


    static std::string MortonToBinary(MortonIndex index, int level){
        std::string str;
        int bits = 1 << ((level * 3) - 1);
        int dim = 0;
        while(bits){
            if(index & bits) str.append("1");
            else str.append("0");
            bits >>= 1;
            // we put a dot each 3 values
            if(++dim == 3){
                str.append(".");
                dim = 0;
            }
        }
        return str;
    }

    /* @brief Compute the index of the cells in neighborhood of a given cell
   * @param OtreeHeight Height of the Octree
   * @param indexes target array to store the MortonIndexes computed
   * @param indexInArray store
   */
    int getNeighborsIndexes(const int OctreeHeight, MortonIndex indexes[26], int indexInArray[26]) const{
        int idxNeig = 0;
        int limite = 1 << (OctreeHeight - 1);
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(this->getX() + idxX,0, limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(this->getY() + idxY,0, limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(this->getZ() + idxZ,0, limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(this->getX() + idxX, this->getY() + idxY, this->getZ() + idxZ);
                        indexes[ idxNeig ] = other.getMortonIndex(OctreeHeight - 1);
                        indexInArray[ idxNeig ] = ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }


    /* @brief Compute the indexes of the neighborhood of the calling cell
   * @param OtreeHeight Height of the Octree
   * @param indexes target array to store the MortonIndexes computed
   */
    int getNeighborsIndexes(const int OctreeHeight, MortonIndex indexes[26]) const{
        int idxNeig = 0;
        int limite = 1 << (OctreeHeight - 1);
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(this->getX() + idxX,0, limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(this->getY() + idxY,0, limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(this->getZ() + idxZ,0, limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(this->getX() + idxX, this->getY() + idxY, this->getZ() + idxZ);
                        indexes[ idxNeig ] = other.getMortonIndex(OctreeHeight - 1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }

    int getInteractionNeighbors(const int inLevel, MortonIndex inNeighbors[/*189+26+1*/216], int inNeighborsPosition[/*189+26+1*/216],
                            const int neighSeparation = 1) const{
        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(this->getX()>>1,this->getY()>>1,this->getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int limite = FMath::pow2(inLevel-1);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(parentCell.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(parentCell.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(parentCell.getZ() + idxZ,0,limite)) continue;

                    // if we are not on the current cell
                    if(neighSeparation<1 || idxX || idxY || idxZ ){
                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        const MortonIndex mortonOther = otherParent.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const int xdiff  = ((otherParent.getX()<<1) | ( (idxCousin>>2) & 1)) - this->getX();
                            const int ydiff  = ((otherParent.getY()<<1) | ( (idxCousin>>1) & 1)) - this->getY();
                            const int zdiff  = ((otherParent.getZ()<<1) | (idxCousin&1)) - this->getZ();

                            // Test if it is a direct neighbor
                            if(FMath::Abs(xdiff) > neighSeparation || FMath::Abs(ydiff) > neighSeparation || FMath::Abs(zdiff) > neighSeparation){
                                // add to neighbors
                                inNeighborsPosition[idxNeighbors] = ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3;
                                inNeighbors[idxNeighbors++] = (mortonOther << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

    int getInteractionNeighbors(const int inLevel, MortonIndex inNeighbors[/*189+26+1*/216], const int neighSeparation = 1) const{
        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(this->getX()>>1,this->getY()>>1,this->getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int limite = FMath::pow2(inLevel-1);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(parentCell.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(parentCell.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(parentCell.getZ() + idxZ,0,limite)) continue;

                    // if we are not on the current cell
                    if(neighSeparation<1 || idxX || idxY || idxZ ){
                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        const MortonIndex mortonOther = otherParent.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const int xdiff  = ((otherParent.getX()<<1) | ( (idxCousin>>2) & 1)) - this->getX();
                            const int ydiff  = ((otherParent.getY()<<1) | ( (idxCousin>>1) & 1)) - this->getY();
                            const int zdiff  = ((otherParent.getZ()<<1) | (idxCousin&1)) - this->getZ();

                            // Test if it is a direct neighbor
                            if(FMath::Abs(xdiff) > neighSeparation || FMath::Abs(ydiff) > neighSeparation || FMath::Abs(zdiff) > neighSeparation){
                                // add to neighbors
                                inNeighbors[idxNeighbors++] = (mortonOther << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

};



#endif //FTREECOORDINATE_HPP


