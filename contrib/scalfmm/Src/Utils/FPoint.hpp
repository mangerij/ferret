// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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
//
#ifndef FPOINT_HPP
#define FPOINT_HPP


// To get memcpy
#include <cstring>
#include <iostream>

#include "FGlobal.hpp"
#include "FMath.hpp"

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * Please read the license
 *
 * This class is a 3D vector. It can be used as a position
 * or as a 3d forces vector etc.
 */
template <class FReal>
class FPoint {
private:
	FReal data[3]; ///< all positions x y z

public:	
    /** @brief Default constructor (sets position to 0/0/0) */
    FPoint<FReal>(){
		data[0] = data[1] = data[2] = FReal(0.0);
	}

    /** @brief Constructor from an array */
    explicit FPoint<FReal>(const FReal inPosition[3]){
		data[0] = inPosition[0];
		data[1] = inPosition[1];
		data[2] = inPosition[2];
	}

    /** @brief Constructor from values */
    explicit FPoint<FReal>(const FReal inX,const FReal inY,const FReal inZ){
		data[0] = inX;
		data[1] = inY;
		data[2] = inZ;
	}

    /** @brief Copy constructor
     * @param other the source class to copy
     */
    FPoint<FReal>(const FPoint<FReal>& other) {
		data[0] = other.data[0];
		data[1] = other.data[1];
		data[2] = other.data[2];
	}

    /** @brief Assignment operator
     * @param other the source class to copy
     */
    FPoint<FReal>(const FPoint<FReal>& other, const FReal addset) {
		data[0] = other.data[0] + addset;
		data[1] = other.data[1] + addset;
		data[2] = other.data[2] + addset;
	}

    /** @brief Copy constructor
     * @param other the source class to copy
     * @return this a reference to the current class
     */
    FPoint<FReal>& operator=(const FPoint<FReal>& other){
        this->data[0] = other.data[0];
        this->data[1] = other.data[1];
        this->data[2] = other.data[2];
        return *this;
    }

    /** @brief Position setter
     * @param other the source class to copy
     * @return this a reference to the current class
     */
    void setPosition(const FReal inX,const FReal inY,const FReal inZ){
        this->data[0] = inX;
        this->data[1] = inY;
        this->data[2] = inZ;
    }

    /** @brief Get x
     * @return this->data[0]
     */
    FReal getX() const{
        return this->data[0];
    }

    /** @brief Get y
     * @return this->data[1]
     */
    FReal getY() const{
        return this->data[1];
    }

    /** @brief Get z
     * @return this->data[2]
     */
    FReal getZ() const{
        return this->data[2];
    }

    /** @brief Set x
     * @param the new x
     */
    void setX(const FReal inX){
        this->data[0] = inX;
    }

    /** @brief Set y
     * @param the new y
     */
    void setY(const FReal inY){
        this->data[1] = inY;
    }

    /** @brief Set z
     * @param the new z
     */
    void setZ(const FReal inZ){
        this->data[2] = inZ;
    }

    /** @brief Add to the x-dimension the inX value
     * @param  inXthe increment in x
     */
    void incX(const FReal inX){
        this->data[0] += inX;
    }

    /** @brief Add to the y-dimension the inY value
     * @param  in<<<<<<y the increment in y
     */
    void incY(const FReal inY){
        this->data[1] += inY;
    }

    /** @brief Add to z-dimension the inZ value
     * @param inZ the increment in z
     */
    void incZ(const FReal inZ){
        this->data[2] += inZ;
    }
    /** @brief Get a pointer on the coordinate of FPoint<FReal>
     * @return the data value array
     */
    FReal * getDataValue(){
        return this->data ;
    }
    /** @brief Get a pointer on the coordinate of FPoint<FReal>
     * @return the data value array
     */
    const FReal *  getDataValue()  const{
        return this->data ;
    }

    /** @brief Compute the distance to the origin
     * @return the norm of the Fpoint
     */
    FReal norm() const {
        return FMath::Sqrt(this->data[0]*this->data[0]+this->data[1]*this->data[1]
                           +this->data[2]*this->data[2]) ;
    }

    /** @brief Compute the distance to the origin
     * @return the square norm of the Fpoint
     */
    FReal norm2() const {
        return (this->data[0]*this->data[0]+this->data[1]*this->data[1]
                +this->data[2]*this->data[2]) ;
    }

    /** @brief Subtracts value from all dimensions
     * 
     * @param inValue the value to substract
     * @return the current object after being subtracted
     */
    FPoint<FReal>& operator-=(const FReal inValue){
        this->data[0] -= inValue;
        this->data[1] -= inValue;
        this->data[2] -= inValue;
        return *this;
    }

    /** @brief Adds value to all dimensions
     * 
     * @param inValue the value to affect
     * @return the current object after being affected
     */
    FPoint<FReal>& operator+=(const FReal inValue){
        this->data[0] += inValue;
        this->data[1] += inValue;
        this->data[2] += inValue;
        return *this;
    }

    /** @brief Subtracts the other vector
     *
     * @param other the value to substract
     * @return the current object after being subtracted
     */
    FPoint<FReal>& operator-=(const FPoint<FReal>& other){
        this->data[0] -= other.data[0];
        this->data[1] -= other.data[1];
        this->data[2] -= other.data[2];
        return *this;
    }

    /** @brief Adds the other vector
     *
     * @param other the value to affect
     * @return the current object after being affected
     */
    FPoint<FReal>& operator+=(const FPoint<FReal>& other){
        this->data[0] += other.data[0];
        this->data[1] += other.data[1];
        this->data[2] += other.data[2];
        return *this;
    }
    /** @brief Divides each dimension by the other position's one
     *
     * @param other the value to affect
     * @return the current object after being affected
     */
    FPoint<FReal>& operator/=(const FPoint<FReal>& other){
        this->data[0] /= other.data[0];
        this->data[1] /= other.data[1];
        this->data[2] /= other.data[2];
        return *this;
    }
    /** @brief Multiplies all dimensions by a value
     *
     * @param other the value to affect
     * @return the current object after being affected
     */
    FPoint<FReal>& operator*=(const FReal value){
        this->data[0] *= value;
        this->data[1] *= value;
        this->data[2] *= value;
        return *this;
    }

    /** @brief Operator F3Position minus FReal
     *
     * This substracts inValue to all dimensions of the inPosition
     * @param inPosition the position to compute
     * @param inValue the value to decrease/substract position
     * @return the resulting position
     */
    friend inline FPoint<FReal> operator-(const FPoint<FReal>& inPosition, const FReal inValue){
        return FPoint<FReal>(inPosition, -inValue);
    }

    /** @brief Operator F3Position plus FReal
     *
     * This affects from inValue all dimensions of the inPosition
     * @param inPosition the position to compute
     * @param inValue the value to increase/affect position
     * @return the resulting position
     */
    friend inline FPoint<FReal> operator+(const FPoint<FReal>& inPosition, const FReal inValue){
        return FPoint<FReal>(inPosition, inValue);
    }

    /** @brief Operator F3Position minus F3Position
     * 
     * This subtracts one from anther
     * @param inPosition the position to reduce
     * @param inOther the position to decrease/substract inPosition
     * @return the resulting position
     */
    friend inline FPoint<FReal> operator-(const FPoint<FReal>& inPosition, const FPoint<FReal>& inOther){
        return FPoint<FReal>(inPosition.data[0] - inOther.data[0], inPosition.data[1] - inOther.data[1], inPosition.data[2] - inOther.data[2]);
    }

    /** @brief Operator F3Position plus F3Position
     *
     * This subtracts one from anther
     * @param inPosition the position to reduce
     * @param inOther the position to increase inPosition
     * @return the resulting position
     */
    friend inline FPoint<FReal> operator+(const FPoint<FReal>& inPosition, const FPoint<FReal>& inOther){
        return FPoint<FReal>(inPosition.data[0] + inOther.data[0], inPosition.data[1] + inOther.data[1], inPosition.data[2] + inOther.data[2]);
    }

    /** @brief Compare two points */
    inline bool operator==(const FPoint<FReal>& pp){
        /* do actual comparison */
        return this->data[0]==pp.data[0] &&   this->data[1]==pp.data[1]&&  this->data[2]==pp.data[2];

    }

    /** @brief Compare two points */
    inline bool operator!=( const FPoint<FReal>& rhs)
        {return !(*this == rhs);}

    /** @brief Operator stream FPoint<FReal> to std::ostream
     * 
     * This can be used to write out a position.
     * @param[in,out] output where to write the position
     * @param[in] inPosition the position to write out
     * @return the output for multiple << operators
     */
    template <class StreamClass>
    friend StreamClass& operator<<(StreamClass& output, const FPoint<FReal>& inPosition){
        output << "(" <<  inPosition.getX() << ", " << inPosition.getY() << ", " << inPosition.getZ() <<")";
        return output;  // for multiple << operators.
    }
    /** @brief Operator stream FPoint<FReal> to std::istream
     *
     * This can be used to write out a position.
     * @param[in,out] input where to write the position
     * @param[out] outPosition the position to write out
     * @return the input for multiple << operators
     */
    template <class StreamClass>
    friend StreamClass& operator>>(StreamClass& input,  FPoint<FReal>& outPosition){
        FReal x,y,z;
        input >> x>> y>> z ;
        outPosition.setPosition(x,y,z);
        return input;  // for multiple << operators.
    }

    /** @brief Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        buffer << data[0] << data[1] << data[2];
    }
    /** @brief Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        buffer >> data[0] >> data[1] >> data[2];
    }
};




#endif //FPOINT_HPP


