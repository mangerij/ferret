// ===================================================================================
// Copyright ScalFmm 2014 INRIA
//
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
#ifndef FGENERATEDISTRIBUTION_HPP
#define FGENERATEDISTRIBUTION_HPP

// @author O. Coulaud

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
//
#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"

/**  return a random number between 0 and 1 */

void initRandom() {
	srand48( static_cast<long int>(time(nullptr))) ;
} ;
template <class FReal>
FReal getRandom() {
	return static_cast<FReal>(drand48());
	//return static_cast<FReal>(rand()/FReal(RAND_MAX));
} ;
//!  \fn   unifRandonPointsOnUnitCube(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the unit cube

//!
//! \param N the number of points uniformly randomly sample on the unit cube
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
template <class FReal>
void unifRandonPointsOnUnitCube(const FSize N , FReal * points) {
	//
	initRandom() ;
	int j = 0;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
		//
        points[j]	  =	getRandom<FReal>()  ;
        points[j+1] =	getRandom<FReal>()  ;
        points[j+2] =	getRandom<FReal>()  ;
		//
	}
};
//!  \fn   unifRandonPointsOnCube(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the cube of length R

//!
//! \param N the number of points uniformly randomly sample on the unit cube
//! \param Lx the the X-length of the  cube
//! \param Ly the the Y-length of the  cube
//! \param Lz the the Z-length of the  cube
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
template <class FReal>
void unifRandonPointsOnCube(const FSize N , const FReal& Lx,  const FReal &Ly,  const FReal& Lz, FReal * points) {
	//
	unifRandonPointsOnUnitCube(N , points) ;
    FSize j =0 ;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
		points[j]	   *= Lx ;
		points[j+1]  *= Ly ;
		points[j+2]  *= Lz ;
	}
};
//!  \fn   unifRandonPointsOnUnitSphere(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the unit sphere

//!
//! \param N the number of points uniformly randomly sample on the unit sphere
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
template <class FReal>
void unifRandonPointsOnUnitSphere(const FSize N , FReal * points) {
	FReal u, v, theta, phi, sinPhi ;
	//
	initRandom() ;
    FSize j = 0 ;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
		//
        u = getRandom<FReal>() ;  v = getRandom<FReal>() ;
        theta  = FMath::FTwoPi<FReal>()*u ;
		phi     = FMath::ACos(2*v-1);
		sinPhi = FMath::Sin(phi);
		//
		points[j]	  =	FMath::Cos(theta)*sinPhi ;
		points[j+1] =	FMath::Sin(theta)*sinPhi ;
		points[j+2] =	2*v-1 ;
		//
	}
};
//!  \fn  nonunifRandonPointsOnElipsoid(const int N , const FReal &a, const FReal &b, const FReal &c, FReal * points)

//! \brief  Generate N points non uniformly distributed on the ellipsoid of  aspect ratio a:b:c

//!
//! \param N the number of points
//! \param a  the x  semi-axe length
//! \param b  the y  semi-axe length
//! \param c  the z  semi-axe length
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
template <class FReal>
void nonunifRandonPointsOnElipsoid(const FSize N , const FReal &a, const FReal &b, const FReal &c, FReal * points) {
	//
	FReal u, v , cosu ;
    FSize j =0 ;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
        u = getRandom<FReal>() ;  v = getRandom<FReal>() ;
        u  = FMath::FPi<FReal>()*u - FMath::FPiDiv2<FReal>();   v   = FMath::FTwoPi<FReal>()*v - FMath::FPi<FReal>();
		cosu = FMath::Cos(u) ;
		points[j]	   = a*cosu*FMath::Cos(v)  ;
		points[j+1]  = b*cosu*FMath::Sin(v)  ;
		points[j+2]  = c*FMath::Sin(u)  ;
	}
};
//!  \fn  nonunifRandonPointsOnElipsoid(const int N , const FReal &a, const FReal &c, FReal * points)

//! \brief  Generate N points uniformly distributed on the ellipsoid of  aspect ratio a:a:c

//!
//! \param N the number of points
//! \param a  the x  semi-axe length
//! \param c  the z  semi-axe length
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
template <class FReal>
void unifRandonPointsOnProlate(const FSize N , const FReal &a, const FReal &c, FReal * points){
	//
	FReal u, w,v ,ksi ;
	FReal e = (a*a*a*a)/(c*c*c*c) ;
	bool isgood = false;
    FSize j =0 , cpt =0 ;
	//
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
		// Select a random point on the prolate
		do {
			cpt++	;
            u = getRandom<FReal>() ;  v = getRandom<FReal>() ;
            u  = 2.0*u - 1.0;   v   = FMath::FTwoPi<FReal>()*v;
			w =FMath::Sqrt(1-u*u) ;
			points[j]	   = a*w*FMath::Cos(v)  ;
			points[j+1]  = a*w*FMath::Sin(v)  ;
			points[j+2]  = c*u ;
			// Accept the position ?
            ksi = a*getRandom<FReal>()  ;
			//			std::cout << "Gradf  "<<  points[j]*points[j] + points[j+1] *points[j+1]  +e*points[j+2] *points[j+2]  << std::endl;
			isgood = (points[j]*points[j] + points[j+1] *points[j+1]  +e*points[j+2] *points[j+2]  < ksi*ksi );
		} while (isgood);
	}
	std::cout.precision(4);
    std::cout << "Total tested points: "<< cpt << " % of rejected points: "<<100*static_cast<FReal>(cpt-N)/static_cast<FReal>(cpt) << " %" <<std::endl;

} ;

//!  \fn  unifRandonPointsOnHyperPara(const int N , const FReal &a, const FReal &b, const FReal &c, FReal * points)

//! \brief  Generate N points uniformly distributed on the hyperbolic paraboloid of  aspect ratio a:b:c

//!
//! \param N the number of points
//! \param a  the x  semi-axe length
//! \param b  the y  semi-axe length
//! \param c  the z  semi-axe length
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
template <class FReal>
void unifRandonPointsOnHyperPara(const FSize N , const FReal &a, const FReal &b, const FReal &c, FReal * points) {
    //
    FReal u, v ;
    FSize j =0 ;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
        u = 2.0*getRandom<FReal>() - 1.0 ;  v = 2.0*getRandom<FReal>() - 1.0 ;
        points[j]    = a*u ;
        points[j+1]  = b*v ;
        points[j+2]  = c*(u*u - v*v)  ;
    }
};


//!  \fn  unifRandonPointsOnSphere(const int N , const FReal R, FReal * points)

//! \brief Generate N points uniformly distributed on the sphere of radius R

//!
//! \param N the number of points uniformly randomly sample on the sphere
//! \param R the radius of the sphere
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
template <class FReal>
void unifRandonPointsOnSphere(const FSize N , const FReal R, FReal * points) {
	//
	unifRandonPointsOnUnitSphere(N , points) ;
    FSize j =0 ;
    for (FSize i = 0 ; i< N ; ++i, j+=4)  {
		points[j]	   *= R ;
		points[j+1]  *= R ;
		points[j+2]  *= R ;
	}
};
//!  \fn void plummerDist(int & cpt, const FReal &R)

//! \brief   Radial Plummer distribution

//!
//! \param cpt : counter to know how many random selections we need to obtain a radius less than R
//! \param R    : Radius of the sphere that contains the particles
//! @return Return the radius according to the Plummer distribution either double type or float type
//!
template <class FReal>
FReal  plummerDist(FSize cpt, const FReal &R) {
	//
	FReal radius ,u ;
	do  {
		//
        u        = FMath::pow (getRandom<FReal>() , 2.0/3.0) ;
		radius = FMath::Sqrt (u/(1.0-u));
		cpt++;
		if(radius  <=R){
			//			std::cout << radius << "    "  <<std::endl;
			return static_cast<FReal>(radius);
		}
	} while (true);
}
//! \fn void unifRandonPlummer(const int N , const FReal R, const FReal M, FReal * points)

//! \brief  Build N points following the Plummer distribution

//! First we construct N points uniformly distributed on the unit sphere. Then the radius in construct according to the Plummr distribution.
//!
//! \param N the number of points following the Plummer distribution
//! \param R the radius of the sphere that contains all the points
//! \param M the total mass of all the particles inside the Sphere or radius R
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
template <class FReal>
void unifRandonPlummer(const FSize N , const FReal R, const FReal M, FReal * points) {
	//
	unifRandonPointsOnUnitSphere(N , points) ;
	//
	FReal r , rm= 0.0;
    //	FReal Coeff =  3.0*M/(4.0*FMath::FPi<FReal>()*R*R*R) ;
	//am1 = 0 ;//1/FMath::pow(1+R*R,2.5);
    FSize cpt = 0 ;
    for (FSize i = 0,j=0 ; i< N ; ++i, j+=4)  {
		// u \in []
		r = plummerDist(cpt,R) ;
		rm = std::max(rm, r);
		points[j]	   *= r ;
		points[j+1]  *= r ;
		points[j+2]  *= r ;
	}

	std::cout << "Total tested points: "<< cpt << " % of rejected points: "
            <<100*static_cast<FReal>(cpt-N)/static_cast<FReal>(cpt) << " %" <<std::endl;

} ;
//
#endif
