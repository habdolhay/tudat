/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include "Tudat/Astrodynamics/Gravitation/asteroidRingAccelerationModel.h"

namespace tudat
{
namespace gravitation
{




//! Compute acceleration from a Ring
//!
//! Based on ring model presented by (Kuchynca et. al 2010)
//! KB and EB are first and second order elliptical integrals with \beta as input evaluated beteen 0 and \pi/2 r_{ring} i
//! a = -\frac{2 G M_{ring} {\pi \alpha (1 - \beta) L^(3/2)} (\alpha r_{body} EB + ((1 - \alpha) KB - EB)r_{ring} i)
//!
Eigen::Vector3d computeAsteroidRingAcceleration(const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration, const double &ringMass)
{
    Eigen::Vector3d acc;

    //Constants (For now determined internally)
    //==========CHANGE THEM LATER=======
    double pi = 3.14159;
    double G = 6.67259e-11;
    double AU = 1.49597870691e11;
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



    //MAKE THESE CLASS INPUTS (in AccelerationSettings.h or Body's properties)
    double RRing = 2.8*AU;
    double MRing = 0.34e-10 * 1.98855e30;
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    //Do rotation arround x-axis by 23.008 degrees (Ring has 0 inclination wrt a pre-defined invraibale plane which is defined has 23.0008 inclination wrt ICRF)
//    Eigen::Matrix3d rotationX;
//    double inclination = 0.40156;
//    rotationX << 1, 0, 0,
//                0, cos(inclination), sin(inclination),
//                0, -sin(inclination), cos(inclination);
//    Eigen::Vector3d transformedPositionOfBodySubjectToAcceleration;
//    transformedPositionOfBodySubjectToAcceleration = rotationX*positionOfBodySubjectToAcceleration;
//    //Same as Range but just in case
//    double transformedRangeOfBodySubjectToAcceleration = transformedPositionOfBodySubjectToAcceleration.norm();
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    //Assume angle I to be equivilant to body's declination in ECLIPJ2000 centered at SSB
    //=========CHANGE LATER=========
    double rangeOfBodySubjectToAcceleration = positionOfBodySubjectToAcceleration.norm();
    double phi = acos(positionOfBodySubjectToAcceleration[2]/rangeOfBodySubjectToAcceleration);
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    Eigen::Vector3d vectorI = rangeOfBodySubjectToAcceleration*sin(phi)*positionOfBodySubjectToAcceleration/(rangeOfBodySubjectToAcceleration*rangeOfBodySubjectToAcceleration);
    Eigen::Vector3d unitVectorI = vectorI/(vectorI.norm());



    double alpha = 2*rangeOfBodySubjectToAcceleration*sin(phi)*RRing/(rangeOfBodySubjectToAcceleration*rangeOfBodySubjectToAcceleration + RRing*RRing);
    double beta  = 2*alpha/(1+alpha);
    double L = rangeOfBodySubjectToAcceleration*rangeOfBodySubjectToAcceleration + RRing*RRing + 2*RRing*rangeOfBodySubjectToAcceleration*sin(phi);

    double KB = boost::math::ellint_1(beta, pi/2.0);
    double EB = boost::math::ellint_2(beta, pi/2.0);

    acc = - 2*G*MRing/(pi*alpha*(1-beta)*pow(L,3/2.0))*(alpha*EB*positionOfBodySubjectToAcceleration + ((1-alpha)*KB - EB)*RRing*unitVectorI);
    return acc;


}



} // namespace gravitation
} // namespace tudat

