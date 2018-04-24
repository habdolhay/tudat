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
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

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
Eigen::Vector3d computeAsteroidRingAcceleration(const Eigen::Vector3d& positionOfBodySubjectToAcceleration, const double &ringMass, std::string& frameOrigin, std::string& frameOrientation)
{
    Eigen::Vector3d acc;


    //Ring radius is set to be constant and fixed to its optimal value(if RRing changes, mass should change accordingly)
    double RRing = 2.8*physical_constants::ASTRONOMICAL_UNIT;
    double MRing = ringMass;


    //Phi is equivilant to body's declination complementary angle in ECLIPJ2000 centered at SSB
    double rangeOfBodySubjectToAcceleration = positionOfBodySubjectToAcceleration.norm();
    double phi = acos(positionOfBodySubjectToAcceleration[2]/rangeOfBodySubjectToAcceleration);

    //VectorI is the projected vector of planet's position onto the asteroid plane = planet's position in ECLIPJ2000 with z=0(simplest way)
    Eigen::Vector3d vectorI = positionOfBodySubjectToAcceleration;
    vectorI[2] = 0;
    Eigen::Vector3d unitVectorI = vectorI/(vectorI.norm());

    double alpha = 2*rangeOfBodySubjectToAcceleration*sin(phi)*RRing/(rangeOfBodySubjectToAcceleration*rangeOfBodySubjectToAcceleration + RRing*RRing);
    double beta  = 2*alpha/(1+alpha);
    double L = rangeOfBodySubjectToAcceleration*rangeOfBodySubjectToAcceleration + RRing*RRing + 2*RRing*rangeOfBodySubjectToAcceleration*sin(phi);

    double KB = boost::math::ellint_1(sqrt(beta), mathematical_constants::PI/2.0);
    double EB = boost::math::ellint_2(sqrt(beta),  mathematical_constants::PI/2.0);

    acc = -2*physical_constants::GRAVITATIONAL_CONSTANT*MRing/(mathematical_constants::PI*alpha*(1-beta)*pow(L,3/2.0))*(alpha*EB*positionOfBodySubjectToAcceleration + ((1-alpha)*KB - EB)*RRing*unitVectorI);
    return acc;


}



} // namespace gravitation
} // namespace tudat

