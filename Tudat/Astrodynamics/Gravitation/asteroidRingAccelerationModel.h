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

#ifndef TUDAT_ASTEROID_RING_MODEL_H
#define TUDAT_ASTEROID_RING_MODEL_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
//#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h"

namespace tudat
{
namespace gravitation
{


Eigen::Vector3d computeAsteroidRingAcceleration(
       const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
       const Eigen::Vector3d& positionOfBodyExertingAcceleration,
       const double& ringMass );


//template< typename StateMatrix = Eigen::Vector3d >
class AsteroidRingGravitationalModel   : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >

{
private:

    //Mass function Typdef
    typedef boost::function< double ( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef boost::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:
        AsteroidRingGravitationalModel(
                const Vector3dReturningFunction positionOfBodySubjectToAccelerationFunction,
                const Vector3dReturningFunction positionOfBodyExertingAccelerationFunction,
                const DoubleReturningFunction ringMassFunction)

            : positionOfBodyExertingAccelerationFunction_(positionOfBodyExertingAccelerationFunction),
              positionOfBodySubjectToAccelerationFunction_(positionOfBodySubjectToAccelerationFunction),
              ringMassFunction_(ringMassFunction)
//              positionOfCentralBodyFunction_(positionOfCentralBodyFunction)


        {
            this->updateMembers( );
        }

        //! Get gravitational acceleration.
        /*!
         * Returns the gravitational acceleration computed using the input parameters provided to the
         * class. This function serves as a wrapper for the computeGravitationalAcceleration()
         * function.
         * \return Computed gravitational acceleration vector.
         */
        Eigen::Vector3d getAcceleration( )
        {
            return computeAsteroidRingAcceleration(
                        positionOfBodyExertingAcceleration_,
                        positionOfBodySubjectToAcceleration_,
                        ringMass_);
        }

        //! Update members.
        /*!
         * Updates class members relevant for computing the central gravitational acceleration. In this
         * case the function simply updates the members in the base class.
         * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
         * \param currentTime Time at which acceleration model is to be updated.
         */
        void updateMembers( const double currentTime = TUDAT_NAN )
        {
            if( !( this->currentTime_ == currentTime ) )
            {

                positionOfBodyExertingAcceleration_ = this->positionOfBodyExertingAccelerationFunction_();
                positionOfBodySubjectToAcceleration_ = this->positionOfBodySubjectToAccelerationFunction_();
                ringMass_ = this->ringMassFunction_();


                currentTime_ = currentTime;

            }
        }

private:

      Vector3dReturningFunction positionOfBodySubjectToAccelerationFunction_;
      Vector3dReturningFunction positionOfBodyExertingAccelerationFunction_;
      DoubleReturningFunction ringMassFunction_;
//      Vector3dReturningFunction positionOfCentralBodyFunction_;


      Eigen::Vector3d positionOfBodyExertingAcceleration_;
      Eigen::Vector3d positionOfBodySubjectToAcceleration_;
      double ringMass_;
//      Eigen::Vector3d positionOfCentralBody_;



    };
}
}

#endif // TUDAT_ASTEROID_RING_MODEL_H
