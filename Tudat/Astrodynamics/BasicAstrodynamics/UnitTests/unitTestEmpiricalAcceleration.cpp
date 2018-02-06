/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace simulation_setup;
using namespace propagators;
using namespace numerical_integrators;
using namespace orbital_element_conversions;
using namespace basic_mathematics;
using namespace basic_astrodynamics;
using namespace gravitation;
using namespace numerical_integrators;
using namespace unit_conversions;
using namespace propagators;


BOOST_AUTO_TEST_SUITE( test_empirical_acceleration )

//! Test if empirical accelerations are computed correctly in propagation.
BOOST_AUTO_TEST_CASE( testEmpiricalAccelerations )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Test three different cases of empirical acceleration values
    for( unsigned testCase = 0; testCase < 3; testCase++ )
    {
        // Set simulation end epoch.
        const double simulationEndEpoch = 2.0 * tudat::physical_constants::JULIAN_DAY;

        // Set numerical integration fixed step size.
        const double fixedStepSize = 15.0;

        // Define body settings for simulation.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
        bodySettings[ "Earth" ] = boost::make_shared< BodySettings >( );
        bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ), "SSB", "J2000" );
        bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

        // Create Earth object
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );


        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );

        // Define empirical acceleration values for current case
        double empiricalAccelerationNorm = 1.0E-8;
        if( testCase == 0 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) ) );
        }
        else if( testCase == 1 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                             Eigen::Vector3d::Zero( ),
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ),
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) ) );
        }
        else if( testCase == 2 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< EmpiricalAccelerationSettings >(
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ),
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ),
                                                             empiricalAccelerationNorm * Eigen::Vector3d::UnitZ( ) ) );
        }
        accelerationMap[  "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );

        // Define dependent variable settings
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        empirical_acceleration, "Asterix", "Earth", 0 ) );

        // Define propagator settings
        TranslationalPropagatorType propagatorType = encke;
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, propagatorType,
                  boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        // Define integrator settings
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( rungeKuttaVariableStepSize, 0.0, fixedStepSize,
                  RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0E-4, 3600.0, 1.0E-14, 1.0E-14 );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

        // Compute multiplier used in tests
        double testMultiplier = 1.0;
        if( testCase == 2 )
        {
            testMultiplier = std::sqrt( 2.0 );
        }

        Eigen::Vector3d vectorRInRswFrame;
        Eigen::Vector3d vectorVInRswFrame;
        Eigen::Vector3d vectorWInRswFrame;

        // Compare computed and expected empirical acceleration and its properties
        for( std::map< double, Eigen::VectorXd >::const_iterator resultIterator = integrationResult.begin( );
             resultIterator != integrationResult.end( ); resultIterator++ )
        {
            Eigen::Vector3d totalEmpiricalAcceleration =
                    dependentVariableResult.at( resultIterator->first ).segment( 0, 3 );

            // Check expected norm of empirical acceleration
            BOOST_CHECK_CLOSE_FRACTION( totalEmpiricalAcceleration.norm( ), testMultiplier * empiricalAccelerationNorm,
                                        8.0 * std::numeric_limits< double >::epsilon( ) );
            Eigen::Matrix3d currentRotationMatrix =
                    reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx( resultIterator->second );
            Eigen::Vector3d totalEmpiricalAccelerationInRswFrame = currentRotationMatrix * totalEmpiricalAcceleration;

            // Compute vector directions in RSW frame
            vectorRInRswFrame = currentRotationMatrix * resultIterator->second.segment( 0, 3 );
            vectorVInRswFrame = currentRotationMatrix * resultIterator->second.segment( 3, 3 );
            vectorWInRswFrame = currentRotationMatrix * ( Eigen::Vector3d( resultIterator->second.segment( 0, 3 ) ).cross(
                                                       Eigen::Vector3d( resultIterator->second.segment( 3, 3 ) ) ) );

            // Check direction of vectors
            BOOST_CHECK_EQUAL( vectorRInRswFrame( 0 ) > 0, true );
            BOOST_CHECK_EQUAL( vectorVInRswFrame( 1 ) > 0, true );
            BOOST_CHECK_EQUAL( vectorWInRswFrame( 2 ) > 0, true );

            // Check direction of R vector in RSW frame
            BOOST_CHECK_EQUAL( vectorRInRswFrame( 0 ) / ( 10.0 * std::numeric_limits< double >::epsilon( ) ) >
                               std::fabs( vectorRInRswFrame( 1 ) ), true );
            BOOST_CHECK_EQUAL( vectorRInRswFrame( 0 ) / ( 10.0 * std::numeric_limits< double >::epsilon( ) ) >
                               std::fabs( vectorRInRswFrame( 2 ) ), true );

            // Check direction of V vector in RSW frame
            BOOST_CHECK_EQUAL( vectorVInRswFrame( 1 ) / 5.0 > std::fabs( vectorVInRswFrame( 0 ) ) , true );
            BOOST_CHECK_EQUAL( vectorVInRswFrame( 1 ) / ( 10.0 * std::numeric_limits< double >::epsilon( ) ) >
                               std::fabs( vectorVInRswFrame( 2 ) ), true );

            // Check direction of W vector in RSW frame
            BOOST_CHECK_EQUAL( vectorWInRswFrame( 2 ) / ( 10.0 * std::numeric_limits< double >::epsilon( ) ) >
                               std::fabs( vectorWInRswFrame( 0 ) ), true );
            BOOST_CHECK_EQUAL( vectorWInRswFrame( 2 ) / ( 10.0 * std::numeric_limits< double >::epsilon( ) ) >
                               std::fabs( vectorWInRswFrame( 1 ) ), true );

            // Compute expected empirical acceleration
            Eigen::VectorXd currentKeplerianState = convertCartesianToKeplerianElements(
                        Eigen::Vector6d( resultIterator->second ), earthGravitationalParameter );
            Eigen::Vector3d expectedAccelerationInRswFrame;
            if( testCase == 0 )
            {
                expectedAccelerationInRswFrame = empiricalAccelerationNorm * Eigen::Vector3d::UnitX( );
            }
            else if( testCase == 1 )
            {
                expectedAccelerationInRswFrame =
                        empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) * std::sin( currentKeplerianState( 5 ) ) +
                        empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) * std::cos( currentKeplerianState( 5 ) );
            }
            else if( testCase == 2 )
            {
                expectedAccelerationInRswFrame =
                        empiricalAccelerationNorm * Eigen::Vector3d::UnitX( ) +
                        empiricalAccelerationNorm * Eigen::Vector3d::UnitY( ) * std::sin( currentKeplerianState( 5 ) ) +
                        empiricalAccelerationNorm * Eigen::Vector3d::UnitZ( ) * std::cos( currentKeplerianState( 5 ) );
            }

            // Compare expected and computed empirical acceleration
            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( expectedAccelerationInRswFrame( i ) - totalEmpiricalAccelerationInRswFrame( i ) ),
                    8.0 * std::numeric_limits< double >::epsilon( ) * empiricalAccelerationNorm );
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )


}

}

