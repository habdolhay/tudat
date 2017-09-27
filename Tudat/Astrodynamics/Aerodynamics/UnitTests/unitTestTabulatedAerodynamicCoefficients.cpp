/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    semiMajorAixscopy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "applicationOutput.h"
#include "Tudat/InputOutput/basicInputOutput.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

BOOST_AUTO_TEST_CASE( testTabulatedDragCoefficient )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace interpolators;
    using namespace input_output;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ ) {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // EARTH
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

    // MOON
    bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );
    double referenceArea = 10.0;

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    // Read CD as semiMajorAixsfunction of altitude
    Eigen::MatrixXd aerodynamicsDataFromFile = readMatrixFromFile(
                tudat_applications::getOutputPath( ) + "tabulatedDragCoefficient.txt" );
    Eigen::VectorXd altitudeInKm = aerodynamicsDataFromFile.col(0);
    std::vector< double > altitudes;
    for ( int i = 0; i < altitudeInKm.rows( ); i++ )
    {
        altitudes.push_back( altitudeInKm( i ) * 1.0E3 );
    }
    Eigen::VectorXd dragCoefficientsFromFile = aerodynamicsDataFromFile.col( 1 );
    std::vector< double > dragCoefficients(
                dragCoefficientsFromFile.data( ), dragCoefficientsFromFile.data( ) + dragCoefficientsFromFile.size( ) );
    std::vector< Eigen::Vector3d > aerodynamicCoefficients;
    for ( unsigned int i = 0; i < dragCoefficients.size( ); i++ )
    {
        aerodynamicCoefficients.push_back( dragCoefficients[ i ] * Eigen::Vector3d::UnitX( ) );
    }

    // Create interpolator
    boost::shared_ptr< InterpolatorSettings > interpolatorSettings =
            boost::make_shared< InterpolatorSettings >( OneDimensionalInterpolatorTypes::linear_interpolator );

    // Tabulated aerodynamic settings
    aerodynamicCoefficientSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                altitudes, aerodynamicCoefficients, referenceArea,
                aerodynamics::altitude_dependent, interpolatorSettings, 1, 1 );

    // Aerodynamics interface
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double meanEarthRadius = 6371.0E3;
    double perigeeAltitude = 150.0E3;
    double apogeeAltitude = 35780.0E3;
    double semiMajorAixs= (apogeeAltitude + perigeeAltitude)/2 + meanEarthRadius;
    double eccentricity = (apogeeAltitude - perigeeAltitude)/(apogeeAltitude + perigeeAltitude + 2*meanEarthRadius);

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = semiMajorAixs;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = eccentricity;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = mathematical_constants::PI / 180.0 * 23.4;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

    // Hybrid termination conditions
    std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;

    // Time limit
    constituentSettings.push_back( boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                                       simulationEndEpoch ) );

    // Altitude limit
    boost::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings =
            boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::altitude_dependent_variable, "Vehicle" ), 100.0E3, 1 );
    constituentSettings.push_back( altitudeTerminationSettings );

    // Stop if ANY of the two is met
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

    // Save dependent variables
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;
    dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                            altitude_dependent_variable, "Vehicle" ) );
    dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                            aerodynamic_moment_coefficients_dependent_variable, "Vehicle" ) );

    boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave, 0 ) ;

    // Translational propagator settings
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > > (
                centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, terminationSettings );


    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > > (
                rungeKutta4, simulationStartEpoch, 300.0 );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CHECK DRAG COEFFICIENTS         ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Retrieve altitudes and drag coefficients from propagation results

    Eigen::Vector3d currentCoefficients;
    double currentAltitude;
    double computedDragCoefficient;
    LinearInterpolatorDoublePointer dragCoefficientInterpolator =
            boost::make_shared<LinearInterpolatorDouble>( LinearInterpolatorDouble( altitudes, dragCoefficients ) );


    for ( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::iterator outputIterator =
          dependentVariableOutput.begin( ); outputIterator != dependentVariableOutput.end( ); ++outputIterator )
    {
        currentAltitude = outputIterator->second( 0 );
        computedDragCoefficient = dragCoefficientInterpolator->interpolate( currentAltitude );
        currentCoefficients = outputIterator->second.segment( 1, 3 );
        BOOST_CHECK_CLOSE_FRACTION( computedDragCoefficient, currentCoefficients( 0 ),
                                    4.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( currentCoefficients( 1 ) ), std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( currentCoefficients( 2 ) ), std::numeric_limits< double >::epsilon( ) );

    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
