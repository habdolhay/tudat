#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/GroundStations/pointingAnglesCalculator.h"
#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;


BOOST_AUTO_TEST_SUITE( test_pointing_angles_calculator )

BOOST_AUTO_TEST_CASE( test_PointingAnglesCalculator )
{
    // Define shape model (sphere)
    boost::shared_ptr< SphericalBodyShapeModel > bodyShape = boost::make_shared< SphericalBodyShapeModel >( 6.371E6 );

    // Define test ground station point.
    double groundStationDistance = 6371.0E3;
    Eigen::Vector3d groundStationPosition;
    groundStationPosition<<groundStationDistance, 0.0, 0.0;

    double degreesToRadians = unit_conversions::convertDegreesToRadians( 1.0 );

    // Test analytically checked azimuth and elevation
    {
        boost::shared_ptr< NominalGroundStationState > stationState = boost::make_shared< NominalGroundStationState >(
                    groundStationPosition, bodyShape );
        boost::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = boost::make_shared< PointingAnglesCalculator >(
                    boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
                    boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

        double testLatitude = 30.0 * degreesToRadians;
        double testLongitude = 0.0;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint;
        testSphericalPoint<<testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        double testAzimuth = pointAnglesCalculator->calculationAzimuthAngle( testCartesianPoint, 0.0 );
        double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

        double expectedAzimuth = 90.0 * degreesToRadians;
        double expectedElevation = 60.0 * degreesToRadians;

        BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 3.0 * std::numeric_limits< double >::epsilon( ) );
    }

    //http://www.movable-type.co.uk/scripts/latlong.html
    {
        boost::shared_ptr< NominalGroundStationState > stationState = boost::make_shared< NominalGroundStationState >(
                    groundStationPosition, bodyShape );
        boost::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = boost::make_shared< PointingAnglesCalculator >(
                    boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
                    boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

        double testLatitude = 21.0 * degreesToRadians;
        double testLongitude = 84.0 * degreesToRadians;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint;
        testSphericalPoint<<testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        double testAzimuth = pointAnglesCalculator->calculationAzimuthAngle( testCartesianPoint, 0.0 );
        double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

        double expectedElevation = mathematical_constants::PI / 2.0 - 9385.0 / 6371.0;
        double expectedAzimuth = mathematical_constants::PI / 2.0 - ( 68.0 + 53.0 / 60.0 + 40.0 / 3600.0 ) * degreesToRadians;

        BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 1.0E-5 );
        BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 1.0E-3 );
    }

    {
        boost::shared_ptr< NominalGroundStationState > stationState = boost::make_shared< NominalGroundStationState >(
                    groundStationPosition, bodyShape );
        boost::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = boost::make_shared< PointingAnglesCalculator >(
                    boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
                    boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

        double testLatitude = -38.0 * degreesToRadians;
        double testLongitude = 234.0 * degreesToRadians;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint;
        testSphericalPoint<<testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        double testAzimuth = pointAnglesCalculator->calculationAzimuthAngle( testCartesianPoint, 0.0 );
        double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

        double expectedElevation = mathematical_constants::PI / 2.0 - 13080.0 / 6371.0;
        double expectedAzimuth = mathematical_constants::PI / 2.0 - ( 225.0 + 59.0 / 60.0 + 56.0 / 3600.0 ) * degreesToRadians;

        BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 1.0E-5 );
        BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 3.0E-2 );

        std::pair< double, double > pointingAngles = pointAnglesCalculator->calculatePointingAngles( testCartesianPoint, 0.0 );

        BOOST_CHECK_CLOSE_FRACTION( pointingAngles.first, testElevation, 3.0E-2 );
        BOOST_CHECK_CLOSE_FRACTION( pointingAngles.second, testAzimuth, 1.0E-5 );
    }

    {

        double poleRightAscension = 56.0 * degreesToRadians;
        double poleDeclination = 45.0 * degreesToRadians;

        Eigen::Quaterniond inertialToBodyFixedFrame = reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    poleDeclination, poleRightAscension, 0.0 );

        groundStationPosition<<1234.0E3, -4539E3, 4298E3;
        boost::shared_ptr< NominalGroundStationState > stationState = boost::make_shared< NominalGroundStationState >(
                    groundStationPosition, bodyShape );

        boost::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = boost::make_shared< PointingAnglesCalculator >(
                    boost::lambda::constant( inertialToBodyFixedFrame ),
                    boost::bind( &NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

        double testLatitude = -38.0 * degreesToRadians;
        double testLongitude = 234.0 * degreesToRadians;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint;
        testSphericalPoint<<testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        Eigen::Vector3d testPointInLocalFrame = pointAnglesCalculator->convertVectorFromInertialToTopocentricFrame( testCartesianPoint, 0.0 );

        Eigen::Vector3d expectedTestPointInLocalFrame = stationState->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) * inertialToBodyFixedFrame *
               testCartesianPoint ;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPointInLocalFrame, expectedTestPointInLocalFrame, std::numeric_limits< double >::epsilon( ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
