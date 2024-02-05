

#include <ns3/core-module.h>
#include <ns3/walker-constellation-helper.h>


using namespace ns3;


int main ( int argc, char *argv[] )
{
    Config::SetDefault("ns3::SatSGP4MobilityModel::StartDateStr", StringValue("2023-03-20 21:24:00"));

    Ptr<WalkerConstellationHelper> constellation = CreateObjectWithAttributes<WalkerConstellationHelper> (
        "Inclination", DoubleValue(74.9),
        "NumOfOrbits", IntegerValue(6),
        "NumOfSats", IntegerValue(9),
        "Altitude", DoubleValue(1200)
    );

    constellation->Initialize();
    constellation->LogInitialPositions("./output/constellation", ".txt");
    
}