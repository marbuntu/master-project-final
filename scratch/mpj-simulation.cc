
#include <ns3/core-module.h>
#include <ns3/log.h>

#include <ns3/walker-constellation-helper.h>
#include <ns3/geocentric-constant-position-mobility-model.h>
#include <ns3/mobility-model.h>
#include <ns3/waypoint-mobility-model.h>
#include <ns3/walker-constellation-helper.h>
#include <ns3/groundstation-helper.h>

#include <ns3/three-gpp-v2v-propagation-loss-model.h>
#include <ns3/channel-condition-model.h>
#include <ns3/three-gpp-spectrum-propagation-loss-model.h>
#include <ns3/three-gpp-propagation-loss-model.h>


#include <ns3/phased-array-model.h>
#include <ns3/uniform-planar-array.h>
#include <ns3/isotropic-antenna-model.h>
#include <ns3/cosine-antenna-model.h>


#include <ns3/simple-net-device.h>


NS_LOG_COMPONENT_DEFINE("MpjSimulation");


using namespace ns3;


/****************************************************************/

static Ptr<ThreeGppPropagationLossModel> m_propLossModel;
static Ptr<ThreeGppChannelConditionModel> m_condModel; //Ptr<ChannelConditionModel> m_condModel;
static Ptr<ThreeGppSpectrumPropagationLossModel> m_specLossModel;

static Ptr<PhasedArrayModel> gsAntenna;
static Ptr<PhasedArrayModel> satAntenna;

std::ofstream result_file;

/****************************************************************/


enum prop_env_scenario {
    NTN_URBAN,
    NTN_SUBURBAN,
    NTN_DENSEURBAN,
    NTN_RURAL
};


struct SimulationParams
{
    Ptr<MobilityModel> sat;
    Ptr<MobilityModel> gs;
    double txPow;
    double noiseFigure;
    Ptr<PhasedArrayModel> txAntenna;
    Ptr<PhasedArrayModel> rxAntenna;
    double frequency;
    double bandwidth;
    double resourceBlockBandwidth;


    SimulationParams (Ptr<MobilityModel> mobSat, Ptr<MobilityModel> mobGs, double pTxPow, double pNoiseFigure,
                        Ptr<PhasedArrayModel> pTxAntenna, Ptr<PhasedArrayModel> pRxAntenna, double pFrequency, 
                        double pBandwidth, double pResourceBockwidth)
    {
        sat = mobSat;
        gs = mobGs;
        txPow = pTxPow;
        noiseFigure = pNoiseFigure;
        txAntenna = pTxAntenna;
        rxAntenna = pRxAntenna;
        frequency = pFrequency;
        bandwidth = pBandwidth;
        resourceBlockBandwidth = pResourceBockwidth;
    }

};


static bool createPropagationEnvironment( double frequency, enum prop_env_scenario scenario )
{
    ObjectFactory propLossModelFactory;
    ObjectFactory chnConditionModelFactory;
    std::string scenario_str = "";

    switch (scenario)
    {
    case NTN_URBAN:
        propLossModelFactory.SetTypeId (ThreeGppNTNUrbanPropagationLossModel::GetTypeId());
        chnConditionModelFactory.SetTypeId (ThreeGppNTNUrbanChannelConditionModel::GetTypeId());
        scenario_str = "NTN-Urban";

        break;

    case NTN_DENSEURBAN:
        propLossModelFactory.SetTypeId (ThreeGppNTNDenseUrbanPropagationLossModel::GetTypeId());
        chnConditionModelFactory.SetTypeId (ThreeGppNTNDenseUrbanChannelConditionModel::GetTypeId());
        scenario_str = "NTN-DenseUrban";
        break;
    
    case NTN_SUBURBAN:
        propLossModelFactory.SetTypeId (ThreeGppNTNSuburbanPropagationLossModel::GetTypeId());


        chnConditionModelFactory.SetTypeId (ThreeGppNTNSuburbanChannelConditionModel::GetTypeId());
        scenario_str = "NTN-Suburban";
        break;

    case NTN_RURAL:
        propLossModelFactory.SetTypeId (ThreeGppNTNRuralPropagationLossModel::GetTypeId());
        m_propLossModel = propLossModelFactory.Create<ThreeGppNTNRuralPropagationLossModel>();
        m_propLossModel->SetAttribute ("Frequency", DoubleValue(frequency));
        m_propLossModel->SetAttribute ("ShadowingEnabled", BooleanValue(true));

        chnConditionModelFactory.SetTypeId (ThreeGppNTNRuralChannelConditionModel::GetTypeId());
        m_condModel = chnConditionModelFactory.Create<ThreeGppNTNRuralChannelConditionModel>();
        scenario_str = "NTN-Rural";
        break;

    default:
        NS_LOG_ERROR("Unknown Scenario Selected!");
        return false;
    }
    

    std::cout << "Selected Scenario: " << scenario_str << "\n";


    // Create Propagation Loss Model


    // Create Spectrum Loss Model
    m_specLossModel = CreateObject<ThreeGppSpectrumPropagationLossModel> ();
    m_specLossModel->SetChannelModelAttribute ("Frequency", DoubleValue(frequency));
    m_specLossModel->SetChannelModelAttribute ("Scenario", StringValue(scenario_str));

    // Create Channel Condition Model
    m_specLossModel->SetChannelModelAttribute ("ChannelConditionModel", PointerValue(m_condModel));
    m_propLossModel->SetChannelConditionModel (m_condModel);

    return true;
}


static bool createAntennaArrays(double frequency, double apertureAntennaEfficiency, double apertureAntennaDiameter )
{
    double apertureAntennaGainDb = 10 * log10( apertureAntennaEfficiency * pow( (M_PI*apertureAntennaDiameter*frequency)/299792458 , 2 ));

    // Groundstation Antenna
    gsAntenna = CreateObjectWithAttributes<UniformPlanarArray> (
        "NumColumns", UintegerValue(2),
        "NumRows", UintegerValue(2),
        "AntennaElement",
        PointerValue (
            CreateObjectWithAttributes <IsotropicAntennaModel> (
                "Gain", DoubleValue (-3)
            )
        )
    );


    // Std Satellite Antenna
    satAntenna = CreateObjectWithAttributes<UniformPlanarArray> (
        "NumColumns", UintegerValue(1),
        "NumRows", UintegerValue(1),
        "AntennaElement",
        PointerValue(
            CreateObjectWithAttributes <CosineAntennaModel> (
                "MaxGain", DoubleValue(apertureAntennaGainDb),
                "VerticalBeamwidth", DoubleValue(120.0),
                "HorizontalBeamwidth", DoubleValue(0.0),
                "Orientation", DoubleValue(180)
            )
        )
    );

    return true;
}

/*
Ptr<SpectrumValue> CreateTxPowerSpectralDensity(double fc, double noiseFigure, double BW, double RB_WIDTH)
{
    unsigned int numRbs = std::floor(BW / RB_WIDTH);
    double f = fc - ( numRbs * RB_WIDTH / 2.0 );
    double noiseFigureDb = noiseFigure;
    double RB_HALF_WIDTH = RB_WIDTH / 2;


    Bands rbs;
    std::vector<int> rbsId;

    for (uint32_t numrb = 0; numrb < numRbs; ++numrb)
    {
        BandInfo rb;
        rb.fl = f;
        f += RB_HALF_WIDTH;
        rb.fc = f;
        f += RB_HALF_WIDTH;
        rb.fh = f;

        rbs.push_back(rb);
        rbsId.push_back(numrb);

    }

    Ptr<SpectrumModel> model = Create<SpectrumModel> (rbs);
    Ptr<SpectrumValue> txPsd = Create<SpectrumValue> (model);

    const double kT_dBm_Hz = -174.0;
    double kT_W_Hz = std::pow ( 10.0, (kT_dBm_Hz - 30) / 10.0 );
    double noiseFigureLin = std::pow( 10.0, noiseFigureDb / 10.0 );

    double noisePwrSpectralDensity = kT_W_Hz * noiseFigureLin;

    for ( auto rbId : rbsId )
    {
        (*txPsd)[rbId] = noisePwrSpectralDensity;
    }

    return txPsd;
}
*/

Ptr<SpectrumValue> CreateTxPowerSpectralDensity(double fc, double pwr, double BW, double RB_WIDTH)
{

    // Number of Bands
    unsigned int numRbs = std::floor(BW / RB_WIDTH);

    // Highest Band Center Frequency
    double f = fc - ( numRbs * RB_WIDTH / 2.0 );

    // Transmit Power
    double powerTx = pwr;

    Bands rbs;
    std::vector<int> rbsId;

    for (uint32_t numrb = 0; numrb < numRbs; numrb++ )
    {
        BandInfo rb;

        // Set lower Boundary 
        rb.fl = f;

        // Set Center Frequncy
        f += RB_WIDTH / 2;
        rb.fc = f;

        // Set upper Boundary
        f += RB_WIDTH / 2;
        rb.fh = f;

        rbs.push_back(rb);
        rbsId.push_back(numrb);
    }

    Ptr<SpectrumModel> model = Create<SpectrumModel> (rbs);
    Ptr<SpectrumValue> txPsd = Create<SpectrumValue> (model);

    double powerTxW = std::pow(10.0, ( powerTx - 30 ) / 10 );
    double txPowerDensity = ( powerTxW / BW );

    for (auto rbId : rbsId )
    {
        (*txPsd)[rbId] = txPowerDensity;
    }

    return txPsd;
}


Ptr<SpectrumValue> CreateNoisePowerSpectralDensity(double fc, double noiseFigure, double BW, double RB_WIDTH)
{
    unsigned int numRbs = std::floor(BW / RB_WIDTH);
    double f = fc - (numRbs * RB_WIDTH / 2.0);
    double noiseFigureDb = noiseFigure; // dB noise Figure

    Bands rbs;
    std::vector<int> rbsId;

    for (uint32_t numrb = 0; numrb < numRbs; ++numrb)
    {
        BandInfo rb;
        rb.fl = f;
        f += RB_WIDTH / 2;
        rb.fc = f;
        f += RB_WIDTH / 2;
        rb.fh = f;

        rbs.push_back(rb);
        rbsId.push_back(numrb);
    }

    Ptr<SpectrumModel> model = Create<SpectrumModel>(rbs);
    Ptr<SpectrumValue> txPsd = Create<SpectrumValue>(model);

    const double kT_dBm_Hz = -174.0;  // dBm/Hz
    double kT_W_Hz = std::pow (10.0, (kT_dBm_Hz - 30) / 10.0); // W/Hz
    double noiseFigureLinear = std::pow (10.0, noiseFigureDb / 10.0);

    double noisePowerSpectralDensity =  kT_W_Hz * noiseFigureLinear;

    for (auto rbId : rbsId)
    {
        (*txPsd)[rbId] = noisePowerSpectralDensity;
    }

    return txPsd;
}

static void DoBeamforming(Ptr<MobilityModel> thisDevMob, Ptr<PhasedArrayModel> thisAntenna, Ptr<MobilityModel> otherDevMob)
{

    PhasedArrayModel::ComplexVector antennaWeights;

    // Get Device Positions
    Vector pos_a = thisDevMob->GetPosition();
    Vector pos_b = otherDevMob->GetPosition();

    // Compute Angles
    Angles complAngle (pos_b, pos_a);

    double angle_h = complAngle.GetAzimuth();
    double angle_v = complAngle.GetInclination();

    //std::cout << "Azimuth: " << angle_h << " Inclination: " << angle_v << std::endl;

    // The total Power is devided equally among the antenna elements
    int totArrayElements = thisAntenna->GetNumberOfElements();
    double power = 1 / sqrt(totArrayElements);


    for (int ind = 0; ind < totArrayElements; ind++)
    {
    
        Vector loc = thisAntenna->GetElementLocation(ind);
        double phase = -2 * M_PI * (sin(angle_v) * cos(angle_h) * loc.x
                        + sin(angle_v) * sin(angle_h) * loc.y
                        + cos(angle_v) * loc.z);

        antennaWeights.push_back( exp( std::complex<double> (0, phase) ) * power );

    }

    thisAntenna->SetBeamformingVector(antennaWeights);

}


double getInclinationAngle(Ptr<const GeocentricConstantPositionMobilityModel> gs, Ptr<const GeocentricConstantPositionMobilityModel> sat)
{
  // NS_ASSERT_MSG(gs.z < 8000, "Altitude of the ground terminal needs to be lower than 8km");
  // NS_ASSERT_MSG(other->GetGeographicPosition().z >= 8000, "Altitude of the HAPS/Satellite needs to be higher than 8km");
  // NS_ASSERT_MSG(gs.z < sat->GetGeographicPosition().z , "Altitude of the argoument node needs to be higher than object node");

  double elev_angle = 0;
  Vector a = gs->GetGeocentricPosition();
  Vector b = sat->GetGeocentricPosition();
  
  double numerator = (a.x*(b.x-a.x)+a.y*(b.y-a.y)+a.z*(b.z-a.z));
  double denominator = sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) * sqrt(pow(b.x-a.x, 2) + pow(b.y-a.y, 2) + pow(b.z-a.z, 2));

  double x = numerator / denominator;

  //This is done to avoid the nan returned by the asin function when numbers are "almost" 1, 
  //for example 1.0000000000000002
  if(x > 1)
  {
    x = 1;
  }

  elev_angle = ((180.0 * M_1_PI) * asin(x)); //asin returns radiants, we convert to degrees

  NS_ASSERT_MSG(!(isnan(elev_angle)), "asin returned a NaN value");

  return elev_angle;
}


static void ComputeSnr(SimulationParams &params)
{

    NS_ASSERT_MSG(params.txAntenna, "params.txAntenna is nullptr!");
    NS_ASSERT_MSG(params.rxAntenna, "params.rxAntenna is nullptr!");

    GeoCoordinate gcoTx (params.sat->GetPosition());
    GeoCoordinate gcoRx (params.gs->GetPosition());

    Ptr<GeocentricConstantPositionMobilityModel> geoTx = CreateObject<GeocentricConstantPositionMobilityModel> ();
    geoTx->SetGeographicPosition(Vector(
       gcoTx.GetLatitude(), gcoTx.GetLongitude(), gcoTx.GetAltitude()
    ));

    Ptr<GeocentricConstantPositionMobilityModel> geoRx = CreateObject<GeocentricConstantPositionMobilityModel> ();
    geoRx->SetGeographicPosition(Vector(
        gcoRx.GetLatitude(), gcoRx.GetLongitude(), gcoRx.GetAltitude()
    ));


    double inclination = getInclinationAngle(geoRx, geoTx);
    double distance = geoRx->GetDistanceFrom(geoTx);


    NodeContainer nodes;
    nodes.Create(2);

    Ptr<SimpleNetDevice> gsDev = CreateObject<SimpleNetDevice>();
    Ptr<SimpleNetDevice> satDev = CreateObject<SimpleNetDevice>();    
    nodes.Get(0)->AddDevice(gsDev);
    gsDev->SetNode(nodes.Get(0));
    nodes.Get(1)->AddDevice(satDev);
    satDev->SetNode(nodes.Get(1));

    nodes.Get(0)->AggregateObject (geoRx);
    nodes.Get(1)->AggregateObject (geoTx);


    //DoBeamforming (geoTx, params.txAntenna, geoRx);
    //DoBeamforming (geoRx, params.rxAntenna, geoTx);
    //m_condModel->SetAttribute()
    
    //if (inclination > 0)
    //    m_condModel->SetLosCondition(ChannelCondition::NLOS);


    Ptr<ChannelCondition> cond = m_condModel->GetChannelCondition(geoTx, geoRx);

    // Create Noise PSD
    Ptr<SpectrumValue> txPsd = CreateTxPowerSpectralDensity(
        params.frequency, params.txPow, params.bandwidth, params.resourceBlockBandwidth
    );
    Ptr<SpectrumValue> rxPsd = txPsd->Copy();

    Ptr<SpectrumValue> noisePsd = CreateNoisePowerSpectralDensity(
        params.frequency, params.noiseFigure, params.bandwidth, params.resourceBlockBandwidth
    );


    // apply the Path Loss
    double propagationGainDb = m_propLossModel->CalcRxPower(0, geoTx, geoRx);
    double propagationGainLin = std::pow( 10.0, (propagationGainDb) / 10.0 );

    // std::cout << propagationGainLin << std::endl;

    *(rxPsd) *= propagationGainLin;

    // Apply fast fading and beamforming gain
    rxPsd = m_specLossModel->CalcRxPowerSpectralDensity (rxPsd, geoTx, geoRx, params.txAntenna, params.rxAntenna);

    // Create Noise PSD
    const double kT_dBm_Hz = -174.0; // dBm/Hz
    double kT_W_Hz = std::pow (10.0, (kT_dBm_Hz - 30) / 10.0);
    double noiseFigureLinear = std::pow (10.0, params.noiseFigure / 10.0);
    double noisePowerSpectralDensity =  kT_W_Hz * noiseFigureLinear;
    //Ptr<SpectrumValue> noisePsd = Create <SpectrumValue> (txPsd->GetSpectrumModel ());

    (*noisePsd) = noisePowerSpectralDensity;

    result_file << Simulator::Now ().GetSeconds () << " " // time [s]
    /*
        << geoTx->GetPosition ().x << " "
        << geoTx->GetPosition ().y << " "
        << geoRx->GetPosition ().x << " "
        << geoRx->GetPosition ().y << " "
    */
        << gcoTx.GetLatitude() << " "
        << gcoTx.GetLongitude() << " "
        << gcoRx.GetLatitude() << " "
        << gcoRx.GetLongitude() << " "
        << inclination << " "
        << distance << " "
        << (cond->GetLosCondition() == ChannelCondition::LOS ? 1 : 0) << " " // channel state
        << Sum(*rxPsd) << " "
        << Sum(*noisePsd) << " "
        << 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd)) << " " // SNR [dB]
        << -propagationGainDb << std::endl; // pathloss [dB]

}




int main ( int argc, char *argv[] )
{

    uint32_t simTime = 60 * 60 * 6; // Seconds
    uint32_t timeRes = 1;    // Seconds

    double frequency = 2e9;
    prop_env_scenario scenario = NTN_RURAL;

    double satEIRPDensity   = 34;   // [ dBW/MHz ]
    double txPower          = 0.0;
    double noiseFigure      = 9.0;  // [ dB ]

    double bandwidth = 30e6;
    double RB_bandwidth = 60e3;
    
    double apertureAntennaEfficiency = 0.6;
    double apertureAntennaDiameter = 2;

    createPropagationEnvironment(frequency, scenario);
    createAntennaArrays(frequency, apertureAntennaEfficiency, apertureAntennaDiameter);

    double apertureAntennaGainDb = 10 * log10( apertureAntennaEfficiency * pow( (M_PI*apertureAntennaDiameter*frequency)/299792458 , 2 ));
    txPower = (satEIRPDensity + 10 * log10( bandwidth/1e6 ) - apertureAntennaGainDb) + 30;
    
    Config::SetDefault("ns3::SatSGP4MobilityModel::StartDateStr", StringValue("2023-03-20 21:24:00"));
    Config::SetDefault ("ns3::ThreeGppChannelModel::UpdatePeriod", TimeValue(MilliSeconds (1))); // update the channel at each iteration
    Config::SetDefault ("ns3::ThreeGppChannelConditionModel::UpdatePeriod", TimeValue(MilliSeconds (0.0))); // do not update the channel condition

    Ptr<GroundstationHelper> gsh = CreateObjectWithAttributes<GroundstationHelper> (
        "Latitude", DoubleValue(53.073635),
        "Longitude", DoubleValue(8.806422)
    );


    Ptr<WalkerConstellationHelper> constellation = CreateObjectWithAttributes<WalkerConstellationHelper> (
        "Inclination", DoubleValue(66.6),
        "Altitude", DoubleValue(850.0),
        "NumOfOrbits", IntegerValue(6),
        "NumOfSats", IntegerValue(9)
    );

    constellation->Initialize();
    constellation->LogInitialPositions("./output/constellation", ".txt");


    Ptr<MobilityModel> satMob; //CreateObject<WaypointMobilityModel> ();
    //satMob->GetObject<WaypointMobilityModel>()->AddWaypoint( Waypoint( Seconds(0), Vector(5.0, 3.0, 10.0)));

    Ptr<MobilityModel> gsMob; //CreateObject<WaypointMobilityModel> ();
    //gsMob->GetObject<WaypointMobilityModel>()->AddWaypoint( Waypoint ( Seconds(0), Vector(0.0, 0.0, 0.0)));

    satMob = constellation->getSatellite(9*3 + 5);
    gsMob = gsh->GetMobilityModel();


    std::cout << "Do Beamforming... " << std::endl;
    DoBeamforming(satMob, satAntenna, gsMob);
    DoBeamforming(gsMob, gsAntenna, satMob);

    std::cout << "Run Simulation - Time (secs): " << simTime << " Resolution (secs): " << timeRes << std::endl; 

    SimulationParams prms(satMob, gsMob, txPower, noiseFigure, satAntenna, gsAntenna, frequency, bandwidth, RB_bandwidth);

    for (long t = 0; t < floor(simTime / timeRes); t++)
    {
        Simulator::Schedule (
            Seconds(timeRes * t),
            &ComputeSnr,
            prms
        );
    }

    if (scenario == NTN_RURAL) {
        result_file.open ("./output/mpj-simulation-output-rural.txt", std::ios::out);
    }
    else if(scenario == NTN_SUBURBAN)
    {
        result_file.open ("./output/mpj-simulation-output-suburban.txt", std::ios::out);
    }
    else if(scenario == NTN_DENSEURBAN)
    {
        result_file.open ("./output/mpj-simulation-output-denseurban.txt", std::ios::out);
    }
    else if (scenario == NTN_URBAN)
    {
        result_file.open ("./output/mpj-simulation-output-urban.txt", std::ios::out);
    }
    else
    {
        std::cerr << "No Valid Scenario Selected!" << std::endl;
    }

    Simulator::Run();
    Simulator::Destroy();

    result_file.close();
}