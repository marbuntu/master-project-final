/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011 CTTC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Mattia Sandri
 */

#include <ns3/log.h>
#include <ns3/double.h>
#include <math.h>
#include <iostream>
#include <ns3/math.h>
#include <cmath>
#include "antenna-model.h"
#include "circular-aperture-antenna-model.h"


namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("CircularApertureAntennaModel");

NS_OBJECT_ENSURE_REGISTERED (CircularApertureAntennaModel);

TypeId
CircularApertureAntennaModel::GetTypeId ()
{
  static TypeId tid = TypeId ("ns3::CircularApertureAntennaModel")
    .SetParent<AntennaModel> ()
    .SetGroupName ("Antenna")
    .AddConstructor<CircularApertureAntennaModel> ()
    .AddAttribute ("AntennaMaxGainDb",
                   "The maximum gain value in dB of the antenna",
                   DoubleValue (1),
                   MakeDoubleAccessor (&CircularApertureAntennaModel::SetMaxGain),
                   MakeDoubleChecker<double> (0.0))
    .AddAttribute ("AntennaCircularApertureRadius",
                   "The radius of the aperture of the antenna, in meters",
                   DoubleValue (0.5),
                   MakeDoubleAccessor (&CircularApertureAntennaModel::SetApertureRadius),
                   MakeDoubleChecker<double> (0.0))
    .AddAttribute ("OperatingFrequency",
                   "The operating frequency of the antenna",
                   DoubleValue (2e9),
                   MakeDoubleAccessor (&CircularApertureAntennaModel::SetOperatingFrequency),
                   MakeDoubleChecker<double> (0.0))
    .AddAttribute ("AntennaInclination",
                   "The inclination angle in rad of the antenna",
                   DoubleValue (0.0),
                   MakeDoubleAccessor (&CircularApertureAntennaModel::SetInclination),
                   MakeDoubleChecker<double> (0.0,M_PI))
  ;
  return tid;
}

CircularApertureAntennaModel::CircularApertureAntennaModel ()
    : AntennaModel ()
{}

CircularApertureAntennaModel::~CircularApertureAntennaModel ()
{}

void
CircularApertureAntennaModel::SetOrientation (Angles a)
{
  NS_LOG_FUNCTION (this << a);
  m_antennaOrientation = a;
}

void
CircularApertureAntennaModel::SetInclination (double theta)
{
  NS_LOG_FUNCTION (this << theta);
  NS_ASSERT_MSG( (0 <= theta && theta <= M_PI) , "Setting invalid inclination(deg): " << RadiansToDegrees(theta));
  m_antennaOrientation = Angles(m_antennaOrientation.GetAzimuth(),theta);
}

void
CircularApertureAntennaModel::SetApertureRadius (double r)
{
  NS_LOG_FUNCTION (this << r);
  m_apertureRadius = r;
}

void
CircularApertureAntennaModel::SetOperatingFrequency (double f)
{
  NS_LOG_FUNCTION (this << f);
  m_operatingFrequency = f;
}

void
CircularApertureAntennaModel::SetMaxGain(double gain)
{
  NS_LOG_FUNCTION (this << gain);
  m_maxGain = gain;
}

double
CircularApertureAntennaModel::GetGainDb (Angles a)
{
  NS_LOG_FUNCTION (this << a);

  double theta = abs( m_antennaOrientation.GetInclination() - a.GetInclination() );  //The angle from the bore sight antenna main beam

  double gain = 0;

  if(theta == 0)
  {
    gain = m_maxGain;
  }
  else //theta =! 0
  {
    double k = ( 2 * M_PI * m_operatingFrequency ) / 299792458;
    double J_1 = std::cyl_bessel_j( 1 , ( k * m_apertureRadius * sin(theta) ) );
    double denominator = k * m_apertureRadius * sin(theta);
    gain = 4 * pow( abs((J_1/denominator)) , 2 );
    gain = 10 * log10(gain) + m_maxGain;
  }

  /*
  if((M_PI_2) <= theta || theta <= (-M_PI_2))
  {
    gain = -200;
  }*/

  return gain;
  }

}
