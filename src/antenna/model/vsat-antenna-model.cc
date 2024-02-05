/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2012 CTTC
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

#include "antenna-model.h"
#include "vsat-antenna-model.h"


namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("VsatAntennaModel");

NS_OBJECT_ENSURE_REGISTERED (VsatAntennaModel);


TypeId 
VsatAntennaModel::GetTypeId ()
{
  static TypeId tid = TypeId ("ns3::VsatAntennaModel")
    .SetParent<AntennaModel> ()
    .SetGroupName("Antenna")
    .AddConstructor<VsatAntennaModel> ()
    .AddAttribute ("Beamwidth",
                   "The 3dB beamwidth (degrees)",
                   DoubleValue (60),
                   MakeDoubleAccessor (&VsatAntennaModel::SetBeamwidth,
                                       &VsatAntennaModel::GetBeamwidth),
                   MakeDoubleChecker<double> (0, 180))
    .AddAttribute ("AzimuthOrientation",
                   "The angle (rad) that expresses the azimuth orientation of the antenna",
                   DoubleValue (0.0),
                   MakeDoubleAccessor (&VsatAntennaModel::SetAzimuth),
                   MakeDoubleChecker<double> (-M_PI_2, M_PI_2))
    .AddAttribute ("InclinationOrientation",
                   "The angle (rad) that expresses the inclination orientation of the antenna. 0 is pointing the sky",
                   DoubleValue (0.0),
                   MakeDoubleAccessor (&VsatAntennaModel::SetInclination),
                   MakeDoubleChecker<double> (0, M_PI))
    .AddAttribute ("MaxGain",
                   "The maximum gain (dB) of the antenna radiation pattern.",
                   DoubleValue (0.0),
                   MakeDoubleAccessor (&VsatAntennaModel::m_maxGain),
                   MakeDoubleChecker<double> ())
  ;
  return tid;
}

void 
VsatAntennaModel::SetBeamwidth (double beamwidthDegrees)
{ 
  NS_LOG_FUNCTION (this << beamwidthDegrees);
  m_beamwidthRadians = DegreesToRadians (beamwidthDegrees);
}

double
VsatAntennaModel::GetBeamwidth () const
{
  return RadiansToDegrees (m_beamwidthRadians);
}

void 
VsatAntennaModel::SetAzimuth (double azimuth)
{
  NS_LOG_FUNCTION (this << azimuth);
  m_orientation = Angles(azimuth,m_orientation.GetInclination());
}

double
VsatAntennaModel::GetAzimuth () const
{
  return m_orientation.GetAzimuth();
}

void 
VsatAntennaModel::SetInclination (double inclination)
{
  NS_LOG_FUNCTION (this << inclination);
  m_orientation = Angles(m_orientation.GetAzimuth(),inclination);
}

double
VsatAntennaModel::GetInclination () const
{
  return m_orientation.GetInclination();
}

double 
VsatAntennaModel::GetGainDb (Angles a)
{
  NS_LOG_FUNCTION (this << a);
    
  double x = sin(m_orientation.GetInclination()) * cos(m_orientation.GetAzimuth());
  double y = sin(m_orientation.GetInclination()) * sin(m_orientation.GetAzimuth());
  double z = cos(m_orientation.GetInclination());
  Vector orientationVec = Vector(x,y,z);

  x = sin(a.GetInclination()) * cos(a.GetAzimuth());
  y = sin(a.GetInclination()) * sin(a.GetAzimuth());
  z = cos(a.GetInclination());
  Vector gainVec = Vector(x,y,z);

  double angleDifference = acos(orientationVec.x * gainVec.x + orientationVec.y * gainVec.y + orientationVec.z * gainVec.z);

  double gainDb = m_maxGain - ( 12 * ( pow (angleDifference / m_beamwidthRadians, 2)) );

  if((M_PI_2) <= angleDifference || angleDifference <= (-M_PI_2))
  {gainDb = -200;}

  NS_LOG_LOGIC ("gain = " << gainDb);
  return 39.7;
}


}

