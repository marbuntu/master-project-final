/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2006,2007 INRIA
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
#include "geocentric-constant-position-mobility-model.h"
#include <math.h>

namespace ns3 {

NS_OBJECT_ENSURE_REGISTERED (GeocentricConstantPositionMobilityModel);

TypeId
GeocentricConstantPositionMobilityModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::GeocentricConstantPositionMobilityModel")
    .SetParent<MobilityModel> ()
    .SetGroupName ("Mobility")
    .AddConstructor<GeocentricConstantPositionMobilityModel> ()
  ;
  return tid;
}

GeocentricConstantPositionMobilityModel::GeocentricConstantPositionMobilityModel ()
{
}
GeocentricConstantPositionMobilityModel::~GeocentricConstantPositionMobilityModel ()
{
}



Vector GeocentricConstantPositionMobilityModel::GetGeographicPosition (void) const
{
  return DoGetGeographicPosition();
}

void GeocentricConstantPositionMobilityModel::SetGeographicPosition (const Vector &position)
{
  DoSetGeographicPosition(position);
}

Vector GeocentricConstantPositionMobilityModel::GetGeocentricPosition (void) const
{
  return DoGetGeocentricPosition();
}

void GeocentricConstantPositionMobilityModel::SetGeocentricPosition (const Vector &position)
{
  DoSetGeocentricPosition(position);
}

double 
GeocentricConstantPositionMobilityModel::GetElevationAngle(Ptr<const GeocentricConstantPositionMobilityModel> other)
{
  return DoGetElevationAngle(other);
}

void GeocentricConstantPositionMobilityModel::SetCoordinateTranslationReferencePoint (const Vector &position)
{
  DoSetCoordinateTranslationReferencePoint(position);
}

Vector 
GeocentricConstantPositionMobilityModel::GetCoordinateTranslationReferencePoint(void) const
{
  return DoGetCoordinateTranslationReferencePoint();
}

/** 
 *  In order to for returned the position to work with the rest of ns-3,
 *  it is converted from geographic to topocentric using the GeograpicPositions method.
 *  The default reference point for conversion is lat=0, lon=0, altitude=0.
 */
Vector
GeocentricConstantPositionMobilityModel::DoGetPosition (void) const
{
  GeographicPositions gp;
  Vector topographicCoordinates = gp.GeographicToTopocentricCoordinates(m_position,m_geographicReferencePoint,GeographicPositions::SPHERE);
  return topographicCoordinates;
}

void
GeocentricConstantPositionMobilityModel::DoSetPosition (const Vector &position)
{
  NS_FATAL_ERROR("Legacy SetPosition has not been implemented for GeocentricConstantPostionMobiltyModel");
  //NEED CONVERSION FROM PLANAR CARTESIAN
  NotifyCourseChange ();
}

double 
GeocentricConstantPositionMobilityModel::DoGetDistanceFrom (Ptr<const GeocentricConstantPositionMobilityModel> other) const
{
  GeographicPositions gp;

  Vector cartesian_coord1 = gp.GeographicToCartesianCoordinates(m_position.x,m_position.y,m_position.z,gp.EarthSpheroidType::SPHERE);
  Vector cartesian_coord2 = other->DoGetGeocentricPosition();

  double distance = sqrt(pow(cartesian_coord1.x - cartesian_coord2.x,2) + pow(cartesian_coord1.y - cartesian_coord2.y,2) + pow(cartesian_coord1.z - cartesian_coord2.z,2));

  return distance;
}

Vector
GeocentricConstantPositionMobilityModel::DoGetGeographicPosition (void) const
{
  return m_position;
}

void
GeocentricConstantPositionMobilityModel::DoSetGeographicPosition (const Vector &position)
{
  m_position = position;
  NotifyCourseChange ();
}

Vector
GeocentricConstantPositionMobilityModel::DoGetGeocentricPosition (void) const
{
  GeographicPositions gp;
  Vector geocentric_coord = gp.GeographicToCartesianCoordinates(m_position.x,m_position.y,m_position.z,gp.EarthSpheroidType::SPHERE);
  return geocentric_coord;
}

void
GeocentricConstantPositionMobilityModel::DoSetGeocentricPosition (const Vector &position)
{
  GeographicPositions gp;
  Vector cartesian_coord = gp.CartesianToGeographicCoordinates(position,gp.SPHERE);
  m_position = cartesian_coord;
  NotifyCourseChange ();
}

double 
GeocentricConstantPositionMobilityModel::DoGetElevationAngle(Ptr<const GeocentricConstantPositionMobilityModel> other)
{
  NS_ASSERT_MSG(m_position.z < 8000, "Altitude of the ground terminal needs to be lower than 8km");
  //NS_ASSERT_MSG(other->GetGeographicPosition().z >= 8000, "Altitude of the HAPS/Satellite needs to be higher than 8km");
  NS_ASSERT_MSG(m_position.z < other->DoGetGeographicPosition().z , "Altitude of the argoument node needs to be higher than object node");

  double elev_angle = 0;
  Vector a = this->DoGetGeocentricPosition();
  Vector b = other->DoGetGeocentricPosition();
  
  double numerator = abs(a.x*(b.x-a.x)+a.y*(b.y-a.y)+a.z*(b.z-a.z));
  double denominator = sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) * sqrt(pow(b.x-a.x, 2) + pow(b.y-a.y, 2) + pow(b.z-a.z, 2));

  double x = numerator/denominator;

  //This is done to avoid the nan returned by the asin function when numbers are "almost" 1, 
  //for example 1.0000000000000002
  if(x > 1)
  {
    x = 1;
  }

  elev_angle = abs((180.0 * M_1_PI) * asin(x)); //asin returns radiants, we convert to degrees

  NS_ASSERT_MSG(!(isnan(elev_angle)), "asin returned a NaN value");

  return elev_angle;
}

void
GeocentricConstantPositionMobilityModel::DoSetCoordinateTranslationReferencePoint(const Vector &refPoint)
{
  m_geographicReferencePoint = refPoint;
}

Vector
GeocentricConstantPositionMobilityModel::DoGetCoordinateTranslationReferencePoint(void) const
{
  return m_geographicReferencePoint;
}

Vector
GeocentricConstantPositionMobilityModel::DoGetVelocity (void) const
{
  return Vector (0.0, 0.0, 0.0);
}

} // namespace ns3
