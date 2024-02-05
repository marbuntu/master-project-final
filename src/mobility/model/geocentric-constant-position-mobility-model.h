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
#ifndef GEOCENTRIC_CONSTANT_POSITION_MOBILITY_MODEL_H
#define GEOCENTRIC_CONSTANT_POSITION_MOBILITY_MODEL_H

#include "mobility-model.h"
#include "geographic-positions.h"

namespace ns3 {

/**
 * \ingroup mobility
 *
 * \brief Mobility model using geocentric euclidean coordinates, as defined in 38.811 chapter 6.3
 */
class GeocentricConstantPositionMobilityModel : public MobilityModel 
{
public:
  /**
   * Register this type with the TypeId system.
   * \return the object TypeId
   */
  static TypeId GetTypeId (void);
  /**
   * Create a position located at coordinates (0,0,0)
   */
  GeocentricConstantPositionMobilityModel ();
  virtual ~GeocentricConstantPositionMobilityModel ();

  virtual double GetElevationAngle(Ptr<const GeocentricConstantPositionMobilityModel> other);
  virtual Vector GetGeographicPosition (void) const;
  virtual void SetGeographicPosition (const Vector &position);
  virtual Vector GetGeocentricPosition (void) const;
  virtual void SetGeocentricPosition (const Vector &position);
  virtual void SetCoordinateTranslationReferencePoint(const Vector &refPoint);
  virtual Vector GetCoordinateTranslationReferencePoint(void) const;

private:
  virtual Vector DoGetPosition (void) const;
  virtual void DoSetPosition (const Vector &position);
  double DoGetDistanceFrom (Ptr<const GeocentricConstantPositionMobilityModel> position) const;
  /**
   * \brief Get the position using geographic (geodetic) coordinates
   * \return Vector containint (latitue, longitute, altitute)
   */
  virtual Vector DoGetGeographicPosition (void) const;
  /**
   * \brief Set the position using geographic coordinates
   * \param position pointer to a Vector containing (latitue, longitute, altitute)
   */
  virtual void DoSetGeographicPosition (const Vector &position);
  /**
   * \brief Get the position using Geocentric Cartesian coordinates
   * \return Vector containint (X, Y, Z)
   */
  virtual Vector DoGetGeocentricPosition (void) const;
  /**
   * \brief Set the position using Geocentric coordinates
   * \param position pointer to a Vector containing (X, Y, Z)
   */
  virtual void DoSetGeocentricPosition (const Vector &position);
  /**
   * \brief Computes elevation angle between a ground terminal and a HAPS/Satellite.
   * After calculating the plane perpendicular to one cartesian position vector,
   * the elevation angle is calculated using https://www.w3schools.blog/angle-between-a-line-and-a-plane.
   * The altitute of the node passed as paramenter needs to be higher.
   * \param other pointer to the HAPS/Satellite mobility model
   * \return the elevation angle in degrees
   */
  virtual double DoGetElevationAngle(Ptr<const GeocentricConstantPositionMobilityModel> other);
  /**
   * \brief Set the reference point for coordinate translation
   * \param refPoint vector containing the geographic reference point
   */
  virtual void DoSetCoordinateTranslationReferencePoint(const Vector &refPoint);
  /**
   * \brief Get the reference point for coordinate translation
   * \return Vector containing geographic reference point
   */
  virtual Vector DoGetCoordinateTranslationReferencePoint(void) const;



  virtual Vector DoGetVelocity (void) const;

  /**
   * the constant Geographic position, in degrees, in the order:
   * latitude
   * longitude
   * altitude
   */
  Vector m_position;
  Vector m_geographicReferencePoint {0,0,0}; //<! This is the point taken as reference when converting from geographic to topographic (aka planar Cartesian)
};

} // namespace ns3

#endif /* GEOCENTRIC_CONSTANT_POSITION_MOBILITY_MODEL_H */