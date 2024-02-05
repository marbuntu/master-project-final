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

#ifndef VSAT_ANTENNA_MODEL_H
#define VSAT_ANTENNA_MODEL_H


#include <ns3/object.h>
#include <ns3/antenna-model.h>

namespace ns3 {

/**
 * \ingroup antenna
 * 
 * \brief  Antenna model based on a VSAT approximation of the main lobe radiation pattern.
 *
 * This class implements the VSAT model as described in some 3GPP document, e.g., R4-092042
 *
 * A similar model appears in:
 *
 * George Calcev and Matt Dillon, "Antenna Tilt Control in CDMA Networks"
 * in Proc. of the 2nd Annual International Wireless Internet Conference (WICON), 2006
 *
 *
 */
class VsatAntennaModel : public AntennaModel
{
public:

  /**
   * \brief Get the type ID.
   * \return The object TypeId.
   */
  static TypeId GetTypeId ();

  // inherited from AntennaModel
  virtual double GetGainDb (Angles a);


  // attribute getters/setters
  /**
   * Set the Beam width
   * \param beamwidthDegrees Beam width in degrees
   */
  void SetBeamwidth (double beamwidthDegrees);
  /**
   * Get the Beam width
   * \return beam width in degrees
   */
  double GetBeamwidth () const;
  /**
   * Set the antenna azimuth orientation
   * \param orientationDegrees antenna azimuth orientation in radiants
   */
  void SetAzimuth (double azimut);
  /**
   * Get the antenna azimuth orientation
   * \return antenna azimuth orientation in radiants
   */
  double GetAzimuth () const;
  /**
   * Set the antenna inclination orientation
   * \param orientationDegrees antenna inclination orientation in radiants
   */
  void SetInclination (double inclination);
  /**
   * Get the antenna inclination orientation
   * \return antenna inclination orientation in radiants
   */
  double GetInclination () const;

private:

  double m_beamwidthRadians; //!< Beam width in radians
  Angles m_orientation {0.0,0.0}; //!< Antenna orientation in radians
  double m_maxGain; //!< Max attenuation
};



} // namespace ns3


#endif // VSAT_ANTENNA_MODEL_H
