// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_SLOWNUCLEONMODEL
#define O2_SLOWNUCLEONMODEL

#include "TObject.h"
class SlowNucleonModel : public TObject {
public:
  SlowNucleonModel() { ; }
  virtual ~SlowNucleonModel() { ; }
  virtual void GetNumberOfSlowNucleons(Int_t /*ncoll*/, Int_t & /*ngp*/,
                                       Int_t & /*ngn*/, Int_t & /*nbp*/,
                                       Int_t & /*nbn*/) const {
    ;
  }
  virtual void GetNumberOfSlowNucleons2(Int_t /*ncoll*/, Int_t & /*ngp*/,
                                        Int_t & /*ngn*/, Int_t & /*nbp*/,
                                        Int_t & /*nbn*/) const {
    ;
  }
  virtual void GetNumberOfSlowNucleons2s(Int_t /*ncoll*/, Int_t & /*ngp*/,
                                         Int_t & /*ngn*/, Int_t & /*nbp*/,
                                         Int_t & /*nbn*/) const {
    ;
  }

protected:
  ClassDef(SlowNucleonModel, 1) // Gray Particle Model
};
#endif
