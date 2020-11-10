#ifndef GENERATORPARAMLIBBASE_H
#define GENERATORPARAMLIBBASE_H
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Baae class for parameterisations used by GeneratorParam
// andreas.morsch@cern.ch

#include <TObject.h>

class TRandom;

class GeneratorParamLibBase : public TObject {
public:
  //
  virtual ~GeneratorParamLibBase() {}
  typedef Double_t (*GenFunc)(const Double_t *, const Double_t *);
  typedef Int_t (*GenFuncIp)(TRandom *);
  virtual GenFunc GetPt(Int_t param, const char *tname) const = 0;
  virtual GenFunc GetY(Int_t param, const char *tname) const = 0;
  virtual GenFuncIp GetIp(Int_t param, const char *tname) const = 0;
  virtual GenFunc GetV2(Int_t, const char *) const { return NoV2; }
  static Double_t NoV2(const Double_t *, const Double_t *) { return 0; }
  ClassDef(GeneratorParamLibBase,
           0) // Library providing y and pT parameterisations
};
#endif
