// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "lgista/RcProbe.hh"

#include "LocalDmpDelayCalc.hh"
#include "ObservedDmpDelayCalc.hh"

namespace lgista {

sta::ArcDelayCalc* makeProbeDelayCalc(sta::StaState* state)
{
  if (useElmoreCeff()) return new ObservedElmoreCeffDelayCalc(state);
  return new ObservedDmpDelayCalc(state);
}

}  // namespace lgista
