// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace sta {
class ArcDelayCalc;
class MinMax;
class Net;
class Parasitic;
class ParasiticNode;
class Pin;
class RiseFall;
class Scene;
class StaState;
class TimingArc;
}

namespace lgista {

using RcProbeId = uint64_t;

// A passive observation point in one complete parasitic network. Bindings
// borrow network/node pointers and MUST NOT outlive that network. Probes add
// neither capacitance nor a logical pin/vertex. The same id may be bound in
// several scene networks; (network,node) identifies this particular binding.
struct RcProbeBinding {
  RcProbeId id = 0;
  const sta::Net* net = nullptr;
  const sta::Parasitic* network = nullptr;
  const sta::ParasiticNode* node = nullptr;
  const sta::Pin* driver = nullptr;
  std::vector<const sta::Pin*> sinks;
};

// Numeric, PtGraph-local reduction result. No borrowed RC node is retained.
// tree_node_idx is used only by the full-tree ElmoreCeff representation.
struct RcProbeMoment {
  RcProbeId id = 0;
  double elmore = 0.0;
  uint32_t tree_node_idx = std::numeric_limits<uint32_t>::max();
};

// Result of exactly one gateDelay call, in the DRIVER library's threshold
// convention. There is intentionally no load-pin threshold adjustment for
// an anonymous observation point. LocalSta must consume these before any
// subsequent gateDelay (including a perturbed sensitivity calculation).
struct RcProbeSample {
  RcProbeId id = 0;
  const sta::Pin* driver = nullptr;
  const sta::TimingArc* arc = nullptr;
  const sta::Scene* scene = nullptr;
  const sta::MinMax* min_max = nullptr;
  const sta::RiseFall* rf = nullptr;
  double wire_delay = std::numeric_limits<double>::quiet_NaN();
  double slew = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
  std::string error;
};

// Cross-cast from ArcDelayCalc to read the last call's private samples.
// Each worker/native calculator owns its own storage; copy() starts empty.
class RcProbeCalculator {
 public:
  virtual ~RcProbeCalculator() = default;
  virtual const std::vector<RcProbeSample>& rcProbeSamples() const = 0;
};

// Explicit opt-in factory. The original makeLocalDelayCalc is unchanged.
// Honors the same cached LRF_USE_ELMORECEFF model selection as reduction.
// Caller owns the returned calculator.
sta::ArcDelayCalc* makeProbeDelayCalc(sta::StaState* state);

}  // namespace lgista
