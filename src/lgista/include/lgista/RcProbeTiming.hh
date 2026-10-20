// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "lgista/RcProbe.hh"
#include "sta/GraphClass.hh"

namespace lgista {

struct RcProbeSinkDelay {
  const sta::Pin* sink = nullptr;
  double wire_delay = std::numeric_limits<double>::quiet_NaN();
};

// A single non-perturbed gate calculation. Sink delays retain their real
// receiver thresholds; probe samples use the driver's output thresholds.
// No value in this record is a max/min over different gate calculations.
struct RcProbeArcContext {
  sta::EdgeId edge_id = 0;
  sta::VertexId driver_id = 0;
  const sta::TimingArc* arc = nullptr;
  int ap_index = -1;
  double gate_delay = std::numeric_limits<double>::quiet_NaN();
  std::vector<RcProbeSample> probes;
  std::vector<RcProbeSinkDelay> sinks;
};

// One matched source path and gate context, never a subtraction of separately
// selected worst paths. Required/slack may be unavailable on unconstrained
// paths; arrival_valid is independent of required_valid. Times are seconds.
// min_max_index: 0=early/min, 1=late/max; rf_index: 0=rise, 1=fall.
struct RcProbeTiming {
  RcProbeId id = 0;
  size_t context_index = 0;
  int ap_index = -1;
  int rf_index = -1;
  int min_max_index = -1;
  uint32_t source_tag_index = 0;
  // Before/after the existing STA TagGroup matching (which can coalesce
  // equivalent propagated tags). Exception states remain part of that match.
  uint32_t propagated_tag_index = 0;
  uint32_t driver_tag_index = 0;
  uint32_t required_tag_index = 0;
  const sta::Pin* required_sink = nullptr;
  double driver_arrival = std::numeric_limits<double>::quiet_NaN();
  double wire_delay = std::numeric_limits<double>::quiet_NaN();
  double arrival = std::numeric_limits<double>::quiet_NaN();
  double required = std::numeric_limits<double>::quiet_NaN();
  double slack = std::numeric_limits<double>::quiet_NaN();
  double slew = std::numeric_limits<double>::quiet_NaN();
  // False for a calculated Liberty arc that carries no legal local data path.
  // True + !arrival_valid means a reachable path is missing an observation.
  bool path_matched = false;
  bool arrival_valid = false;
  bool required_valid = false;
  std::string error;
};

}  // namespace lgista
