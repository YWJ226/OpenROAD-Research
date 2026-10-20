// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <memory>

#include "dcalc/DmpCeff.hh"
#include "LocalElmoreCeffDelayCalc.hh"
#include "lgista/RcProbe.hh"

namespace lgista {

// Explicit probe-enabled twin of the ordinary DMP/Elmore calculator. The
// base gateDelay owns the DMP excitation and the real-load calculation.
// Passive probes are sampled with that SAME excitation before the next call.
class ObservedDmpDelayCalc : public sta::DmpCeffDelayCalc,
                             public RcProbeCalculator {
 public:
  explicit ObservedDmpDelayCalc(sta::StaState* state);
  sta::ArcDelayCalc* copy() override;
  void copyState(const sta::StaState* state) override;
  std::string_view name() const override { return "lgista_observed_dmp"; }
  sta::ArcDcalcResult gateDelay(
      const sta::Pin* driver, const sta::TimingArc* arc,
      const sta::Slew& in_slew, float load_cap,
      const sta::Parasitic* parasitic,
      const sta::LoadPinIndexMap& loads, const sta::Scene* scene,
      const sta::MinMax* min_max) override;
  sta::ArcDcalcResult inputPortDelay(
      const sta::Pin* port, float slew, const sta::RiseFall* rf,
      const sta::Parasitic* parasitic,
      const sta::LoadPinIndexMap& loads, const sta::Scene* scene,
      const sta::MinMax* min_max) override;
  const std::vector<RcProbeSample>& rcProbeSamples() const override
  { return samples_; }

 protected:
  void loadDelaySlew(const sta::Pin* load, double driver_slew,
                    const sta::RiseFall* rf,
                    const sta::LibertyLibrary* driver_library,
                    const sta::Parasitic* parasitic,
                    double& delay, double& slew) override;

 private:
  std::unique_ptr<sta::ArcDelayCalc> input_calculator_;
  std::vector<RcProbeSample> samples_;
};

// The same passive interface for the opt-in full-tree ElmoreCeff model.
// The parent computes all real loads unchanged; probes read its cached
// moments with the final, just-computed driver slew.
class ObservedElmoreCeffDelayCalc : public LocalElmoreCeffDelayCalc,
                                    public RcProbeCalculator {
 public:
  explicit ObservedElmoreCeffDelayCalc(sta::StaState* state);
  sta::ArcDelayCalc* copy() override;
  std::string_view name() const override { return "lgista_observed_elmore_ceff"; }
  sta::ArcDcalcResult gateDelay(
      const sta::Pin* driver, const sta::TimingArc* arc,
      const sta::Slew& in_slew, float load_cap,
      const sta::Parasitic* parasitic,
      const sta::LoadPinIndexMap& loads, const sta::Scene* scene,
      const sta::MinMax* min_max) override;
  sta::ArcDcalcResult inputPortDelay(
      const sta::Pin* port, float slew, const sta::RiseFall* rf,
      const sta::Parasitic* parasitic,
      const sta::LoadPinIndexMap& loads, const sta::Scene* scene,
      const sta::MinMax* min_max) override;
  const std::vector<RcProbeSample>& rcProbeSamples() const override
  { return samples_; }

 private:
  std::vector<RcProbeSample> samples_;
};

}  // namespace lgista
