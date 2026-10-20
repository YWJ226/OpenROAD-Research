// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#include "ObservedDmpDelayCalc.hh"

#include <cmath>
#include <utility>
#ifdef LGISTA_RC_PROBE_TESTING
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#endif

#include "LocalDmpDelayCalc.hh"
#include "PtElmoreCeff.hh"
#include "PtPiElmore.hh"
#include "dcalc/DmpDelayCalc.hh"
#include "parasitics/ConcreteParasiticsPvt.hh"
#include "sta/Delay.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "sta/Parasitics.hh"
#include "sta/Scene.hh"
#include "sta/TimingArc.hh"

namespace lgista {
namespace {
template <class T>
const T* localModel(const sta::Parasitic* parasitic)
{
  // Parasitic itself is not polymorphic; all local networks/models use
  // the ConcreteParasitic implementation, just as the ordinary EC calc.
  return parasitic ? dynamic_cast<const T*>(
      static_cast<const sta::ConcreteParasitic*>(parasitic)) : nullptr;
}

RcProbeSample sampleFor(RcProbeId id, const sta::Pin* driver,
                       const sta::TimingArc* arc, const sta::Scene* scene,
                       const sta::MinMax* min_max)
{
  RcProbeSample sample;
  sample.id = id;
  sample.driver = driver;
  sample.arc = arc;
  sample.scene = scene;
  sample.min_max = min_max;
  sample.rf = arc->toEdge()->asRiseFall();
  return sample;
}

void validateSample(RcProbeSample& sample)
{
  sample.valid = std::isfinite(sample.wire_delay)
      && std::isfinite(sample.slew)
      && sample.wire_delay >= 0.0 && sample.slew >= 0.0;
  if (!sample.valid) sample.error = "invalid_rc_probe_delay_or_slew";
}

#ifdef LGISTA_RC_PROBE_TESTING
// Built only into isolated validation binaries. An ordinary factory result
// is an independent reference for all genuine loads, including all four
// statistical Delay fields. A second test pass aliases each real load as a
// passive probe, then converts to that real pin's thresholds for comparison.
struct ProbeCheckCounts {
  std::atomic<unsigned long long> gates{0}, loads{0}, aliases{0}, failures{0};
  std::atomic<unsigned> ap_rf_mask{0};
};
ProbeCheckCounts probe_checks[2];

void printProbeChecks()
{
  for (int model = 0; model < 2; ++model) {
    const auto& c = probe_checks[model];
    std::fprintf(stderr,
        "[lgista_rc_probe_test] model=%s checked_gate_count=%llu real_pin_count=%llu sink_probe_alias_count=%llu ap_rf_mask=%u failures=%llu status=%s\n",
        model == 0 ? "dmp" : "elmore_ceff", c.gates.load(), c.loads.load(),
        c.aliases.load(), c.ap_rf_mask.load(), c.failures.load(),
        c.failures.load() == 0 ? "pass" : "fail");
  }
}

[[noreturn]] void failProbeTest(int model, const char* message)
{
  probe_checks[model].failures.fetch_add(1, std::memory_order_relaxed);
  throw std::runtime_error(message);
}

bool exactDelay(const sta::Delay& a, const sta::Delay& b)
{
  return a.mean() == b.mean() && a.meanShift() == b.meanShift()
      && a.stdDev2() == b.stdDev2() && a.skewness() == b.skewness();
}

sta::ArcDcalcResult checkOriginalResult(
    int model, sta::StaState* state, sta::ArcDcalcResult& observed,
    const sta::Pin* driver, const sta::TimingArc* arc, const sta::Slew& in_slew,
    float load_cap, const sta::Parasitic* parasitic,
    const sta::LoadPinIndexMap& loads, const sta::Scene* scene,
    const sta::MinMax* min_max)
{
  static std::once_flag once;
  std::call_once(once, [] { std::atexit(printProbeChecks); });
  std::unique_ptr<sta::ArcDelayCalc> original(makeLocalDelayCalc(state));
  auto reference = original->gateDelay(
      driver, arc, in_slew, load_cap, parasitic, loads, scene, min_max);
  if (!exactDelay(observed.gateDelay(), reference.gateDelay())
      || !exactDelay(observed.drvrSlew(), reference.drvrSlew()))
    failProbeTest(model, "rc_probe_test_changed_real_gate_result");
  for (const auto& [load, index] : loads) {
    (void) load;
    if (!exactDelay(observed.wireDelay(index), reference.wireDelay(index))
        || !exactDelay(observed.loadSlew(index), reference.loadSlew(index)))
      failProbeTest(model, "rc_probe_test_changed_real_load_result");
  }
  probe_checks[model].gates.fetch_add(1, std::memory_order_relaxed);
  probe_checks[model].loads.fetch_add(loads.size(), std::memory_order_relaxed);
  probe_checks[model].ap_rf_mask.fetch_or(
      1U << (min_max->index() * 2 + arc->toEdge()->asRiseFall()->index()),
      std::memory_order_relaxed);
  return reference;
}

void checkLoadAlias(int model, double delay, double slew,
                    const sta::ArcDcalcResult& reference, size_t index)
{
  // ArcDcalcResult stores float means, while DMP's helper returns double.
  if (static_cast<float>(delay) != reference.wireDelay(index).mean()
      || static_cast<float>(slew) != reference.loadSlew(index).mean())
    failProbeTest(model, "rc_probe_test_real_sink_alias_mismatch");
  probe_checks[model].aliases.fetch_add(1, std::memory_order_relaxed);
}
#endif
}  // namespace

ObservedDmpDelayCalc::ObservedDmpDelayCalc(sta::StaState* state)
    : sta::DmpCeffDelayCalc(state),
      input_calculator_(sta::makeDmpCeffElmoreDelayCalc(state))
{
}

sta::ArcDelayCalc* ObservedDmpDelayCalc::copy()
{
  // Construct fresh algorithm state; never copy a previous gate's waveform.
  return new ObservedDmpDelayCalc(this);
}

void ObservedDmpDelayCalc::copyState(const sta::StaState* state)
{
  sta::DmpCeffDelayCalc::copyState(state);
  input_calculator_->copyState(state);
  samples_.clear();
}

sta::ArcDcalcResult ObservedDmpDelayCalc::gateDelay(
    const sta::Pin* driver, const sta::TimingArc* arc,
    const sta::Slew& in_slew, float load_cap,
    const sta::Parasitic* parasitic, const sta::LoadPinIndexMap& loads,
    const sta::Scene* scene, const sta::MinMax* min_max)
{
  samples_.clear();
  const auto* pi = localModel<PtPiElmore>(parasitic);
  // Match the exact base-class branch that initializes a fresh DMP model.
  // The fallback branch must NEVER sample a previous call's dmp_alg_.
  const bool supported = arc->gateTableModel(scene, min_max) && pi;
  auto result = sta::DmpCeffDelayCalc::gateDelay(
      driver, arc, in_slew, load_cap, parasitic, loads, scene, min_max);
  if (!pi) return result;
  for (const auto& probe : pi->rcProbes()) {
    auto sample = sampleFor(probe.id, driver, arc, scene, min_max);
    if (!supported) {
      sample.error = "unsupported_dmp_rc_probe_model";
    } else if (!std::isfinite(probe.elmore) || probe.elmore < 0.0) {
      sample.error = "invalid_rc_probe_elmore";
    } else {
      // Both protected implementations ignore the pin argument. Use a
      // real driver anyway, and deliberately do NOT thresholdAdjust: an
      // anonymous probe uses the driver's library threshold convention.
      sample.wire_delay = probe.elmore;
      sample.slew = sta::delayAsFloat(result.drvrSlew());
      loadDelaySlewElmore(driver, probe.elmore, sample.wire_delay, sample.slew);
      validateSample(sample);
    }
    samples_.push_back(std::move(sample));
  }
#ifdef LGISTA_RC_PROBE_TESTING
  if (supported && !pi->rcProbes().empty()) {
    const auto reference = checkOriginalResult(
        0, this, result, driver, arc, in_slew, load_cap, parasitic,
        loads, scene, min_max);
    const auto* rf = arc->toEdge()->asRiseFall();
    const auto* library = arc->from()->libertyCell()->libertyLibrary();
    for (const auto& [load, index] : loads) {
      float elmore = 0.0f;
      bool exists = false;
      pi->findElmore(load, elmore, exists);
      if (!exists) failProbeTest(0, "rc_probe_test_real_sink_has_no_elmore");
      double delay = elmore;
      double slew = sta::delayAsFloat(result.drvrSlew());
      // Exactly the same protected helper used by the passive samples.
      loadDelaySlewElmore(driver, elmore, delay, slew);
      thresholdAdjust(load, library, rf, delay, slew);
      checkLoadAlias(0, delay, slew, reference, index);
    }
  }
#endif
  return result;
}

void ObservedDmpDelayCalc::loadDelaySlew(
    const sta::Pin* load, double driver_slew, const sta::RiseFall* rf,
    const sta::LibertyLibrary* driver_library, const sta::Parasitic* parasitic,
    double& delay, double& slew)
{
  // Equivalent to OpenSTA DmpCeffElmoreDelayCalc's real-load path. Keeping
  // this override in lgista exposes the protected observation helper while
  // retaining the existing DMP excitation and threshold conversion.
  delay = 0.0;
  slew = driver_slew;
  float elmore = 0.0f;
  bool exists = false;
  if (parasitic) parasitics_->findElmore(parasitic, load, elmore, exists);
  if (exists) loadDelaySlewElmore(load, elmore, delay, slew);
  thresholdAdjust(load, driver_library, rf, delay, slew);
}

sta::ArcDcalcResult ObservedDmpDelayCalc::inputPortDelay(
    const sta::Pin* port, float slew, const sta::RiseFall* rf,
    const sta::Parasitic* parasitic, const sta::LoadPinIndexMap& loads,
    const sta::Scene* scene, const sta::MinMax* min_max)
{
  samples_.clear();
  // HBT optimization currently requires a combinational cell driver.
  // Other input ports in its local graph retain the original model.
  return input_calculator_->inputPortDelay(
      port, slew, rf, parasitic, loads, scene, min_max);
}

ObservedElmoreCeffDelayCalc::ObservedElmoreCeffDelayCalc(sta::StaState* state)
    : LocalElmoreCeffDelayCalc(state)
{
}

sta::ArcDelayCalc* ObservedElmoreCeffDelayCalc::copy()
{
  return new ObservedElmoreCeffDelayCalc(this);
}

sta::ArcDcalcResult ObservedElmoreCeffDelayCalc::gateDelay(
    const sta::Pin* driver, const sta::TimingArc* arc,
    const sta::Slew& in_slew, float load_cap,
    const sta::Parasitic* parasitic, const sta::LoadPinIndexMap& loads,
    const sta::Scene* scene, const sta::MinMax* min_max)
{
  samples_.clear();
  auto result = LocalElmoreCeffDelayCalc::gateDelay(
      driver, arc, in_slew, load_cap, parasitic, loads, scene, min_max);
  const auto* ec = localModel<PtElmoreCeff>(parasitic);
  if (!ec) return result;
  // Preserve the same float operations used for genuine loads by the
  // ordinary EC calculator, including the final sqrt rounding.
  const float driver_slew = sta::delayAsFloat(result.drvrSlew());
  const float driver_slew_sq = driver_slew * driver_slew;
  const float factor = dcalcSlewFactor();
  const float factor_sq = factor * factor;
  for (const auto& probe : ec->rcProbes()) {
    auto sample = sampleFor(probe.id, driver, arc, scene, min_max);
    if (!arc->gateModel(scene, min_max)) {
      sample.error = "unsupported_elmore_ceff_rc_probe_model";
    } else if (probe.tree_node_idx >= ec->tree().size()) {
      sample.error = "invalid_rc_probe_tree_index";
    } else {
      const auto& node = ec->tree()[probe.tree_node_idx];
      sample.wire_delay = node.delay;
      sample.slew = std::sqrt(driver_slew_sq + factor_sq * node.impulse_sq);
      validateSample(sample);
    }
    samples_.push_back(std::move(sample));
  }
#ifdef LGISTA_RC_PROBE_TESTING
  if (arc->gateModel(scene, min_max) && !ec->rcProbes().empty()) {
    const auto reference = checkOriginalResult(
        1, this, result, driver, arc, in_slew, load_cap, parasitic,
        loads, scene, min_max);
    const auto* rf = arc->toEdge()->asRiseFall();
    const auto* library = arc->to()->libertyLibrary();
    for (const auto& [load, index] : loads) {
      const auto* record = ec->findLoadByPin(load);
      if (!record || record->tree_node_idx >= ec->tree().size())
        failProbeTest(1, "rc_probe_test_real_sink_has_no_tree_node");
      const auto& node = ec->tree()[record->tree_node_idx];
      double delay = node.delay;
      double slew = std::sqrt(driver_slew_sq + factor_sq * node.impulse_sq);
      thresholdAdjust(load, library, rf, delay, slew);
      checkLoadAlias(1, delay, slew, reference, index);
    }
  }
#endif
  return result;
}

sta::ArcDcalcResult ObservedElmoreCeffDelayCalc::inputPortDelay(
    const sta::Pin* port, float slew, const sta::RiseFall* rf,
    const sta::Parasitic* parasitic, const sta::LoadPinIndexMap& loads,
    const sta::Scene* scene, const sta::MinMax* min_max)
{
  samples_.clear();
  return LocalElmoreCeffDelayCalc::inputPortDelay(
      port, slew, rf, parasitic, loads, scene, min_max);
}

}  // namespace lgista
