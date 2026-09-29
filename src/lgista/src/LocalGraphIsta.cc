// SPDX-License-Identifier: BSD-3-Clause
//
// LocalGraphIsta facade implementation. See LocalGraphIsta.hh.

#include "lgista/LocalGraphIsta.hh"

#include <chrono>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "sta/FuncExpr.hh"
#include "sta/Liberty.hh"
#include "sta/Mode.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/TimingRole.hh"
#include "sta/SearchPred.hh"
#include "sta/Variables.hh"
#include <unordered_set>
#include "search/Sim.hh"

#include "db_sta/dbSta.hh"
#include "sta/Network.hh"
#include "sta/Sta.hh"

#include "LocalSta.hh"
#include "PtGraph.hh"
#include "LocalDmpDelayCalc.hh"
#include "IstaWriteback.hh"

namespace lgista {

static inline double now_ms() {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

LocalGraphIsta::LocalGraphIsta(sta::dbSta* sta) : sta_(sta)
{
  // Construction is deferred to prepare(): LocalSta::copyState and the local
  // delay calculator need a built/levelized global timing graph, which only
  // exists after updateTiming(). Building them here (before any timing run)
  // dereferences a not-yet-built graph and faults.
}

LocalGraphIsta::~LocalGraphIsta()
{
  delete arc_delay_calc_;
  delete local_sta_;
}

void LocalGraphIsta::prepare()
{
  if (prepared_)
    return;
  // Bring the global timing graph fully up to date (ensures graph build +
  // levelization + arrivals + requireds) so the local graph can read upstream
  // slews / arrivals from it. updateTiming(true) is robust even when the
  // design has no clocks/constraints (findRequireds alone can fault then).
  sta_->updateTiming(/*full=*/true);
  local_sta_ = new LocalSta(sta_);
  local_sta_->copyState(sta_);
  arc_delay_calc_ = lgista::makeLocalDelayCalc(sta_);
  prepared_ = true;
}

PtGraph* LocalGraphIsta::analyzeInstance(sta::Instance* inst, bool print_report)
{
  if (inst == nullptr)
    return nullptr;
  if (!prepared_)
    prepare();

  PtGraph* pt_graph = local_sta_->makePtGraph(inst, /*update_timing_first=*/true);
  if (pt_graph == nullptr)
    return nullptr;

  local_sta_->findLocalDelays(pt_graph, arc_delay_calc_);
  local_sta_->findLocalArrivals(pt_graph);
  local_sta_->findLocalRequireds(pt_graph);

  if (print_report)
    local_sta_->printLocalTiming(pt_graph);

  return pt_graph;
}

PtGraph* LocalGraphIsta::analyzeInstance(const std::string& inst_name,
                                         bool print_report)
{
  sta::Network* network = sta_->network();
  sta::Instance* inst = network->findInstance(inst_name.c_str());
  if (inst == nullptr) {
    std::printf("LocalGraphIsta: instance '%s' not found\n", inst_name.c_str());
    return nullptr;
  }
  return analyzeInstance(inst, print_report);
}

void LocalGraphIsta::commitToGlobal(PtGraph* pt_graph)
{
  if (pt_graph != nullptr)
    updateTimingFromPtGraph(pt_graph);
}

void LocalGraphIsta::initializeReadOnly()
{
  ensureInit();
}

void LocalGraphIsta::ensureInit()
{
  if (prepared_)
    return;
  local_sta_ = new LocalSta(sta_);
  local_sta_->copyState(sta_);
  arc_delay_calc_ = lgista::makeLocalDelayCalc(sta_);
  prepared_ = true;
}

bool LocalGraphIsta::analyzeAndCommit(sta::Instance* inst, int depth)
{
  if (inst == nullptr)
    return false;
  ensureInit();

  PtGraph pt_graph(sta_);
  if (!local_sta_->makePtGraph(&pt_graph, inst, -1, depth))
    return false;

  local_sta_->findLocalDelays(&pt_graph, arc_delay_calc_);
  local_sta_->findLocalArrivals(&pt_graph);
  local_sta_->findLocalRequireds(&pt_graph);

  updateTimingFromPtGraph(&pt_graph);
  return true;
}

std::unique_ptr<PtGraph> LocalGraphIsta::analyzeOnly(sta::Instance* inst, int depth)
{
  return analyzeOnly(inst, depth, nullptr);
}

std::unique_ptr<PtGraph> LocalGraphIsta::analyzeOnly(sta::Instance* inst, int depth,
                                                     AnalyzeTimings* tm)
{
  if (inst == nullptr)
    return nullptr;
  ensureInit();

  double t0 = tm ? now_ms() : 0.0;
  auto pt_graph = std::make_unique<PtGraph>(sta_);
  if (!local_sta_->makePtGraph(pt_graph.get(), inst, -1, depth))
    return nullptr;
  if (tm) { tm->make_ptgraph_ms += now_ms() - t0; t0 = now_ms(); }

  // Reduce parasitics into PtGraph-local PtPiElmore so the delay calc uses
  // them at priority #1 (over the global scene). For nets present in the
  // per-engine local container (trial candidate parasitics) this yields the
  // candidate RC; nets absent there fall back to the scene network. Without
  // this, localParasiticLoad would read the scene PiElmore (current geometry)
  // and ignore any trial parasitics built into the local container.
  local_sta_->localParasitics()->recomputePtParasitics(pt_graph.get());
  if (tm) { tm->recompute_para_ms += now_ms() - t0; t0 = now_ms(); }

  local_sta_->findLocalDelays(pt_graph.get(), arc_delay_calc_);
  if (tm) { tm->delays_ms += now_ms() - t0; t0 = now_ms(); }
  local_sta_->findLocalArrivals(pt_graph.get());
  if (tm) { tm->arrivals_ms += now_ms() - t0; t0 = now_ms(); }
  local_sta_->findLocalRequireds(pt_graph.get());
  if (tm) { tm->requireds_ms += now_ms() - t0; }

  // NOTE: deliberately NO updateTimingFromPtGraph() — global graph stays
  // frozen; caller reads local values from the returned PtGraph.
  return pt_graph;
}

namespace {

using TrialPortMap
    = std::unordered_map<const sta::Pin*, sta::LibertyPort*>;

// Evaluate Liberty conditions without calling Sim::evalExpr, whose BDD cache
// is shared. Exhaustive Boolean evaluation preserves identities such as
// A | !A which a three-valued recursive evaluator would leave unknown.
sta::LogicValue trialCondition(
    const sta::FuncExpr* expression,
    const std::unordered_map<const sta::LibertyPort*, const sta::Pin*>& pins,
    const sta::Sim* sim)
{
  if (!expression)
    return sta::LogicValue::one;
  std::unordered_map<const sta::LibertyPort*, bool> values;
  std::vector<const sta::LibertyPort*> unknown;
  for (const auto* port : expression->ports()) {
    const auto it = pins.find(port);
    if (it == pins.end())
      throw std::runtime_error("analyzeTrial: unmapped Liberty condition port");
    const auto value = sim->simValue(it->second);
    if (value == sta::LogicValue::zero || value == sta::LogicValue::one)
      values[port] = value == sta::LogicValue::one;
    else if (value == sta::LogicValue::unknown)
      unknown.push_back(port);
    else
      throw std::runtime_error("analyzeTrial: edge-valued case analysis unsupported");
  }
  if (unknown.size() > 16)
    throw std::runtime_error("analyzeTrial: Liberty condition exceeds 16 unknown ports");
  const auto evaluate = [&](const auto& self, const sta::FuncExpr* expr) -> bool {
    using Op = sta::FuncExpr::Op;
    switch (expr->op()) {
      case Op::port: return values.at(expr->port());
      case Op::zero: return false;
      case Op::one: return true;
      case Op::not_: return !self(self, expr->left());
      case Op::and_: return self(self, expr->left()) && self(self, expr->right());
      case Op::or_: return self(self, expr->left()) || self(self, expr->right());
      case Op::xor_: return self(self, expr->left()) != self(self, expr->right());
    }
    throw std::runtime_error("analyzeTrial: unsupported Liberty condition operator");
  };
  unsigned outcomes = 0;
  for (size_t assignment = 0; assignment < (size_t{1} << unknown.size()); ++assignment) {
    for (size_t i = 0; i < unknown.size(); ++i)
      values[unknown[i]] = (assignment & (size_t{1} << i)) != 0;
    outcomes |= evaluate(evaluate, expression) ? 2 : 1;
    if (outcomes == 3)
      return sta::LogicValue::unknown;
  }
  return outcomes == 2 ? sta::LogicValue::one : sta::LogicValue::zero;
}

void buildMappedTrial(sta::dbSta* sta,
                      sta::Instance* instance,
                      sta::LibertyCell* target,
                      const TrialPortMap& overrides,
                      PtGraph* graph)
{
  auto* network = sta->network();
  auto* source = network->libertyCell(instance);
  if (!target || !source || target->hasSequentials() || source->hasSequentials())
    throw std::runtime_error("analyzeTrial: requires combinational source and target cells");
  if (overrides.empty())
    throw std::runtime_error("analyzeTrial: explicit pin overrides required");
  std::unordered_map<const sta::LibertyPort*, const sta::Pin*> target_pins;
  std::unordered_map<const sta::LibertyPort*, sta::VertexId> target_vertices;
  for (const auto& [pin, port] : overrides) {
    auto* source_port = pin ? network->libertyPort(pin) : nullptr;
    if (!pin || network->instance(pin) != instance || !source_port || !port
        || port->libertyCell() != target)
      throw std::runtime_error("analyzeTrial: invalid mapped port ownership");
    if (source_port->direction() != port->direction()
        || (!port->direction()->isInput() && !port->direction()->isOutput()))
      throw std::runtime_error("analyzeTrial: unsupported mapped port direction");
    if (!target_pins.emplace(port, pin).second)
      throw std::runtime_error("analyzeTrial: non-bijective mapped ports");
  }
  std::unique_ptr<sta::InstancePinIterator> pins(network->pinIterator(instance));
  while (pins->hasNext()) {
    const auto* pin = pins->next();
    auto* port = network->libertyPort(pin);
    if (port && !port->isPwrGnd() && !overrides.count(pin))
      throw std::runtime_error("analyzeTrial: incomplete source pin map: "
                               + std::string(port->name()));
  }
  sta::LibertyCellPortBitIterator target_ports(target);
  while (target_ports.hasNext()) {
    auto* port = target_ports.next();
    if (!port->isPwrGnd() && !target_pins.count(port))
      throw std::runtime_error("analyzeTrial: incomplete target pin map: "
                               + std::string(port->name()));
  }
  for (auto& vertex : graph->ptVertices()) {
    const auto found = overrides.find(vertex.pin());
    if (found != overrides.end()) {
      vertex.setLibertyPort(found->second);
      target_vertices.emplace(found->second, vertex.objectIdx());
    }
  }
  graph->setMappedCellTrial(overrides);
  graph->setRefGate(target);
  auto* scene = graph->scene();
  const auto* sim = scene->mode()->sim();
  const auto* sdc = scene->sdc();
  std::vector<std::pair<sta::TimingArcSet*, bool>> enabled;
  bool have_search_arc = false;
  bool identity_ports = source == target;
  bool annotate_conditions = false;
  for (const auto& [pin, port] : overrides) {
    identity_ports &= network->libertyPort(pin) == port;
    annotate_conditions |= sim->simValue(pin) != sta::LogicValue::unknown;
  }
  // Native Sim annotates conditional edges only for instances whose pin
  // values changed during constant propagation. A fresh replacement with all
  // unknown signal pins retains both conditional paths, even if a condition
  // is syntactically a tautology/contradiction. Identity trials instead read
  // the existing native annotations, including their incremental history.
  for (auto* arcs : target->timingArcSets()) {
    if (arcs->role() != sta::TimingRole::combinational()
        || !arcs->modeName().empty() || !arcs->from() || !arcs->to()
        || arcs->arcCount() == 0)
      throw std::runtime_error("analyzeTrial: unsupported target timing arc role or mode");
    if (!target_pins.count(arcs->from()) || !target_pins.count(arcs->to()))
      throw std::runtime_error("analyzeTrial: target arc has an unmapped endpoint");
    const auto* from_pin = target_pins.at(arcs->from());
    const auto* to_pin = target_pins.at(arcs->to());
    if (sdc->isDisabled(instance, from_pin, to_pin, arcs->role())
        || sdc->isDisabled(arcs)
        || sdc->disabledLibPorts()->count(arcs->from())
        || sdc->disabledLibPorts()->count(arcs->to())
        || (!sta->variables()->condDefaultArcsEnabled() && arcs->isCondDefault()))
      continue;
    bool search_enabled = !sim->isConstant(from_pin) && !sim->isConstant(to_pin);
    if (identity_ports) {
      auto* to_vertex = sta->graph()->pinDrvrVertex(to_pin);
      bool found = false;
      sta::VertexInEdgeIterator native_edges(to_vertex, sta->graph());
      while (native_edges.hasNext()) {
        auto* edge = native_edges.next();
        if (edge->timingArcSet() == arcs && edge->from(sta->graph())->pin() == from_pin) {
          search_enabled &= !sim->isDisabledCond(edge);
          found = true;
          break;
        }
      }
      if (!found)
        throw std::runtime_error("analyzeTrial: identity timing arc missing from native graph");
    } else if (annotate_conditions) {
      if (arcs->cond()) {
        if (trialCondition(arcs->cond(), target_pins, sim) == sta::LogicValue::zero)
          search_enabled = false;
      } else {
        for (auto* other : target->timingArcSets(arcs->from(), arcs->to())) {
          if (other->cond()
              && trialCondition(other->cond(), target_pins, sim) == sta::LogicValue::one) {
            search_enabled = false;
            break;
          }
        }
      }
    }
    if (!target_vertices.count(arcs->from()) || !target_vertices.count(arcs->to()))
      throw std::runtime_error("analyzeTrial: target arc endpoint absent from local graph");
    if (graph->ptVertex(target_vertices.at(arcs->from())).level()
        >= graph->ptVertex(target_vertices.at(arcs->to())).level())
      throw std::runtime_error("analyzeTrial: target arc changes frozen topological order");
    for (const auto* arc : arcs->arcs()) {
      if (!arc->fromEdge()->asRiseFall() || !arc->toEdge()->asRiseFall())
        throw std::runtime_error("analyzeTrial: unsupported transition in target arc");
      for (const auto* min_max : sta::MinMax::range()) {
        if (!arc->gateModel(scene, min_max))
          throw std::runtime_error("analyzeTrial: missing target min/max gate model");
      }
    }
    enabled.emplace_back(arcs, search_enabled);
    have_search_arc |= search_enabled;
  }
  if (!have_search_arc)
    throw std::runtime_error("analyzeTrial: no enabled target timing arcs");
  const size_t count = graph->edgeCount();
  for (size_t i = 1; i < count; ++i) {
    const auto& edge = graph->edge(i);
    if (edge.type() != PtEdgeType::Sentinel && !edge.isWire()
        && overrides.count(graph->ptVertex(edge.ptFromId()).pin())
        && overrides.count(graph->ptVertex(edge.ptToId()).pin()))
      graph->deleteEdge(i);
  }
  for (const auto& [arcs, search_enabled] : enabled) {
    const auto id = graph->makeVirtualEdge(target_vertices.at(arcs->from()),
                                           target_vertices.at(arcs->to()), arcs, false);
    graph->edge(id).setTrialSearchEnabled(search_enabled);
  }

  // OpenSTA DcalcPred deliberately includes Sim-disabled conditional models
  // when aggregating slew. makePtGraph uses SearchPred and omitted them, so
  // restore those local neighboring-cell edges for delay calculation only.
  // The arrival/search mask still follows the native frozen SearchPred.
  sta::SearchPred1 search(sta);
  std::unordered_set<const sta::Edge*> present;
  for (size_t id = 1; id < graph->edgeCount(); ++id) {
    auto& edge = graph->edge(id);
    if (edge.type() != PtEdgeType::Sentinel && edge.hasBase())
      present.insert(edge.edge());
  }
  for (const auto& to : graph->ptVertices()) {
    if (!to.vertex() || !to.isDriver() || overrides.count(to.pin())) continue;
    sta::VertexInEdgeIterator edges(to.vertex(), sta->graph());
    while (edges.hasNext()) {
      auto* edge = edges.next();
      if (edge->isWire() || edge->role() != sta::TimingRole::combinational()
          || edge->isDisabledLoop() || sdc->isDisabledConstraint(edge)
          || sdc->isDisabledCondDefault(edge)
          || (edge->isBidirectInstPath() && !sta->variables()->bidirectInstPathsEnabled()))
        continue;
      auto* from_vertex = edge->from(sta->graph());
      const auto* from_pin = from_vertex->pin();
      const auto* net = network->net(from_pin);
      if (sdc->isDisabledConstraint(from_pin)
          || (net && (network->isPower(net) || network->isGround(net)))) continue;
      if (present.count(edge)) continue;
      auto* from = graph->ptVertex(from_vertex);
      if (!from) {
        // Downstream inputs outside the selected local cone remain its frozen
        // boundary. Upstream RefDriver slew, however, needs every input model.
        if (to.type() == PtVertexType::RefDriver)
          throw std::runtime_error("analyzeTrial: missing upstream conditional delay input: "
                                   + network->pathName(from_pin));
        continue;
      }
      if (sim->isConstant(from_vertex) || sim->isConstant(to.vertex()))
        throw std::runtime_error("analyzeTrial: unsupported constant conditional delay boundary");
      const auto id = graph->makeTrialBaseEdge(edge, from->objectIdx(), to.objectIdx());
      graph->edge(id).setTrialSearchEnabled(
          search.searchFrom(from_vertex, scene->mode())
          && search.searchThru(edge, scene->mode())
          && search.searchTo(to.vertex(), scene->mode()));
      present.insert(edge);
    }
  }
  if (!graph->topoSortVertices())
    throw std::runtime_error("analyzeTrial: mapped local graph is cyclic");
}

}  // namespace

std::unique_ptr<PtGraph> LocalGraphIsta::analyzeTrial(
    sta::Instance* inst,
    int depth,
    sta::LibertyCell* target,
    const std::unordered_map<const sta::Pin*, sta::LibertyPort*>& overrides)
{
  if (!inst)
    return nullptr;
  ensureInit();
  auto graph = std::make_unique<PtGraph>(sta_);
  if (!local_sta_->makePtGraph(graph.get(), inst, -1, depth))
    return nullptr;
  buildMappedTrial(sta_, inst, target, overrides, graph.get());
  // Port replacement must precede reduction: the target input capacitances
  // are loads of the upstream drivers, even though DB connectivity is frozen.
  local_sta_->localParasitics()->recomputePtParasitics(graph.get());
  auto* scene = graph->scene();
  for (const auto* min_max : sta::MinMax::range()) {
    graph->setAnalysisPoint(scene, min_max);
    local_sta_->findLocalDelays(graph.get(), arc_delay_calc_);
  }
  for (const auto* min_max : sta::MinMax::range()) {
    graph->setAnalysisPoint(scene, min_max);
    local_sta_->findLocalArrivals(graph.get());
    local_sta_->findLocalRequireds(graph.get());
  }
  graph->setAnalysisPoint(scene, sta::MinMax::max());
  return graph;
}

sta::Parasitics* LocalGraphIsta::localParasiticsContainer()
{
  return local_sta_ ? local_sta_->localParasitics() : nullptr;
}

size_t LocalGraphIsta::analyzeAll(bool print_report)
{
  if (!prepared_)
    prepare();

  sta::Network* network = sta_->network();
  size_t built = 0;
  sta::LeafInstanceIterator* iter = network->leafInstanceIterator();
  while (iter->hasNext()) {
    sta::Instance* inst = iter->next();
    if (analyzeInstance(inst, print_report) != nullptr)
      ++built;
  }
  delete iter;
  return built;
}

}  // namespace lgista
