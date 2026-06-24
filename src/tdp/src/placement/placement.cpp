// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "placement.h"
#include "MIS3D.h"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "sta/Path.hh"
#include "../eva/eva.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Pin3D.h"
#include "../infrastructure/Common3D.h"

#include <iostream>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <limits>
#include <queue>
#include <stack>
#include <list>
#include <vector>
#include <string>

using utl::TDP;

namespace tdp {

namespace {
struct TarjanState {
    int index_counter = 0;
    std::unordered_map<Node3D*, int> index;
    std::unordered_map<Node3D*, int> lowlink;
    std::unordered_map<Node3D*, bool> on_stack;
    std::stack<Node3D*> stk;
    std::vector<std::vector<Node3D*>> sccs;
};

void tarjanDFS(
    Node3D* v,
    TarjanState& state,
    const std::unordered_map<Node3D*, std::unordered_set<Node3D*>>& forward_of)
{
    state.index[v] = state.lowlink[v] = state.index_counter++;
    state.stk.push(v);
    state.on_stack[v] = true;

    auto it = forward_of.find(v);
    if (it != forward_of.end()) {
        for (auto* w : it->second) {
            if (state.index.find(w) == state.index.end()) {
                tarjanDFS(w, state, forward_of);
                state.lowlink[v] = std::min(state.lowlink[v], state.lowlink[w]);
            } else if (state.on_stack[w]) {
                state.lowlink[v] = std::min(state.lowlink[v], state.index[w]);
            }
        }
    }

    if (state.lowlink[v] == state.index[v]) {
        std::vector<Node3D*> scc;
        Node3D* w;
        do {
            w = state.stk.top();
            state.stk.pop();
            state.on_stack[w] = false;
            scc.push_back(w);
        } while (w != v);
        state.sccs.push_back(std::move(scc));
    }
}
} // anonymous namespace


Placement3D::Placement3D(Timingdp* timingdp)
    : timingdp_(timingdp),
      sta_(timingdp->getSta()),
      logger_(timingdp ? timingdp->getLogger() : nullptr),
      timing_engine_(timingdp->gettimingengine())
{
    if (timingdp_) {
        arch_ = timingdp_->getArchitecture3D();
        network_ = timingdp_->getNetwork3D();
    }
}

Placement3D::~Placement3D()
{
    // 清理资源
}

void Placement3D::buildTraverseOrder(Timingdp* timingdp) {
    traverse_list.clear();
    if (!timingdp) {
        logger_->error(TDP, 1300, "Error: Timingdp is null");
        return;
    }

    auto* network_3d = timingdp->getNetwork3D();
    if (!network_3d) {
        logger_->error(TDP, 1301, "Error: Network3D is null");
        return;
    }

    auto all_cells = network_3d->getCells();

    // ================================================================
    // Phase 1: 计算每个 cell 的出度 + 构建反向邻接表（谁驱动我）
    // ================================================================
    std::unordered_map<Node3D*, int> out_degree;
    std::unordered_map<Node3D*, std::unordered_set<Node3D*>> drivers_of;

    for (auto* cell : all_cells) {
        if (!cell)
            continue;

        std::unordered_set<Node3D*> downstream;

        if (cell->isTerminal()) {
            // ===== 终端 =====
            bool is_pi = false;
            for (auto* pin : cell->getPins()) {
                if (pin && isPrimaryInput(pin)) { is_pi = true; break; }
            }

            if (is_pi) {
                // PI: 找下游 cell
                for (auto* pin : cell->getPins()) {
                    if (!pin) continue;
                    Net3D* net = pin->getNet();
                    if (!net) continue;
                    if (net->db_net && timingdp->isIdealClockNet(net->db_net)) continue;
                    for (auto* net_pin : net->getPins()) {
                        if (net_pin == pin) continue;
                        Node3D* sink = net_pin->getNode();
                        if (!sink || sink == cell) continue;
                        if (downstream.insert(sink).second) {
                            drivers_of[sink].insert(cell);
                        }
                    }
                }
            }
            // PO: 不处理

        } else if (isSequentialCell(cell)) {
            // 寄存器：Q 输出视为伪 PI，下游 sink 记录 FF 为驱动者
            for (auto* pin : cell->getPins()) {
                if (!pin || pin->getDirection() != Pin3D::Dir_OUT)
                    continue;

                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && timingdp->isIdealClockNet(net->db_net)) continue;

                for (auto* net_pin : net->getPins()) {
                    if (net_pin == pin) continue;
                    if (net_pin->getDirection() != Pin3D::Dir_IN && !isPrimaryOutput(net_pin)) continue;

                    Node3D* sink = net_pin->getNode();
                    if (!sink || sink == cell) continue;

                    drivers_of[sink].insert(cell);
                }
            }
        } else {
            // ===== 普通组合逻辑 cell：Dir_OUT → net → sink =====
            for (auto* pin : cell->getPins()) {
                if (!pin || pin->getDirection() != Pin3D::Dir_OUT)
                    continue;

                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && timingdp->isIdealClockNet(net->db_net)) continue;

                for (auto* net_pin : net->getPins()) {
                    if (net_pin == pin) continue;
                    if (net_pin->getDirection() != Pin3D::Dir_IN && !isPrimaryOutput(net_pin)) continue;

                    Node3D* sink = net_pin->getNode();
                    if (!sink || sink == cell) continue;

                    if (downstream.insert(sink).second) {
                        drivers_of[sink].insert(cell);
                    }
                }
            }
        }

        out_degree[cell] = static_cast<int>(downstream.size());
    }

    // ================================================================
    // Phase 2: Kahn 算法 —— 出度 = 0 的 cell 先入队
    // ================================================================
    std::queue<Node3D*> q;
    std::unordered_set<Node3D*> visited;
    std::unordered_set<Node3D*> remaining;


    for (auto* cell : all_cells) {
        if (out_degree[cell] == 0) {
            q.push(cell);
        }
        if (out_degree[cell] > 0 && !isSequentialCell(cell)) {
            if(cell->isTerminal()) continue;
            remaining.insert(cell);
        }
    }

    int cycle_break_count = 0;

    while (!q.empty()) {
        Node3D* cell = q.front();
        q.pop();

        traverse_list.push_back(cell);
        visited.insert(cell);
        remaining.erase(cell); 

        auto it = drivers_of.find(cell);
        if (it != drivers_of.end()) {
            for (auto* driver : it->second) {
                if (--out_degree[driver] == 0 && !visited.count(driver)) {
                    if (!isSequentialCell(driver)) {
                        q.push(driver);
                    }
                }
            }
        }

        if (q.empty() && !remaining.empty()) {
            // 1. 在 remaining 子图上构建正向邻接表
            std::unordered_map<Node3D*, std::unordered_set<Node3D*>> forward_of;
            for (auto* sink : remaining) {
                auto dit = drivers_of.find(sink);
                if (dit != drivers_of.end()) {
                    for (auto* driver : dit->second) {
                        if (remaining.count(driver)) {
                            forward_of[driver].insert(sink);
                        }
                    }
                }
            }

            // 2. Tarjan SCC
            TarjanState state;
            for (auto* c : remaining) {
                if (state.index.find(c) == state.index.end()) {
                    tarjanDFS(c, state, forward_of);
                }
            }

            // 3. 对每个真环 (size > 1)，切断一次
            for (auto& scc : state.sccs) {
                if (scc.size() <= 1) continue;

                Node3D* cut_cell = nullptr;
                int min_deg = std::numeric_limits<int>::max();

                for (auto* c : scc) {
                    if (out_degree[c] == 1) {
                        cut_cell = c;
                        break;
                    }
                    if (out_degree[c] < min_deg) {
                        min_deg = out_degree[c];
                    }
                }
                if (!cut_cell) {
                    for (auto* c : scc) {
                        if (out_degree[c] == min_deg) {
                            cut_cell = c;
                            break;
                        }
                    }
                }
                if (!cut_cell) continue;

                auto fwd_it = forward_of.find(cut_cell);
                if (fwd_it != forward_of.end()) {
                    for (auto* target : fwd_it->second) {
                        if (remaining.count(target)) {
                            drivers_of[target].erase(cut_cell);
                            break;
                        }
                    }
                }

                out_degree[cut_cell] = 0;
                q.push(cut_cell);
                remaining.erase(cut_cell);

                cycle_break_count++; 

                // logger_->warn(TDP, 1951,
                //     "Cycle of size {} broken at cell {}",
                //     scc.size(), cut_cell->name);
            }

            // 4. 处理非环的孤立 cell（下游全是 sequential）
            if (q.empty() && !remaining.empty()) {
                for (auto& scc : state.sccs) {
                    if (scc.size() != 1) continue;
                    auto* c = scc[0];
                    if (!remaining.count(c)) continue;
                    out_degree[c] = 0;
                    q.push(c);
                    remaining.erase(c);
                }
            }
        }
    }

    
    // ================================================================
    // Debug 输出
    // ================================================================
    int total_cells = network_->getCells().size();
    int traversed = traverse_list.size();
    std::cout << "Traverse list: " << traversed 
          << " / Total cells: " << total_cells 
          << " / Remaining: " << remaining.size()
          << " / Cycles broken: " << cycle_break_count << std::endl;
    


    // int count = 0;
    // std::cout << "Last 10 cells in traverse_list:" << std::endl;
    // for (auto it = traverse_list.rbegin();
    //      it != traverse_list.rend() && count < 10; ++it, ++count) {
    //     std::cout << "  " << (*it)->name
    //               << " (layer: " << (*it)->layer << ")" << std::endl;
    // }

    if (traverse_list.empty()) {
        logger_->warn(TDP, 1303, "WARNING: Traverse list is empty!");
    } else {
        logger_->info(TDP, 1304, "Reverse topological order built, {} cells", traverse_list.size());
    }
}


bool Placement3D::isInputPin(const Pin3D* pin) {
    return pin->getDirection() == Pin3D::Dir_IN;
}

// 检查引脚是否是输出引脚
bool Placement3D::isOutputPin(const Pin3D* pin) {
    return pin->getDirection() == Pin3D::Dir_OUT;
}

// 检查引脚是否是主输入（边界输入引脚）
bool Placement3D::isPrimaryInput(const Pin3D* pin) {
    // 边界引脚且方向为输入
    return pin->getDbBTerm() != nullptr && 
           pin->getDirection() == Pin3D::Dir_IN;
}

// 检查引脚是否是主输出（边界输出引脚）
bool Placement3D::isPrimaryOutput(const Pin3D* pin) {
    // 边界引脚且方向为输出
    return pin->getDbBTerm() != nullptr && 
           pin->getDirection() == Pin3D::Dir_OUT;
}

bool Placement3D::isPrimaryInput(const Node3D* cell) const {
    if (!cell || !cell->isTerminal()) return false;
    if (!cell->db_bterm_) return false;
    return cell->db_bterm_->getIoType() == odb::dbIoType::INPUT;
}

bool Placement3D::isPrimaryOutput(const Node3D* cell) const {
    if (!cell || !cell->isTerminal()) return false;
    if (!cell->db_bterm_) return false;
    return cell->db_bterm_->getIoType() == odb::dbIoType::OUTPUT;
}

// 判断一个 cell 是否是时序元件（寄存器/FF/latch）
bool Placement3D::isSequentialCell(const Node3D* cell) const{
    if (!cell || cell->isTerminal()) return false;
    if (cell->db_inst && cell->db_inst->getMaster()) {
        return cell->db_inst->getMaster()->isSequential();
    }
    return false;
}

//find cell row 
Architecture3D::Row3D* Placement3D::findCellRow(Node3D* cell, 
                                                const std::vector<Architecture3D::Row3D*>& rows) {
    
    if (!cell || rows.empty()) {
        return nullptr;
    }
    // 获取cell的位置和尺寸
    Pointdp3D cell_loc = cell->location;
    DBU cell_width = cell->getWidth().v;
    DBU cell_height = cell->getHeight().v;
    
    // 遍历所有行，查找包含cell的行
    for (auto* row : rows) {
        if (!row) continue;
        
        // 获取行的边界
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        
        // 检查cell是否在行内（考虑cell的尺寸）
        bool in_y = (cell_loc.y.v >= row_bottom) && 
                    (cell_loc.y.v + cell_height <= row_top);
        
        bool in_x = (cell_loc.x.v >= row_left) && 
                    (cell_loc.x.v + cell_width <= row_right);
        
        if (in_y && in_x) {
            return row;  // 找到cell所在的行
        }
    }
    
   std::cout << "WARNING: cell row is empty!" << std::endl;
    
    return nullptr;  // 没找到任何匹配的行
}


void Placement3D::rundetailedplacement(){


        LRParams params;
        params.lambda_wirelength = 1.0;
        params.lambda_timing = 1.0;
        params.lambda_alignment = 0.5;
        params.max_iterations = 10;
        params.convergence_threshold = 0.001;
        params.window_y = 5;
        params.window_x = 20;
        

        LSP3D tdp_lsp(timingdp_, params);
        DetailedMis3D tdp_mis(timingdp_, &tdp_lsp);
        buildTraverseOrder(timingdp_);
        verifyTraverseOrder();
        tdp_lsp.traverse_list = traverse_list;

        // 调用 initWeight 函数
        tdp_lsp.place();

        // tdp_mis.run();
        //test file

        
}
// void placeall (){


// }

void Placement3D::verifyTraverseOrder() const
{
    if (!network_) {
        if (logger_)
            logger_->warn(TDP, 9100, "verifyTraverseOrder: network_ is null");
        return;
    }

    // 1. 建位置索引
    std::unordered_map<const Node3D*, int> pos;
    int idx = 0;
    for (auto* cell : traverse_list) {
        pos[cell] = idx++;
    }

    // 2. 统计
    int total_arcs      = 0;
    int violations      = 0;
    int skipped_fixed   = 0;
    int skipped_missing = 0;

    int cur_pos = 0;
    for (auto* cell : traverse_list) {
        if (!cell || cell->is_fixed ) {
            cur_pos++;
            continue;
        }

        for (auto* pin : cell->getPins()) {
            if (!pin) continue;
            if (pin->getDirection() != Pin3D::Dir_IN) continue;

            Net3D* net = pin->getNet();
            if (!net) continue;

            // 在 net 上找 driver pin（输出引脚）
            Pin3D* driver_pin = nullptr;
            for (auto* net_pin : net->getPins()) {
                if (!net_pin || net_pin == pin) continue;
                if (net_pin->getDirection() == Pin3D::Dir_OUT) {
                    driver_pin = net_pin;
                    break;
                }
            }
            if (!driver_pin) continue;

            Node3D* driver_cell = driver_pin->getNode();
            if (!driver_cell || driver_cell == cell) continue;

            if (driver_cell->isTerminal() || driver_cell->is_fixed || isSequentialCell(driver_cell)) {
                skipped_fixed++;
                continue;
            }

            total_arcs++;

            auto it = pos.find(driver_cell);
            if (it == pos.end()) {
                skipped_missing++;
                continue;
            }

            int driver_pos = it->second;

            // 逆拓扑序：消费者(load) 应在 生产者(driver) 之前
            if (cur_pos >= driver_pos) {
                violations++;
                // std::cout << "[REVERSE-TOPO VIOLATION #" << violations << "] "
                //           << cell->name << " (pos=" << cur_pos
                //           << ") >= driver "
                //           << driver_cell->name
                //           << " (pos=" << driver_pos << ")\n";
            }
        }
        cur_pos++;
    }

    std::cout << "\n=== Traverse Order Verification ===\n"
              << "  Total cells:         " << traverse_list.size() << "\n"
              << "  Total driver arcs:   " << total_arcs << "\n"
              << "  Violations:          " << violations << "\n"
              << "  Fixed/terminal skip: " << skipped_fixed << "\n"
              << "  Driver not in list:  " << skipped_missing << "\n"
              << "  Result:              "
              << (violations == 0 ? "PASS" : "FAIL")
              << std::endl;
}

void Placement3D::refreshPinTiming(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;
    auto cell_pins = cell->getPins();
    /////test
    //  auto pctDiff = [](double v1, double v2) -> double {
    //     double denom = std::max(std::abs(v1), std::abs(v2));
    //     const double eps = 1e-15;
    //     return (denom < eps) ? 0.0 : std::abs(v1 - v2) / denom * 100.0;
    // };
    // if (cell->name == "_647_"){
    //     logger_->info(TDP, 1933, " sta im");
    // }
    ////test
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        auto timing_info = timing_engine_->getPinTimingInfo(sta_pin);
        if (timing_info.is_valid) {
            pin->setEarlyArrival(timing_info.early_arrival);
            pin->setLateArrival(timing_info.late_arrival);
            pin->setEarlyRequired(timing_info.early_required);
            pin->setLateRequired(timing_info.late_required);
            pin->setEarlySlack(timing_info.early_slack);
            pin->setLateSlack(timing_info.late_slack);
        }
        //test////////////
        // auto timing_info2 = timing_engine_->getPinTimingInfo(sta_pin);
        // double max_diff = 0.0;
        // // max_diff = std::max(max_diff, pctDiff(timing_info.early_arrival, full.early_arrival));
        // // max_diff = std::max(max_diff, pctDiff(timing_info.late_arrival,  full.late_arrival));
        // // max_diff = std::max(max_diff, pctDiff(timing_info2.early_required, full.early_required));
        // // max_diff = std::max(max_diff, pctDiff(bypass.late_required,  full.late_required));
        // max_diff = std::max(max_diff, pctDiff(timing_info.early_slack,   timing_info2.early_slack));
        // max_diff = std::max(max_diff, pctDiff(timing_info.late_slack,    timing_info2.late_slack));
        // if (max_diff > 3.0) {
        //     logger_->info(TDP, 1943, "all sta im");
        // }
        //test////////////
    }
}

void Placement3D::refreshPinTiming2(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;
    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        auto timing_info = timing_engine_->getPinTimingInfoDirect(sta_pin);
        if (timing_info.is_valid) {
            pin->setEarlyArrival(timing_info.early_arrival);
            pin->setLateArrival(timing_info.late_arrival);
            pin->setEarlyRequired(timing_info.early_required);
            pin->setLateRequired(timing_info.late_required);
            pin->setEarlySlack(timing_info.early_slack);
            pin->setLateSlack(timing_info.late_slack);
        }
    }
}

void Placement3D::compareLocalVsFullSTA(Node3D* cell) {
    if (!cell || !timing_engine_) return;
    sta::dbNetwork* network = timing_engine_->getNetwork();
    if (!network) return;

    struct PinTimingPair {
        std::string name;
        TimingEngine3D::PinTimingInfo local;
        TimingEngine3D::PinTimingInfo full;
    };
    std::vector<PinTimingPair> results;

    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        sta::Pin* sta_pin = network->findPin(pin->getName().c_str());
        if (!sta_pin) continue;
        PinTimingPair p;
        p.name = pin->getName();
        p.local = timing_engine_->getPinTimingInfoDirect(sta_pin);
        results.push_back(p);
    }

    for (auto& p : results) {
        sta::Pin* sta_pin = network->findPin(p.name.c_str());
        if (!sta_pin) continue;
        p.full = timing_engine_->getPinTimingInfo(sta_pin);
    }

    auto pctDiff = [](double a, double b) -> double {
        double denom = std::max(std::abs(a), std::abs(b));
        if (denom < 1e-18) return 0.0;
        return std::abs(a - b) / denom * 100.0;
    };

    double max_arr_pct = 0.0, max_req_pct = 0.0, max_slk_pct = 0.0;
    std::string worst_pin;

    for (const auto& p : results) {
        if (!p.local.is_valid || !p.full.is_valid) continue;
        double arr_pct = pctDiff(p.local.late_arrival, p.full.late_arrival);
        double req_pct = pctDiff(p.local.late_required, p.full.late_required);
        double slk_pct = pctDiff(p.local.late_slack, p.full.late_slack);
        double slk_abs = std::abs(p.local.late_slack - p.full.late_slack);

        if (slk_pct > max_slk_pct) {
            max_slk_pct = slk_pct;
            worst_pin = p.name;
        }
        max_arr_pct = std::max(max_arr_pct, arr_pct);
        max_req_pct = std::max(max_req_pct, req_pct);

        if (arr_pct > 1.0 || req_pct > 1.0) {
            logger_->info(TDP, 890,
                "  pin={} late_arr: local={:.4e} full={:.4e} diff={:.1f}%"
                "  late_req: local={:.4e} full={:.4e} diff={:.1f}%"
                "  late_slk: local={:.4e} full={:.4e} diff={:.1f}% abs={:.4e}",
                p.name,
                p.local.late_arrival, p.full.late_arrival, arr_pct,
                p.local.late_required, p.full.late_required, req_pct,
                p.local.late_slack, p.full.late_slack, slk_pct, slk_abs);
        }
    }

    // logger_->info(TDP, 891,
    //     "cell={} max_arr_diff={:.2f}% max_req_diff={:.2f}% max_slk_diff={:.2f}% worst_pin={}",
    //     cell->name, max_arr_pct, max_req_pct, max_slk_pct, worst_pin);
}

// ============================================================================
// sta::Pin → Node3D 映射辅助函数
// ============================================================================
static Node3D* staPinToNode3D(sta::dbNetwork* db_network,
                               sta::Pin* sta_pin,
                               Network3D* network_3d)
{
    odb::dbITerm* iterm = nullptr;
    odb::dbBTerm* bterm = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    db_network->staToDb(sta_pin, iterm, bterm, moditerm);
    
    if (iterm) {
        odb::dbInst* inst = iterm->getInst();
        if (inst) return network_3d->getCell(inst);
    }
    if (bterm) {
        return network_3d->getNode(bterm);
    }
    return nullptr;
}

// ============================================================================
// 构建关键路径优先遍历顺序
//   取 top-k 违反约束路径，收集沿途 cell（排除终端 PI/PO），
//   路径回溯顺序即为逆拓扑序（PO→PI），剩余 cell 保持原顺序追加。
// ============================================================================
std::vector<Node3D*> LSP3D::buildCriticalPathPriorityOrder(
    TimingEngine3D* timing_engine,
    int top_k)
{
    std::vector<Node3D*> priority;
    std::unordered_set<Node3D*> in_priority;

    if (!timing_engine || !network_) return priority;

    sta::dbSta* sta = timing_engine->getSta();
    sta::dbNetwork* db_network = sta->getDbNetwork();

    auto paths = timing_engine->getViolatingPaths();

    int path_count = 0;
    for (sta::Path* p : paths) {
        if (path_count >= top_k) break;

        sta::Path* step = p;
        while (step) {
            sta::Pin* sta_pin = step->pin(sta);
            if (sta_pin) {
                Node3D* cell = staPinToNode3D(db_network, sta_pin, network_);
                if (cell
                    && !cell->isTerminal()
                    && !cell->is_fixed
                    && !in_priority.count(cell)) {
                    in_priority.insert(cell);
                    priority.push_back(cell);
                }
            }
            step = step->prevPath();
        }
        path_count++;
    }

    // for (auto* cell : traverse_list) {
    //     if (cell
    //         && !cell->is_fixed
    //         && !in_priority.count(cell)) {
    //         priority.push_back(cell);
    //     }
    // }

    return priority;
}

} // namespace tdp
