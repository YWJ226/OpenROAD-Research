// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "placement.h"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "../eva/eva.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"
#include "../infrastructure/Node3D.h"
#include "../infrastructure/Net3D.h"
#include "../infrastructure/Pin3D.h"
#include "../infrastructure/Common3D.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits> 
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <vector>
#include <omp.h>

using utl::TDP;

// ============================================================================
// 层 RC 常量（Ω/m, F/m）
// ============================================================================
// bottom: metal layer 9,  upper: metal layer 18
struct LayerRC { double r_per_m; double c_per_m; };
static constexpr LayerRC kLayerRC[2] = {
    // bottom (layer 9):   R=1.4386e-03 Ω/m, C=7.96529e-02 F/m
    {1.4386e-03, 7.96529e-02},
    // upper  (layer 18):  R=3.5875e-02 Ω/m, C=1.51771e-01 F/m
    {3.5875e-02, 1.51771e-01}
};
static int getLayerIndex(const std::string& layer) {
    return (layer == "bottom") ? 0 : 1;
}

namespace tdp {

LSP3D::LSP3D(Timingdp* timingdp, const LRParams& params)
    : Placement3D(timingdp), params_(params)
{
    if (timingdp_) {
        arch_ = timingdp_->getArchitecture3D();
        network_ = timingdp_->getNetwork3D();
    }
    // Initialize local timing engine for parallel candidate evaluation
    auto* sg = new SteinerTreeGenerator(timingdp_);
    local_timing_ = std::make_unique<LocalTimingPropagator>(
        sg, timing_engine_, sta_);
}

LSP3D::~LSP3D()
{
    // 清理资源
}

// ============================================================================
// 主要布局方法
// ============================================================================

void LSP3D::place()
{
    if (!timingdp_ || !arch_ || !network_) {
        logger_->error(TDP, 1305, "Error: Database, architecture, or network not initialized");
        return;
    }
    
    initWeight(timingdp_,timing_engine_);
    // 获取遍历顺序（反向拓扑顺序）
    logger_ -> info(TDP, 1934, "Starting discrete local search...");
    

    for (int iter = 0; iter < params_.max_iterations; ++iter) {
        ///test
        sta_->updateTiming(false);
        sta_->findRequireds();

        auto stats = timing_engine_->getTimingStats();
        logger_->info(TDP, 1294, "Iter {}: Early WNS={:.6f}ns TNS={:.6f}ns | Late WNS={:.f}ns TNS={:.6f}ns",
            iter,
            stats.early_wns * 1e9, stats.early_tns * 1e9,
            stats.late_wns * 1e9, stats.late_tns * 1e9);        
        updatePinWeights();
        logger_->info(TDP, 1292, " update pin weight over" );   
        double L_lambda = 0.0;
        for (auto* cell : network_->getCells()) {
            L_lambda += calculateCellCost(cell, 1.0, 1.0);
        }
        logger_->info(TDP, 1293, "L(λ) = {}", L_lambda );   

        discreteLocalSearchAll(timing_engine_); 
        // 每轮结束后重建所有 Steiner 树，更新寄生参数
        SteinerTreeGenerator stn_gen(timingdp_);
        stn_gen.createSteninerRoute();
            
    }
    
}

void LSP3D::updatePinWeights()
{
    if (!network_) {
        logger_->error(TDP, 1306, "Error: Network not initialized");
        return;
    }
    
    // 获取所有单元
    // const auto& all_cells = network_->getCells();
    
    int updated_pins = 0;
    // 遍历所有单元
    for (auto* cell : traverse_list) {
        if (!cell) continue;
        refreshPinTiming2(cell);
        
        // 获取单元的所有引脚
        auto cell_pins = cell->getPins();

        if(isSequentialCell(cell)){


            for (auto* pin : cell_pins) {
                if (!pin || pin->getDirection() != Pin3D::Dir_IN) continue;

                Net3D* net = pin->getNet();
                if (!net) continue;
                if (net->db_net && timingdp_->isIdealClockNet(net->db_net)) continue;

                double current_late_weight = pin->getLateWeight();
                double current_early_weight = pin->getEarlyWeight();
                
                // 获取时序信息
                double late_arrival = pin->getLateArrival();
                double late_required = pin->getLateRequired();
                double early_arrival = pin->getEarlyArrival();
                double early_required = pin->getEarlyRequired();
                
                // 避免除零错误
                if (late_required != 0.0 && early_arrival != 0.0) {

                    if(early_required < 0.0){
                            early_required = early_arrival / 10;
                    }

                    double new_late_weight = current_late_weight * (late_arrival / late_required);
                    double new_early_weight = current_early_weight * (early_required / early_arrival); 
                    if(new_late_weight < 0.0  || new_early_weight < 0.0){
                        logger_ -> info(TDP, 1321, "DDX weight flase");
                    }                  
                    pin->setLateWeight(new_late_weight);
                    pin->setEarlyWeight(new_early_weight);
                }
            }
            continue;

        }
        else{

            if (cell->getBTerm() != nullptr){
            for (auto* pin : cell_pins) {
                if (!pin) continue;
                updated_pins++;
                // 只更新输出引脚
                // if (pin->getDbBTerm() != nullptr && pin->getDirection() == Pin3D::Dir_OUT) {
                if (pin->getDbBTerm() != nullptr && pin->getDirection() == Pin3D::Dir_OUT) {
                //     std::cout << "this is outport: " << pin->getName() << std::endl;}
                // if (pin->getDirection() == Pin3D::Dir_OUT) {
                    // 获取当前权重
                    double current_late_weight = pin->getLateWeight();
                    double current_early_weight = pin->getEarlyWeight();
                    
                    // 获取时序信息
                    double late_arrival = pin->getLateArrival();
                    double late_required = pin->getLateRequired();
                    double early_arrival = pin->getEarlyArrival();
                    double early_required = pin->getEarlyRequired();
                    
                    // 避免除零错误
                    if (late_required != 0.0 && early_arrival != 0.0) {

                        if(early_required < 0.0){
                            early_required = early_arrival / 10;
                        }

                        double new_late_weight = current_late_weight * (late_arrival / late_required);
                        double new_early_weight = current_early_weight * (early_required / early_arrival);                   
                        if(new_late_weight < 0.0 || new_early_weight <0.0){
                            logger_ -> info(TDP, 1322, "PO weight flase");
                        }                            
                        pin->setLateWeight(new_late_weight);
                        pin->setEarlyWeight(new_early_weight);
                        

                    }
                }
            }
            }
            

        }
        
 
        balanceCellWeights(cell);
        distributeCellWeights(cell);
    }
    
    std::cout << "Updated weights for " << updated_pins << " output pins" << std::endl;
}

/**
 * @brief 平衡单元内部输入引脚权重（基于公式 28）
 * 
 * 对每个输入引脚 i，找到上游驱动引脚 driver，对每个输出引脚 j：
 *   λ_i,j^L = λ_i,j^L × (a_i^L + d_i,j^L) / a_j^L
 *   λ_i,j^E = λ_i,j^E × a_i^E / (a_i^E + d_i,j^E)
 * 对每路输出 j 的结果取算术平均，更新输入引脚权重。
 */
void LSP3D::balanceCellWeights(Node3D* cell) {
    if (!cell || !network_ || !timing_engine_) {
        return;
    }

    auto input_pins  = getCellInputPins(cell);
    auto output_pins = getCellOutputPins(cell);
    if (input_pins.empty() || output_pins.empty()) {
        return;
    }

    for (auto* in_pin : input_pins) {
        if (!in_pin) continue;

        // Pin3D* driver = findDriverPin(in_pin->getNet());
        // if (!driver) continue;

        double a_i_L = in_pin->getLateArrival();
        double a_i_E = in_pin->getEarlyArrival();

        double sum_late  = 0.0;
        double sum_early = 0.0;
        int    cnt       = 0;

        for (auto* out_pin : output_pins) {
            if (!out_pin) continue;

            auto arc = timing_engine_->getArcDelay(in_pin->getName(),
                                                    out_pin->getName());
            if (!arc.is_valid) continue;

            double a_j_L = out_pin->getLateArrival();
            double a_j_E = out_pin->getEarlyArrival();

            if (a_j_L != 0.0) {
                sum_late += in_pin->getLateWeight()
                          * (a_i_L + arc.late_delay) / a_j_L;
            }
            if (a_i_E + arc.early_delay != 0.0) {
                sum_early += in_pin->getEarlyWeight()
                           * a_j_E / (a_i_E + arc.early_delay);
            }
            cnt++;
        }

        if (cnt > 0) {
            if(sum_late < 0.0 || sum_early < 0.0){
                logger_ -> info(TDP, 1323, "weight flase");
            }    
            in_pin->setLateWeight(sum_late / cnt);
            in_pin->setEarlyWeight(sum_early / cnt);
        }
    }
}

/**
 * @brief 在单元内部重新分布权重，确保输入总权重和输出总权重平衡
 * @param cell 要重新分布权重的单元
 */
void LSP3D::distributeCellWeights(Node3D* cell) {
    if (!cell || !network_) {
        return;
    }
    
    // 初始化权重总和
    double in_late_weight = 0.0;
    double in_early_weight = 0.0;
    double out_late_weight = 0.0;
    double out_early_weight = 0.0;
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    // 1. 计算输入引脚和输出引脚的总权重
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_IN) {
            Net3D* net = pin->getNet();
            if (!net) continue;
            if (net->db_net && timingdp_->isIdealClockNet(net->db_net)) {
                continue;
            }
            in_late_weight += pin->getLateWeight();
            in_early_weight += pin->getEarlyWeight();
        } else if (pin->getDirection() == Pin3D::Dir_OUT) {
            // 输出引脚
            // out_late_weight += pin->getLateWeight();
            // out_early_weight += pin->getEarlyWeight();
            Net3D* net = pin->getNet();
                if (net) {
                    auto net_pins = net->getPins();
                    for (auto* sink_pin : net_pins) {
                        if (sink_pin == pin) continue;
                        out_late_weight += sink_pin->getLateWeight();
                        out_early_weight += sink_pin->getEarlyWeight();
                    }
                }
        }
    }
    // 2. 重新分布输入引脚权重
    if (in_late_weight > 0.0 && in_early_weight > 0.0) {
        // 计算权重调整比例
        double late_ratio = (out_late_weight > 0.0) ? (out_late_weight / in_late_weight) : 1.0;
        double early_ratio = (out_early_weight > 0.0) ? (out_early_weight / in_early_weight) : 1.0;
        
        //std::cout << "  Adjustment ratios: late=" << late_ratio << ", early=" << early_ratio << std::endl;
        
        // 调整所有输入引脚权重
        for (auto* pin : cell_pins) {
            if (!pin || pin->getDirection() != Pin3D::Dir_IN) {
                continue;
            }
            
            // 保存原始权重用于调试
            double old_late = pin->getLateWeight();
            double old_early = pin->getEarlyWeight();
            
            // 应用调整
            double new_late = old_late * late_ratio;
            double new_early = old_early * early_ratio;
            
            pin->setLateWeight(new_late);
            pin->setEarlyWeight(new_early);

        }
    }
}
// ============================================================================
// 窗口寻找
// ============================================================================
void LSP3D::computeTimingVectors(
    Node3D* cell_j,
    Node3D* connected_cell,
    bool is_fan_in,
    double late_lm,
    double early_lm,
    TimingVector& late_vector,
    TimingVector& early_vector) {
    
    if (!cell_j || !connected_cell) {
        return;
    }
    
    // 获取单元中心位置
    Pointdp3D center_j = getCellCenter(cell_j);
    Pointdp3D center_conn = getCellCenter(connected_cell);
    
    // 计算方向向量（从cell_j指向connected_cell）
    double dx = center_conn.x.v - center_j.x.v;
    double dy = center_conn.y.v - center_j.y.v;
    
    // 计算距离
    double distance = std::sqrt(dx * dx + dy * dy);
    if (distance < 1e-10) {
        return;  // 避免除零
    }
    
    // 单位向量
    double unit_x = dx / distance;
    double unit_y = dy / distance;
    
    if (is_fan_in) {
        // 扇入单元：
        // - 晚时序向量：指向互连单元（吸引力）
        // - 早时序向量：背向互连单元（排斥力）
        late_vector = TimingVector(unit_x * late_lm, unit_y * late_lm, late_lm);
        early_vector = TimingVector(-unit_x * early_lm, -unit_y * early_lm, early_lm);
    } else {
        // 扇出单元：
        // - 晚时序向量：指向互连单元（吸引力）
        // - 早时序向量：背向互连单元（排斥力）
        late_vector = TimingVector(unit_x * late_lm, unit_y * late_lm, late_lm);
        early_vector = TimingVector(-unit_x * early_lm, -unit_y * early_lm, early_lm);
    }
}

/**
 * @brief 获取单元中心位置
 */
Pointdp3D LSP3D::getCellCenter(Node3D* cell) const {
    if (!cell) {
        return Pointdp3D(0, 0, 0);
    }
    
    DBU center_x = cell->location.x.v + cell->getWidth().v / 2;
    DBU center_y = cell->location.y.v + cell->getHeight().v / 2;
    
    return Pointdp3D(DbuX{center_x}, DbuY{center_y}, cell->location.z);
}

/**
 * @brief 计算合成向量和方向分量
 * @param cell 当前单元
 * @param late_components 输出的晚时序方向分量
 * @param early_components 输出的早时序方向分量
 */
bool LSP3D::computeDirectionComponents(
    Node3D* cell,
    DirectionComponents& late_components,
    DirectionComponents& early_components) {
    
    if (!cell || !network_) {
        return false;
    }
    
    // 初始化合成向量
    TimingVector late_synthetic(0.0, 0.0, 0.0);
    TimingVector early_synthetic(0.0, 0.0, 0.0);
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        Net3D* net = pin->getNet();
        if (!net) continue;
        
        // 获取引脚方向
        int dir_int = pin->getDirection();
        Pin3D::Direction dir = static_cast<Pin3D::Direction>(dir_int);
        
        
        if (dir == Pin3D::Dir_IN) {
            // 输入引脚：查找驱动单元
            // 获取LM值（这里需要从引脚或网络获取）
            double late_lm = pin->getLateWeight();  // 使用权重作为LM的代理
            double early_lm = pin->getEarlyWeight();
            Pin3D* driver_pin = findDriverPin(net, pin);
            if (!driver_pin) continue;
            
            Node3D* driver_cell = driver_pin->getNode();
            if (!driver_cell) continue;
            
            // 计算向量
            TimingVector late_vec, early_vec;
            computeTimingVectors(cell, driver_cell, true, late_lm, early_lm, 
                                late_vec, early_vec);
            
            // 累加到合成向量
            late_synthetic.x += late_vec.x;
            late_synthetic.y += late_vec.y;
            early_synthetic.x += early_vec.x;
            early_synthetic.y += early_vec.y;
            
        } else if (dir == Pin3D::Dir_OUT) {
            // 输出引脚：查找所有接收单元
            auto net_pins = net->getPins();
            
            for (auto* sink_pin : net_pins) {
                if (sink_pin == pin) continue;  // 跳过驱动引脚
                // 获取LM值（这里需要从引脚或网络获取）
                double late_lm = sink_pin->getLateWeight();  // 使用权重作为LM的代理
                double early_lm = sink_pin->getEarlyWeight();
                
                Node3D* sink_cell = sink_pin->getNode();
                if (!sink_cell) continue;
                
                // 计算向量
                TimingVector late_vec, early_vec;
                computeTimingVectors(cell, sink_cell, false, late_lm, early_lm,
                                    late_vec, early_vec);
                
                // 累加到合成向量
                late_synthetic.x += late_vec.x;
                late_synthetic.y += late_vec.y;
                early_synthetic.x += early_vec.x;
                early_synthetic.y += early_vec.y;
            }
        }
    }
    
    // 分解为方向分量
    decomposeToComponents(late_synthetic, late_components);
    decomposeToComponents(early_synthetic, early_components);
    
    return true;
}

/**
 * @brief 将向量分解为方向分量
 */
void LSP3D::decomposeToComponents(
    const TimingVector& vector,
    DirectionComponents& components) {
    
    // 水平分量
    if (vector.x > 0) {
        components.right += std::abs(vector.x);
    } else {
        components.left +=  std::abs(vector.x);
    }
    
    // 垂直分量
    if (vector.y > 0) {
        components.up += std::abs(vector.y);
    } else {
        components.down +=  std::abs(vector.y);
    }
}

/**
 * @brief 确定搜索窗口位置
 * @param cell 当前单元
 * @param window_width 搜索窗口宽度
 * @param window_height 搜索窗口高度
 * @return 搜索窗口
 */
LSP3D::SearchWindow LSP3D::determineSearchWindow(
    Node3D* cell,
    DBU window_width,
    DBU window_height) {
    
    if (!cell) {
        return SearchWindow();
    }
    
    // 计算方向分量
    DirectionComponents late_components, early_components;
    if (!computeDirectionComponents(cell, late_components, early_components)) {
        // 如果计算失败，返回以单元为中心的默认窗口
        // std::cerr << "Error: No find cell windows " << std::endl;
        DBU x0 = cell->location.x.v - window_width / 2;
        DBU y0 = cell->location.y.v - window_height / 2;
        return SearchWindow(x0, y0, window_width, window_height);
    }
    
    // 合并晚时序和早时序分量
    DirectionComponents total_components;
    total_components.left = late_components.left + early_components.left;
    total_components.right = late_components.right + early_components.right;
    total_components.up = late_components.up + early_components.up;
    total_components.down = late_components.down + early_components.down;
    
    // 获取单元左下角位置
    DBU xj = cell->location.x.v;
    DBU yj = cell->location.y.v;
    
    // 初始窗口位置（单元位于窗口中心）
    DBU x0 = xj - window_width / 2 ;
    DBU y0 = yj - window_height / 2 ;
    
    // 计算调整后的窗口位置
    DBU x0_adjusted = x0;
    DBU y0_adjusted = y0;
    
    // 计算水平调整
    double total_horizontal = total_components.left + total_components.right;
    if (total_horizontal > 1e-10) {
        double left_ratio = total_components.left / total_horizontal;
        x0_adjusted = xj - left_ratio * window_width;
    }
    
    // 计算垂直调整
    double total_vertical = total_components.down + total_components.up;
    if (total_vertical > 1e-10) {
        double down_ratio = total_components.down / total_vertical;
        y0_adjusted = yj - down_ratio * window_height;
    }
    
    return SearchWindow(x0_adjusted, y0_adjusted, window_width, window_height);
}

/**
 * @brief 在搜索窗口覆盖的多个行中生成候选位置（考虑空位）
 * @param cell 当前单元
 * @param window 搜索窗口
 * @param num_candidates 候选位置数量
 * @param step_x X方向步长（站点间距）
 * @return 候选位置列表
 */
std::vector<Pointdp3D> LSP3D::generateCandidatePositions(
    Node3D* cell,
    const SearchWindow& window,
    int num_candidates,
    DBU step_x) {
    
    std::vector<Pointdp3D> candidates;
    
    if (!cell || !arch_ || !network_ || window.width <= 0 || num_candidates <= 0) {
        return candidates;
    }
    
    // 获取单元所在的层
    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    const auto& all_rows = arch_->getRows(layer_type);
    
    if (all_rows.empty()) {
        return candidates;
    }
    
    // 1. 找到搜索窗口覆盖的行
    std::vector<Architecture3D::Row3D*> covered_rows;
    
    for (auto* row : all_rows) {
        if (!row) continue;
        
        // 检查行是否与搜索窗口重叠
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;
        
        // Y方向重叠检查
        bool y_overlap = (row_bottom <= window.y + window.height) && 
                        (row_top >= window.y);
        
        // // X方向重叠检查
        // bool x_overlap = (row_left <= window.x + window.width) && 
        //                 (row_right >= window.x);
        
        if (y_overlap) {
            covered_rows.push_back(row);
        }
    }
    
    if (covered_rows.empty()) {
        return candidates;
    }
    
    // 2. 获取站点间距（使用第一个行的参数）
    if (step_x == 0) {
        step_x = covered_rows[0]->getSiteSpacing().v;
        if (step_x == 0) {
            step_x = covered_rows[0]->getSiteWidth().v;
        }
        if (step_x == 0) {
            step_x = 100;  // 默认值
        }
    }
    
    // 3. 获取单元尺寸
    DBU cell_width = cell->getWidth().v;
    DBU cell_height = cell->getHeight().v;
    
    // 4. 为每个覆盖的行生成候选位置
    for (auto* row : covered_rows) {
        if (candidates.size() >= static_cast<size_t>(num_candidates)) {
            break;
        }
        
        // 获取行参数
        DBU row_bottom = row->getBottom().v;
        DBU row_top = row->getTop().v;
        DBU row_left = row->getLeft().v;
        DBU row_right = row->getRight().v;
        
        // 检查单元高度是否适合行高度
        if (cell_height > (row_top - row_bottom)) {
            continue;  // 单元太高，不适合这个行
        }
        
        // 计算行与搜索窗口的交集
        DBU intersect_left = std::max(window.x, row_left);
        DBU intersect_right = std::min(window.x + window.width, row_right);
        DBU intersect_width = intersect_right - intersect_left;
        
        if (intersect_width <= 0) {
            continue;  // 没有交集
        }

        DBU site_spacing = row->getSiteSpacing().v;
        DBU site_width = row->getSiteWidth().v;
        
        // 5. 找到行中的可用空位
        std::vector<std::pair<DBU, DBU>> available_slots = 
            findAvailableSlotsInRowAligned(row, cell, intersect_left, intersect_right, site_width, site_spacing);
        
        // 6. 在可用空位中生成候选位置
        for (const auto& slot : available_slots) {
            DBU slot_start = slot.first;
            DBU slot_end = slot.second;
            DBU slot_width = slot_end - slot_start;
            
            if (slot_width < cell_width) {
                continue;  // 空位太小，放不下单元
            }
            
            // 计算空位内可能的候选位置数量
            int max_positions_in_slot = static_cast<int>((slot_width - cell_width) / cell_width + 1 );
            
            // 在空位内生成候选位置
            for (int i = 0; i < max_positions_in_slot; ++i) {
                DBU x = slot_start + i * step_x * (cell_width/step_x);
                
                // 确保单元完全在空位内
                if (x + cell_width > slot_end) {
                    break;
                }
                
                // 确保位置在搜索窗口内
                if (x < window.x || x + cell_width > window.x + window.width) {
                    continue;
                }
                
                // 计算Y坐标（单元在行底部）
                DBU y = row_bottom;
                
                // 创建候选位置
                candidates.emplace_back(DbuX{x}, DbuY{y}, cell->location.z);
                
                // 如果达到所需数量，提前退出
                if (candidates.size() >= static_cast<size_t>(num_candidates)) {
                    return candidates;
                }
            }
        }
    }
    
    return candidates;
}

std::vector<Pointdp3D> LSP3D::generateCandidatePositions2(
    Node3D* cell,
    const SearchWindow& window,
    int num_candidates)
{
    std::vector<Pointdp3D> candidates;
    if (num_candidates <= 0 || !cell || !arch_) return candidates;

    DBU W = window.width;
    DBU H = window.height;

    int num_x = std::max(1, (int)std::sqrt(num_candidates * (double)W / H));
    int num_y = std::max(1, num_candidates / num_x);

    DBU step_x = W / num_x;
    DBU step_y = H / num_y;

    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    const auto& rows = arch_->getRows(layer_type);
    if (rows.empty()) return candidates;

    for (int j = 0; j < num_y && (int)candidates.size() < num_candidates; j++) {
        DBU y = window.y + j * step_y;

        Architecture3D::Row3D* nearest_row = nullptr;
        DBU min_dist = std::numeric_limits<DBU>::max();
        for (auto* row : rows) {
            if (!row) continue;
            DBU row_bottom = row->getBottom().v;
            DBU dist = std::abs(y - row_bottom);
            if (dist < min_dist) {
                min_dist = dist;
                nearest_row = row;
            }
        }
        if (!nearest_row) continue;
        DBU snapped_y = nearest_row->getBottom().v;

        DBU site_spacing = nearest_row->getSiteSpacing().v;
        DBU row_left = nearest_row->getLeft().v;
        DBU row_right = nearest_row->getRight().v;

        for (int i = 0; i < num_x && (int)candidates.size() < num_candidates; i++) {
            DBU x = window.x + i * step_x;

            x = alignToSite(x, site_spacing, /*round_up=*/false);
            x = std::max(row_left, x);
            x = std::min(row_right, x);

            candidates.emplace_back(DbuX{x}, DbuY{snapped_y}, DbuZ{0});
        }
    }
    return candidates;
}

/**
 * @brief 时序驱动的搜索窗口定位（主函数）
 * @param cell 要优化的单元
 * @param num_candidates 候选位置数量
 * @param window_width 搜索窗口宽度（如果为0则使用默认值）
 * @param window_height 搜索窗口高度（如果为0则使用默认值）
 * @return 候选位置列表
 */
std::vector<Pointdp3D> LSP3D::timingDrivenWindowSearch(
    Node3D* cell,
    int num_candidates,
    int window_x,
    int window_y){
    
    std::vector<Pointdp3D> candidates;
    DBU window_width, window_height;

    if (!cell || !arch_ || num_candidates <= 0) {
        return candidates;
    }
    
    // 确定搜索窗口尺寸
    if (window_x == 0 || window_y == 0) {
        // 使用默认窗口尺寸
        window_width = cell->getWidth().v * 5;
        window_height = cell->getHeight().v * 3;
    } else{
        window_width = cell->getWidth().v * window_x;
        window_height = cell->getHeight().v * window_y;
    }
    
    // 确定搜索窗口位置
    SearchWindow window = determineSearchWindow(cell, window_width, window_height);
    
    // 确保窗口在合法范围内
    window = clampSearchWindow(window, cell);
    
    // 生成候选位置
    candidates = generateCandidatePositions(cell, window, num_candidates);
    // candidates = generateCandidatePositions2(cell, window, num_candidates);
    
    // 为每个候选位置设置正确的层信息
    std::string cell_layer = cell->layer;
    DBU z_offset = (cell_layer == "bottom") ? 0 : 1000;  // 简化：底层Z=0，上层Z=1000
    
    for (auto& candidate : candidates) {
        candidate.z = DbuZ{z_offset};
    }
    
    return candidates;
}

/**
 * @brief 确保搜索窗口在合法范围内
 */
LSP3D::SearchWindow LSP3D::clampSearchWindow(
    const SearchWindow& window,
    Node3D* cell) {
    
    if (!cell || !arch_) {
        return window;
    }
    
    SearchWindow clamped = window;
    std::string cell_layer = cell->layer;
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(cell_layer);
    
    // 获取层的边界
    DbuX layer_min_x = arch_->getMinX(layer_type);
    DbuX layer_max_x = arch_->getMaxX(layer_type);
    DbuY layer_min_y = arch_->getMinY(layer_type);
    DbuY layer_max_y = arch_->getMaxY(layer_type);
    
    // 调整X方向
    if (clamped.x < layer_min_x.v) {
        clamped.x = layer_min_x.v;
    }
    
    if (clamped.x + clamped.width > layer_max_x.v) {
        clamped.x = layer_max_x.v - clamped.width;
        if (clamped.x < layer_min_x.v) {
            clamped.x = layer_min_x.v;
            clamped.width = layer_max_x.v - layer_min_x.v;
        }
    }
    
    // 调整Y方向
    if (clamped.y < layer_min_y.v) {
        clamped.y = layer_min_y.v;
    }
    
    if (clamped.y + clamped.height > layer_max_y.v) {
        clamped.y = layer_max_y.v - clamped.height;
        if (clamped.y < layer_min_y.v) {
            clamped.y = layer_min_y.v;
            clamped.height = layer_max_y.v - layer_min_y.v;
        }
    }
    
    return clamped;
}


/**
 * @brief 在行中找到可用的空位
 * @param row 要检查的行
 * @param exclude_cell 要排除的单元（当前单元）
 * @param range_start 检查范围的起始X坐标
 * @param range_end 检查范围的结束X坐标
 * @return 可用空位列表（每个空位是[start, end]对）
 */
std::vector<std::pair<DBU, DBU>> LSP3D::findAvailableSlotsInRow(
    Architecture3D::Row3D* row,
    Node3D* exclude_cell,
    DBU range_start,
    DBU range_end) {
    
    std::vector<std::pair<DBU, DBU>> available_slots;
    
    if (!row || !network_) {
        return available_slots;
    }
    
    // 获取行参数
    DBU row_left = row->getLeft().v;
    DBU row_right = row->getRight().v;
    DBU row_bottom = row->getBottom().v;
    DBU row_top = row->getTop().v;
    
    // 调整检查范围，确保在行内
    range_start = std::max(range_start, row_left);
    range_end = std::min(range_end, row_right);
    
    if (range_start >= range_end) {
        return available_slots;
    }
    
    // 1. 获取行中所有单元（排除当前单元）
    std::vector<Node3D*> row_cells;

    // 从network_中获取所有单元
    const auto& all_cells = network_->getCells();
    
    // 过滤出在当前行内的单元
    for (auto* cell : all_cells) {
        if (!cell || cell == exclude_cell) {
            continue;
        }
        
        // 检查单元是否在当前行内
        DBU cell_bottom = cell->location.y.v;
        DBU cell_top = cell_bottom + cell->getHeight().v;
        
        // Y方向检查：单元是否在当前行内
        bool y_in_row = (cell_bottom >= row_bottom) && (cell_top <= row_top);
        
        if (!y_in_row) {
            continue;  // 单元不在当前行内
        }
        
        // 检查单元是否在检查范围内（X方向）
        DBU cell_left = cell->location.x.v;
        DBU cell_right = cell_left + cell->getWidth().v;
        
        // 检查是否有重叠
        bool x_overlap = (cell_right >= range_start) && (cell_left <= range_end);
        
        if (x_overlap) {
            row_cells.push_back(cell);
        }
    }
    
    // 2. 按X坐标排序单元
    std::sort(row_cells.begin(), row_cells.end(),
              [](Node3D* a, Node3D* b) {
                  return a->location.x.v < b->location.x.v;
              });
    
    // 3. 找到占用区域
    std::vector<std::pair<DBU, DBU>> occupied_ranges;
    
    for (auto* cell : row_cells) {
        DBU cell_left = cell->location.x.v;
        DBU cell_right = cell_left + cell->getWidth().v;
        
        // 只考虑在检查范围内的单元
        if (cell_right <= range_start || cell_left >= range_end) {
            continue;
        }
        
        // 调整到检查范围内
        cell_left = std::max(cell_left, range_start);
        cell_right = std::min(cell_right, range_end);
        
        if (cell_left < cell_right) {
            occupied_ranges.emplace_back(cell_left, cell_right);
        }
    }
    
    // 4. 合并重叠的占用区域
    std::vector<std::pair<DBU, DBU>> merged_ranges;
    
    if (!occupied_ranges.empty()) {
        std::sort(occupied_ranges.begin(), occupied_ranges.end());
        
        DBU current_start = occupied_ranges[0].first;
        DBU current_end = occupied_ranges[0].second;
        
        for (size_t i = 1; i < occupied_ranges.size(); ++i) {
            if (occupied_ranges[i].first <= current_end) {
                // 重叠，合并
                current_end = std::max(current_end, occupied_ranges[i].second);
            } else {
                // 不重叠，保存当前区域
                merged_ranges.emplace_back(current_start, current_end);
                current_start = occupied_ranges[i].first;
                current_end = occupied_ranges[i].second;
            }
        }
        
        // 保存最后一个区域
        merged_ranges.emplace_back(current_start, current_end);
    }
    
    // 5. 从占用区域推导出可用空位
    DBU current_pos = range_start;
    
    for (const auto& occupied : merged_ranges) {
        if (current_pos < occupied.first) {
            // 找到可用空位
            available_slots.emplace_back(current_pos, occupied.first);
        }
        current_pos = std::max(current_pos, occupied.second);
    }
    
    // 检查最后一个区域之后的空位
    if (current_pos < range_end) {
        available_slots.emplace_back(current_pos, range_end);
    }
    
    return available_slots;
}

/**
 * @brief 在行中找到可用的空位（考虑站点对齐）
 * @param row 要检查的行
 * @param exclude_cell 要排除的单元（当前单元）
 * @param range_start 检查范围的起始X坐标
 * @param range_end 检查范围的结束X坐标
 * @param site_width 站点宽度
 * @param site_spacing 站点间距
 * @return 可用空位列表（每个空位是[start, end]对，已对齐到站点）
 */
std::vector<std::pair<DBU, DBU>> LSP3D::findAvailableSlotsInRowAligned(
    Architecture3D::Row3D* row,
    Node3D* exclude_cell,
    DBU range_start,
    DBU range_end,
    DBU site_width,
    DBU site_spacing) {
    
    std::vector<std::pair<DBU, DBU>> available_slots;
    
    if (!row || !network_) {
        return available_slots;
    }
    
    // 获取站点参数
    if (site_spacing == 0) {
        site_spacing = row->getSiteSpacing().v;
    }
    if (site_width == 0) {
        site_width = row->getSiteWidth().v;
    }
    
    if (site_spacing == 0) {
        site_spacing = 100;  // 默认值
    }
    
    // 1. 先找到原始的空位
    auto raw_slots = findAvailableSlotsInRow(row, exclude_cell, range_start, range_end);
    
    // 2. 对每个空位进行站点对齐
    for (const auto& raw_slot : raw_slots) {
        DBU slot_start = raw_slot.first;
        DBU slot_end = raw_slot.second;
        // DBU slot_width = slot_end - slot_start;
        
        // 对齐到站点边界
        DBU aligned_start = alignToSite(slot_start, site_spacing, true);  // 向上对齐
        DBU aligned_end = alignToSite(slot_end, site_spacing, false);     // 向下对齐
        
        // 确保对齐后的空位有效
        if (aligned_start < aligned_end) {
            available_slots.emplace_back(aligned_start, aligned_end);
        }
    }
    
    return available_slots;
}

/**
 * @brief 将坐标对齐到站点边界
 * @param coord 要对齐的坐标
 * @param site_spacing 站点间距
 * @param round_up 是否向上对齐（true=向上，false=向下）
 * @return 对齐后的坐标
 */
DBU LSP3D::alignToSite(DBU coord, DBU site_spacing, bool round_up) const {
    if (site_spacing == 0) {
        return coord;
    }
    
    DBU remainder = coord % site_spacing;
    
    if (remainder == 0) {
        return coord;  // 已经在站点边界上
    }
    
    if (round_up) {
        return coord + (site_spacing - remainder);
    } else {
        return coord - remainder;
    }
}

// ============================================================================
// 翻转方法
// ============================================================================

/**
 * @brief 翻转单元（使用MX方向）
 * @param cell 要翻转的单元
 */
void LSP3D::flipCell(Node3D* cell) const {
    if (!cell) {
        return;
    }
    
    // 获取当前方向
    odb::dbOrientType current_orient = cell->getOrient();
    
    // 根据当前方向计算翻转后的方向
    odb::dbOrientType new_orient = current_orient;
    
    // 应用MX翻转（关于X轴镜像）
    switch (current_orient) {
        case odb::dbOrientType::R0:
            new_orient = odb::dbOrientType::MX;
            break;
        case odb::dbOrientType::R90:
            new_orient = odb::dbOrientType::MXR90;
            break;
        case odb::dbOrientType::R180:
            new_orient = odb::dbOrientType::R0;  // MX + R180 = R0
            break;
        case odb::dbOrientType::R270:
            new_orient = odb::dbOrientType::MYR90;  // MX + R270 = MYR90
            break;
        case odb::dbOrientType::MX:
            new_orient = odb::dbOrientType::R0;  // 两次MX = R0
            break;
        case odb::dbOrientType::MY:
            new_orient = odb::dbOrientType::R180;  // MX + MY = R180
            break;
        case odb::dbOrientType::MYR90:
            new_orient = odb::dbOrientType::R270;  // MX + MYR90 = R270
            break;
        case odb::dbOrientType::MXR90:
            new_orient = odb::dbOrientType::R90;  // MX + MXR90 = R90
            break;
        default:
            new_orient = odb::dbOrientType::MX;  // 默认使用MX
            break;
    }
    
    // 设置新方向
    cell->setOrient(new_orient);
    
    // std::cout << "  Cell " << cell->name << " flipped from " 
    //           << orientToString(current_orient) << " to " 
    //           << orientToString(new_orient) << std::endl;
}

/**
 * @brief 将方向转换为字符串
 */
std::string LSP3D::orientToString(const odb::dbOrientType& orient) const {
    switch (orient) {
        case odb::dbOrientType::R0: return "R0";
        case odb::dbOrientType::R90: return "R90";
        case odb::dbOrientType::R180: return "R180";
        case odb::dbOrientType::R270: return "R270";
        case odb::dbOrientType::MY: return "MY";
        case odb::dbOrientType::MYR90: return "MYR90";
        case odb::dbOrientType::MX: return "MX";
        case odb::dbOrientType::MXR90: return "MXR90";
        default: return "UNKNOWN";
    }
}



// ============================================================================
// 成本计算方法
// ============================================================================
double LSP3D::calculateCellCost(Node3D* cell, double alpha, double omega,
                               TimingEngine3D* timing_engine) const {
    if (!cell || !network_) {
        return 0.0;
    }
    
    //计算线长成本
    // double wirelength_cost = calculateCellHPWL(cell);
    
    //计算时序成本
    double timing_cost = calculateCellTimingCost(cell, alpha, omega, timing_engine_);
    
    //组合成本（可以根据需要调整权重）
    //double total_cost = wirelength_cost + timing_cost;
    
    return timing_cost;
}

std::vector<Pin3D*> LSP3D::getCellInputPins(Node3D* cell) const {
    std::vector<Pin3D*> input_pins;
    
    if (!cell || !network_) {
        return input_pins;
    }
    
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_IN) {
            input_pins.push_back(pin);
        }
    }
    
    return input_pins;
}

std::vector<Pin3D*> LSP3D::getCellOutputPins(Node3D* cell) const {
    std::vector<Pin3D*> output_pins;
    
    if (!cell || !network_) {
        return output_pins;
    }
    
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_OUT) {
            output_pins.push_back(pin);
        }
    }
    
    return output_pins;
}

void LSP3D::cellArcDelays(Node3D* cell) const {
    if (!network_ || !timing_engine_) return;
    if (!cell || cell->is_fixed) return;
    auto output_pins = getCellOutputPins(cell);
    auto input_pin = getCellInputPins(cell);
    if (input_pin.empty()) return;

    for (auto* out_pin : output_pins) {
        double min_early = 0.0 , max_late = 0.0;
        for (const auto& in_name : input_pin) {
            auto info = timing_engine_->getArcDelay(in_name->getName(), out_pin->getName());
            if (info.is_valid) {
                min_early = std::min(min_early, info.early_delay);
                max_late  = std::max(max_late,  info.late_delay);
            }
        }
        out_pin->setCachedArcDelay(min_early, max_late);
    }

    for (auto* in_pin : input_pin) {
        double min_early = 0.0 , max_late = 0.0;
        for (const auto& out_pin : output_pins) {
            auto info = timing_engine_->getArcDelay(in_pin->getName(), out_pin->getName());
            if (info.is_valid) {
                min_early = std::min(min_early, info.early_delay);
                max_late  = std::max(max_late,  info.late_delay);
            }
        }
        in_pin->setCachedArcDelay(min_early, max_late);
    }
    // std::cout << "Cached arc delays for " << cached << " output pins" << std::endl;
}
/**
 * @brief 计算单元的时序成本（包含单元内部弧延迟）
 * @param cell 要计算时序成本的单元
 * @param alpha 晚时序权重
 * @param omega 早时序权重
 * @param timing_engine 时序引擎
 * @return 时序成本
 */
double LSP3D::calculateCellTimingCost(Node3D* cell, double alpha, double omega, 
                                     TimingEngine3D* timing_engine) const {
    if (!cell || !network_) {
        return 0.0;
    }
    
    double cost_late = 0.0;
    double cost_early = 0.0;
    
    cellArcDelays(cell);
    // //对每个sink cell的arc进行更新
    std::unordered_set<Node3D*> visited;
    auto out_pins = getCellOutputPins(cell);
    for (auto* out_pin : out_pins) {
        Net3D* net = out_pin->getNet();
        if (!net) continue;
        for (auto* sink_pin : net->getPins()) {
            if (sink_pin == out_pin) continue;
            Node3D* sink_cell = sink_pin->getNode();
            if (sink_cell && visited.insert(sink_cell).second) {
                cellArcDelays(sink_cell);
            }
        }
    }
    
    // 获取单元的所有引脚
    auto cell_pins = cell->getPins();
    
    for (auto* pin : cell_pins) {
        if (!pin) continue;
        
        // 获取引脚方向
        int dir_int = pin->getDirection();
        Pin3D::Direction dir = static_cast<Pin3D::Direction>(dir_int);

        Net3D* net = pin->getNet();
        if (net && net->db_net && timingdp_->isIdealClockNet(net->db_net)) {
            continue;
        }
        
        if (dir == Pin3D::Dir_IN) {
            // 输入引脚：计算从驱动引脚到当前引脚的延迟
            // Net3D* net = pin->getNet();
            if (!net) continue;
            
            // 找到驱动引脚（source）
            Pin3D* driver_pin = findDriverPin(net, pin);
            if (!driver_pin) continue;

            // 获取单元的内部弧延迟
            double arc_delay_early = 0.0;
            double arc_delay_late = 0.0;
            
            if (pin->isArcCached()) {
                arc_delay_early = pin->getCachedArcEarly();
                arc_delay_late  = pin->getCachedArcLate();
            }
            
            // 计算延迟（线网延迟 + 驱动单元的内部弧延迟）
            double delay_early = pin->getEarlyArrival() - driver_pin->getEarlyArrival() + arc_delay_early;
            double delay_late = pin->getLateArrival() - driver_pin->getLateArrival() + arc_delay_late;
            // double delay_early = pin->getEarlyArrival() - driver_pin->getEarlyArrival();
            // double delay_late = pin->getLateArrival() - driver_pin->getLateArrival();
            // 获取权重
            double early_weight = pin->getEarlyWeight();
            double late_weight = pin->getLateWeight();
            
            cost_early += delay_early * early_weight;
            cost_late += delay_late * late_weight;
            
        } else if (dir == Pin3D::Dir_OUT) {
            // 输出引脚：计算到所有接收引脚的延迟
            // Net3D* net = pin->getNet();
            if (!net) continue;
            
            auto net_pins = net->getPins();
            
            for (auto* sink_pin : net_pins) {
                if (sink_pin == pin) continue;  // 跳过驱动引脚
                
                // 计算线网延迟
                double net_delay_early = sink_pin->getEarlyArrival() - pin->getEarlyArrival();
                double net_delay_late = sink_pin->getLateArrival() - pin->getLateArrival();
                
                // 获取接收单元的内部弧延迟
                double arc_delay_early = 0.0;
                double arc_delay_late = 0.0;

                if(!isSequentialCell(sink_pin->getNode())){
                    arc_delay_early = sink_pin->getCachedArcEarly();
                    arc_delay_late = sink_pin->getCachedArcLate();
                }

            
                // 总延迟 = 线网延迟 + 单元内部弧延迟
                double total_delay_early = net_delay_early + arc_delay_early;
                double total_delay_late = net_delay_late + arc_delay_late;
                
                // 获取权重
                double early_weight = sink_pin->getEarlyWeight();
                double late_weight = sink_pin->getLateWeight();
                
                cost_early += total_delay_early * early_weight;
                cost_late += total_delay_late * late_weight;
                
            }
        }
    }
    
    // 组合时序成本
    double total_timing_cost = alpha * cost_late - omega * cost_early;
    
    return total_timing_cost;
}


// ============================================================================
// LR算法特定方法
// ============================================================================
//local search
/**
 * @brief 执行离散局部搜索（针对单个单元）
 * @param cell 要优化的单元
 * @param search_windows 搜索窗口数量（默认33）
 * @param alpha 晚时序权重（默认5.0）
 * @param omega 早时序权重（默认1.0）
 * @return 是否找到了更好的位置
 */
bool LSP3D::discreteLocalSearch(Node3D* cell, int search_windows, double alpha, double omega,
                               TimingEngine3D* timing_engine) {
    if (!cell || cell->is_fixed || !arch_ || !network_) {
        return false;
    }
    double time_sta = 0.0;
    double time_fresh = 0.0;
    // 获取当前单元信息
    Pointdp3D current_location = cell->location;
    std::string current_layer = cell->layer;
    odb::dbOrientType current_orient = cell->getOrient();  // 保存当前方向
    //生成候选位置
    std::vector<Pointdp3D> candidates = timingDrivenWindowSearch(cell, search_windows, params_.window_x, params_.window_y);
    if (candidates.empty()) {
        // std::cerr << "Warning: No candidate positions found for cell " << cell->name << std::endl;
        return false;
    }
    // 2. 计算当前位置的成本
    // ========= 解析模型模式 =========
    // computeReferenceDelays(cell);
    // double current_cost = calculateCellCostAnalytical(cell, alpha, omega);
    // ========= STA 模式（已注释） =========
    SteinerTreeGenerator steiner_generator(timingdp_);
    // steiner_generator.writeSpef("tdp_initial");
    double current_cost = calculateCellCost(cell, alpha, omega, timing_engine_);
    // 获取行信息
    Architecture3D::LayerType layer_type = Architecture3D::stringToLayer(current_layer);
    // 获取当前层的所有行
    const auto& rows = arch_->getRows(layer_type);
    if (rows.empty()) {
        std::cerr << "Error: No rows found for layer " << current_layer << std::endl;
        return false;
    }
    Architecture3D::Row3D* row = findCellRow(cell, rows);
    if (!row) {
        std::cerr << "Error: Cannot find row for cell " << cell->name 
                << " at (" << current_location.x.v << ", " << current_location.y.v 
                << ") on layer " << current_layer << std::endl;
        return false;
    }
    
    // 3. 评估每个候选位置
    double best_cost = current_cost;
    Pointdp3D best_location = current_location;
    bool found_better = false;
    odb::dbOrientType candidate_orient = current_orient;
    odb::dbOrientType best_orient = current_orient;
    odb::dbInst* inst = cell->db_inst;

    // ========= 解析模型模式：不修改路由，不需要 global_router =========
    auto* global_router = timingdp_->getGlobalRouter();
    auto& global_routes = global_router->getRoutes();
    std::vector<odb::dbNet*> affected_nets = cell->connected_nets;
    // ========= 解析模型模式：不修改路由，不需要保存/恢复 =========
    std::map<odb::dbNet*, grt::GRoute> saved_routes;
    for (auto* db_net : affected_nets) {
        if (!db_net) continue;
        auto it = global_routes.find(db_net);
        if (it != global_routes.end())
            saved_routes[db_net] = it->second;
    }

    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& candidate = candidates[i];
        
        // 临时移动单元到候选位置
        cell->location = candidate;
        bool need = (std::abs(((candidate.y.v - row->getBottom().v) /  row->getHeight().v )% 2) == 1);
        
        
        if (!need) {
            // 保存当前方向
            odb::dbOrientType saved_orient = cell->getOrient();
            
            // 执行翻转
            flipCell(cell);
            candidate_orient = cell->getOrient();

            // 恢复方向（用于后续恢复）
            cell->setOrient(saved_orient);

        }
        // 计算候选位置的成本
        timingdp_->setInstanceLocation3D(inst, cell->location.x.v, cell->location.y.v, cell->location.z.v);

        // ========= STA 模式（已注释） =========
        steiner_generator.updateSteinerTrees(cell);
        // steiner_generator.writeSpef("LSP_output");
        auto t1 = std::chrono::high_resolution_clock::now();
        timingdp_->executeTimingDrivenI(cell, 2);
        auto t2 = std::chrono::high_resolution_clock::now();
        refreshAffectedPinTiming(cell);
        compareLocalVsFullSTA(cell);
        auto t3 = std::chrono::high_resolution_clock::now();
        double sta_ms    = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double refresh_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
        time_sta += sta_ms;
        time_fresh += refresh_ms;
        double candidate_cost = calculateCellCost(cell, alpha, omega, timing_engine_);
        
        // 检查是否找到更好的位置
        if (candidate_cost < best_cost) {
            best_cost = candidate_cost;
            best_location = candidate;
            best_orient = candidate_orient;
            found_better = true;
            
            // std::cout << "  -> Found better position! New cost: " << best_cost 
            //           << " (improvement: " << (current_cost - best_cost) << ")" << std::endl;
        }
        // Restore original routes for next candidate
        // ========= STA 模式（已注释） =========
        for (const auto& [db_net, route] : saved_routes)
            global_routes[db_net] = route;
    }
    // logger_ ->info(TDP, 881, "sta_inter time{}", time_sta);
    // logger_ ->info(TDP, 882, "fresh_inter time{}", time_fresh);
    // 4. 恢复单元位置
    cell->location = current_location;
    // 恢复方向
    cell->setOrient(current_orient);

    // 5. 如果找到更好的位置，应用它
    if (found_better) {
        cell->location = best_location;
        cell->setOrient(best_orient);
        timingdp_->setInstanceLocation3D(inst, best_location.x.v, best_location.y.v, best_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        // sta_->updateTiming(false);
        timingdp_->executeTimingDrivenI(cell, 1);
        refreshAffectedPinTiming(cell);
        return true;
    }
    else{
        timingdp_->setInstanceLocation3D(inst, current_location.x.v,
                                  current_location.y.v,
                                  current_location.z.v);
        steiner_generator.updateSteinerTrees(cell);
        // sta_->updateTiming(false);
    }
    
    return false;
}

// ============================================================================
// prepareSnapshot — Phase 0: collect boundary conditions from STA graph
// ============================================================================

void LSP3D::prepareSnapshot(Node3D* cell, LocalTimingSnapshot& out)
{
    if (!local_timing_) return;
    local_timing_->buildSnapshot(cell, out);
}

// ============================================================================
// discreteLocalSearchParallel — 3-phase parallel candidate evaluation
// ============================================================================

bool LSP3D::discreteLocalSearchParallel(Node3D* cell, int search_windows,
                                        double alpha, double omega,
                                        TimingEngine3D* timing_engine)
{
    if (!cell || cell->is_fixed || !arch_ || !network_ || !local_timing_)
        return false;

    // --- save current state ---
    Pointdp3D current_location = cell->location;
    odb::dbOrientType current_orient = cell->getOrient();

    // --- generate candidates (unchanged) ---
    std::vector<Pointdp3D> candidates = timingDrivenWindowSearch(
        cell, search_windows, params_.window_x, params_.window_y);
    if (candidates.empty())
        return false;

    // --- compute current cost (runs STA once) ---
    SteinerTreeGenerator steiner_gen(timingdp_);
    double current_cost = calculateCellCost(cell, alpha, omega, timing_engine_);

    // --- get row / inst ---
    Architecture3D::LayerType layer_type =
        Architecture3D::stringToLayer(cell->layer);
    const auto& rows = arch_->getRows(layer_type);
    Architecture3D::Row3D* row = findCellRow(cell, rows);
    if (!row) return false;
    odb::dbInst* inst = cell->db_inst;

    // ====================================================================
    // Phase 0: snapshot boundary conditions
    // ====================================================================
    LocalTimingSnapshot snap;
    prepareSnapshot(cell, snap);

    // ====================================================================
    // Phase 1: parallel candidate evaluation (no STA, no global writes)
    // ====================================================================
    int num_cand = static_cast<int>(candidates.size());
    std::vector<double> costs(num_cand, std::numeric_limits<double>::max());
    std::vector<int> cand_need_flip(num_cand, 0);

    // Pre-compute flip info
    for (int i = 0; i < num_cand; ++i) {
        bool need = (std::abs(((candidates[i].y.v - row->getBottom().v)
                               / row->getHeight().v) % 2) == 1);
        cand_need_flip[i] = need ? 1 : 0;
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_cand; ++i) {
        costs[i] = local_timing_->evaluateCost(
            snap, candidates[i], alpha, omega);
    }

    // ====================================================================
    // Phase 2: apply best candidate with full STA verification
    // ====================================================================
    int best_idx = 0;
    double best_cost = current_cost;
    for (int i = 0; i < num_cand; ++i) {
        if (costs[i] < best_cost) {
            best_cost = costs[i];
            best_idx = i;
        }
    }

    if (best_cost < current_cost) {
        // Apply best position
        cell->location = candidates[best_idx];
        if (cand_need_flip[best_idx])
            flipCell(cell);

        timingdp_->setInstanceLocation3D(
            inst, cell->location.x.v, cell->location.y.v, cell->location.z.v);
        steiner_gen.updateSteinerTrees(cell);
        timingdp_->executeTimingDrivenI(cell);
        refreshAffectedPinTiming(cell);
        return true;
    } else {
        // Restore
        cell->location = current_location;
        cell->setOrient(current_orient);
        timingdp_->setInstanceLocation3D(
            inst, current_location.x.v, current_location.y.v,
            current_location.z.v);
        steiner_gen.updateSteinerTrees(cell);
    }

    return false;
}

/**
 * @brief 执行全局离散局部搜索（遍历所有单元）
 */
void LSP3D::discreteLocalSearchAll(TimingEngine3D* timing_engine) {
    if (!network_) {
        std::cerr << "Error: Network not initialized" << std::endl;
        return;
    }
    
    int moved_cells = 0;
    int total_cells = 0;
    
    auto priority_order = buildCriticalPathPriorityOrder(timing_engine, 200);

    int total_movable = std::count_if(priority_order.begin(), priority_order.end(),
        [](Node3D* c) { return c && !c->is_fixed; });
    
    const int bar_width = 40;
    int report_step = std::max(1, total_movable / 1000);

    double totaltime = 0.0;
    
    // 遍历所有单元（关键路径优先）
    for (auto* cell : priority_order) {
        if (!cell || cell->is_fixed || isSequentialCell(cell)) {
            // 固定单元保持原位
            continue;
        }
        
        total_cells++;
        auto t4 = std::chrono::high_resolution_clock::now();
        // 执行局部搜索 - parallel version
        // if (discreteLocalSearchParallel(cell, 22, 5.0, 1.0, timing_engine_)) {
        if (discreteLocalSearch(cell, 20, 5.0, 1.0,timing_engine_)) {
            moved_cells++;
        }
        auto t5 = std::chrono::high_resolution_clock::now();

        double dis_ms    = std::chrono::duration<double, std::milli>(t5 - t4).count();
        // logger_ ->info (TDP, 1899, "one inter time: {}", dis_ms);

        totaltime += dis_ms;

        // if (total_cells == total_movable/2) {
        //     compareLocalVsFullSTA(cell);
        // }
        // else if (total_cells == total_movable -1) {
        //     compareLocalVsFullSTA(cell);
        // }
        
        if (total_cells % report_step == 0 || total_cells == total_movable) {
            double pct = 100.0 * total_cells / total_movable;
            int pos = bar_width * total_cells / total_movable;
            std::cout << "\r进度: [" << std::string(pos, '#')
                      << std::string(bar_width - pos, '.') << "] "
                      << std::fixed << std::setprecision(1)
                      << std::setw(5) << pct << "% ("
                      << total_cells << "/" << total_movable << ")"
                      << std::flush;
        }
    }
    
    std::cout << std::endl;
    std::cout << "Discrete local search completed. Moved " 
              << moved_cells << " out of " << total_cells << " movable cells" << std::endl;
    logger_ ->info (TDP, 1991, "total time: {}", totaltime);
}

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * 查找网络的驱动引脚
 */
Pin3D* LSP3D::findDriverPin(Net3D* net, Pin3D* exclude_pin) const {
    if (!net) return nullptr;
    
    auto net_pins = net->getPins();
    
    for (auto* pin : net_pins) {
        if (pin == exclude_pin) continue;
        
        if (pin->getDirection() == Pin3D::Dir_OUT) {
            return pin;  // 找到输出引脚（驱动）
        }
    }
    
    return nullptr;
}

// ============================================================================
// 时序缓存刷新
// ============================================================================
void LSP3D::refreshAffectedPinTiming(Node3D* cell) {
    refreshPinTiming2(cell);
    for (auto* pin : cell->getPins()) {
        if (pin->getDirection() != Pin3D::Dir_OUT) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;
        for (auto* net_pin : net->getPins()) {
            Node3D* pin_cell = net_pin->getNode();
            if (pin_cell && pin_cell != cell) {
                refreshPinTiming2(pin_cell);
            }
        }
    }
}

// ============================================================================
// 解析模型辅助方法
// ============================================================================

double LSP3D::dbuToMeters(DBU dbu) const {
    int dbu_per_micron = timingdp_->getBlock()->getDbUnitsPerMicron();
    return static_cast<double>(dbu) / dbu_per_micron * 1e-6;
}

double LSP3D::getPinCapacitance(Pin3D* pin) const {
    if (!pin) return 0.0;
    sta::Pin* sta_pin = sta_->getDbNetwork()->findPin(pin->getName());
    if (!sta_pin) return 0.0;
    sta::LibertyPort* lport = sta_->network()->libertyPort(sta_pin);
    if (!lport) return 0.0;
    return static_cast<double>(lport->capacitance());  // 法拉
}

// ============================================================================
// 解析模型 cost 计算（论文 Algorithm 1）
// ============================================================================


// ============================================================================
// 捕获 STA 参考值（每 iteration/每次接受后调一次）
// ============================================================================

void LSP3D::computeReferenceDelays(Node3D* cell) {
    if (!cell || !timing_engine_) return;

    refreshAffectedPinTiming(cell);

    constexpr double KD = 0.69;
    int l_idx = getLayerIndex(cell->layer);
    double R = kLayerRC[l_idx].r_per_m;
    double C = kLayerRC[l_idx].c_per_m;

    for (auto* pin : cell->getPins()) {
        if (!pin) continue;
        pin->ref_late_delay = pin->getLateArrival();
        pin->ref_early_delay = pin->getEarlyArrival();
    }

    // fanout sink pin 参考值
    auto output_pins = getCellOutputPins(cell);
    for (auto* pin : output_pins) {
        if (!pin) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;

        for (auto* sink_pin : net->getPins()) {
            if (sink_pin == pin) continue;
            Node3D* sink_cell = sink_pin->getNode();
            if (!sink_cell) continue;

            sink_pin->ref_late_delay = sink_pin->getLateArrival();
            sink_pin->ref_early_delay = sink_pin->getEarlyArrival();

            double l = std::abs(static_cast<double>(
                (cell->location.x.v - sink_cell->location.x.v) +
                (cell->location.y.v - sink_cell->location.y.v)));
            double l_m = dbuToMeters(static_cast<DBU>(l));
            double C_pin = getPinCapacitance(sink_pin);
            sink_pin->ref_elmore_delay = KD * R * l_m * (C * l_m / 2.0 + C_pin);
        }
    }

    // fanin driver pin 的 HPWL 参考值
    auto input_pins = getCellInputPins(cell);
    for (auto* pin : input_pins) {
        if (!pin) continue;
        Net3D* net = pin->getNet();
        if (!net) continue;
        Pin3D* driver = findDriverPin(net, pin);
        Node3D* drv_cell = driver ? driver->getNode() : nullptr;
        if (!drv_cell) continue;

        double l = std::abs(static_cast<double>(
            (cell->location.x.v - drv_cell->location.x.v) +
            (cell->location.y.v - drv_cell->location.y.v)));
        pin->ref_hpwl = l;
    }
}

// ============================================================================
// 网络相关辅助方法
// ============================================================================
void LSP3D::initWeight(Timingdp* timingdp, TimingEngine3D* timing_engine)
{

    auto* network_3d = timingdp->getNetwork3D();
        if (!network_3d) {
            std::cout << "Error: Network3D not available" << std::endl;
        }else{

            auto& nets = timingdp->getNets();

            for (auto* db_net : nets) {
                auto* net = network_3d->getNet(db_net);
                if(!net){
                std::cout << "  db_net pointer: " << db_net 
                  << ", name: " << db_net->getName() 
                  << ", signal type: " << (db_net->getSigType().getString()) << std::endl;
                    continue;
                }
                auto pins = net->getPins();
                
                for (auto* pin : pins) {
                    
                    // 设置默认权重为1.0
                    pin->setEarlyWeight(1.0);
                    pin->setLateWeight(1.0);
                    pin->setEarlyWeightini(1.0);
                    pin->setLateWeightini(1.0);
                    
                    // 获取引脚名称
                    std::string pin_name = pin->getName();
                    
                    // 尝试从timing_engine获取时序信息
                    try {
                        // 获取引脚时序信息
                        auto timing_info = timing_engine_->getPinTimingInfo(pin_name);
                        
                        if (timing_info.is_valid) {
                            // 设置真实的时序信息
                            pin->setEarlySlack(timing_info.early_slack);
                            pin->setLateSlack(timing_info.late_slack);
                            pin->setEarlyArrival(timing_info.early_arrival);
                            pin->setLateArrival(timing_info.late_arrival);
                            pin->setEarlyRequired(timing_info.early_required);
                            pin->setLateRequired(timing_info.late_required);
                        } else {
                            // 如果无法获取时序信息，使用默认值
                            pin->setEarlySlack(0.0);
                            pin->setLateSlack(0.0);
                            pin->setEarlyArrival(0.0);
                            pin->setLateArrival(0.0);
                            pin->setEarlyRequired(0.0);
                            pin->setLateRequired(0.0);
                            std::cout << "Error timing info "<< std::endl;
                        }
                    } catch (const std::exception& e) {
                        // 如果获取失败，使用默认值
                        // std::cerr << "Warning: Failed to get timing info for pin " 
                                // << pin_name << ": " << e.what() << std::endl;
                        pin->setEarlySlack(0.0);
                        pin->setLateSlack(0.0);
                        pin->setEarlyArrival(0.0);
                        pin->setLateArrival(0.0);
                        pin->setEarlyRequired(0.0);
                        pin->setLateRequired(0.0);
                    }
                }
            }
        }

}

} // namespace tdp
