// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#pragma once

#include "tdp/Timingdp.h"
#include "utl/Logger.h"
#include "sta/Sta.hh"
#include "../eva/TimingEngine3D.h"
#include "../eva/LocalTimingPropagator.h"
#include "../infrastructure/Architecture3D.h"
#include "../infrastructure/Network3D.h"

namespace tdp {

struct LRParams {
    double lambda_wirelength = 1.0;
    double lambda_timing = 1.0;
    double lambda_alignment = 1.0;
    int max_iterations = 20;
    double convergence_threshold = 1e-3;
    int    window_y = 4;
    int    window_x = 4;
};
/**
 * @class Placement3D
 * @brief Base class for 3D placement algorithms
 * 
 * Following DPL design pattern, algorithms should be implemented
 * in separate classes that use Timingdp for data access.
 */
class Placement3D
{
public:
    Placement3D(Timingdp* timingdp);
    virtual ~Placement3D();
    
    static Architecture3D::Row3D* findCellRow(Node3D* cell, 
                                            const std::vector<Architecture3D::Row3D*>& rows);
    // 倒序遍历
    void buildTraverseOrder(Timingdp* timingdp);
    void rundetailedplacement();
    void verifyTraverseOrder() const;

    // 刷新单个单元的时序缓存（从 STA 拉取最新值）
    void refreshPinTiming(Node3D* cell);
    void refreshPinTiming2(Node3D* cell);
    void compareLocalVsFullSTA(Node3D* cell);
    bool isSequentialCell(const Node3D* cell) const;
    bool isPrimaryInput(const Node3D* cell) const;
    bool isPrimaryOutput(const Node3D* cell) const;

protected:
    Timingdp* timingdp_ = nullptr;
    sta::dbSta* sta_ = nullptr;
    utl::Logger* logger_ = nullptr;
    Architecture3D* arch_ = nullptr;
    Network3D* network_ = nullptr;
    TimingEngine3D* timing_engine_ = nullptr;
    std::list<Node3D*> traverse_list;

    bool isInputPin(const Pin3D* pin);
    bool isOutputPin(const Pin3D* pin);
    bool isPrimaryInput(const Pin3D* pin);
    bool isPrimaryOutput(const Pin3D* pin);
};

/**
 * @class LSP3D
 * @brief Timing-driven 3D placement algorithm using LR method
 * 
 * This algorithm:
 * 1. Considers both upper and lower layers simultaneously
 * 2. Moves cells based on cost sensitivity minimization
 * 3. Uses LR (Lagrangian Relaxation) method
 * 4. Handles inter-layer timing optimization
 */
class LSP3D : public Placement3D
{
public:

    
    LSP3D(Timingdp* timingdp, const LRParams& params = LRParams());
    ~LSP3D() override;
    
    // 主要布局方法
    void place();

    // 窗口寻找相关方法
    std::vector<Pointdp3D> timingDrivenWindowSearch(
        Node3D* cell,
        int num_candidates = 33,
        int window_width = 0,
        int window_height = 0);

    // 时序权重初始化
    void initWeight(Timingdp* timingdp, TimingEngine3D* timing_engine);


    // 成本计算方法
    double calculateCellCost(Node3D* cell, double alpha, double omega,
                             TimingEngine3D* timing_engine = nullptr) const;
    
    double calculateCellTimingCost(Node3D* cell, double alpha, double omega,
                                   TimingEngine3D* timing_engine) const;
    
    // 局部搜索方法
    bool discreteLocalSearch(Node3D* cell, int search_windows, double alpha, double omega,
                             TimingEngine3D* timing_engine = nullptr);
    
    void discreteLocalSearchAll(TimingEngine3D* timing_engine = nullptr);

    // 构建关键路径优先遍历顺序（top-k 违反约束路径上的 cell 优先）
    std::vector<Node3D*> buildCriticalPathPriorityOrder(
        TimingEngine3D* timing_engine,
        int top_k);

    // Parallel local search: Phase 0 + Phase 1 (OpenMP) + Phase 2
    void prepareSnapshot(Node3D* cell, LocalTimingSnapshot& out);
    bool discreteLocalSearchParallel(Node3D* cell, int search_windows,
                                     double alpha, double omega,
                                     TimingEngine3D* timing_engine = nullptr);

    // 辅助函数
    Pin3D* findDriverPin(Net3D* net, Pin3D* exclude_pin) const;
    std::vector<Pin3D*> getCellInputPins(Node3D* cell) const;
    void flipCell(Node3D* cell) const;
    std::string orientToString(const odb::dbOrientType& orient) const;
    // 权重更新
    void updatePinWeights();
    void balanceCellWeights(Node3D* cell);
    void distributeCellWeights(Node3D* cell);
    // cell Arc delay 计算
    void cellArcDelays(Node3D* cell) const;
    // 辅助函数
    std::vector<Pin3D*> getCellOutputPins(Node3D* cell) const;

    // 刷新受影响的 pin 时序缓存（cell + 其 net 上所有 pin）
    void refreshAffectedPinTiming(Node3D* cell);

    
    // 解析模型：为 cell 捕获 STA 参考值（每 iteration/每次接受后调一次）
    void computeReferenceDelays(Node3D* cell);

    // 解析模型辅助：DBU → 米
    double dbuToMeters(DBU dbu) const;
    // 解析模型辅助：获取 pin 的 Liberty 输入电容（法拉）
    double getPinCapacitance(Pin3D* pin) const;


    
private:
    LRParams params_;
    int current_iteration_ = 0;
    double current_cost_ = 0.0;
    double previous_cost_ = 0.0;

    std::unique_ptr<LocalTimingPropagator> local_timing_;
    // Timingdp* timingdp_ = nullptr;
    // sta::dbSta* sta_ = nullptr;
    // utl::Logger* logger_;
    
    // 敏感度数据
    struct CellSensitivity {
        double timing_sensitivity = 0.0;
        double wirelength_sensitivity = 0.0;
        double alignment_sensitivity = 0.0;
        double total_sensitivity = 0.0;
    };
    
    std::unordered_map<Node3D*, CellSensitivity> sensitivities_;
    // 辅助数据结构
    struct TimingVector {
        double x;
        double y;
        double magnitude;
        TimingVector() : x(0.0), y(0.0), magnitude(0.0) {}
        TimingVector(double x_val, double y_val, double mag) 
            : x(x_val), y(y_val), magnitude(mag) {}
    };
    
    struct DirectionComponents {
        double left;
        double right;
        double up;
        double down;
        DirectionComponents() : left(0.0), right(0.0), up(0.0), down(0.0) {}
    };
    
    struct SearchWindow {
        DBU x;
        DBU y;
        DBU width;
        DBU height;
        SearchWindow() : x(0), y(0), width(0), height(0) {}
        SearchWindow(DBU x_val, DBU y_val, DBU w, DBU h)
            : x(x_val), y(y_val), width(w), height(h) {}
    };
    
    // 私有辅助方法
    void computeTimingVectors(
        Node3D* cell_j,
        Node3D* connected_cell,
        bool is_fan_in,
        double late_lm,
        double early_lm,
        TimingVector& late_vector,
        TimingVector& early_vector);
    
    Pointdp3D getCellCenter(Node3D* cell) const;
    
    bool computeDirectionComponents(
        Node3D* cell,
        DirectionComponents& late_components,
        DirectionComponents& early_components);
    
    void decomposeToComponents(
        const TimingVector& vector,
        DirectionComponents& components);
    
    SearchWindow determineSearchWindow(
        Node3D* cell,
        DBU window_width,
        DBU window_height);
    
    std::vector<Pointdp3D> generateCandidatePositions(
        const SearchWindow& window,
        int num_candidates,
        DBU step_x = 0,
        DBU step_y = 0);
    
    SearchWindow clampSearchWindow(
        const SearchWindow& window,
        Node3D* cell);
    
    // 候选位置生成相关方法
    std::vector<Pointdp3D> generateCandidatePositions(
        Node3D* cell,
        const SearchWindow& window,
        int num_candidates,
        DBU step_x = 0);

    std::vector<Pointdp3D> generateCandidatePositions2(
        Node3D* cell,
        const SearchWindow& window,
        int num_candidates);
    
    std::vector<std::pair<DBU, DBU>> findAvailableSlotsInRow(
        Architecture3D::Row3D* row,
        Node3D* exclude_cell,
        DBU range_start,
        DBU range_end);
    
    std::vector<std::pair<DBU, DBU>> findAvailableSlotsInRowAligned(
        Architecture3D::Row3D* row,
        Node3D* exclude_cell,
        DBU range_start,
        DBU range_end,
        DBU site_width = 0,
        DBU site_spacing = 0);
    
    DBU alignToSite(DBU coord, DBU site_spacing, bool round_up) const;
};

} // namespace tdp