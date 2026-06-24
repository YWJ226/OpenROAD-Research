// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "HBTDatabaseInterface.h"
#include <set>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "sta/Network.hh"

// 包含必要的3D-TDP头文件以获取完整类型定义
#include "tdp/Timingdp.h"
#include "infrastructure/Coordinates.h"
#include "infrastructure/Network3D.h"
#include "infrastructure/Net3D.h"
#include "infrastructure/Node3D.h"

namespace hbta {

// ============================================================================
// HBTDatabaseInterface 实现
// ============================================================================

HBTDatabaseInterface::HBTDatabaseInterface() {
    // 初始化统计信息
    extraction_stats_ = ExtractionStats();
}

bool HBTDatabaseInterface::initialize(tdp::Timingdp* timingdp,
                                     const std::string& hbt_array_json, 
                                     int dbu_per_micron) {
    if (!timingdp) {
        std::cerr << "ERROR: Timingdp pointer is null" << std::endl;
        return false;
    }
    
    timingdp_ = timingdp;
    hbt_array_json_path_ = hbt_array_json;
    dbu_per_micron_ = dbu_per_micron;
    
    // 初始化HBT分配器
    if (!allocator_.initialize(hbt_array_json, dbu_per_micron)) {
        std::cerr << "ERROR: Failed to initialize HBT allocator" << std::endl;
        return false;
    }
    
    // 重置统计信息
    extraction_stats_ = ExtractionStats();
    cell_layer_cache_.clear();
    
    std::cout << "HBTDatabaseInterface initialized successfully" << std::endl;
    std::cout << "  - HBT array JSON: " << hbt_array_json << std::endl;
    std::cout << "  - DBU per micron: " << dbu_per_micron << std::endl;
    
    return true;
}

std::vector<InterLayerConnection> HBTDatabaseInterface::extractInterLayerConnections() {
    std::vector<InterLayerConnection> connections;
    
    if (!timingdp_) {
        std::cerr << "ERROR: Timingdp not initialized" << std::endl;
        return connections;
    }
    
    // 获取3D网络
    auto* network_3d = timingdp_->getNetwork3D();
    if (!network_3d) {
        std::cerr << "ERROR: Network3D not available" << std::endl;
        return connections;
    }
    
    // 重置统计信息
    extraction_stats_ = ExtractionStats();
    
    // 获取数据库中的所有网络
    auto& db_nets = timingdp_->getNets();
    extraction_stats_.total_nets = db_nets.size();
    
    std::cout << "Extracting inter-layer connections from " << db_nets.size() << " nets..." << std::endl;
    
    // 遍历所有数据库网络
    for (auto* db_net : db_nets) {
        if (!db_net) continue;
        
        // 跳过电源网络
        if (db_net->getSigType().isSupply()) {
            continue;
        }
        // 获取对应的Net3D对象
        auto* net_3d = network_3d->getNet(db_net);
        if (!net_3d) {
            // 如果没有对应的Net3D，跳过
            continue;
        }
        
        // 检查是否是跨层网络
        if (isInterLayerNetwork(net_3d)) {
            extraction_stats_.inter_layer_nets++;
            
            // 转换为InterLayerConnection
            auto connection = convertNetToConnection(net_3d);
            if (!connection.net_name.empty()) {
                connections.push_back(connection);
            } else {
                extraction_stats_.failed_conversions++;
            }
        }
    }
    
    std::cout << "Found " << connections.size() << " inter-layer connections" << std::endl;
    return connections;
}

InterLayerConnection HBTDatabaseInterface::convertNetToConnection(tdp::Net3D* net) {
    InterLayerConnection connection;
    
    if (!net || !net->db_net) {
        return connection;
    }
    
    connection.net_name = net->name;
    
    // 通过 STA 获取 driver pin
    auto* sta = timingdp_->getSta();
    auto* db_network = sta->getDbNetwork();
    auto* network = sta->network();
    
    sta::Net* sta_net = db_network->dbToSta(net->db_net);
    sta::PinSet* drivers = network->drivers(sta_net);
    if (!drivers || drivers->empty()) {
        return connection;
    }
    const sta::Pin* drvr_pin = *drivers->begin();
    
    // STA driver pin → db pin
    odb::dbITerm* drvr_iterm = nullptr;
    odb::dbBTerm* drvr_bterm = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    db_network->staToDb(drvr_pin, drvr_iterm, drvr_bterm, moditerm);
    
    // 通过 Network3D 映射 O(1) 查找 driver cell
    auto* network_3d = timingdp_->getNetwork3D();
    tdp::Node3D* driver_cell = nullptr;
    if (drvr_iterm) {
        driver_cell = network_3d->getCell(drvr_iterm->getInst());
    } else if (drvr_bterm) {
        driver_cell = network_3d->getNode(drvr_bterm);
    }
    
    if (!driver_cell) {
        return connection;
    }
    
    // driver 信息（使用真实 layer）
    connection.driver_cell_name = driver_cell->name;
    connection.driver_location = convertToHBTAPoint(driver_cell->location);
    connection.driver_layer = driver_cell->layer;
    
    // 其余 cell 均为 load
    for (auto* cell : net->connected_cells) {
        if (!cell || cell == driver_cell) continue;
        std::string layer = getCellLayer(cell);
        if (layer == "bottom") extraction_stats_.bottom_layer_cells++;
        else if (layer == "upper") extraction_stats_.upper_layer_cells++;
        connection.load_cell_names.push_back(cell->name);
        connection.load_locations.push_back(convertToHBTAPoint(cell->location));
        connection.load_layers.push_back(layer);
    }
    
    // 计算建议的HBT位置
    connection.suggested_hbt_location = calculateSuggestedHBTLocation(net);
    
    return connection;
}

Point HBTDatabaseInterface::calculateSuggestedHBTLocation(tdp::Net3D* net) {
    
    return Point(0.0, 0.0);
}

bool HBTDatabaseInterface::isInterLayerNetwork(tdp::Net3D* net) const {
    if (!net || net->connected_cells.empty()) {
        return false;
    }
    // 只要出现两种及以上层标签，就视为跨层网络。
    std::set<std::string> layers;
    for (auto* cell : net->connected_cells) {
        if (!cell) continue;
        layers.insert(getCellLayer(cell));
    }
    
    return layers.size() > 1;
}

std::string HBTDatabaseInterface::getCellLayer(tdp::Node3D* cell) const {
    if (!cell) {
        return "unknown";
    }
    
    // 先查缓存，避免在大网表中重复字符串处理。
    auto it = cell_layer_cache_.find(cell->name);
    if (it != cell_layer_cache_.end()) {
        return it->second;
    }
    
    // 直接从 Node3D 读取层标签。
    std::string layer = cell->layer;
    
    // 标准化层名称
    if (layer.empty()) {
        layer = "unknown";
    }
    
    // 更新缓存（使用const_cast，因为缓存应该是mutable的）
    // 在实际实现中，应该将cell_layer_cache_声明为mutable
    cell_layer_cache_[cell->name] = layer;
    
    return layer;
}

Point HBTDatabaseInterface::convertToHBTAPoint(const tdp::Pointdp3D& point) const {
    // 旧实现：将DBU坐标转换为微米坐标
    // double x_um = static_cast<double>(point.x.v) / dbu_per_micron_;
    // double y_um = static_cast<double>(point.y.v) / dbu_per_micron_;
    // return Point(x_um, y_um);

    // 新实现：统一使用 DBU，避免 HBTArray(DBU) 与网络坐标(um) 混用。
    double x_dbu = static_cast<double>(point.x.v);
    double y_dbu = static_cast<double>(point.y.v);
    
    return Point(x_dbu, y_dbu);
}

tdp::Point HBTDatabaseInterface::convertToTDPPoint(const Point& point) const {
    return tdp::Point(tdp::DbuX(convertToTDPDBU(point.x)),
                      tdp::DbuY(convertToTDPDBU(point.y)));
}

int64_t HBTDatabaseInterface::convertToTDPDBU(double value) const {
    return static_cast<int64_t>(std::llround(value));
}

void HBTDatabaseInterface::syncHBTAssignmentsToNetwork3D(
    const std::vector<InterLayerConnection>& connections) {
    if (timingdp_ == nullptr) {
        return;
    }

    auto* network_3d = timingdp_->getNetwork3D();
    if (network_3d == nullptr) {
        return;
    }

    for (const auto& connection : connections) {
        if (auto* net_3d = network_3d->getNet(connection.net_name)) {
            net_3d->hbt_connections.clear();
        }
    }

    const tdp::DBU pitch = static_cast<tdp::DBU>(
        convertToTDPDBU(allocator_.getHBTArray().getPitchDBU()));

    for (const auto& connection : connections) {
        const auto* hbt = allocator_.getHBTForNet(connection.net_name);
        if (hbt == nullptr) {
            continue;
        }

        auto* net_3d = network_3d->getNet(connection.net_name);
        if (net_3d == nullptr) {
            continue;
        }

        tdp::Net3D::HBTConnection hbt_conn;
        hbt_conn.location = convertToTDPPoint(hbt->location);
        hbt_conn.pitch = pitch;
        hbt_conn.aligned = false;

        auto set_endpoint = [&](const std::string& cell_name,
                                const std::string& layer) {
            auto* cell = network_3d->getCell(cell_name);
            if (cell == nullptr) {
                return;
            }
            if (layer == "bottom" && hbt_conn.bottom_cell == nullptr) {
                hbt_conn.bottom_cell = cell;
            } else if (layer == "upper" && hbt_conn.upper_cell == nullptr) {
                hbt_conn.upper_cell = cell;
            }
        };

        set_endpoint(connection.driver_cell_name, connection.driver_layer);
        for (size_t i = 0; i < connection.load_cell_names.size(); ++i) {
            const std::string layer =
                i < connection.load_layers.size()
                    ? connection.load_layers[i]
                    : std::string();
            set_endpoint(connection.load_cell_names[i], layer);
        }

        net_3d->addHBTConnection(hbt_conn);
        net_3d->updateHBTAlignment();
    }
}

void HBTDatabaseInterface::reset() {
    timingdp_ = nullptr;
    allocator_.resetAllAllocations();
    extraction_stats_ = ExtractionStats();
    cell_layer_cache_.clear();
}

// ============================================================================
// ExtractionStats 实现
// ============================================================================

void HBTDatabaseInterface::ExtractionStats::print() const {
    std::cout << "=== HBT Network Extraction Statistics ===" << std::endl;
    std::cout << "  Total nets analyzed: " << total_nets << std::endl;
    std::cout << "  Inter-layer nets: " << inter_layer_nets << std::endl;
    std::cout << "  Bottom layer cells: " << bottom_layer_cells << std::endl;
    std::cout << "  Upper layer cells: " << upper_layer_cells << std::endl;
    std::cout << "  Failed conversions: " << failed_conversions << std::endl;
    
    if (total_nets > 0) {
        double inter_layer_ratio = static_cast<double>(inter_layer_nets) / total_nets * 100.0;
        std::cout << "  Inter-layer ratio: " << inter_layer_ratio << "%" << std::endl;
    }
}

} // namespace hbta
