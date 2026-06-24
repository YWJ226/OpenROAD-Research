// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "eva.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <set>

#include "db_sta/SpefWriter.hh"

#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "stt/SteinerTreeBuilder.h"
#include "utl/Logger.h"
#include "sta/Network.hh"
#include "sta/MinMax.hh"
#include "sta/Sdc.hh"
#include "est/SteinerTree.h" 
#include "grt/GlobalRouter.h"
#include "grt/GRoute.h"
#include "grt/RoutePt.h"
#include "sta/Parasitics.hh"
#include "sta/ParasiticsClass.hh"
#include "sta/ArcDelayCalc.hh"
#include "sta/Units.hh"



#include "tdp/Timingdp.h"
#include "../infrastructure/Network3D.h"

using utl::TDP;
using est::SteinerTree;
using est::PinLoc;
using sta::NetConnectedPinIterator;
using sta::PinSet;


namespace tdp {

// ============================================================================
// SteinerTreeGenerator 实现
// ============================================================================

SteinerTreeGenerator::SteinerTreeGenerator(Timingdp* timingdp)
    : tdp_db_(timingdp),
      stt_builder_(std::make_unique<stt::SteinerTreeBuilder>(tdp_db_->getDatabase(), logger_)),
      logger_(timingdp ? timingdp->getLogger() : nullptr),
      sta_(timingdp->getSta()),
      global_router_(timingdp->getGlobalRouter()),
      network_(timingdp->getSta()->network()),
      db_network_(timingdp->getSta()->getDbNetwork()),
      block_(timingdp->getBlock()),
      tech_(timingdp->getDatabase()->getTech()),
      min_max_(sta::MinMax::max()),
      est_(timingdp->getEstimateParasitics())
{
    // 初始化STT构建器
    if (tdp_db_ && tdp_db_->getDatabase()) {
        stt_builder_ = std::make_unique<stt::SteinerTreeBuilder>(
            tdp_db_->getDatabase(), logger_);
    } else {
        logger_->error(utl::TDP, 100, "Failed to initialize SteinerTreeGenerator: invalid database");
    }
    
    // 初始化统计信息
    stats_ = SteinerTreeStats{};
}

SteinerTreeGenerator::~SteinerTreeGenerator() = default;

// ============================================================================
// 主要接口：为Net3D生成斯坦纳树
// ============================================================================
void SteinerTreeGenerator::createSteninerRoute() {
    initSceneParasitics();

    sta::LibertyLibrary* default_lib = db_network_->defaultLibertyLibrary();
    db_network_->sta::Network::clear();
    db_network_->setDefaultLibertyLibrary(default_lib);

    auto& nets = tdp_db_->getNets();
    generateSteinerTrees(nets);

}

void SteinerTreeGenerator::generateSteinerTree(Net3D* net_3d) {
    if (!net_3d) {
        logger_->error(utl::TDP, 101, "Cannot generate Steiner tree for null net");
        return;
    }
    // 确保成员变量已初始化
    if (!network_) {
        network_ = tdp_db_->getSta()->network();
    }
    if (!db_network_) {
        db_network_ = tdp_db_->getSta()->getDbNetwork();
    }

    // 更新统计信息
    stats_.num_nets_processed++;
    

    // 1. 收集引脚位置
    auto dbnet = net_3d->db_net;
    if (!dbnet) {
    logger_->error(TDP, 1201, 
                  "Net3D {} does not have associated dbNet", 
                  net_3d->name);
    return;
    }
    sta::Net* cur_net = db_network_->dbToSta(dbnet);

    if (dbnet->getITerms().size() + dbnet->getBTerms().size() < 2) {
        return;
    }

    if(db_network_->staToDb(cur_net)->isSpecial()){
        return;
    }

    if (isPadNet(cur_net)) {
        makePadParasitic(cur_net);
        return;
    }
    
    PinCoordinates pin_locs = estimateWireSteiner(cur_net ,net_3d);
    int drvr_index = pin_locs.drvr_idx;
    std::vector<odb::Point> all_locs;
    all_locs.reserve(pin_locs.pinCount());
    for (size_t i = 0; i < pin_locs.pinCount(); i++) {
        all_locs.emplace_back(pin_locs.x_coords[i], pin_locs.y_coords[i]);
    }

    //判断是否跨层
    grt::GRoute result_route;
    if (isInterLayerNetwork(net_3d)){
        std::vector<odb::Point> bottom_pins;
        std::vector<odb::Point> upper_pins;
        std::vector<int> bottom_indices;
        std::vector<int> upper_indices;

        int drvr_index = pin_locs.drvr_idx;
        int drvr_layer = -1;
    
        for (size_t i = 0; i < pin_locs.pinCount(); i++) {
            int layer = pin_locs.layer_coords[i];
            odb::Point point(pin_locs.x_coords[i], pin_locs.y_coords[i]);
            
            if (layer == 0) {  // bottom层
                bottom_pins.push_back(point);
                bottom_indices.push_back(i);
                if (i == drvr_index) drvr_layer = 0;
            } else if (layer == 1) {  // upper层
                upper_pins.push_back(point);
                upper_indices.push_back(i);
                if (i == drvr_index) drvr_layer = 1;
            }
        }

        // 3. 添加HBT位置
        std::vector<odb::Point> hbt_locations;
        for (const auto& hbt_conn : net_3d->hbt_connections) {
            hbt_locations.push_back(odb::Point(hbt_conn.location.x.v, hbt_conn.location.y.v));
        }
        
        // 4. 为每层生成斯坦纳树
        stt::Tree bottom_tree, upper_tree;
        
        // 生成bottom层斯坦纳树
        if (!bottom_pins.empty()) {
            std::vector<odb::Point> bottom_points = bottom_pins;
            bottom_points.insert(bottom_points.end(), hbt_locations.begin(), hbt_locations.end());
            
            std::vector<int> bottom_x, bottom_y;
            for (const auto& point : bottom_points) {
                bottom_x.push_back(point.x());
                bottom_y.push_back(point.y());
            }
            
            // 确定驱动引脚索引：如果驱动在bottom层，使用原始索引；否则使用HBT索引
            int bottom_drvr_idx = (drvr_layer == 0) ? 
                std::distance(bottom_indices.begin(), 
                            std::find(bottom_indices.begin(), bottom_indices.end(), drvr_index)) : 
                bottom_pins.size();  // HBT作为驱动
                
            bottom_tree = callFluteAlgorithm(bottom_x, bottom_y, bottom_drvr_idx);
        }
        
        // 生成upper层斯坦纳树
        if (!upper_pins.empty()) {
            std::vector<odb::Point> upper_points = upper_pins;
            upper_points.insert(upper_points.end(), hbt_locations.begin(), hbt_locations.end());
            
            std::vector<int> upper_x, upper_y;
            for (const auto& point : upper_points) {
                upper_x.push_back(point.x());
                upper_y.push_back(point.y());
            }
            
            // 确定驱动引脚索引：如果驱动在upper层，使用原始索引；否则使用HBT索引
            int upper_drvr_idx = (drvr_layer == 1) ? 
                std::distance(upper_indices.begin(), 
                            std::find(upper_indices.begin(), upper_indices.end(), drvr_index)) : 
                upper_pins.size();  // HBT作为驱动
                
            upper_tree = callFluteAlgorithm(upper_x, upper_y, upper_drvr_idx);
        }
        
        //return {bottom_tree, upper_tree};
        int bottom_layer = 9;
        int upper_layer = 18;
        grt::GRoute bottom_route = convertSteinerTreeToGRoute(bottom_tree, bottom_layer);
        grt::GRoute upper_route = convertSteinerTreeToGRoute(upper_tree, upper_layer);
        grt::GRoute combined_route;
        
        // 合并bottom和upper层的route
        combined_route.insert(combined_route.end(), bottom_route.begin(), bottom_route.end());
        combined_route.insert(combined_route.end(), upper_route.begin(), upper_route.end());
        
        // 添加HBT作为via连接
        for (const auto& hbt_conn : net_3d->hbt_connections) {
            odb::Point hbt_loc(hbt_conn.location.x.v, hbt_conn.location.y.v);
            grt::GSegment hbt_via;
            hbt_via.init_x = hbt_loc.x();
            hbt_via.init_y = hbt_loc.y();
            hbt_via.init_layer = bottom_layer;
            hbt_via.final_x = hbt_loc.x();
            hbt_via.final_y = hbt_loc.y();
            hbt_via.final_layer = upper_layer;
            hbt_via.is_jumper = false;
            combined_route.push_back(hbt_via);
        }
        if (global_router_ && !combined_route.empty() && dbnet) {
            for (size_t i = 0; i < pin_locs.pinCount(); i++) {
                int conn = getPinConnLayer(pin_locs.pins[i]);
                int rt_layer = pin_locs.layer_coords[i] == 0
                    ? bottom_layer : upper_layer;
                if (conn > 0 && rt_layer > 0 && conn != rt_layer) {
                    grt::GSegment via;
                    via.init_x = all_locs[i].x();
                    via.init_y = all_locs[i].y();
                    via.init_layer = conn;
                    via.final_x = all_locs[i].x();
                    via.final_y = all_locs[i].y();
                    via.final_layer = rt_layer;
                    via.is_jumper = false;
                    combined_route.push_back(via);
                }
            }
        }
        result_route = std::move(combined_route);

    }else{
        stats_.total_pins += all_locs.size();

        // 3. 提取X、Y坐标
        std::vector<int> x_coords, y_coords;
        x_coords.reserve(all_locs.size());
        y_coords.reserve(all_locs.size());
        
        for (const auto& point : all_locs) {
            x_coords.push_back(point.x());
            y_coords.push_back(point.y());
        }
        
        // 4. 调用FLUTE算法生成斯坦纳树
        int layer = -1;
        stt::Tree tree = callFluteAlgorithm(x_coords, y_coords, drvr_index);
        if (pin_locs.layer_coords[0] == 0){
            layer = 9;
        }else {
            layer = 18;
        }
        grt::GRoute route = convertSteinerTreeToGRoute(tree, layer);
        
        // 6. 更新统计信息
        stats_.total_steiner_points += (tree.branchCount() - all_locs.size());
        stats_.total_wirelength += tree.length;

        if (global_router_ && !route.empty() && dbnet) {
            for (size_t i = 0; i < pin_locs.pinCount(); i++) {
                int conn = getPinConnLayer(pin_locs.pins[i]);
                if (conn > 0 && conn != layer) {
                    grt::GSegment via;
                    via.init_x = all_locs[i].x();
                    via.init_y = all_locs[i].y();
                    via.init_layer = conn;
                    via.final_x = all_locs[i].x();
                    via.final_y = all_locs[i].y();
                    via.final_layer = layer;
                    via.is_jumper = false;
                    route.push_back(via);
                }
            }
        }

        result_route = std::move(route);

    }
    
    if (global_router_ && !result_route.empty() && dbnet) {
        estimateTdpWireParasitic(dbnet, result_route, pin_locs, all_locs);
        global_router_->setRoute(dbnet, result_route);
    }
}

void SteinerTreeGenerator::updateSteinerTrees(Node3D* cell){

    if (!cell || !global_router_) return;

    // sta::LibertyLibrary* default_lib = db_network_->defaultLibertyLibrary();
    // db_network_->sta::Network::clear();
    // db_network_->setDefaultLibertyLibrary(default_lib);

    Network3D* network_3d = tdp_db_->getNetwork3D();
    for (auto* db_net : cell->connected_nets) {
        if (!db_net) continue;
        bool all_modes_ideal_clock = true;
        sta::Net* sta_net = db_network_->dbToSta(db_net);
        PinSet* drivers = network_->drivers(sta_net);
        if (drivers && !drivers->empty()) {
            const sta::Pin* drvr_pin = *drivers->begin();
            for (sta::Mode* mode : sta_->modes()) {
                if (!sta_->isIdealClock(drvr_pin, mode)) {
                all_modes_ideal_clock = false;
                break;
                }
            }
        }
        if (all_modes_ideal_clock) {
            continue;
        }

        Net3D* net3d = network_3d->getNet(db_net);
        if (net3d) {
            generateSteinerTree(net3d);
            auto& routes = global_router_->getRoutes();
            auto it = routes.find(db_net);
            if (it != routes.end() && !it->second.empty()) {
                sta_->delaysInvalidFromFanin(sta_net);
            }
        }
    }
}

void SteinerTreeGenerator::generateSteinerTrees(std::vector<odb::dbNet*>& nets) {
    
    Network3D* network_3d = tdp_db_ -> getNetwork3D();
    // std::vector<stt::Tree> trees;
    // trees.reserve(nets.size());
     // 获取GlobalRouter
    // bool has_global_router = (global_router_ != nullptr);
    
    for (auto* db_net : nets){

        bool all_modes_ideal_clock = true;
        sta::Net* sta_net = db_network_->dbToSta(db_net);
        PinSet* drivers = network_->drivers(sta_net);
        if (drivers && !drivers->empty()) {
            const sta::Pin* drvr_pin = *drivers->begin();
            for (sta::Mode* mode : sta_->modes()) {
                if (!sta_->isIdealClock(drvr_pin, mode)) {
                all_modes_ideal_clock = false;
                break;
                }
            }
        }
        if (all_modes_ideal_clock) {
            continue;
        }

        Net3D* net = network_3d->getNet(db_net);
        if (net) {
            generateSteinerTree(net);

        }
    }
}

// ============================================================================
// 辅助方法实现
// ============================================================================
std::vector<odb::Point> SteinerTreeGenerator::addHBTLocations(
    Net3D* net_3d, const std::vector<odb::Point>& pin_locs) {
    
    std::vector<odb::Point> all_points = pin_locs;
    
    if (!net_3d) {
        return all_points;
    }
    
    // 添加HBT位置
    for (const auto& hbt_conn : net_3d->hbt_connections) {
        all_points.push_back(odb::Point(hbt_conn.location.x.v, hbt_conn.location.y.v));
    }
////////debug
    // 先复制一份并排序
    std::vector<odb::Point> sorted_points = all_points;
    std::sort(sorted_points.begin(), sorted_points.end(),
        [](const odb::Point& a, const odb::Point& b) {
            return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
        });

    // 检查排序后的向量中是否有相邻重复
    bool has_duplicates = std::adjacent_find(sorted_points.begin(), sorted_points.end(),
        [](const odb::Point& a, const odb::Point& b) {
            return a.x() == b.x() && a.y() == b.y();
        }) != sorted_points.end();

    if (has_duplicates) {
        logger_->warn(TDP, 113, "Found duplicate points in net (after sorting)");
    }
    
    return all_points;
}

stt::Tree SteinerTreeGenerator::callFluteAlgorithm(
    const std::vector<int>& x_coords,
    const std::vector<int>& y_coords,
    int                     drvr_index) {
    
    if (x_coords.size() != y_coords.size() || x_coords.size() < 2) {
        logger_->warn(TDP, 104, 
                      "Invalid coordinates for FLUTE algorithm: {} points",
                      x_coords.size());
        return stt::Tree{};
    }
    
    try {
        // 调用STT的FLUTE算法，这里使用默认的alpha=0.0，表示纯FLUTE算法
        return stt_builder_->makeSteinerTree(x_coords, y_coords, drvr_index);
        
    } catch (const std::exception& e) {
        logger_->error(TDP, 105, 
                      "FLUTE algorithm failed: {}", e.what());
        return stt::Tree{};
    }
}


void SteinerTreeGenerator::initSceneParasitics()
{
    sta_->deleteParasitics();
    for (sta::Scene* scene : sta_->scenes()) {
        sta::Parasitics* parasitics
            = sta_->makeConcreteParasitics(scene->name(), "");
        scene->setParasitics(parasitics, sta::MinMaxAll::minMax());
    }
}

bool SteinerTreeGenerator::isInterLayerNetwork(Net3D* net) const {
    if (!net) {
        return false;
    }
    
    // 检查网络连接的引脚是否来自不同层
    const auto& pins = net->getPins();
    if (pins.empty()) {
        return false;
    }
    
    std::string first_layer = "";
    for (const auto* pin : pins) {
        if (pin && pin->getNode()) {
            std::string pin_layer = pin->getNode()->layer;
            if (first_layer.empty()) {
                first_layer = pin_layer;
            } else if (pin_layer != first_layer) {
                return true;
            }
        }
    }
    
    return false;
}

PinCoordinates SteinerTreeGenerator::estimateWireSteiner(const sta::Net* net, const Net3D* net_3d)
{
    PinCoordinates result;
  PinSet* drivers = network_->drivers(net);
  if (drivers && !drivers->empty()) {
    const sta::Pin* drvr_pin = *drivers->begin();
   result = estimateWireSteiner(drvr_pin, net, net_3d);
  }
  return result;
}

PinCoordinates SteinerTreeGenerator::estimateWireSteiner(const sta::Pin* drvr_pin,
                                               const sta::Net* net, const Net3D* net_3d)
{
    PinCoordinates result;
  if (!network_->isPower(net) && !network_->isGround(net)) {
      result = collectPinLocations(drvr_pin, net_3d);
    }
    return result;
}

static void connectedPins(const sta::Net* net,
                          sta::Network* network,
                          sta::dbNetwork* db_network,
                          // Return value.
                          std::vector<PinLoc>& pins);

static void connectedPins(const sta::Net* net,
                          sta::Network* network,
                          sta::dbNetwork* db_network,
                          // Return value.
                          std::vector<PinLoc>& pins)
{
  NetConnectedPinIterator* pin_iter = network->connectedPinIterator(net);
  while (pin_iter->hasNext()) {
    const sta::Pin* pin = pin_iter->next();
    odb::dbITerm* iterm;
    odb::dbBTerm* bterm;
    odb::dbModITerm* moditerm;
    db_network->staToDb(pin, iterm, bterm, moditerm);
    //
    // only accumuate the flat pins (in hierarchical mode we may
    // hit moditerms/modbterms).
    //
    if (iterm || bterm) {
      odb::Point loc = db_network->location(pin);
      pins.push_back({pin, loc});
    }
  }
  delete pin_iter;
}

PinCoordinates SteinerTreeGenerator::collectPinLocations(const sta::Pin* drvr_pin, const Net3D* net_3d)
{
    PinCoordinates result;

    odb::dbNet* db_net = db_network_->findFlatDbNet(drvr_pin);
    
    if (!db_net) {
        logger_->error(TDP, 107, "Cannot find dbNet for driver pin");
        return result; 
    }

    sta::Net* net = db_network_->dbToSta(db_net);

    SteinerTree* tree = new SteinerTree(drvr_pin, db_network_, logger_);
    std::vector<PinLoc>& pinlocs = tree->pinlocs();
    // Find all the connected pins
    connectedPins(net, network_, db_network_, pinlocs);
    // Sort pins by location because connectedPins order is not deterministic.
    std::stable_sort(pinlocs.begin(), pinlocs.end(), [=](const PinLoc& pin1, const PinLoc& pin2) {
        return pin1.loc.getX() < pin2.loc.getX()
                || (pin1.loc.getX() == pin2.loc.getX()
                    && pin1.loc.getY() < pin2.loc.getY());
    });
    int pin_count = pinlocs.size();
    if (pin_count < 2) {
        logger_->warn(TDP, 108, "Net has only {} pin(s)", pin_count);
        delete tree; // 记得释放内存
        return result; // 返回空结果
    }
    
    // 准备坐标向量
    result.x_coords.reserve(pin_count);
    result.y_coords.reserve(pin_count);
    result.layer_coords.reserve(pin_count);
    result.pins.reserve(pin_count);
    
    bool is_placed = true;
    result.drvr_idx = -1;
    
    for (int i = 0; i < pin_count; i++) {
        const PinLoc& pinloc = pinlocs[i];
        
        // 记录驱动引脚索引
        if (pinloc.pin == drvr_pin) {
            result.drvr_idx = i;
        }
        
        result.layer_coords.push_back(findPinLayer(pinloc, net_3d));
        // // 收集坐标，snap 到路由网格
        // odb::Point grid_loc = global_router_
        //                           ? global_router_->getPositionOnGrid(pinloc.loc)
        //                           : pinloc.loc;
        // 收集坐标，snap 到路由网格
        odb::Point grid_loc = pinloc.loc;
        result.x_coords.push_back(grid_loc.x());
        result.y_coords.push_back(grid_loc.y());
        result.pins.push_back(pinloc.pin);
        
        // 检查引脚是否已放置
        is_placed &= db_network_->isPlaced(pinloc.pin);
        
        // 如果需要，添加到树的映射中
        tree->locAddPin(pinloc.loc, pinloc.pin);
    }
    
    // 清理内存
    delete tree;
    
    if (!is_placed) {
        logger_->warn(utl::TDP, 109, "Some pins in net are not placed");
    }
    
    if (result.drvr_idx == -1) {
        logger_->error(utl::TDP, 110, "Driver pin not found in pin list");
    }
    
    return result;
}

////////////////////
// 引脚分层
////////////////////
int SteinerTreeGenerator::findPinLayer(const PinLoc& pinloc, const Net3D* net_3d){

    const auto& pins = net_3d->getPins();
    for (const auto* pin : pins) {
        if (pin) {

            Pin3Dlocation coord = getPin3DLocation(pin);
            int layer = -2;  
            if (pin->getNode()) {
                std::string node_str = pin->getNode()->layer;
                if (node_str == "bottom") {
                    layer = 0;
                } else if (node_str == "upper") {
                    layer = 1;
                }
            }
            if ( layer == -2 ){
                logger_ -> error(TDP , 1168, "layer error");
            }

            if (coord.x == pinloc.loc.x() && coord.y == pinloc.loc.y()) {
                return layer;
            }
        }
    }
    return -2;
}

int SteinerTreeGenerator::getPinConnLayer(const sta::Pin* sta_pin) const
{
    odb::dbITerm* iterm = nullptr;
    odb::dbBTerm* bterm = nullptr;
    odb::dbModITerm* moditerm = nullptr;
    db_network_->staToDb(sta_pin, iterm, bterm, moditerm);

    int highest = -1;

    if (iterm) {
        odb::dbMTerm* mterm = iterm->getMTerm();
        for (odb::dbMPin* mpin : mterm->getMPins()) {
            for (odb::dbBox* box : mpin->getGeometry()) {
                odb::dbTechLayer* tech_layer = box->getTechLayer();
                if (tech_layer
                    && tech_layer->getType() == odb::dbTechLayerType::ROUTING) {
                    highest = std::max(
                        highest, tech_layer->getRoutingLevel());
                }
            }
        }
    }

    if (bterm) {
        for (odb::dbBPin* bpin : bterm->getBPins()) {
            for (odb::dbBox* box : bpin->getBoxes()) {
                odb::dbTechLayer* tech_layer = box->getTechLayer();
                if (tech_layer
                    && tech_layer->getType() == odb::dbTechLayerType::ROUTING) {
                    highest = std::max(
                        highest, tech_layer->getRoutingLevel());
                }
            }
        }
    }

    return highest;
}

Pin3Dlocation SteinerTreeGenerator::getPin3DLocation(const Pin3D* pin) const {
    if (!pin || !tdp_db_) {
        return Pin3Dlocation(0, 0, 0);
    }
    
    // 将Pin3D转换为sta::Pin
    odb::dbBTerm* Bterm = pin->getDbBTerm();
    odb::dbITerm* Iterm = pin->getDbITerm();

    sta::Pin* sta_pin = nullptr;

    if (Bterm) {
        // 边界引脚
        sta_pin = db_network_->dbToSta(Bterm);
    } else if (Iterm) {
        // 内部引脚
        sta_pin = db_network_->dbToSta(Iterm);
    } else {
        return Pin3Dlocation(0, 0, 0);
        logger_->error(TDP, 1901, "Error pin location");
    }
    
    // 获取OpenSTA坐标
    odb::Point sta_loc = db_network_->location(sta_pin);
    
    // 获取Z坐标（从Pin3D或根据层判断）
    int z = pin->getOffsetZ().v;
    if (pin->getNode()) {
        // 根据层设置Z：bottom=0, upper=1000
        z = (pin->getNode()->layer == "upper") ? 1000 : 0;
    }
    
    return Pin3Dlocation(sta_loc.x(), sta_loc.y(), z);
}

////////////////////
//转化布线
///////////////////
grt::GRoute SteinerTreeGenerator::convertSteinerTreeToGRoute(
    const stt::Tree& tree, 
    int routing_layer) {
    
    grt::GRoute route;
    
    if (tree.branch.empty()) {
        logger_->error(TDP, 120, "Empty Steiner tree");
        return route;
    }
    
    // 遍历所有branch
    for (int i = 0; i < tree.branchCount(); i++) {
        const auto& branch = tree.branch[i];
        
        // 获取起点坐标
        int start_x = branch.x;
        int start_y = branch.y;
        
        // 获取终点坐标（通过邻居索引）
        int end_idx = branch.n;
        if (end_idx == i) continue;
        if (end_idx < 0 || end_idx >= tree.branchCount()) {
            logger_->warn(TDP, 121, "Invalid neighbor index in branch {}: {}", i, end_idx);
            continue;
        }
        
        const auto& end_branch = tree.branch[end_idx];
        int end_x = end_branch.x;
        int end_y = end_branch.y;
        
        // 创建GSegment（直接字段赋值，绕过构造函数的 min/max bug）
        grt::GSegment segment;
        segment.init_x = start_x;
        segment.init_y = start_y;
        segment.init_layer = routing_layer;
        segment.final_x = end_x;
        segment.final_y = end_y;
        segment.final_layer = routing_layer;
        segment.is_jumper = false;

        
        route.push_back(segment);
    }
    
    // logger_->info(TDP, 122, "Converted Steiner tree to GRoute with {} segments", route.size());
    return route;
}

// ============================================================================
// generateLocalRoute — like generateSteinerTree but without global writes
// ============================================================================

grt::GRoute SteinerTreeGenerator::generateLocalRoute(
    Net3D* net_3d,
    const Pin3D* moved_pin,
    int new_x, int new_y)
{
    grt::GRoute empty_route;

    if (!net_3d) return empty_route;

    if (!network_)  network_  = tdp_db_->getSta()->network();
    if (!db_network_) db_network_ = tdp_db_->getSta()->getDbNetwork();

    auto dbnet = net_3d->db_net;
    if (!dbnet) return empty_route;

    sta::Net* cur_net = db_network_->dbToSta(dbnet);
    if (!cur_net) return empty_route;

    if (dbnet->getITerms().size() + dbnet->getBTerms().size() < 2)
        return empty_route;
    if (db_network_->staToDb(cur_net)->isSpecial())
        return empty_route;
    if (isPadNet(cur_net))
        return empty_route;

    // Collect pin locations (with override for moved_pin)
    PinCoordinates pin_locs = estimateWireSteiner(cur_net, net_3d);
    int drvr_index = pin_locs.drvr_idx;

    // Override moved_pin position
    if (moved_pin) {
        sta::Pin* moved_sta = db_network_->findPin(
            moved_pin->getName().c_str());
        for (size_t i = 0; i < pin_locs.pinCount(); ++i) {
            if (pin_locs.pins[i] == moved_sta) {
                pin_locs.x_coords[i] = new_x;
                pin_locs.y_coords[i] = new_y;
                break;
            }
        }
    }

    std::vector<odb::Point> all_locs;
    all_locs.reserve(pin_locs.pinCount());
    for (size_t i = 0; i < pin_locs.pinCount(); i++)
        all_locs.emplace_back(pin_locs.x_coords[i], pin_locs.y_coords[i]);

    grt::GRoute result_route;

    if (isInterLayerNetwork(net_3d)) {
        // Simplified inter-layer handling — same as generateSteinerTree
        std::vector<odb::Point> bottom_pins, upper_pins;
        std::vector<int> bottom_indices, upper_indices;
        int drvr_layer = -1;

        for (size_t i = 0; i < pin_locs.pinCount(); i++) {
            int layer = pin_locs.layer_coords[i];
            odb::Point pt(pin_locs.x_coords[i], pin_locs.y_coords[i]);
            if (layer == 0) {
                bottom_pins.push_back(pt);
                bottom_indices.push_back(static_cast<int>(i));
                if (static_cast<int>(i) == drvr_index) drvr_layer = 0;
            } else if (layer == 1) {
                upper_pins.push_back(pt);
                upper_indices.push_back(static_cast<int>(i));
                if (static_cast<int>(i) == drvr_index) drvr_layer = 1;
            }
        }

        std::vector<odb::Point> hbt_locations;
        for (const auto& hbt : net_3d->hbt_connections)
            hbt_locations.emplace_back(hbt.location.x.v, hbt.location.y.v);

        stt::Tree bottom_tree, upper_tree;

        if (!bottom_pins.empty()) {
            auto bp = bottom_pins;
            bp.insert(bp.end(), hbt_locations.begin(), hbt_locations.end());
            std::vector<int> bx, by;
            for (auto& p : bp) { bx.push_back(p.x()); by.push_back(p.y()); }
            int bi = (drvr_layer == 0) ? static_cast<int>(std::distance(
                bottom_indices.begin(),
                std::find(bottom_indices.begin(), bottom_indices.end(), drvr_index)))
                : static_cast<int>(bottom_pins.size());
            bottom_tree = callFluteAlgorithm(bx, by, bi);
        }
        if (!upper_pins.empty()) {
            auto up = upper_pins;
            up.insert(up.end(), hbt_locations.begin(), hbt_locations.end());
            std::vector<int> ux, uy;
            for (auto& p : up) { ux.push_back(p.x()); uy.push_back(p.y()); }
            int ui = (drvr_layer == 1) ? static_cast<int>(std::distance(
                upper_indices.begin(),
                std::find(upper_indices.begin(), upper_indices.end(), drvr_index)))
                : static_cast<int>(upper_pins.size());
            upper_tree = callFluteAlgorithm(ux, uy, ui);
        }

        int bl = 9, ul = 18;
        auto br = convertSteinerTreeToGRoute(bottom_tree, bl);
        auto ur = convertSteinerTreeToGRoute(upper_tree, ul);
        result_route.insert(result_route.end(), br.begin(), br.end());
        result_route.insert(result_route.end(), ur.begin(), ur.end());

        for (const auto& hbt : net_3d->hbt_connections) {
            grt::GSegment hbt_via;
            hbt_via.init_x = hbt.location.x.v;
            hbt_via.init_y = hbt.location.y.v;
            hbt_via.init_layer = bl;
            hbt_via.final_x = hbt.location.x.v;
            hbt_via.final_y = hbt.location.y.v;
            hbt_via.final_layer = ul;
            hbt_via.is_jumper = false;
            result_route.push_back(hbt_via);
        }
    } else {
        std::vector<int> x_coords, y_coords;
        for (const auto& pt : all_locs) {
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
        }
        int layer = (pin_locs.layer_coords.empty()
                     || pin_locs.layer_coords[0] == 0) ? 9 : 18;
        stt::Tree tree = callFluteAlgorithm(x_coords, y_coords, drvr_index);
        result_route = convertSteinerTreeToGRoute(tree, layer);
    }

    // Add pin connection layer-crossing segments (same as original)
    int base_layer = (pin_locs.layer_coords.empty()
                      || pin_locs.layer_coords[0] == 0) ? 9 : 18;
    for (size_t i = 0; i < pin_locs.pinCount(); i++) {
        int conn = getPinConnLayer(pin_locs.pins[i]);
        if (conn > 0 && conn != base_layer) {
            grt::GSegment via;
            via.init_x = all_locs[i].x();
            via.init_y = all_locs[i].y();
            via.init_layer = conn;
            via.final_x = all_locs[i].x();
            via.final_y = all_locs[i].y();
            via.final_layer = base_layer;
            via.is_jumper = false;
            result_route.push_back(via);
        }
    }

    return result_route;
}

// ============================================================================
// getSegmentRC — public wrapper around tdpLayerRC
// ============================================================================

void SteinerTreeGenerator::getSegmentRC(
    int wire_length_dbu,
    int init_x, int init_y,
    int final_x, int final_y,
    int layer_id,
    float& res, float& cap) const
{
    const sta::Scene* corner = nullptr;
    if (tdp_db_ && tdp_db_->getSta()
        && !tdp_db_->getSta()->scenes().empty())
        corner = tdp_db_->getSta()->scenes()[0];
    tdpLayerRC(wire_length_dbu, init_x, init_y, final_x, final_y,
               layer_id, corner, res, cap);
}

void SteinerTreeGenerator::estimateTdpWireParasitic(
    odb::dbNet* db_net, grt::GRoute& route,
    const PinCoordinates& pin_locs,
    const std::vector<odb::Point>& all_locs)
{
    sta::Net* sta_net = db_network_->dbToSta(db_net);
    // PinSet* drivers = network_->drivers(sta_net);
    // const sta::Pin* drvr_pin
    //     = (drivers && !drivers->empty()) ? *drivers->begin() : nullptr;

    // // ---- 调试打印 net 整体信息----
    // std::cout << "=== Net: " << db_net->getName()
    //           << "  segments=" << route.size()
    //           << "  pins=" << pin_locs.pinCount()
    //           << std::endl;
    // if (drvr_pin)
    //     std::cout << "  driver: " << network_->name(drvr_pin) << std::endl;
    // for (size_t i = 0; i < pin_locs.pinCount(); i++) {
    //     if (pin_locs.pins[i] != drvr_pin) {
    //         std::cout << "  load: " << network_->name(pin_locs.pins[i])
    //                   << " layer=" << pin_locs.layer_coords[i]
    //                   << " at (" << all_locs[i].x() << "," << all_locs[i].y() << ")"
    //                   << std::endl;
    //     }
    // }        

    for (sta::Scene* corner : sta_->scenes()) {
        std::map<grt::RoutePt, sta::ParasiticNode*> node_map;
        sta::Parasitics* parasitics = corner->parasitics(min_max_);
        sta::Parasitic* parasitic
            = parasitics->makeParasiticNetwork(sta_net, false);

        std::map<grt::RoutePt, std::vector<const sta::Pin*>> pin_at_routept;
        for (size_t i = 0; i < pin_locs.pinCount(); i++) {
            int conn_layer = getPinConnLayer(pin_locs.pins[i]);
            odb::Point grid_pt = all_locs[i];
            pin_at_routept[grt::RoutePt(grid_pt.x(), grid_pt.y(), conn_layer)]
                .push_back(pin_locs.pins[i]);
        }

        size_t resistor_id = 1;
        std::set<const sta::Pin*> connected_pins;

        for (grt::GSegment& segment : route) {

            // // ----segment 调试 ----
            // std::cout << "  seg[" << &segment - route.data() << "] "
            //           << "init(" << segment.init_x << "," << segment.init_y << ",L" << segment.init_layer << ")";
            // {
            //     auto it = pin_at_routept.find(grt::RoutePt(segment.init_x, segment.init_y, segment.init_layer));
            //     if (it != pin_at_routept.end())
            //         for (auto* p : it->second)
            //             std::cout << " [pin:" << network_->name(p) << "]";
            // }
            // std::cout << " -> final(" << segment.final_x << "," << segment.final_y << ",L" << segment.final_layer << ")";
            // {
            //     auto it = pin_at_routept.find(grt::RoutePt(segment.final_x, segment.final_y, segment.final_layer));
            //     if (it != pin_at_routept.end())
            //         for (auto* p : it->second)
            //             std::cout << " [pin:" << network_->name(p) << "]";
            // }
            // std::cout << std::endl;

            sta::ParasiticNode* n1 = ensureTdpParasiticNode(
                parasitics, parasitic, segment.init_x, segment.init_y,
                segment.init_layer, node_map, sta_net);
            sta::ParasiticNode* n2 = ensureTdpParasiticNode(
                parasitics, parasitic, segment.final_x, segment.final_y,
                segment.final_layer, node_map, sta_net);

            if (segment.init_layer != segment.final_layer) {
                if ((segment.init_layer == 18 && segment.final_layer == 9) ||
                    (segment.init_layer == 9 && segment.final_layer == 18)){
                        parasitics->incrCap(n1, 0.0);
                        parasitics->makeResistor(
                            parasitic, resistor_id++, 6.3, n1, n2);
                        parasitics->incrCap(n2, 0.0);
                    }else{
                        parasitics->makeResistor(
                            parasitic, resistor_id++, 7.2, n1, n2);
                }
                        
            } else {
                float res, cap;
                tdpLayerRC(segment.length(),
                           segment.init_x, segment.init_y,
                           segment.final_x, segment.final_y,
                           segment.init_layer, corner, res, cap);
                parasitics->incrCap(n1, cap / 2.0);
                parasitics->makeResistor(
                    parasitic, resistor_id++, res, n1, n2);
                parasitics->incrCap(n2, cap / 2.0);
            }

            auto it1 = pin_at_routept.find(
                grt::RoutePt(segment.init_x, segment.init_y, segment.init_layer));
            if (it1 != pin_at_routept.end()) {
                for (const sta::Pin* pin : it1->second) {
                    sta::ParasiticNode* pin_node
                        = parasitics->ensureParasiticNode(
                            parasitic, pin, network_);
                    if (connected_pins.find(pin) == connected_pins.end()) {
                        parasitics->makeResistor(
                            parasitic, resistor_id++, 0.001, pin_node, n1);
                        connected_pins.insert(pin);
                    }
                }
            }

            auto it2 = pin_at_routept.find(
                grt::RoutePt(segment.final_x, segment.final_y, segment.final_layer));
            if (it2 != pin_at_routept.end()) {
                for (const sta::Pin* pin : it2->second) {
                    sta::ParasiticNode* pin_node
                        = parasitics->ensureParasiticNode(
                            parasitic, pin, network_);
                    if (connected_pins.find(pin) == connected_pins.end()) {
                        parasitics->makeResistor(
                            parasitic, resistor_id++, 0.001, pin_node, n2);
                        connected_pins.insert(pin);
                    }
                }
            }
        }

        sta_->arcDelayCalc()->reduceParasitic(
            parasitic, sta_net, corner, sta::MinMaxAll::all());

        // for (const sta::RiseFall* rf : {sta::RiseFall::rise(), sta::RiseFall::fall()}) {
        //     sta::Parasitic* pi = drvr_pin
        //                             ? parasitics->findPiElmore(drvr_pin, rf, min_max_)
        //                             : nullptr;
        //     float c2 = 0.0f;
        //     float rpi = 0.0f;
        //     float c1 = 0.0f;

        //     if (pi != nullptr && parasitics->isPiModel(pi)) {
        //         parasitics->piModel(pi, c2, rpi, c1);
        //     }

        //     std::cout << "pi=" << (pi ? "found" : "null")
        //     << " rf="
        //     << (rf == sta::RiseFall::rise() ? "rise" : "fall")
        //     << " pi_rpi=" << rpi
        //     << " pi_c1=" << c1
        //     << " pi_c2=" << c2
        //     << std::endl;
        // }

        //         // parasitics->deleteParasiticNetwork(sta_net);
    }
}

sta::ParasiticNode* SteinerTreeGenerator::ensureTdpParasiticNode(
    sta::Parasitics* parasitics, sta::Parasitic* parasitic,
    int x, int y, int layer,
    std::map<grt::RoutePt, sta::ParasiticNode*>& node_map,
    sta::Net* net) const
{
    grt::RoutePt pt(x, y, layer);
    auto& node = node_map[pt];
    if (!node) {
        node = parasitics->ensureParasiticNode(
            parasitic, net, node_map.size(), network_);
    }
    return node;
}

double SteinerTreeGenerator::dbuToMeters(int dbu) const
{
    return (double)dbu / (tech_->getDbUnitsPerMicron() * 1E+6);
}

void SteinerTreeGenerator::tdpLayerRC(int wire_length_dbu,
                                       int init_x, int init_y,
                                       int final_x, int final_y,
                                       int layer_id,
                                       const sta::Scene* corner,
                                       float& res, float& cap) const
{
    double len = dbuToMeters(wire_length_dbu);

    double dx = wire_length_dbu > 0
        ? dbuToMeters(std::abs(final_x - init_x)) / len
        : 0.5;
    double dy = 1.0 - dx;

    double r_per_meter = 0.0;
    double c_per_meter = 0.0;

    if (est_) {
        if(layer_id == 18){
        c_per_meter = dx * est_->wireClkHCapacitance(corner)
                    + dy * est_->wireClkVCapacitance(corner);
        r_per_meter = dx * est_->wireClkHResistance(corner)
                    + dy * est_->wireClkVResistance(corner);
        }else{
        r_per_meter = dx * est_->wireSignalHResistance(corner)
                    + dy * est_->wireSignalVResistance(corner);
        c_per_meter = dx * est_->wireSignalHCapacitance(corner)
                    + dy * est_->wireSignalVCapacitance(corner);
        }
    }

    if (r_per_meter == 0.0 && c_per_meter == 0.0) {
        odb::dbTechLayer* layer = tech_->findRoutingLayer(layer_id);
        double w = block_->dbuToMicrons(layer->getWidth());
        if (layer->getResistance() != 0) {
            r_per_meter = layer->getResistance() / w * 1E+6;
        }
        c_per_meter = (w * layer->getCapacitance()
                     + 2 * layer->getEdgeCapacitance())
                     * 1E+6 * 1E-12;
    }

    res = static_cast<float>(len * r_per_meter);
    cap = static_cast<float>(len * c_per_meter);
}

bool SteinerTreeGenerator::isPadNet(const sta::Net* net) const
{
    const sta::Pin *pin1, *pin2;
    net2Pins(net, pin1, pin2);
    return pin1 && pin2
           && ((network_->isTopLevelPort(pin1) && isPadPin(pin2))
               || (network_->isTopLevelPort(pin2) && isPadPin(pin1)));
}

void SteinerTreeGenerator::net2Pins(const sta::Net* net,
                                    const sta::Pin*& pin1,
                                    const sta::Pin*& pin2) const
{
    pin1 = nullptr;
    pin2 = nullptr;
    NetConnectedPinIterator* pin_iter = network_->connectedPinIterator(net);
    if (pin_iter->hasNext())
        pin1 = pin_iter->next();
    if (pin_iter->hasNext())
        pin2 = pin_iter->next();
    delete pin_iter;
}

bool SteinerTreeGenerator::isPadPin(const sta::Pin* pin) const
{
    sta::Instance* inst = network_->instance(pin);
    return inst && !network_->isTopInstance(inst) && isPad(inst);
}

bool SteinerTreeGenerator::isPad(const sta::Instance* inst) const
{
    odb::dbInst* db_inst;
    odb::dbModInst* mod_inst;
    db_network_->staToDb(inst, db_inst, mod_inst);
    if (mod_inst)
        return false;
    const auto type = db_inst->getMaster()->getType().getValue();
    switch (type) {
        case odb::dbMasterType::CORE:
        case odb::dbMasterType::CORE_ANTENNACELL:
        case odb::dbMasterType::CORE_FEEDTHRU:
        case odb::dbMasterType::CORE_TIEHIGH:
        case odb::dbMasterType::CORE_TIELOW:
        case odb::dbMasterType::CORE_WELLTAP:
        case odb::dbMasterType::ENDCAP:
        case odb::dbMasterType::ENDCAP_PRE:
        case odb::dbMasterType::ENDCAP_POST:
        case odb::dbMasterType::CORE_SPACER:
        case odb::dbMasterType::BLOCK:
        case odb::dbMasterType::BLOCK_BLACKBOX:
        case odb::dbMasterType::BLOCK_SOFT:
        case odb::dbMasterType::ENDCAP_TOPLEFT:
        case odb::dbMasterType::ENDCAP_TOPRIGHT:
        case odb::dbMasterType::ENDCAP_BOTTOMLEFT:
        case odb::dbMasterType::ENDCAP_BOTTOMRIGHT:
        case odb::dbMasterType::ENDCAP_LEF58_BOTTOMEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_TOPEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_RIGHTEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_LEFTEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_RIGHTBOTTOMEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_LEFTBOTTOMEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_RIGHTTOPEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_LEFTTOPEDGE:
        case odb::dbMasterType::ENDCAP_LEF58_RIGHTBOTTOMCORNER:
        case odb::dbMasterType::ENDCAP_LEF58_LEFTBOTTOMCORNER:
        case odb::dbMasterType::ENDCAP_LEF58_RIGHTTOPCORNER:
        case odb::dbMasterType::ENDCAP_LEF58_LEFTTOPCORNER:
        case odb::dbMasterType::COVER:
        case odb::dbMasterType::RING:
            return false;
        case odb::dbMasterType::COVER_BUMP:
        case odb::dbMasterType::PAD:
        case odb::dbMasterType::PAD_AREAIO:
        case odb::dbMasterType::PAD_INPUT:
        case odb::dbMasterType::PAD_OUTPUT:
        case odb::dbMasterType::PAD_INOUT:
        case odb::dbMasterType::PAD_POWER:
        case odb::dbMasterType::PAD_SPACER:
            return true;
    }
    return false;
}

void SteinerTreeGenerator::makePadParasitic(const sta::Net* net)
{
    const sta::Pin *pin1, *pin2;
    net2Pins(net, pin1, pin2);
    for (sta::Scene* corner : sta_->scenes()) {
        sta::Parasitics* parasitics = corner->parasitics(min_max_);
        sta::Parasitic* parasitic
            = parasitics->makeParasiticNetwork(net, false);
        sta::ParasiticNode* n1
            = parasitics->ensureParasiticNode(parasitic, pin1, network_);
        sta::ParasiticNode* n2
            = parasitics->ensureParasiticNode(parasitic, pin2, network_);
        parasitics->makeResistor(parasitic, 1, 0.001, n1, n2);
        sta_->arcDelayCalc()->reduceParasitic(
            parasitic, net, corner, sta::MinMaxAll::all());
        parasitics->deleteParasiticNetwork(net);
    }
}
///debug Parasitics
void SteinerTreeGenerator::writeSpef(const std::string& base_filename)
{
    std::vector<std::ofstream> spef_files;
    std::map<sta::Scene*, std::ostream*> spef_streams;
    for (sta::Scene* corner : sta_->scenes()) {
        std::string filename = base_filename + "_" + corner->name() + ".spef";
        spef_files.emplace_back(filename);
        if (spef_files.back().is_open()) {
            spef_streams[corner] = &spef_files.back();
            logger_->info(TDP, 200, "Writing SPEF for {} to {}", corner->name(), filename);
        } else {
            logger_->warn(TDP, 201, "Cannot open SPEF file: {}", filename);
        }
    }

    if (spef_streams.empty()) {
        logger_->error(TDP, 202, "No SPEF files could be opened.");
        return;
    }

    sta::SpefWriter spef_writer(logger_, sta_, spef_streams);
    spef_writer.writeHeader();
    spef_writer.writePorts();

    int net_count = 0;
    sta::NetIterator* net_iter = network_->netIterator(network_->topInstance());
    while (net_iter->hasNext()) {
        sta::Net* net = net_iter->next();
        if (network_->isPower(net) || network_->isGround(net)) {
            continue;
        }

        for (sta::Scene* corner : sta_->scenes()) {
            sta::Parasitics* parasitics = corner->parasitics(min_max_);
            sta::Parasitic* parasitic = parasitics->findParasiticNetwork(net);
            if (parasitic) {
                spef_writer.writeNet(corner, net, parasitic, parasitics);
                net_count++;
            }
        }
    }
    delete net_iter;

    logger_->info(TDP, 203, "SPEF written: {} nets to {}_*.spef", net_count, base_filename);
}

void SteinerTreeGenerator::clearAllParasitics()
{
    for (sta::Scene* corner : sta_->scenes()) {
        sta::NetIterator* iter = network_->netIterator(network_->topInstance());
        while (iter->hasNext()) {
            sta::Net* net = iter->next();
            corner->parasitics(min_max_)->deleteParasiticNetwork(net);
        }
        delete iter;
    }
}

} // namespace tdp
