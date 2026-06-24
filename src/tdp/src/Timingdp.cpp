// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, TDP Authors

#include "tdp/Timingdp.h"
#include "infrastructure/Architecture3D.h"
#include "infrastructure/Network3D.h"
#include "infrastructure/Common3D.h"
#include "infrastructure/Net3D.h"
#include "infrastructure/SegmentManager3D.h"
#include "infrastructure/Segment3D.h"
#include "infrastructure/Node3D.h"
#include "eva/eva.h"
#include "eva/TimingEngine3D.h"
#include "hbta/HBTACostModel.h"
#include "placement/placement.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include "grt/GRoute.h"
#include "grt/GlobalRouter.h"
#include "sta/Sta.hh"
#include "sta/Search.hh"
#include "sta/Graph.hh"
#include "rsz/Resizer.hh"
#include "odb/db.h"
#include "odb/defout.h"
#include "db_sta/dbSta.hh"
#include "utl/Logger.h"
#include "lgista/LocalGraphIsta.hh"

using utl::TDP;

namespace tdp {

// ============================================================================
// Timingdp实现（使用新的转换系统）
// ============================================================================

Timingdp::Timingdp(odb::dbDatabase* db,
                   utl::Logger* logger, 
                   sta::dbSta* sta, 
                   rsz::Resizer* resizer, 
                   grt::GlobalRouter* global_router,
                   dpl::Opendp* opendp) 
    : db_(db), sta_(sta), logger_(logger), 
      rs_(resizer), global_router_(global_router),
      estimate_parasitics_(resizer->getEstimateParasitics()), 
      opendp_(opendp),
      segment_manager_(nullptr) {
    // 初始化3D架构和网络（类似DPL的Opendp）
    timing_engine_ = std::make_unique<TimingEngine3D>(sta_, this);
    arch_3d_ = std::make_unique<Architecture3D>();
    network_3d_ = std::make_unique<Network3D>();
}

void Timingdp::setSta(sta::dbSta* sta) {
    sta_ = sta;
}

Timingdp::~Timingdp() {
    if (segment_manager_) {
        delete segment_manager_;
        segment_manager_ = nullptr;
    }
    destroyDatabase();
}

// ============================================================================
// 数据库管理（类似DPL的Opendp）
// ============================================================================
void Timingdp::detailedPlacement(const int max_displacement_x,
                               const int max_displacement_y,
                               const std::string& report_file_name,
                               int margin_um)
{
  importDb();
//   adjustNodesOrient();
//   for (const auto& node : network_3d_->getNode()) {
//     if (node->getType() == Node::CELL && !node->isFixed()) {
//       node->setPlaced(false);
//     }
//   }

    if (max_displacement_x == 0 || max_displacement_y == 0) {
    // defaults
    max_displacement_x_ = 500;
    max_displacement_y_ = 100;
    } else {
    max_displacement_x_ = max_displacement_x;
    max_displacement_y_ = max_displacement_y;
    }

    logger_->info(TDP, 6, "import right");

    std::string hbt_array_json = "simple/hbt_array.json";
    if (!initializeHBTAllocator(hbt_array_json, 1000)) {
        std::cerr << "WARNING: Failed to initialize HBT allocator, but continuing..." << std::endl;
    } else {
        std::cout << "✓ HBT allocator initialized successfully" << std::endl;
        
        // 执行HBT分配
        auto allocation_result = allocateHBTsForDesign();
        
        std::cout << "   - Total inter-layer nets: " << allocation_result.total_nets << std::endl;
        std::cout << "   - Successful allocations: " << allocation_result.successful_allocations << std::endl;
        std::cout << "   - Failed allocations: " << allocation_result.failed_allocations << std::endl;
        
        if (allocation_result.total_nets > 0) {
            std::cout << "   - Success rate: " 
                        << (static_cast<double>(allocation_result.successful_allocations) / allocation_result.total_nets * 100.0)
                        << "%" << std::endl;
            std::cout << "   - Total estimated wirelength: " << allocation_result.total_wirelength << " microns" << std::endl;
        }
    }
    rs_->resizeSlackPreamble();

    // 先初始化 GRT 的 grid 和 db_net_map_
    int min_layer, max_layer;
    global_router_->getMinMaxLayer(min_layer, max_layer);
    global_router_->initFastRoute(min_layer, max_layer);
    SteinerTreeGenerator steiner_generator(this);
    steiner_generator.createSteninerRoute();
    executeTimingDriven(true);
    auto stats = timing_engine_->getTimingStats();
    // runTimingDrivenHBTReallocation(last_hbt_connections_);
    
    

    Placement3D placetdp(this);
    placetdp.rundetailedplacement();
    steiner_generator.createSteninerRoute();
    sta_->updateTiming(true);  
    auto stats2 = timing_engine_->getTimingStats();
    logger_->info(TDP, 1007,
        "Early: WNS={:.20f}ns TNS={:.6f}ns -> WNS={:.6f}ns (Δ={:.6f}) TNS={:.6f}ns (Δ={:.6f})",
        stats.early_wns * 1e9, stats.early_tns * 1e9,
        stats2.early_wns * 1e9, (stats2.early_wns - stats.early_wns)/stats.early_wns,
        stats2.early_tns * 1e9 , (stats2.early_tns - stats.early_tns)/stats.early_tns);
    logger_->info(TDP, 1008,
        "Late : WNS={:.6f}ns TNS={:.6f}ns -> WNS={:.6f}ns (Δ={:.6f}) TNS={:.6f}ns (Δ={:.6f})",
        stats.late_wns * 1e9, stats.late_tns * 1e9,
        stats2.late_wns * 1e9, (stats2.late_wns - stats.late_wns)/stats.late_wns,
        stats2.late_tns * 1e9, (stats2.late_tns - stats.late_tns)/stats.late_tns);

    odb::DefOut def_before(logger_);
    def_before.setVersion(odb::DefOut::DEF_5_8);
    def_before.writeBlock(block_, "before_legalization.def");

    legalizeTier("bottom", margin_um, max_displacement_x_, max_displacement_y_);
    legalizeTier("upper", margin_um, max_displacement_x_, max_displacement_y_);

    odb::DefOut def_after(logger_);
    def_after.setVersion(odb::DefOut::DEF_5_8);
    def_after.writeBlock(block_, "after_legalization.def");
    steiner_generator.createSteninerRoute();
    sta_->updateTiming(true);  
    auto stats3 = timing_engine_->getTimingStats();
    logger_->info(TDP, 1009,
        "Late : WNS={:.20f} TNS={:.20f} -> WNS={:.20f} (Δ={:.6f}) TNS={:.20f} (Δ={:.6f})",
        stats.late_wns, stats.late_tns,
        stats3.late_wns, (stats3.late_wns - stats.late_wns)/stats.late_wns,
        stats3.late_tns, (stats3.late_tns - stats.late_tns)/stats.late_tns);
}

void Timingdp::legalizeTier(const std::string& tier,
                            int margin_um,
                            int max_displacement_x,
                            int max_displacement_y)
{
    std::string site_name = (tier == "bottom")
        ? "FreePDK45_38x28_10R_NP_162NW_34O"
        : "asap7sc7p5t";
    odb::dbSite* target_site = nullptr;
    for (auto* lib : db_->getLibs()) {
        target_site = lib->findSite(site_name.c_str());
        if (target_site) break;
    }
    if (!target_site) {
        logger_->error(TDP, 700, "Cannot find site '{}'", site_name);
    }

    for (auto* row : block_->getRows()) {
        odb::dbRow::destroy(row);
    }

    int margin_dbu = block_->micronsToDbu(margin_um);
    odb::Rect die = block_->getDieArea();
    odb::Rect core(die.xMin() + margin_dbu,
                   die.yMin() + margin_dbu,
                   die.xMax() - margin_dbu,
                   die.yMax() - margin_dbu);

    sta::dbNetwork* network = sta_->getDbNetwork();
    ifp::InitFloorplan floorplan(block_, logger_, network);
    floorplan.makeRows(core, target_site);

    std::string other_tier = (tier == "bottom") ? "upper" : "bottom";
    std::vector<std::pair<odb::dbMaster*, odb::dbMasterType>> saved_types;
    for (auto* lib : db_->getLibs()) {
        for (auto* master : lib->getMasters()) {
            std::string name = master->getName();
            if (name.find(other_tier) != std::string::npos) {
                saved_types.emplace_back(master, master->getType());
                master->setType(odb::dbMasterType::COVER);
            }
        }
    }

    opendp_->detailedPlacement(0, 0, "", false);

    for (auto& [master, type] : saved_types) {
        master->setType(type);
    }
}

bool Timingdp::destroyDatabase() {
    importClear();  // 清除转换数据（类似DPL的importClear）
    
    // if (db_) {
    //     odb::dbDatabase::destroy(db_);
    //     db_ = nullptr;
    //     chip_ = nullptr;
    //     block_ = nullptr;
    // }
    
    //logger_->info("Destroyed 3D placement database");
    return true;
}

// ============================================================================
// 3D PDK参数提取（仿照DPL从数据库获取）
// ============================================================================

bool Timingdp::extractPDKParametersFromDatabase() {
    if (!block_) {
        //logger_->error("No block available for PDK parameter extraction");
        return false;
    }
    
    try {
        logger_->info(TDP, 13, "Extracting PDK parameters from database...");
        
        // 清空现有配置
        pdk_config_.clear();
        
        // 分析所有row，提取层信息（仿照DPL的dbToOpendp）
        std::map<DBU, LayerInfo> height_to_layer_info;
        
        for (odb::dbRow* row : block_->getRows()) {
            if (!row->getSite()) {
                //logger_->warn("Row has no site, skipping");
                continue;
            }
            
            odb::dbSite* site = row->getSite();
            DBU height = site->getHeight();
            DBU width = site->getWidth();
            std::string tech_lib = site->getLib()->getName();
            
            // 记录不同高度的层（仿照DPL处理混合工艺）
            if (height_to_layer_info.find(height) == height_to_layer_info.end()) {
                LayerInfo layer_info;
                layer_info.height = height;
                layer_info.site_width = width;
                layer_info.site_name = site->getName();
                layer_info.tech_lib = tech_lib;
                layer_info.row_count = 1;
                
                // 这里我们根据工艺库名称是否匹配用户设置的层工艺库来判断
                if (!bottom_tech_lib_.empty() && tech_lib.find(bottom_tech_lib_) != std::string::npos) {
                    layer_info.layer_type = "bottom";
                } else if (!upper_tech_lib_.empty() && tech_lib.find(upper_tech_lib_) != std::string::npos) {
                    layer_info.layer_type = "upper";
                } else {
                    // 默认：如果无法确定，标记为unknown
                    layer_info.layer_type = "unknown";
                }
                
                height_to_layer_info[height] = layer_info;
            } else {
                height_to_layer_info[height].row_count++;
            }
            
            //logger_->debugCheck("Row site: " + site->getName() + 
            //              ", width=" + std::to_string(width) + 
            //              ", height=" + std::to_string(height) +
            //              ", tech_lib=" + tech_lib);
        }
        
        // 将层信息存储到PDK配置中
        int layer_index = 0;
        for (const auto& [height, layer_info] : height_to_layer_info) {
            PDKLayerConfig layer_config;
            layer_config.layer_id = layer_index++;
            layer_config.layer_name = layer_info.layer_type + "_layer_" + std::to_string(layer_index);
            layer_config.height = layer_info.height;
            layer_config.site_width = layer_info.site_width;
            layer_config.site_name = layer_info.site_name;
            layer_config.tech_lib = layer_info.tech_lib;
            layer_config.row_count = layer_info.row_count;
            
            pdk_config_.layers.push_back(layer_config);
            
            logger_->info(TDP, 8, "Detected layer: " + layer_config.layer_name + 
                          " (height=" + std::to_string(layer_config.height) + 
                          ", width=" + std::to_string(layer_config.site_width) + 
                          ", tech_lib=" + layer_config.tech_lib +
                          ", rows=" + std::to_string(layer_config.row_count) + ")");
        }
        
        // 从数据库提取HBT间距和层间间距（仿照DPL从LEF文件获取）
        extract3DSpacingParametersFromDatabase();
        
        // 更新3D架构
        if (arch_3d_ && !pdk_config_.layers.empty()) {
            arch_3d_->setHBTPitch(pdk_config_.hbt_pitch);
            arch_3d_->setInterLayerSpacing(pdk_config_.inter_layer_spacing);
            
            // 设置层特定参数
            for (const auto& layer_config : pdk_config_.layers) {
                Architecture3D::LayerType layer_type = 
                    (layer_config.layer_name.find("bottom") != std::string::npos) ? 
                    Architecture3D::LAYER_BOTTOM : Architecture3D::LAYER_UPPER;
                
                arch_3d_->setSiteWidth(layer_type, DbuX{layer_config.site_width});
                // 设置行高度（通过Row3D类）
                // arch_3d_->setHeight(DbuY{layer_config.height}); // 这个方法不存在
            }
        }
        
        // logger_->info("Successfully extracted PDK parameters: " + 
        //               std::to_string(pdk_config_.layers.size()) + " layers detected");
        return true;
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to extract PDK parameters: " + std::string(e.what()));
        return false;
    }
}

// ============================================================================
// 数据访问
// ============================================================================

int64_t Timingdp::getDieWidth() const {
    if (!block_) return 0;
    
    odb::Rect die_area = block_->getDieArea();
    return die_area.xMax() - die_area.xMin();
}

int64_t Timingdp::getDieHeight() const {
    if (!block_) return 0;
    
    odb::Rect die_area = block_->getDieArea();
    return die_area.yMax() - die_area.yMin();
}

std::vector<odb::dbInst*>& Timingdp::getMovableInstances() {
    return movable_insts_;
}

std::vector<odb::dbInst*>& Timingdp::getFixedInstances() {
    return fixed_insts_;
}

std::vector<odb::dbInst*>& Timingdp::getAllInstances() {
    return all_insts_;
}

std::vector<odb::dbNet*>& Timingdp::getNets() {
    return nets_;
}

// ============================================================================
// 3D特定操作
// ============================================================================


bool Timingdp::setInstanceLocation3D(odb::dbInst* inst, int64_t x, int64_t y, int64_t z) {
    if (!inst) return false;
    
    // 设置XY位置（Z位置目前只用于内部表示）
    inst->setLocation(x, y);
    
    // 更新3D网络中的位置
    if (network_3d_) {
        auto* node = network_3d_->getCell(inst);
        if (node) {
            // 确保Z坐标是0或1（bottom层=0, upper层=1）
            DBU clamped_z = (z > 0) ? 1 : 0;
            node->location = Pointdp3D(DbuX{x}, DbuY{y}, DbuZ{clamped_z});
        }
    }
    
    return true;
}


// ============================================================================
// 日志管理
// ============================================================================

utl::Logger* Timingdp::getLogger() const {
    return logger_;
}
// bool Timingdp::addHBTConnection(const std::string& bottom_inst, const std::string& upper_inst) {
//     // TODO: 实现HBT连接添加
//     logger_->warn("HBT connection addition not yet implemented");
//     return true;
// }

// std::vector<std::pair<std::string, std::string>> Timingdp::getHBTConnections() const {
//     // TODO: 实现HBT连接获取
//     return {};
// }

// ============================================================================
// 验证和调试
// ============================================================================

bool Timingdp::check3DPlacementLegality() const {
    if (!arch_3d_ || !network_3d_) {
        return false;
    }
    
    return network_3d_->checkPlacementLegality(arch_3d_.get());
}

/////

bool Timingdp::setupLayerArchitecture() {
    if (!arch_3d_ || !block_) {
        return false;
    }
    
    // 首先提取PDK参数（仿照DPL）
    if (!extractPDKParametersFromDatabase()) {
        // logger_->error("Failed to extract PDK parameters for layer architecture");
        return false;
    }
    
    // 根据提取的参数设置层架构
    for (const auto& layer_config : pdk_config_.layers) {
        // 确定层类型
        Architecture3D::LayerType layer_type;
        if (layer_config.layer_name.find("bottom") != std::string::npos) {
            layer_type = Architecture3D::LAYER_BOTTOM;
        } else if (layer_config.layer_name.find("upper") != std::string::npos) {
            layer_type = Architecture3D::LAYER_UPPER;
        } else {
            // 未知层类型，根据高度推断
            layer_type = (layer_config.height <= 150) ? 
                         Architecture3D::LAYER_BOTTOM : Architecture3D::LAYER_UPPER;
        }
        
        // 添加层（如果尚未添加）
        if (!arch_3d_->hasLayer(layer_type)) {
            arch_3d_->addLayer(layer_type, layer_config.layer_name);
        }
        
        // 设置层参数（仿照DPL从数据库获取）
        arch_3d_->setSiteWidth(layer_type, DbuX{layer_config.site_width});
        // arch_3d_->setHeight(DbuY{layer_config.height}); // 这个方法不存在
        
        logger_->info(TDP, 9, "Configured layer " + layer_config.layer_name + 
                      ": width=" + std::to_string(layer_config.site_width) + 
                      ", height=" + std::to_string(layer_config.height));
    }
    
    // 确保至少有两个层（对于3D设计）
    if (!arch_3d_->hasLayer(Architecture3D::LAYER_BOTTOM)) {
        arch_3d_->addLayer(Architecture3D::LAYER_BOTTOM, "default_bottom");
        //logger_->warn("Added default bottom layer");
    }
    
    if (!arch_3d_->hasLayer(Architecture3D::LAYER_UPPER)) {
        arch_3d_->addLayer(Architecture3D::LAYER_UPPER, "default_upper");
        //logger_->warn("Added default upper layer");
    }
    
    return true;
}

void Timingdp::cacheDesignData() {
    logger_ ->info (TDP, 14, "Importing database to cacheDesignDat ......");
    cacheInstances();
    cacheNets();
}

void Timingdp::cacheInstances() {
    if (!block_) return;
    
    movable_insts_.clear();
    fixed_insts_.clear();
    all_insts_.clear();
    
    for (auto* inst : block_->getInsts()) {
        all_insts_.push_back(inst);
        
        if (inst->isFixed()) {
            fixed_insts_.push_back(inst);
        } else {
            movable_insts_.push_back(inst);
        }
    }
}

void Timingdp::cacheNets() {
    if (!block_) return;
    
    nets_.clear();
    
    for (auto* net : block_->getNets()) {
        nets_.push_back(net);
    }
}

void Timingdp::cleanup() {
    movable_insts_.clear();
    fixed_insts_.clear();
    all_insts_.clear();
    nets_.clear();
}

// ============================================================================
// 3D间距参数提取（仿照DPL从LEF文件获取）
// ============================================================================

void Timingdp::extract3DSpacingParametersFromDatabase() {
    if (!db_) {
        //logger_->warn("Database or chip not available for 3D spacing extraction");
        return;
    }
    
    try {
        //logger_->info("Extracting 3D spacing parameters from database...");
        
        // 重置间距参数
        pdk_config_.hbt_pitch = hbt_pitch_;  // 使用用户设置的HBT间距
        pdk_config_.inter_layer_spacing = inter_layer_spacing_;  // 使用用户设置的层间间距
        
        // 获取所有技术库
        std::vector<odb::dbLib*> libs;
        for (auto lib : db_->getLibs()) {
            libs.push_back(lib);
        }
        
        // 查找hb_layer层（从LEF文件中）- 仅用于信息记录，不覆盖用户设置
        bool found_hb_layer = false;
        for (odb::dbLib* lib : libs) {
            odb::dbTech* tech = lib->getTech();
            if (!tech) continue;
            
            odb::dbTechLayer* hb_layer = tech->findLayer("hb_layer");
            if (hb_layer) {
                DBU lef_hbt_pitch = hb_layer->getSpacing();
                found_hb_layer = true;
                //logger_->info("Found hb_layer in LEF: spacing=" + std::to_string(lef_hbt_pitch) + " DBU");
                
                // 如果用户没有设置HBT间距，使用LEF文件中的值
                if (hbt_pitch_ == 0) {
                    pdk_config_.hbt_pitch = lef_hbt_pitch;
                    hbt_pitch_ = lef_hbt_pitch;  // 更新用户变量
                    //logger_->info("Using HBT pitch from LEF: " + std::to_string(lef_hbt_pitch) + " DBU");
                } else {
                    //logger_->info("Using user-specified HBT pitch: " + std::to_string(hbt_pitch_) + " DBU (LEF value: " + std::to_string(lef_hbt_pitch) + " DBU)");
                }
                break;
            }
        }
        
        // 如果没有找到hb_layer，使用用户设置的值
        if (!found_hb_layer) {
            if (hbt_pitch_ == 0) {
                //logger_->warn("hb_layer not found in LEF files and no user HBT pitch specified. Using default 0.");
            } else {
                //logger_->info("Using user-specified HBT pitch: " + std::to_string(hbt_pitch_) + " DBU");
            }
        }
        
        // logger_->info("3D spacing parameters: HBT pitch=" + std::to_string(pdk_config_.hbt_pitch) +
        //               " DBU, inter-layer spacing=" + std::to_string(pdk_config_.inter_layer_spacing) + " DBU");
        
    } catch (const std::exception& e) {
        //logger_->error("Failed to extract 3D spacing parameters: " + std::string(e.what()));
    }
}

// ============================================================================
// HBT分配接口实现（新增）
// ============================================================================

bool Timingdp::initializeHBTAllocator(const std::string& hbt_array_json, int dbu_per_micron) {
    //logger_->info("Initializing HBT allocator with JSON: " + hbt_array_json);

    // 读取 HBT array，并创建 HBTA allocator。
    hbt_interface_ = std::make_unique<hbta::HBTDatabaseInterface>();
    
    // 这里只初始化接口；真正分配在 allocateHBTsForDesign() 中触发。
    if (!hbt_interface_->initialize(this, hbt_array_json, dbu_per_micron)) {
        //logger_->error("Failed to initialize HBT database interface");
        hbt_interface_.reset();
        return false;
    }
    
    //logger_->info("HBT allocator initialized successfully");
    return true;
}

std::vector<hbta::InterLayerConnection> Timingdp::extractInterLayerConnections() {

    if (!hbt_interface_) {
        return {};
    }
    
    // 将需要 HBT 的跨层 net 转成 HBTA 使用的连接描述。
    auto connections = hbt_interface_->extractInterLayerConnections();
    
    // 打印提取统计，便于检查跨层 net 是否正确识别。
    auto stats = hbt_interface_->getExtractionStats();
    stats.print();
    
    return connections;
}

hbta::HBTAllocator* Timingdp::getHBTAllocator() {
    if (!hbt_interface_) {
        return nullptr;
    }

    // allocator 的所有权仍在 hbt_interface_ 中。
    return &hbt_interface_->getAllocator();
}

const hbta::HBTAllocator* Timingdp::getHBTAllocator() const {
    if (!hbt_interface_) {
        return nullptr;
    }

    // const 版本用于只读访问 allocator 状态。
    return &hbt_interface_->getAllocator();
}

hbta::HBTAllocator::BatchNetAllocationResult Timingdp::allocateHBTsForDesign() {
    hbta::HBTAllocator::BatchNetAllocationResult result;
    
    // 获取HBT分配器
    auto* hbt_allocator = getHBTAllocator();
    if (!hbt_allocator) {
        //logger_->error("HBT allocator not available");
        result.total_nets = 0;
        result.failed_allocations = 0;
        return result;
    }
    
    // 提取跨层 net，并缓存给后续 timing-driven 重匹配复用。
    auto connections = extractInterLayerConnections();
    last_hbt_connections_ = connections;
    
    if (connections.empty()) {
        result.total_nets = 0;
        result.successful_allocations = 0;
        result.failed_allocations = 0;
        return result;
    }

    // 初始分配：HBTA 内部使用四叉树分区 + 局部 Hungarian。
    result = hbt_allocator->allocateHBTsForNets(connections);
    
    // 获取分配器统计信息
    auto stats = hbt_allocator->getStatistics();

    // 同步回 Network3D，后续 routing/STA 才能看到 HBT 连接。
    if (result.successful_allocations > 0) {
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
    }

    return result;
}

std::unordered_map<std::string, double> Timingdp::collectHBTNetSlacks(
    const std::vector<hbta::InterLayerConnection>& connections) {
    std::unordered_map<std::string, double> slack_by_net;
    if (!rs_) {
        return slack_by_net;
    }

    std::unordered_map<std::string, odb::dbNet*> db_net_by_name;
    db_net_by_name.reserve(nets_.size());
    for (auto* db_net : nets_) {
        if (db_net != nullptr) {
            db_net_by_name.emplace(db_net->getName(), db_net);
        }
    }

    for (const auto& connection : connections) {
        auto net_it = db_net_by_name.find(connection.net_name);
        if (net_it == db_net_by_name.end()) {
            continue;
        }

        auto slack = rs_->resizeNetSlack(net_it->second);
        if (!slack.has_value() || !std::isfinite(static_cast<double>(*slack))) {
            continue;
        }

        slack_by_net[connection.net_name] = static_cast<double>(*slack);

        if (network_3d_) {
            if (auto* net_3d = network_3d_->getNet(connection.net_name)) {
                net_3d->updateTimingInfo(static_cast<double>(*slack));
            }
        }
    }

    return slack_by_net;
}

bool Timingdp::runTimingDrivenHBTReallocation(
    const std::vector<hbta::InterLayerConnection>& connections) {
    auto* allocator = getHBTAllocator();
    if (!allocator || connections.empty()) {
        return false;
    }

    const auto old_stats = timing_engine_->getTimingStats();

    auto slack_by_net = collectHBTNetSlacks(connections);
    if (slack_by_net.empty()) {
        logger_->warn(TDP,
                      1070,
                      "Timing-driven HBT reallocation skipped: no STA slack was available for inter-layer nets.");
        return false;
    }

    hbta::HBTAllocator::TimingDrivenReallocationConfig config;
    config.critical_slack_threshold = 0.0;
    config.critical_nets_percent = 100.0;
    config.local_bonding_box_radius_in_hbt_pitch = 2.0;

    if (allocator->hasPendingTimingDrivenReallocation()) {
        allocator->rollbackPendingTimingDrivenReallocation();
    }

    const auto problem = allocator->prepareTimingDrivenReallocationProblem(
        connections, slack_by_net, config);

    if (problem.empty()) {
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
        logger_->info(TDP,
                      1071,
                      "Timing-driven HBT reallocation produced no candidate.");
        return false;
    }

    const auto baseline_assignment = allocator->captureAssignmentSnapshot();
    const double failed_sta_cost = std::numeric_limits<double>::max() / 16.0;
    std::vector<std::vector<double>> sta_cost_matrix(
        problem.connections.size(),
        std::vector<double>(problem.candidate_hbts.size(), failed_sta_cost));

    std::unordered_map<std::string, double> criticality_by_net;
    criticality_by_net.reserve(problem.connections.size());
    for (const auto& connection : problem.connections) {
        double criticality = 0.0;
        auto slack_it = slack_by_net.find(connection.net_name);
        if (slack_it != slack_by_net.end()) {
            // 与 HBTA critical net 筛选使用同一套 criticality 公式。
            criticality = hbta::HBTACostModel::criticalityFromSlack(
                slack_it->second,
                problem.worst_slack,
                config.critical_slack_threshold,
                config.normalization_epsilon);
        }
        criticality_by_net[connection.net_name] = criticality;
    }

    auto restore_baseline_assignment = [&]() {
        allocator->restoreAssignmentSnapshot(baseline_assignment);
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
    };

    for (size_t net_idx = 0; net_idx < problem.connections.size(); ++net_idx) {
        bool row_has_valid_sta_cost = false;
        const auto& candidate_connection = problem.connections[net_idx];

        for (size_t hbt_idx = 0; hbt_idx < problem.candidate_hbts.size(); ++hbt_idx) {
            restore_baseline_assignment();

            const auto& candidate_hbt = problem.candidate_hbts[hbt_idx];
            if (!allocator->forceAssignHBTToNetForEvaluation(
                    candidate_connection, candidate_hbt.hbt_id)) {
                continue;
            }

            hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
            SteinerTreeGenerator trial_steiner_generator(this);
            trial_steiner_generator.createSteninerRoute();
            if (!executeTimingDriven(true)) {
                continue;
            }

            const auto trial_slack_by_net = collectHBTNetSlacks(connections);
            auto slack_it = trial_slack_by_net.find(candidate_connection.net_name);
            if (slack_it == trial_slack_by_net.end()
                || !std::isfinite(slack_it->second)) {
                continue;
            }

            auto criticality_it = criticality_by_net.find(candidate_connection.net_name);
            const double criticality = criticality_it == criticality_by_net.end()
                                           ? 0.0
                                           : criticality_it->second;

            // Hungarian 最小化 cost，因此加权后的 STA slack 越大，cost 越小。
            sta_cost_matrix[net_idx][hbt_idx] = -criticality * slack_it->second;
            row_has_valid_sta_cost = true;
        }

        if (!row_has_valid_sta_cost) {
            restore_baseline_assignment();
            SteinerTreeGenerator restore_steiner_generator(this);
            restore_steiner_generator.createSteninerRoute();
            executeTimingDriven(true);
            logger_->warn(TDP,
                          1075,
                          "Timing-driven HBT reallocation skipped: no valid STA candidate cost for net {}.",
                          candidate_connection.net_name);
            return false;
        }
    }

    restore_baseline_assignment();

    auto result = allocator->commitTimingDrivenReallocationWithCostMatrix(
        problem, sta_cost_matrix);

    if (!allocator->hasPendingTimingDrivenReallocation()
        || result.total_nets == 0
        || !result.allSuccessful()) {
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
        SteinerTreeGenerator restore_steiner_generator(this);
        restore_steiner_generator.createSteninerRoute();
        executeTimingDriven(true);
        logger_->info(TDP,
                      1076,
                      "Timing-driven HBT reallocation produced no valid assignment.");
        return false;
    }

    hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
    SteinerTreeGenerator steiner_generator(this);
    steiner_generator.createSteninerRoute();
    if (!executeTimingDriven(true)) {
        allocator->rollbackPendingTimingDrivenReallocation();
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
        steiner_generator.createSteninerRoute();
        executeTimingDriven(true);
        logger_->warn(TDP,
                      1072,
                      "Timing-driven HBT reallocation rolled back: tentative timing update failed.");
        return false;
    }

    const auto new_stats = timing_engine_->getTimingStats();
    hbta::HBTAllocator::TimingDrivenAcceptanceConfig accept_config;
    const bool accept_reallocation =
        hbta::HBTAllocator::shouldAcceptTimingDrivenReallocation(
            old_stats.late_wns,
            old_stats.late_tns,
            new_stats.late_wns,
            new_stats.late_tns,
            accept_config);

    if (!accept_reallocation) {
        allocator->rollbackPendingTimingDrivenReallocation();
        hbt_interface_->syncHBTAssignmentsToNetwork3D(connections);
        steiner_generator.createSteninerRoute();
        executeTimingDriven(true);
        logger_->info(TDP,
                      1073,
                      "Timing-driven HBT reallocation rejected: WNS {:.20f}->{:.20f}, TNS {:.20f}->{:.20f}.",
                      old_stats.late_wns,
                      new_stats.late_wns,
                      old_stats.late_tns,
                      new_stats.late_tns);
        return false;
    }

    allocator->acceptPendingTimingDrivenReallocation();

    logger_->info(TDP,
                  1074,
                  "Timing-driven HBT reallocation accepted for {} nets: WNS {:.20f}->{:.20f}, TNS {:.20f}->{:.20f}.",
                  result.successful_allocations,
                  old_stats.late_wns,
                  new_stats.late_wns,
                  old_stats.late_tns,
                  new_stats.late_tns);
    return true;
}

/////////////////////
///////Segment///////
////////////////////

void Timingdp::initializeSegmentManager() {
    if (arch_3d_) {
        // 动态分配段管理器
        if (segment_manager_) {
            delete segment_manager_;  // 清理旧的（如果存在）
        }
        logger_->info(TDP, 15, "Segment manager initialized");
        segment_manager_ = new SegmentManager3D(arch_3d_.get());
        segment_manager_->initializeSegments();
        //logger_->info("Segment manager initialized");
    } else {
        //logger_->error("Cannot initialize segment manager: Architecture3D not available");
    }
}

// 实现从 .h 文件移出的方法
std::vector<Segment3D*> Timingdp::getCellSegments(const Node3D* cell) const {
    return segment_manager_ ? segment_manager_->getCellSegments(cell) : std::vector<Segment3D*>();
}

bool Timingdp::addCellToSegment(Node3D* cell, int seg_id) {
    return segment_manager_ ? segment_manager_->addCellToSegment(cell, seg_id) : false;
}

bool Timingdp::removeCellFromSegment(const Node3D* cell, int seg_id) {
    return segment_manager_ ? segment_manager_->removeCellFromSegment(cell, seg_id) : false;
}

void Timingdp::initializeSegments() {
    if (segment_manager_) {
        segment_manager_->initializeSegments();
    }
}

void Timingdp::createnewrow(){
    //设置层架构（仿照DPL的createArchitecture）
    if (!setupLayerArchitecture()) {
        logger_->error(TDP, 591, "Failed to setup layer architecture");
    }
    
    //仿照DPL的dbToOpendp::createArchitecture）
    if (!convertRows()) {
        logger_->error(TDP, 592, "Failed to convert rows");
    }

}

void Timingdp::reportLegalizationStats() const
{
  logger_->report("Placement Analysis");
  logger_->report("---------------------------------");
//   logger_->report("total displacement   {:10.1f} u",
//                   block_->dbuToMicrons(displacement_sum_));
//   logger_->metric("design__instance__displacement__total",
//                   block_->dbuToMicrons(displacement_sum_));
//   logger_->report("average displacement {:10.1f} u",
//                   block_->dbuToMicrons(displacement_avg_));
//   logger_->metric("design__instance__displacement__mean",
//                   block_->dbuToMicrons(displacement_avg_));
//   logger_->report("max displacement     {:10.1f} u",
//                   block_->dbuToMicrons(displacement_max_));
//   logger_->metric("design__instance__displacement__max",
//                   block_->dbuToMicrons(displacement_max_));
//   logger_->report("original HPWL        {:10.1f} u",
//                   block_->dbuToMicrons(hpwl_before_));
//   odb::WireLengthEvaluator eval(block_);
//   const double hpwl_legal = eval.hpwl();
//   logger_->report("legalized HPWL       {:10.1f} u",
//                   block_->dbuToMicrons(hpwl_legal));
//   logger_->metric("route__wirelength__estimated",
//                   block_->dbuToMicrons(hpwl_legal));
//   const int hpwl_delta
//       = (hpwl_before_ == 0.0)
//             ? 0.0
//             : round((hpwl_legal - hpwl_before_) / hpwl_before_ * 100);
//   logger_->report("delta HPWL           {:10} %", hpwl_delta);
  logger_->report("");
}


//////////////
/////test sta
/////////////
bool Timingdp::executeTimingDriven(bool run_journal_restore)
{
    // rs_->findResizeSlackstdp(run_journal_restore);
    sta_->ensureGraph();
    // 更新时序
    sta_->updateTiming(true);
    sta_->findRequireds();

    return true;
}
void Timingdp::executeTimingDrivenI(Node3D* cell, int depth)
{
    if (!cell || !cell->db_inst) return;
    sta::Instance* sta_inst = sta_->getDbNetwork()->dbToSta(cell->db_inst);
    if (!sta_inst) return;

    if (!ista_)
        ista_ = std::make_unique<lgista::LocalGraphIsta>(sta_);

    static bool compare_env = (std::getenv("LGISTA_COMPARE") != nullptr);
    if (compare_env) {
        struct PinSlew { float rise; float fall; };
        std::map<std::string, PinSlew> old_slew;
        sta::DcalcAPIndex ap0 = 0;
        {
            sta_->updateTiming(false);
            sta::InstancePinIterator* pit = sta_->network()->pinIterator(sta_inst);
            while (pit->hasNext()) {
                sta::Pin* pin = pit->next();
                sta::Vertex *v, *bidir;
                sta_->graph()->pinVertices(pin, v, bidir);
                if (v) {
                    std::string name = sta_->network()->name(pin);
                    old_slew[name].rise = sta::delayAsFloat(
                        sta_->graph()->slew(v, sta::RiseFall::rise(), ap0));
                    old_slew[name].fall = sta::delayAsFloat(
                        sta_->graph()->slew(v, sta::RiseFall::fall(), ap0));
                }
            }
            delete pit;
        }

        ista_->analyzeAndCommit(sta_inst, depth);

        double max_slew_delta = 0;
        int pin_count = 0;
        sta::InstancePinIterator* pit2 = sta_->network()->pinIterator(sta_inst);
        while (pit2->hasNext()) {
            sta::Pin* pin = pit2->next();
            sta::Vertex *v, *bidir;
            sta_->graph()->pinVertices(pin, v, bidir);
            if (!v) continue;
            std::string name = sta_->network()->name(pin);
            auto it = old_slew.find(name);
            if (it == old_slew.end()) continue;
            double new_r = sta::delayAsFloat(sta_->graph()->slew(v, sta::RiseFall::rise(), ap0));
            double new_f = sta::delayAsFloat(sta_->graph()->slew(v, sta::RiseFall::fall(), ap0));
            max_slew_delta = std::max(max_slew_delta,
                std::max(std::abs(new_r - (double)it->second.rise),
                         std::abs(new_f - (double)it->second.fall)));
            ++pin_count;
        }
        delete pit2;
        printf("[LGISTA_COMPARE] cell=%s pins=%d max_slew_delta=%.4e\n",
               sta_->network()->name(sta_inst), pin_count, max_slew_delta);
        fflush(stdout);
        return;
    }

    ista_->analyzeAndCommit(sta_inst, depth);
}

bool Timingdp::isIdealClockNet(odb::dbNet* db_net) const {
    if (!db_net || !sta_) {
        return false;
    }
    sta::dbNetwork* db_network = sta_->getDbNetwork();
    sta::Net* sta_net = db_network->dbToSta(db_net);
    if (!sta_net) {
        return false;
    }
    sta::Network* network = sta_->network();
    sta::PinSet* drivers = network->drivers(sta_net);
    if (drivers && !drivers->empty()) {
        const sta::Pin* drvr_pin = *drivers->begin();
        for (sta::Mode* mode : sta_->modes()) {
            if (!sta_->isIdealClock(drvr_pin, mode)) {
                return false;
            }
        }
    }
    return true;
}


} // namespace tdp

