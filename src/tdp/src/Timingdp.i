// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

// clang-format off
%{

#include "ord/OpenRoad.hh"
#include "tdp/Timingdp.h"
#include "utl/Logger.h"

%}

%include "../../Exception.i"

%import <std_vector.i>
%import "dbtypes.i"

%inline %{

namespace tdp {

void
timing_detailed_placement_cmd(int max_displacment_x,
                       int max_displacment_y,
                       int m,
                       const char* report_file_name){
  tdp::Timingdp *timingdp = ord::OpenRoad::openRoad()->getTimingdp();
  timingdp->detailedPlacement(max_displacment_x, max_displacment_y, std::string(report_file_name), m);
}

void
report_legalization_stats()
{
  tdp::Timingdp *timingdp = ord::OpenRoad::openRoad()->getTimingdp();
  timingdp->reportLegalizationStats();
}

void
create_new_row()
{
  tdp::Timingdp *timingdp = ord::OpenRoad::openRoad()->getTimingdp();
  timingdp->createnewrow();
}

void
design_set()
{
  tdp::Timingdp *timingdp = ord::OpenRoad::openRoad()->getTimingdp();
  timingdp->designset();
}

// void
// check_placement_cmd(bool verbose, const char* report_file_name)
// {
//   dpl::Opendp *opendp = ord::OpenRoad::openRoad()->getOpendp();
//   opendp->checkPlacement(verbose, std::string(report_file_name));
// }

// void
// set_debug_cmd(float min_displacement,
//               const odb::dbInst* debug_instance,
//               int jump_moves,
//               bool iterative_placement)
// {
//   dpl::Opendp* opendp = ord::OpenRoad::openRoad()->getOpendp();
//   opendp->setJumpMoves(jump_moves);
//   opendp->setIterativePlacement(iterative_placement);
//   if (dpl::Graphics::guiActive()) {
//       std::unique_ptr<DplObserver> graphics = std::make_unique<dpl::Graphics>(
//           opendp, min_displacement, debug_instance);
//       opendp->setDebug(graphics);
//   }
// }

} // namespace

%} // inline
