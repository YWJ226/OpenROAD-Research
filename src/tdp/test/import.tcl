source helpers.tcl
set test_name tdp01
read_lef  ./tech/asap7_nangate45_2A6M10M.lef
read_lef  ./tech/asap7sc7p5t_28_R_1x_220121a.upper.processed.lef
read_lef  ./tech/NangateOpenCellLibrary.macro.mod.bottom.processed.lef

read_liberty ./tech/NangateOpenCellLibrary_typical.bottom.processed.lib
read_liberty ./tech/asap7sc7p5t_AO_RVT_FF_nldm_211120.upper.processed.lib
read_liberty ./tech/asap7sc7p5t_INVBUF_RVT_FF_nldm_220122.upper.processed.lib
read_liberty ./tech/asap7sc7p5t_OA_RVT_FF_nldm_211120.upper.processed.lib
read_liberty ./tech/asap7sc7p5t_SEQ_RVT_FF_nldm_220123.upper.processed.lib
read_liberty ./tech/asap7sc7p5t_SIMPLE_RVT_FF_nldm_211120.upper.processed.lib


read_def  ./simple/gcd/gcd_3D.lg.def
read_verilog ./simple/gcd/gcd_3D.lg.v
read_sdc ./simple/gcd/3_place.sdc

source ./simple/setRC_2A6M10M.tcl

set_routing_layers -signal M2-M3_add -clock M2-M3_add
estimate_parasitics -placement
report_tns
report_checks -to [all_outputs] -path_delay max -digits 6 

timing_detailed_placement  -max_displacement {500 100} -m 2

# estimate_parasitics -placement
# report_tns
exit
