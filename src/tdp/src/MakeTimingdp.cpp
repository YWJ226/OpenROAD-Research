// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "tdp/MakeTimingdp.h"

#include "tcl.h"
#include "utl/decode.h"

extern "C" {
extern int Tdp_Init(Tcl_Interp* interp);
}

namespace tdp {

// Tcl files encoded into strings.
extern const char* tdp_tcl_inits[];

void initTimingdp(Tcl_Interp* tcl_interp)
{
  // Define swig TCL commands.
  Tdp_Init(tcl_interp);
  // Eval encoded sta TCL sources.
  utl::evalTclInit(tcl_interp, tdp::tdp_tcl_inits);
}

}  // namespace tdp
