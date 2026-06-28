// Simple multi-stage combinational fixture for lgista tests (Nangate45).
//
//   in -> u1(INV) -> n1 -> u2(BUF) -> n2 -+-> u3(INV) -> n3 -+
//                                         |                  |-> u5(NAND2) -> out
//                                         +-> u4(INV) -> n4 -+
//
// Topological levels: u1 < u2 < {u3,u4} < u5.
// Has a fanout (n2 -> u3,u4) and a reconvergence (n3,n4 -> u5), so write-back
// in topological order exercises cross-stage propagation.
module chain (in, out);
  input in;
  output out;
  wire n1, n2, n3, n4;
  INV_X1   u1 (.A(in), .ZN(n1));
  BUF_X1   u2 (.A(n1), .Z(n2));
  INV_X1   u3 (.A(n2), .ZN(n3));
  INV_X1   u4 (.A(n2), .ZN(n4));
  NAND2_X1 u5 (.A1(n3), .A2(n4), .ZN(out));
endmodule
