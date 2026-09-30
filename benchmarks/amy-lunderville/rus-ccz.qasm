/*
 * Figure 2b of Amy & Lunderville (POPL 2025): the repeat-until-success loop of
 * rus.qasm with each Toffoli written as a CCZ conjugated by Hadamards on the target.
 * Reconstructed from the figure; not part of the Feynman artifact.
 */
include "stdgates.inc";

gate ccz p, q, r {
    t p;
    t q;
    t r;
    cx p, q;
    cx q, r;
    cx r, p;
    tdg p;
    tdg q;
    t r;
    cx q, p;
    tdg p;
    cx q, r;
    cx r, p;
    cx p, q;
}

qubit psi;
qubit[2] anc;
bit[2] flags = "11";

reset psi;
h psi;
t psi;

while(int[2](flags) != 0) {
  reset anc[0];
  reset anc[1];
  h anc[0];
  h anc[1];
  h psi;
  ccz anc[0], anc[1], psi;
  h psi;
  s psi;
  h psi;
  ccz anc[0], anc[1], psi;
  h psi;
  z psi;
  h anc[0];
  h anc[1];
  measure anc[0:1] -> flags[0:1];
}

tdg psi;
h psi;
