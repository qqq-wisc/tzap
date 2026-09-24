OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
ccz q[0],q[1],q[2];
ccz q[2],q[0],q[1];
