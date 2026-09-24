OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
cz q[0],q[1];
ccx q[0],q[1],q[2];
ccz q[0],q[1],q[2];
