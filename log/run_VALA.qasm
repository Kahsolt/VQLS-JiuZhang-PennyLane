OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
ry(-0.047208991809220814) q[0];
ry(0.7808882554867012) q[1];
cx q[0],q[1];
ry(0.509390326201012) q[0];
