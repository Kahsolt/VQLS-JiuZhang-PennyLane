OPENQASM 2.0;
include "qelib1.inc";
qreg q0[2];
ry(-0.04720265722146309) q0[0];
ry(0.7808892648873235) q0[1];
cx q0[0],q0[1];
ry(0.5093863183720745) q0[0];
