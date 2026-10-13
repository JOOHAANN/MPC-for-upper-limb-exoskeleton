clc;
clear;
close all;

m = 5;

M = 45.8542898388227;

L1 = 0.286467647976371;

L2 = 0.798896309352907;

phi = 1.04710175828713;

I = m * L2^2 + (4 * M * L1^2) / 3;

c = 10;

g = 9.81;

load('Torque_real_1.mat');

t = 0:0.001:0.001 * length(Torque_real_1);

data = [t(1:end - 1)', Torque_real_1];
