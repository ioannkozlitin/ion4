clear; clc; close all;

Raizer;
load xe_res;

xe = flip(flip(xe, 1), 2);

A = abs(xe - xe_Saha);

max_deviation = max(max(A))