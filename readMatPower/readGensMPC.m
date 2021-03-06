function [nGen, Vnom, Vmin, Vmax, V0, Pd, Qd, Pav, Sinj, A, B, C, D, PF] = readGensMPC(testCase, nBuses);
% Buses
V0 = testCase.bus(1,8);
Vnom = testCase.bus(2:end,8);
Vmin = testCase.bus(2:end,13);
Vmax = testCase.bus(2:end,12);
% Generators
nGen = size(testCase.gen,1); % number of generators
baseMVA = testCase.baseMVA; % power rating
% Demand 
Pd = testCase.bus(2:end,3)/baseMVA; % active power w/o SLACK BUS
Qd = testCase.bus(2:end,4)/baseMVA; % reactive power w/o SLACK BUS
%% Defining Pav and Sinj
Pav = zeros(nBuses-1,1); % solar PV vector
Sinj = zeros(nBuses-1,1); % inverter capacity vector 
%% Active and Reactive Power Costs 
a = 5; b = 10; c = 2; d = 4;
% Active Power Quadratic costs is defined as follows
A = 0.5 * eye(nBuses-1) * a; % x*A*x' 
B = b * ones(nBuses-1,1);
C = 0.5 * eye(nBuses-1) * c;
D = d * ones(nBuses-1,1);
%%% Define parameters 
PF = 0.8; % selected power factor
end