function [YBus, ZBus, Ymn, nBuses, nLines, nB] = readLinesMPC(testCase);
%% mpc data
% Lines
nBuses = size(testCase.bus,1); % number of buses 
nLines = size(testCase.branch,1); % number of lines
% Admittance matrix
YBusSlack = makeYbus(testCase); % Admittance matrix
YBus = YBusSlack(2:end, 2:end); % Admittance matrix: NO SLACK BUS
ZBus = inv(YBus); % Impedance matrix: NO SLACK BUS
Ymn = 2*[real(conj(YBusSlack)) zeros(nBuses,nBuses); zeros(nBuses,nBuses) real(conj(YBusSlack))];
nB = size(YBus,1); % w/o SLACK BUS
%figure; imagesc (abs (inv (YBus))) %% -> matrix visualisastion
end