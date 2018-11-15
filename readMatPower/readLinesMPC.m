function [YBus, ZBus, Ysc, Aa, Ymn, Imax, nBuses, nLines, nB] = readLinesMPC(testCase)
%% mpc data
% Lines
nBuses = size(testCase.bus,1); % number of buses 
nLines = size(testCase.branch,1); % number of lines
baseMVA = testCase.baseMVA; % power rating
% Admittance matrix
YBusSlack = makeYbus(testCase); % Admittance matrix
YBus = YBusSlack(2:end, 2:end); % Admittance matrix: NO SLACK BUS
ZBus = inv(YBus); % Impedance matrix: NO SLACK BUS
Ymn = 2*[real(conj(YBusSlack)) zeros(nBuses,nBuses); zeros(nBuses,nBuses) real(conj(YBusSlack))];
R = testCase.branch(:,3); % resistance
X = testCase.branch(:,4); % reactance
Ysc = 1./ (R - 1j*X);
nB = size(YBus,1); % w/o SLACK BUS
%%
Imax = testCase.branch(:,6)/baseMVA; % current limits
i2e = testCase.bus(:, 1);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(testCase.bus, 1))';
Cf = sparse(1:nLines, e2i(testCase.branch(:,1)), testCase.branch(:,11), nLines, nBuses);
Ct = sparse(1:nLines, e2i(testCase.branch(:,2)), testCase.branch(:,11), nLines, nBuses);
Aa = spdiags(ones(nBuses,1), 0, nLines, nLines)* Cf - Ct; 
%figure; imagesc (abs (inv (YBus))) %% -> matrix visualisastion
end