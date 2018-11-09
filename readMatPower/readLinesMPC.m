function [linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Yline, Y, lineMaxFlow, OriginBusLoc] = readLinesMPC(mpc);
%% convert branch impedances from Ohms to p.u.
Vbase = mpc.bus(1,10)*10^3; % in Volts 11,000
Sbase = mpc.baseMVA * 1e6; % in VA 10e6
% Lines
nBuses = size(mpc.bus,1); % number of buses 
nLines = size(mpc.branch,1); % number of lines
%
OriginBusLoc = mpc.bus(:,1); % number of buses 
OriginLinesFrom = mpc.branch(:,1);
OriginLinesTo = mpc.branch(:,2);
[tf, linesFrom]=ismember(OriginLinesFrom,OriginBusLoc,'rows');
[tf, linesTo]=ismember(OriginLinesTo,OriginBusLoc,'rows');
R = mpc.branch(:,3); % resistance
j = sqrt(-1);
%delta_R = 1e-4; 
%R(R < delta_R) = delta_R; % adding small resistance to every tansformer with zero resistance
X = mpc.branch(:,4); % reactance
Z = R + 1j * X;
Y = 1./Z;
Yline = 1 ./ (R - 1j * X); % line admittance bus
lineMaxFlow = mpc.branch(:,6); % rateA, (MVA long term rating)
% Line incidence matrix for SDP
linesMatFrom = sparse(1:nLines, linesFrom, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
linesMatTo = sparse(1:nLines, linesTo, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
end