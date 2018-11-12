function [solar, loadHH, loadTotal, timestamp, T, penetration, nBuses] = dataInput(testCase, solarCap)
%% inputs
load('loadProfile.mat', 'all30min'); load('solarProfile.mat', 'solarHH'); % load input files
nBuses = size(testCase.bus,1)-1;
timestamp = table2array(all30min(1:48,1));
solar = table2array(solarHH); % normalised solar output
solar = solar*solarCap; % multiply by capacity increase coeficient
T = length(solar); % number of periods
if nBuses < 25 
    loadtable = all30min(1:48,[2 9:9+nBuses-2]);
    loadHH = table2array(loadtable)*2; % household load
    loadTotal = sum(loadHH,2); % aggregated load
else
    loadHH = testCase.bus(:,3); % household load
    loadTotal = sum(loadHH,2); % aggregated loads
end
% Peak load and peak solar 
loadP = max(sum(loadHH,1));
solarP = sum(testCase.gen(2:end,9));
penetration = solarP*solarCap/loadP;
nBuses = size(testCase.bus,1);
end