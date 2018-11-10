function [solarHH, loadHH, loadTotal, timestamp, T] = dataInput(testCase);
%% inputs
load('loadProfile.mat', 'all30min'); load('solarProfile.mat') % load input files
nBuses = size(testCase.bus,1)-1;
loadtable = all30min(1:48,[2 9:9+nBuses-2]);
timestamp = table2array(all30min(1:48,1));
loadHH = table2array(loadtable)*2; % household load
solarHH = table2array(solarHH); % normalised solar output
loadTotal = sum(loadHH,2); % aggregated load 
T = length(solarHH); % number of periods
end