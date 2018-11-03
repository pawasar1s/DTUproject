clear all; clc; close all; 
addpath(genpath('matpower6.0')); addpath(genpath('IEEECases')); addpath(genpath('gurobi')); 
%% load profiles 
load('workspace_30min_XA10.mat', 'all30min')
load('solarProfile.mat')
loadtable = all30min(1:48,[2 9:19]);
timestamp = table2array(all30min(1:48,1));
loadHH = table2array(loadtable)*2;
solarHH = table2array(solarHH);
%% load profiles
% figure(1)
% plot(loadHH)
% xlabel('Periods') 
% ylabel('Household Loads [kW]') 
% xlim([1 48]);
% set(gcf,'color','w'); 
% grid on
% figure(10)
% plot(sum(loadHH,2))
% xlabel('Periods') 
% ylabel('Household Loads [kW]') 
%% ACOPF
mpc = IEEE_18BUS_PV; % AC OPF
ACOPF_V = complex(zeros(size(mpc.bus,1),size(loadHH,1))); % Bus Voltage Vector 
ACOPF_I2R = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Loss Vector
idxHH = find(mpc.bus(:,2) == 2);
for t = 1:size(loadHH,1)
    mpc = IEEE_18BUS_PV; % AC OPF
    mpc.gen(2:end,10) = mpc.gen(2:end,10)*solarHH(t);
    mpc.gen(2:end,9) = mpc.gen(2:end,9)*solarHH(t);
    mpc.gen(2:end,4) = mpc.gen(2:end,4)*solarHH(t);
    mpc.gen(2:end,5) = mpc.gen(2:end,5)*solarHH(t);
    mpc.bus(idxHH,3) = loadHH(t);
    mpc.bus(idxHH,4) = loadHH(t)*0.6;
%loadScen = 1;%[1 1.5 2.0 2.5]; 
%changePd = 6; % No of Bus with demand change
% B9_ACOPF_V = complex(zeros(size(mpc.bus,1),size(loadHH,1))); % Bus Voltage Vector 
% B9_ACOPF_I2R = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Loss Vector
% B18_ACOPF_V = complex(zeros(size(mpc.bus,1),size(loadHH,1))); % Bus Voltage Vector 
% B18_ACOPF_I2R = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Loss Vector
    ACOPF_I = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Current
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); 
%mpc.bus(changePd,3) = mpc.bus(changePd,3)*loadScen(l); 
    ACOPF = runopf(mpc,mpopt); % 
    ACOPF_V(:,t) = ACOPF.bus(:,8); 
    ACOPF_I2R(:,t) = get_losses(ACOPF); 
% B9_ACOPF_V(:,t) = ACOPF.bus(9,8); 
% B9_ACOPF_I2R(:,t) = get_losses(ACOPF); 
% B18_ACOPF_V(:,t) = ACOPF.bus(18,8); 
% B18_ACOPF_I2R(:,t) = get_losses(ACOPF); 
%ACOPF_I(:,l) = linesCurrent(mpc, ACOPF_V); 
    clear mpc
end
%% 
figure(2)
pol0 = animatedline('Marker','*','Color','r');
pol1 = animatedline('Color','r');
pol2 = animatedline('Color','b');
pol3 = animatedline('Color','k');
axis([1 size(loadHH,1) 0.94 1.08])
t = 0:minutes(30):hours(23.5);
tt = 1:48;
set(gcf,'color','w'); 
UB = ones(size(solarHH,1),1)*1.05;
for k = 1:length(tt)
    addpoints(pol0,tt(k),UB(k));
    addpoints(pol1,tt(k),ACOPF_V(3,k));
    addpoints(pol2,tt(k),ACOPF_V(9,k));
    addpoints(pol3,tt(k),ACOPF_V(18,k));
    %xticks(0:0.5:24); 
    drawnow
    grid on
    xlabel('Time (Hour)')  
    title('Voltage profile over day')
    pause(.15)
end
figure(3)
plot(t,ACOPF_V(3,:),'r')
hold on
plot(t,ACOPF_V(9,:),'b')
plot(t,ACOPF_V(18,:),'k')
plot(t,ones(size(solarHH,1),1)*1.05,'Marker','*','Color','r')
%plot(t,ones(size(solarHH,1),1)*0.95,'Marker','*','Color','r')
xlabel('Time (Hour)') 
ylabel('Voltage [p.u.]') 
%xlim([1 size(loadHH,1)]); 
ylim([0.94 1.08]) 
%xticks(0:2:size(loadHH,1)); 
hold off
title('Voltage profile over day')
legend({'Pole 3','Pole 9','Pole 18'},'Location','northwest')
set(gcf,'color','w'); 
grid on

%% Branch Data
mpc = IEEE_18BUS_PV; % AC OPF
addpath(genpath('cvx')); addpath(genpath('readMatPower')); 
rmpath('cvx/lib/narginchk_')
[genMatrix,nGen, genLoc, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = readGensMPC(mpc);
[linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Ybus, Yline, lineMaxFlow, OriginBusLoc] = readLinesMPC(mpc);
%% Network Vizualisation
%network2 = networkViz(mpc,linesFrom, linesTo, Z);
%% save voltages 
if exist('mpc.bus_name','var') == 0
    BusName = mpc.bus_name;
    VoltageTable = table(ACOPF_V, 'rownames', BusName); % compare Guggilam and ACOPF
    VoltageTable.Properties.VariableNames = {'ACOPF_V'};
    %writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',true);
    %
    LineLossTable = table(ACOPF_I2R, 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
    LineLossTable.Properties.VariableNames = {'ACOPF_I2R' };
    %writetable(LineLossTable, 'LineLossTable', 'WriteRowNames',true);
    %
    %ParTable = table(Pav, round(Pc,4), Pav-round(Pc,4), round(Pd,4), round(Qc,4), round(Qd,4), round(Sinj,4), 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
    %ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Pd', 'Qc','Qd','Sinj'}
    %writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',true);
end
%% plot voltages  
f4 = figure(4);
movegui(f4,'northwest');
%plot(1:nBuses,voltageVec) 
%hold on 
%plot([3 6 9 12 15 18],voltageVec([3 6 9 12 15 18]),'b*')  % poles
%hold on 
plot(1:size(ACOPF.bus,1),ACOPF_V(:,25),'r--')
hold on
plot(1:size(ACOPF.bus,1),ACOPF_V(:,40),'b--')
plot([3 6 9 12 15 18],ACOPF_V([3 6 9 12 15 18],25),'k*') % poles
plot([5 7 8 11 19],ACOPF_V([5 7 8 11 19],25),'ro') % big solar
plot([3 6 9 12 15 18],ACOPF_V([3 6 9 12 15 18],40),'k*') % poles
plot([5 7 8 11 19],ACOPF_V([5 7 8 11 19],40),'ro') % big solar
xlabel('bus') 
ylabel('Voltage [p.u.]') 
xlim([1 size(ACOPF.bus,1)]); ylim([0.95 1.07]) 
xticks(0:1:size(ACOPF.bus,1)); 
legend('V @ 12.30pm','V @ 8.00pm','Poles','Large Solar','Location','SouthWest')
title('Voltage Profile')
set(gcf,'color','w'); 
grid on
grid minor
%% line losses
f5 = figure(5);
movegui(f5,'southwest');
%plot(1:size(ACOPF.branch,1),ACOPF_I2R(:,[25 40]),'+')
%hold on
plot([1 4 7 10 13 16],ACOPF_I2R([1 4 7 10 13 16],[25 40]),'o') 
hold on
plot([2 3 5 6 8 9 11 12 14 15 17 18],ACOPF_I2R([2 3 5 6 8 9 11 12 14 15 17 18],[25 40]),'*') 
xlabel('Lines') 
ylabel('Line losses [kW]') 
legend('Pole-to-pole @ 12.30pm','Pole-to-pole @ 8.00pm','Drop Lines @ 12.30pm',...
    'Drop Lines @ 8.0pm','Location','NorthEast')
title('Line Losses')
xlim([0 nLines]);% ylim([0.9 1.05]) 
%xticks(1:1:size(ACOPF.branch,1)); 
set(gcf,'color','w'); 
grid on
grid minor