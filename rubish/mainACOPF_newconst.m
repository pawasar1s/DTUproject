clear; clc; close all; 
% mpver % - MAtpower version and installed solvers
load('loadProfile.mat', 'all30min'); load('solarProfile.mat') % load input files
addpath(genpath('matpower7.0b1')); addpath(genpath('IEEECases')); addpath(genpath('gurobi')); 
addpath(genpath('cvx')); addpath(genpath('readMatPower')); 
rmpath('cvx/lib/narginchk_') % remove this function to avoid a potential name conflict
% inputs
loadtable = all30min(1:48,[2 9:19]);
timestamp = table2array(all30min(1:48,1));
loadHH = table2array(loadtable)*2;
solarHH = table2array(solarHH);
loadTotal = sum(loadHH,2); % Total power demand 
T = length(solarHH); % number of periods
% plot graphs ???
plotting = 0; % YES = 1, NO = 0
% ACOPF
mpc = IEEE_18BUS_PV; % AC OPF
% case_info(mpc) % info about the Matpower case file
%mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); 
mpopt = mpoption('model','AC', 'pf.alg', 'ISUM', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT');
[genMatrix,nGen, genLoc, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = readGensMPC(mpc);
[linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Ybus, Yline, lineMaxFlow, OriginBusLoc] = readLinesMPC(mpc);
%%
om = opf_model(mpc);
mpc = toggle_OID(mpc, 'ON')
%%
I = speye(nGen); %% identity matrix
Ar = [I I];
PF = tan(acos(0.8));
QinjMax = PF*mpc.gen(:,3);
QabsMax = -PF*mpc.gen(:,3);
%Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;
%% add them to the model
%om.add_var('R', ng, [], Rmin, Rmax);
om.add_lin_constraint('QinjLim', Ar, QabsMax, QinjMax, {'Qg'});

%%
g_fcn = @(x)power_flow_fcn(x, mpc, Pf, Qf, Pt, Qt);
hess_fcn = @(x, lam)power_flow_hess(x, lam, Pf, Qf, Pt, Qt);
om.add_nln_constraint('InvLimit', nBuses, 0, g_fcn, hess_fcn);
% 0 - for inequality constraints, 1 - for equality. 
%%
%om.add_lin_constraint('PF_low', speye(nGen), 1, -1,  {'Pg', 'Qg'})

%%
ACOPF_V = complex(zeros(T,size(mpc.bus,1))); % Bus Voltage Vector 
ACOPF_I2R = complex(zeros(T,size(mpc.bus,1)-1)); % Line Loss Vector
ACOPF_f = zeros(T,1); % Objective Function Vector 
%ACOPF_I = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Current
ACOPF_PgQg = complex(zeros(T,2));
idxHH = find(mpc.bus(:,2) == 2);
for t = 1 : T
    mpc = IEEE_18BUS_PV; % AC OPF
    mpc.gen(2:end,10) = mpc.gen(2:end,10)*solarHH(t);
    mpc.gen(2:end,9) = mpc.gen(2:end,9)*solarHH(t);
    mpc.gen(2:end,4) = mpc.gen(2:end,4)*solarHH(t);
    mpc.gen(2:end,5) = mpc.gen(2:end,5)*solarHH(t);
    mpc.bus(idxHH,3) = loadHH(t,:);
    mpc.bus(idxHH,4) = loadHH(t,:)*0.6; 
    %loadScen = 1;%[1 1.5 2.0 2.5]; 
    %changePd = 6; % No of Bus with demand change
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); 
    %mpc.bus(changePd,3) = mpc.bus(changePd,3)*loadScen(l); 
    ACOPF = runopf(mpc,mpopt); % 
    ACOPF_V(t,:) = ACOPF.bus(:,8); 
    ACOPF_I2R(t,:) = get_losses(ACOPF); 
    ACOPF_f(t) = ACOPF.f;
    ACOPF_PgQg(t,:) = ACOPF.gen(1,[2 3]);
    %ACOPF_I(:,l) = linesCurrent(mpc, ACOPF_V); 
    clear mpc 
end
% ACOPF.om % shows the list of variables and constraints 
%% Branch Data
mpc = IEEE_18BUS_PV; % AC OPF
[genMatrix,nGen, genLoc, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = readGensMPC(mpc);
[linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Ybus, Yline, lineMaxFlow, OriginBusLoc] = readLinesMPC(mpc);
% Network Vizualisation
%network2 = networkViz(mpc,linesFrom, linesTo, Z);
%% save voltages 
if exist('mpc.bus_name','var') == 0
%     BusName = mpc.bus_name;
%     VoltageTable = table(ACOPF_V, 'rownames', BusName); % compare Guggilam and ACOPF
%     VoltageTable.Properties.VariableNames = {'ACOPF_V'};
%     %writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',true);
%     %
%     LineLossTable = table(ACOPF_I2R, 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
%     LineLossTable.Properties.VariableNames = {'ACOPF_I2R' };
    %writetable(LineLossTable, 'LineLossTable', 'WriteRowNames',true);
    %
    %ParTable = table(Pav, round(Pc,4), Pav-round(Pc,4), round(Pd,4), round(Qc,4), round(Qd,4), round(Sinj,4), 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
    %ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Pd', 'Qc','Qd','Sinj'}
    %writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',true);
end
%% Save a new case file
% savecase('OID.mat', mpc, 1)
% mpc = loadcase(casefile)
%% GRAPHS
if plotting == 1
    % ============  LOAD PROFILES
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
    % Voltage profile for 3 poles over day
    UB = ones(size(solarHH,1),1)*1.05;
    t = 0:minutes(30):hours(23.5);
    tt = 1:48;
    % figure(2) ============  ANIMATION
    % pol0 = animatedline('Marker','*','Color','r');
    % pol1 = animatedline('Color','r');
    % pol2 = animatedline('Color','b');
    % pol3 = animatedline('Color','k');
    % axis([1 T 0.94 1.08])
    % set(gcf,'color','w');
    % for k = 1:length(tt)
    %     addpoints(pol0,tt(k),UB(k));
    %     addpoints(pol1,tt(k),ACOPF_V(k,3));
    %     addpoints(pol2,tt(k),ACOPF_V(k,9));
    %     addpoints(pol3,tt(k),ACOPF_V(k,18));
    %     %xticks(0:0.5:24);
    %     drawnow
    %     grid on
    %     xlabel('Time (Hour)')
    %     title('Voltage profile over day')
    %     pause(.15)
    % end
    figure(3) % ============  VOLTAGE OVER DAY FOR 3 BUSES
    plot(t,ACOPF_V(:,3),'r')
    hold on
    plot(t,ACOPF_V(:,9),'b')
    plot(t,ACOPF_V(:,18),'k')
    plot(t,ones(T,1)*1.05,'Marker','*','Color','r')
    %plot(t,ones(size(solarHH,1),1)*0.95,'Marker','*','Color','r')
    xlabel('Time (Hour)')
    ylabel('Voltage [p.u.]')
    ylim([0.94 1.08])
    %xticks(0:2:size(loadHH,1));
    hold off
    title('Voltage profile over day')
    legend({'Pole 3','Pole 9','Pole 18'},'Location','northwest')
    set(gcf,'color','w');
    grid on
    % ============  FEEDER VOLTAGE PROFILE
    f4 = figure(4);
    movegui(f4,'northwest');
    %plot(1:nBuses,voltageVec)
    %hold on
    %plot([3 6 9 12 15 18],voltageVec([3 6 9 12 15 18]),'b*')  % poles
    %hold on
    plot(1:nBuses,ACOPF_V(26,:),'r--')
    hold on
    plot(1:nBuses,ACOPF_V(40,:),'b--')
    plot([3 6 9 12 15 18],ACOPF_V(26,[3 6 9 12 15 18]),'k*') % poles
    plot([5 7 8 11 19],ACOPF_V(26,[5 7 8 11 19]),'ro') % big solar
    plot([3 6 9 12 15 18],ACOPF_V(40,[3 6 9 12 15 18]),'k*') % poles
    plot([5 7 8 11 19],ACOPF_V(40,[5 7 8 11 19]),'ro') % big solar
    xlabel('bus')
    ylabel('Voltage [p.u.]')
    xlim([1 nBuses]); ylim([0.95 1.07])
    xticks(0:1:nBuses);
    legend('V @ 1.00pm','V @ 8.00pm','Poles','Large Solar','Location','SouthWest')
    title('Voltage Profile')
    set(gcf,'color','w');
    grid on
    grid minor
    % ============  LINE LOSSES
    f5 = figure(5);
    movegui(f5,'southwest');
    %plot(1:size(ACOPF.branch,1),ACOPF_I2R(:,[25 40]),'+')
    %hold on
    plot([1 4 7 10 13 16], ACOPF_I2R([26 40],[1 4 7 10 13 16]),'o')
    hold on
    plot([2 3 5 6 8 9 11 12 14 15 17 18], ACOPF_I2R([26 40],[2 3 5 6 8 9 11 12 14 15 17 18]),'*')
    xlabel('Lines')
    ylabel('Line losses [kW]')
    legend('Pole-to-pole @ 1.00pm','Pole-to-pole @ 8.00pm','Drop Lines @ 12.30pm',...
        'Drop Lines @ 8.0pm','Location','NorthEast')
    title('Line Losses')
    xlim([0 nLines]);% ylim([0.9 1.05])
    %xticks(1:1:size(ACOPF.branch,1));
    set(gcf,'color','w');
    grid on
    grid minor
    % ============  ACTIVE/REACTIVE POWER PURCHASE
    f6 = figure(6);
    movegui(f6,'northeast');
    plot(t, ACOPF_PgQg)
    hold on
    plot(t, loadTotal, 'g')
    plot(t, abs(sum(ACOPF_I2R,2)), 'k','LineWidth',2)
    % yyaxis right
    % plot(t, ACOPF_f,'*')
    hold off
    xlabel('Time [Hours]')
    ylabel('Grid Purchases [kW]')
    legend('Active Grid Power (P)','Reactive Grid Power (Q)','Total Load','Line Losses','Location','SouthWest')
    title('Grid Purchases and Cost')
    grid on
    set(gcf,'color','w');
end
%% 
