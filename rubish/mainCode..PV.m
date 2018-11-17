clear; clc; close all; 
addpath(genpath('matpower6.0')) 
addpath(genpath('IEEECases')); addpath(genpath('readMatPower'));
addpath(genpath('C:\Users\plus0002\OneDrive\PhD\DTUproject\cvx')) 
% checkcode mainCode -cyc % -> Applies McCabe Complexity (should keep below 10)
mpc = IEEE_9BUS; 
%% AC OPF
loadScen = 1;% 1.5 2.0 5.0]; 
ACOPF_V = complex(zeros(size(mpc.bus,1),length(loadScen))); % Bus Voltage Vector 
ACOPF_I2R = complex(zeros(size(mpc.bus,1)-1,length(loadScen))); % Line Loss Vector
ACOPF_I = complex(zeros(size(mpc.bus,1)-1,length(loadScen))); % Line Current
changePd = 3; % No of Bus with demand change
for l = 1 : length(loadScen) 
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); 
    mpc.bus(changePd,3) = mpc.bus(changePd,3)*loadScen(l); 
    ACOPF = runopf(mpc,mpopt); % 
    ACOPF_V(:,l) = ACOPF.bus(:,8); 
    ACOPF_I2R(:,l) = get_losses(ACOPF); 
    ACOPF_I(:,l) = linesCurrent(mpc, ACOPF_V); 
end
%% Defining Network Topology
[genMatrix,nGen, genLoc, PMin, PMax, QMin, QMax, nBuses, busLoc, ~, ~, Pd, Qd] = readGensMPC(mpc);
[linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Ybus, Yline, lineMaxFlow, OriginBusLoc] = readLinesMPC(mpc);
% adds bus names if given in MatPower file 
if exist('mpc.bus_name','var') == 1 
    BusName = mpc.bus_name; 
end 
[ZBus, YBus, ZBusSlack, YBusSlack, Ymn] = readYBusMPC(linesFrom, linesTo, Z, Ybus, nLines, nBuses); 
%% Active and Reactive Power Costs 
a = 5; b = 10; c = 2; d = 4;
% Active Power Quadratic costs is defined as follows
A = eye(nBuses-1) * a; % x*A*x' 
B = b * ones(nBuses-1,1);
C = eye(nBuses-1) * c;
D = d * ones(nBuses-1,1);
%%% Define parameters 
PF = 0.8; % selected power factor
V0 = 1; % slack bus voltage
Vnom = ones(nBuses-1,1)*1; % nominal voltage vector [p.u.]
Vmin = ones(nBuses-1,1)*0.95; % min voltage [p.u.]
Vmax = ones(nBuses-1,1)*1.05; % max voltage [p.u.]
Pd = Pd(2:end); % w/o slack bus
Qd = Qd(2:end); % w/o slack bus 
Nnode = size(YBus,1);
%% Optimisation problem
voltageVec = []; % creating empty vector for results
lossVec = []; % lone loss vector
lossVec_old = []; % lone loss vector
Sinj = zeros(nBuses-1,1); % inverter capacity vector 
Pav = zeros(nBuses-1,1); % solar PV vector
%Pav(mpc.gen(2:end,1)-1) = mpc.gen(2:end,10); % solar PV output (minGen from MatPower [kW]
%Pav = [ 0 0 0 0.05]';
%Sinj = [ 0 0 0 0.05]';
%Sinj(mpc.gen(2:end,1)-1) = mpc.gen(2:end,10)*1.1; % inverter capacity vector 
%changePd = 5; % between 2-10
%Sinj(changePd-1) = Pd(changePd-1)*1.1; % solar PV output [MW]
%Pav(changePd-1) = Pd(changePd-1)*1; % solar PV output [MW]
%loadCoef = [1 1.5 2 2.5];
%solarCoef = [1 2.5 5 10];
for k = 1% : length(loadCoef)
    %Pd(changePd) = Pd(changePd)*loadCoef(k);
    %Qd(changePd) = Qd(changePd)*loadCoef(k);
    %Pav(changePd-1) = Pav(changePd-1)*solarCoef(k);
    %Sinj(changePd-1) = Sinj(changePd-1)*solarCoef(k);
% START ====================================================
cvx_clear % clears any model being constructed
cvx_begin 
    variable V(nBuses-1,1) complex; % Voltage vector [p.u.]
    variable Pc(nBuses-1,1); % Curtailed PV vector [MW]
    variable Qc(nBuses-1,1); % Reative power absorbed/produced by inverter [MW] 
    % OBJECTIVE FUNCTION ====================================
    %Obj1 = 0;
    Obj1_Vec = [real(V0); real(V); imag(V0); imag(V)]; 
    Obj1 = 0.1 * quad_form(Obj1_Vec, Ymn); %Eq.1
    %Obj2 = 0;
    Obj2 = 100 * quad_form(Pc, A) + quad_form(Qc, C) + B' * Pc + D' * abs(Qc); %Eq.2
    %Obj3 = 0;
    Obj3_Vec = eye(Nnode) - (1/Nnode).*ones(Nnode,1)*ones(1,Nnode);
    Obj3 = norm(Obj3_Vec*diag(V),2) 
    minimize(Obj1 + Obj2 + Obj3) % Eq.4
    % CONSTRAINTS ===========================================
       subject to 
       % Solar PV constraints
       0 <= Pc <= Pav; %Eq.9
       (Qc).^2 <= (Sinj).^2-(Pav-Pc).^2; %Eq.10
       abs(Qc) <= tan(acos(PF))*(Pav-Pc); %Eq.11
       % Voltage
       real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); %Eq.5
       imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); %Eq.5
       %real(V(1)) == Vnom + real(ZBus(1,:))*(Pav - Pc - Pd) + imag(ZBus(1,:))*(Qc - Qd); %Eq.5
       %real(V(2:end)) == real (V(1:end-1)) + real(ZBus(2:end,:))*(Pav - Pc - Pd) + imag(ZBus(2:end,:))*(Qc - Qd); %Eq.5
%       real(V(2:end)) == Vnom(2:end) + real(ZBus(2:end,:))*(Pav - Pc - Pd) + imag(ZBus(2:end,:))*(Qc - Qd); %Eq.5
       %imag(V(1)) == imag(ZBus(1,:))*(Pav - Pc - Pd) - real(ZBus(1,:))*(Qc - Qd); %Eq.6
       %imag(V(2:end)) == imag (V(1:end-1)) + imag(ZBus(2:end,:))*(Pav - Pc - Pd) - real(ZBus(2:end,:))*(Qc - Qd); %Eq.6
%       imag(V(2:end)) == imag(ZBus(2:end,:))*(Pav - Pc - Pd) - real(ZBus(2:end,:))*(Qc - Qd); %Eq.6
       % Min & Max voltage magnitude limit 
       Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; % Eq.7 & 8
cvx_end
% ==========================================================
%Pd(changePd) = mpc.bus(changePd,3)/100;
% Voltage vector
Vfull = [V0; V];
voltageVec = [voltageVec Vfull];
% Line Losses
[lineLoss, Vdrop] = linesLoss(mpc, Vfull);
lossVec = [lossVec, lineLoss];
% Current 

% intermLos = [];
% %intermCurrent = []
% for m = 1 : nBuses-1
%    lineLoss_old = real((conj(YBusSlack(m,m+1))))*((real(Vfull(m))+real(Vfull(m+1)))^2+(imag(Vfull(m))+imag(Vfull(m+1)))^2); 
%    %lineCurrent = YBusSlack(m,m+1)*Vfull(m)
%    % Need to double check Loss implementing I^2 Z 
%    intermLos = [intermLos lineLoss_old];
%    %intermCurrent = [intermCurrent lineCurrent];
% end
% intermLos = intermLos*mpc.Vbase^2/mpc.Sbase; % Sbase =10e6 [VA] and Vbase= 11,000 [Volts]
% intermLos = round(intermLos'/1000,4); % from kW to MW
% lossVec_old = [lossVec_old intermLos]; 
end
%% Results
if exist('mpc.bus_name','var') == 1
    VoltageTable = table(voltageVec, 'rownames', BusName) % compare Guggilam and ACOPF
    VoltageTable.Properties.VariableNames = {'V Guggilam' };
    %writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',true);
    %
    LineLossTable = table(lossVec, ACOPF_Losses, 'rownames', BusName(2:end)) % compare Guggilam and ACOPF
    LineLossTable.Properties.VariableNames = {'Line losses Gug' };
    %writetable(LineLossTable, 'LineLossTable', 'WriteRowNames',true);
    %
    ParTable = table(Pav, round(Pc,4), Pav-round(Pc,4), round(Pd,4), round(Qc,4), round(Qd,4), round(Sinj,4), 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
    ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Pd', 'Qc','Qd','Sinj'};
    %writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',true);
else 
    VoltageTable = table(voltageVec) % compare Guggilam and ACOPF
    VoltageTable.Properties.VariableNames = {'Guggilam'};
    %writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',false);
    %
    LineLossTable = table(lossVec) % compare Guggilam and ACOPF
    LineLossTable.Properties.VariableNames = {'Guggilam'};
    %writetable(LineLossTable, 'LineLossTable', 'WriteRowNames',true);
    %
    ParTable = table(Pav, round(Pc,4), Pav-round(Pc,4), round(Pd,4), round(Qc,4), round(Qd,4),round(Sinj,4)); % compare Guggilam and ACOPF
    ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Pd', 'Qc','Qd','Sinj'};
    %writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',false);
end
%% plot voltages 
figure(1) 
plot(1:nBuses,voltageVec) 
%hold on 
%plot([3 6 9 12 15 18],voltageVec([3 6 9 12 15 18]),'b*')  % poles
hold on 
plot(1:nBuses,ACOPF_V,'--')
%hold on
%plot([3 6 9 12 15 18],ACOPFVoltVec([3 6 9 12 15 18]),'r*') % poles
xlabel('bus') 
ylabel('Voltage [p.u.]') 
xlim([1 nBuses]); ylim([0.93 1.05]) 
xticks(0:1:nBuses); 
legend('Guggilam','ACOPF')
%legend('P=1','P=1.5','P=2','P=2.5','AC=1','AC=1.5','AC=2','AC=2.5')
set(gcf,'color','w'); 
grid on
grid minor
% plot line losses
figure(12)
plot(1:nBuses-1,lossVec)
hold on
plot(1:nBuses-1,ACOPF_V,'--')
xlabel('Lines') 
ylabel('Line losses [MW]') 
%xlim([1 nBuses]); ylim([0.9 1.05]) 
xticks(1:1:nBuses); 
legend('P=1','P=1.5','P=2','P=2.5','AC=1','AC=1.5','AC=2','AC=2.5')
set(gcf,'color','w'); 
grid on
% plot line current 
%% plot P&Q
% figure(2)
% xCenter = 0;
% yCenter = 0;
% theta = 0 : 0.01 : 2*pi;
% radius = 1;
% Sx = radius * cos(theta) + xCenter;
% Sy = radius * sin(theta) + yCenter;
% M = [imag(V(1)) imag(V(end)); ...
%      real(V(1)) real(V(end))];
% M1 = [-0.6 0.6; ...
%      0.8 0.8];
% plotv(M,'-')
% hold on 
% plotv(M1,'--')
% hold on
% plot(Sx, Sy,'b','LineWidth',2);
% xlim([-0.9 0.9]); ylim([0 1]) 
% xlabel('Reactive Power') 
% ylabel('Active Power') 
% legend('Grid','Last Bus')
% set(gcf,'color','w'); 