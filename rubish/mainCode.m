clear; clc; 
addpath(genpath('matpower7.0b1')) 
addpath(genpath('IEEECases')); addpath(genpath('readMatPower'));
addpath(genpath('C:\Users\plus0002\OneDrive\PhD\DTUproject\cvx')) 
% checkcode mainCode -cyc % -> Applies McCabe Complexity (should keep below 10)
%  figure; imagesc (abs (inv (YBus)))
per = 12;
mpc = IEEE_18BUS; 
[solarHH, loadHH, loadTotal, timestamp, T] = dataInput();
idxHH = find(mpc.bus(:,2) == 2);
mpc.bus(idxHH,3) = loadHH(per,:)'; % Pg
mpc.bus(idxHH,4) = loadHH(per,:)'*0.6; % Qg
%% AC OPF
loadScen = [1]; 
changePd = 6; % No of Bus with demand change
ACOPF_V = complex(zeros(size(mpc.bus,1),length(loadScen))); % Bus Voltage Vector 
ACOPF_I2R = complex(zeros(size(mpc.bus,1)-1,length(loadScen))); % Line Loss Vector
ACOPF_I = complex(zeros(size(mpc.bus,1)-1,length(loadScen))); % Line Current
for l = 1 : length(loadScen) 
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); 
    mpc.bus(changePd,3) = mpc.bus(changePd,3)*loadScen(l); 
    ACOPF = runopf(mpc,mpopt); % 
    ACOPF_V(:,l) = ACOPF.bus(:,8); 
    ACOPF_I2R(:,l) = get_losses(ACOPF); 
    %ACOPF_I(:,l) = linesCurrent(mpc, ACOPF_V); 
end
%% Defining Network Topology
[genMatrix,nGen, genLoc, PMin, PMax, QMin, QMax, nBuses, busLoc, ~, ~, Pd, Qd] = readGensMPC_old(mpc);
[linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, Z, Y, Yline, lineMaxFlow, OriginBusLoc] = readLinesMPC_old(mpc);
% adds bus names if given in MatPower file 
if exist('mpc.bus_name','var') == 1 
    BusName = mpc.bus_name; 
end 
%[ZBus, YBus, YBusSlack, Ymn] = readYBusMPC(linesFrom, linesTo, Z, Y, nLines, nBuses); 
YBusSlack = makeYbus(mpc);
YBus = YBusSlack(2:end, 2:end);
ZBus = inv(YBus);
Ymn = 2*[real(conj(YBusSlack)) zeros(nBuses,nBuses); zeros(nBuses,nBuses) real(conj(YBusSlack))];
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
Vmin = ones(nBuses-1,1)*0.90; % min voltage [p.u.]
Vmax = ones(nBuses-1,1)*1.05; % max voltage [p.u.]
Pd = Pd(2:end); % w/o slack bus
Qd = Qd(2:end); % w/o slack bus 
Nnode = size(YBus,1);
%% Optimisation problem - No Solar
Sinj = zeros(nBuses-1,1); % inverter capacity vector 
Pav = zeros(nBuses-1,1); % solar PV vector
idxPV = find(mpc.bus(:,2) == 2); % idx for buses with PV systems
Pav(mpc.gen(2:end,1)-1) = mpc.gen(2:end,9); % solar PV output (minGen from MatPower [kW]
Sinj(mpc.gen(2:end,1)-1) = mpc.gen(2:end,9); % solar PV output (minGen from MatPower [kW]
voltageVec = complex(zeros(size(mpc.bus,1),length(loadScen))); % empty vector for voltage scenarios
lossVec = complex(zeros(size(mpc.bus,1)-1,length(loadScen))); % lone loss vector
changePd = changePd; % between 2-10
loadCoef = [1];
for k = 1 : length(loadCoef)
    Pd(changePd) = Pd(changePd)*loadCoef(k);
    Qd(changePd) = Qd(changePd)*loadCoef(k)/2;
% START ====================================================
cvx_clear % clears any model being constructed
cvx_begin 
    variable V(nBuses-1,1) complex; % Voltage vector [p.u.]
    variable Pc(nBuses-1,1); % Curtailed PV vector [MW]
    variable Qc(nBuses-1,1); % Reative power absorbed/produced by inverter [MW] 
    % OBJECTIVE FUNCTION ====================================
    Obj1_Vec = [real(V0); real(V); imag(V0); imag(V)]; 
    Obj1 = 0.5 * quad_form(Obj1_Vec, Ymn); %Eq.1
    Obj2 = 0.5 * quad_form(Pc, A) + quad_form(Qc, C) + B' * Pc + D' * abs(Qc); %Eq.2
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
voltageVec(:,k) = Vfull;
% Line Losses
[lineLoss, Vdrop] = linesLoss(mpc, Vfull);
lossVec(:,k) = lineLoss;
% Current 
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
xlim([1 nBuses]); ylim([0.90 1.05]) 
xticks(0:1:nBuses); 
%legend('P1','P1.5','P2','P2.5','AC1','AC1.5','AC2','AC2.5')
set(gcf,'color','w'); 
grid on; grid minor
% plot line losses
figure(12)
plot(1:nBuses-1,lossVec)
legend('P=1','P=1.5')%,'P=2','P=2.5','AC=1','AC=1.5','AC=2','AC=2.5')
set(gcf,'color','w'); 
hold on
plot(1:nBuses-1,ACOPF_I2R,'--')
xlabel('Lines') 
ylabel('Line losses [kW]') 
xticks(1:1:nBuses); 
grid on; grid minor
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