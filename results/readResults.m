clear all;
load('PcQcTable_noCont.mat');
load('PcQcTable_GuG.mat');
load('PcQcTable_VAR.mat');
load('OtherRes_noCont.mat');
load('OtherRes_GuG.mat');
load('OtherRes_VAR.mat');
load('VoltageTable_noCont.mat');
load('VoltageTable_GuG.mat');
load('VoltageTable_VAR.mat');
%%
solarCap = [1 2 2.5]; % capacity increase
testCase = IEEE_18BUS_PV; 
% Load Data
[solar, loadHH, loadTotal, timestamp, T, penetration, nBuses] = dataInput(testCase, solarCap(3));
%% Reactive Power
figure(21);bar([PcQcTAble_Gug.store_Gug_Qc(:,3) PcQcTAble_VAR.store_VAR_Qc(:,3) PcQcTAble_noCont.store_Gug_Qc(:,3)],'DisplayName','PcQcTAble_OID.store_Gug_Pg')
legend({'OID', 'Volt/VAR', 'No Control'})
%ylim([0 0.05])
xlabel('Time')
ylabel('Reactive Power [kW]')
title('Inverter Q over the day')
grid on
set(gcf,'color','w');  
%% Active Power
figure(22);bar([PcQcTAble_Gug.store_Gug_Pc(:,3) PcQcTAble_VAR.store_VAR_Pc(:,3) PcQcTAble_noCont.store_Gug_Pc(:,3)],'DisplayName','PcQcTAble_OID.store_Gug_Pg')
legend({'OID', 'Volt/VAR', 'No Control'})
%ylim([0 0.05])
xlabel('Time')
ylabel('Reactive Power [kVAR]')
title('Inverter P curtailment over the day')
grid on
set(gcf,'color','w');  
figure(23);bar([PcQcTAble_Gug.store_Gug_Qg(:,3) PcQcTAble_VAR.store_VAR_Qg(:,3) PcQcTAble_noCont.store_Gug_Qg(:,3)],'DisplayName','PcQcTAble_OID.store_Gug_Pg')
legend({'OID', 'Volt/VAR', 'No Control'})
%%ylim([0 0.05])
xlabel('Time')
ylabel('Reactive Power [kVAR]')
title('Reactive power from grid')
grid on
set(gcf,'color','w');  
% Active Power
figure(24);bar([PcQcTAble_Gug.store_Gug_Pg(:,3) PcQcTAble_VAR.store_VAR_Pg(:,3) PcQcTAble_noCont.store_Gug_Pg(:,3)],'DisplayName','PcQcTAble_OID.store_Gug_Pg')
legend({'OID', 'Volt/VAR', 'No Control'})
%ylim([0 0.05])
xlabel('Time')
ylabel('Avtive Power [kW]')
title('Active Power from grid')
grid on
set(gcf,'color','w'); 
% voltage profile
figure(25);plot(1:48, [VoltageTable_Gug.store_Gug_V18(:,3) VoltageTable_VAR.store_VAR_V18(:,3) VoltageTable_noCont.store_Gug_V18(:,3)],'DisplayName','VoltageTable_OID.store_Gug_Pg')
legend({'OID', 'Volt/VAR', 'No Control'})
%ylim([0 0.05])
xlabel('Time')
ylabel('Voltage Magnitude [p.u.]')
title('Voltage Profile')
grid on
set(gcf,'color','w');  
%% Voltage duration curve 
f32 = figure(32)
ret = duration_plot([VoltageTable_Gug.store_Gug_V18(:) VoltageTable_VAR.store_VAR_V18(:) VoltageTable_noCont.store_Gug_V18(:)], f32, 'x_label','% of time', 'y_label','Voltage Mangnitude [p.u.]', 'Legend', {'OID','Volt/VAR', 'No Control'}, 'title', 'Voltage Duration Curve')
%%
%plot(VoltageTable_VAR.store_VAR_V18,'*'); legend({'PV=1','PV=2','PV=3'})
%% ============  Reative Power vs Voltage Magnitude (Volt/Var)
f33 = figure(33);
movegui(f33,'northeast');        
scatter(VoltageTable_Gug.store_Gug_V18(:), abs(PcQcTAble_Gug.store_Gug_Qc(:)),'b*', 'DisplayName','VoltageTable_OID.store_Gug_V18')
hold on 
scatter(VoltageTable_VAR.store_VAR_V18(:), abs(PcQcTAble_VAR.store_VAR_Qc(:)),'r*', 'DisplayName','VoltageTable_VAR.store_Gug_V18')
hold on
scatter(VoltageTable_noCont.store_Gug_V18(:), abs(PcQcTAble_noCont.store_Gug_Qc(:)),'k*', 'DisplayName','VoltageTable_OID.store_Gug_V18')
xlim([0.93 1.1])
xlabel('Voltage Magnitude [p.u.]')
ylabel('Reactive Power [kVAR]')
title('Volt/VAR control')
legend({'OID', 'Volt/VAR', 'No Control'})
grid on
set(gcf,'color','w');    
% ============  Active Power vs Voltage Magnitude (Volt/Var)
f34 = figure(34);
movegui(f34,'southeast');        
scatter(VoltageTable_Gug.store_Gug_V18(:), abs(PcQcTAble_Gug.store_Gug_Pc(:)),'b*', 'DisplayName','VoltageTable_OID.store_Gug_V18')
hold on 
scatter(VoltageTable_VAR.store_VAR_V18(:), abs(PcQcTAble_VAR.store_VAR_Pc(:)),'r*', 'DisplayName','VoltageTable_VAR.store_Gug_V18')
hold on
scatter(VoltageTable_noCont.store_Gug_V18(:), abs(PcQcTAble_noCont.store_Gug_Pc(:)),'k*', 'DisplayName','VoltageTable_OID.store_Gug_V18')
xlim([0.93 1.1])
xlabel('Voltage Magnitude [p.u.]')
ylabel('Active Power [kW]')
title('PV output curtailment with Volt/VAR control')
legend({'OID', 'Volt/VAR', 'No Control'})
grid on
set(gcf,'color','w');
%% Line Losses 
figure(35);bar([OtherRes_Gug.store_Gug_I2R(:,3) OtherRes_VAR.store_VAR_I2R(:,3) OtherRes_noCont.store_Gug_I2R(:,3)])
legend({'OID', 'Volt/VAR', 'No Control'})
xlabel('Time')
ylabel('Active Power Losses [kW]')
title('Line Losses')
grid on
set(gcf,'color','w');  
figure(36);bar([OtherRes_Gug.store_Gug_I2R(:,3)./loadTotal*100 OtherRes_VAR.store_VAR_I2R(:,3)./loadTotal*100 OtherRes_noCont.store_Gug_I2R(:,3)./loadTotal*100])
legend({'OID', 'Volt/VAR', 'No Control'})
xlabel('Time')
ylabel('Active Power Losses [%]')
title('Line Losses (% of Load)')
grid on
set(gcf,'color','w');  
figure(37);plot([loadTotal sum(testCase.gen(2:end,9)).*solar])
legend({'load', 'solar'})
xlabel('Time')
ylabel('Active Power [kW]')
title('Load and solar PV output')
grid on
set(gcf,'color','w');  
%% Costs
tariff = 0.3 % $/kWh

