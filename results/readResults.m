%clear all;
%load('PcQcTable_noCont.mat');
load('PcQcTable_OIDnew.mat');
load('PcQcTable_VARnew.mat');
%load('OtherRes_noCont.mat');
load('OtherRes_OIDnew.mat');
load('OtherRes_VARnew.mat');
%load('VoltageTable_noCont.mat');
load('VoltageTable_OIDnew.mat');
load('VoltageTable_VARnew.mat');
load('Penetration_OIDnew.mat');
load('Penetration_VARnew.mat');
%% Penetration graph OID
figure(1);semilogy([Penetration_OIDnew.store_Gug_Penet],store_Gug_PVGen)
hold on
semilogy(Penetration_OIDnew.store_Gug_Penet,sum(PcQcTAble_OIDnew.store_Gug_Pc,1))
legend({'PV_{av}','PV_{c}'}, 'Location', 'northwest')
ylim([1 10^4])
xlabel('Penetration')
ylabel('Active Power [kW]')
title('PV hosting capacity chart')
grid on
set(gcf,'color','w'); 
%% Penetration graph VOlt/VAR
figure(2);semilogy([Penetration_VARnew.store_VAR_Penet],Penetration_VARnew.store_VAR_PVGen)
hold on
semilogy(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Pc,1))
legend({'PV_{av}','PV_{c}'}, 'Location', 'northwest')
ylim([1 10^4])
xlabel('Penetration')
ylabel('Active Power [kW]')
title('PV hosting capacity chart')
grid on
set(gcf,'color','w'); 
%% OID: Qc and PC
figure(3)
plot(Penetration_OIDnew.store_Gug_Penet,sum(PcQcTAble_OIDnew.store_Gug_Qc,1), 'k')
hold on
plot(Penetration_OIDnew.store_Gug_Penet,sum(PcQcTAble_OIDnew.store_Gug_Pc,1),'b')
plot(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Qc,1),'r')
plot(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Pc,1),'r--')
legend({'OID: Q_{c}','OID: P_{c}','V/V: Q_{c}','V/V: P_{c}'}, 'Location', 'northwest')
%ylim([1 10^4])
xlabel('Penetration')
ylabel('Active/Reactive power')
title('Inverter operation chart')
grid on
set(gcf,'color','w'); 
% %% Volt/VAR: Qc and PC
% figure(4)
% plot(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Qc,1))
% hold on
% plot(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Pc,1))
% legend({'Q_{c}','PV_{c}'}, 'Location', 'northwest')
% %ylim([1 10^4])
% xlabel('Penetration')
% ylabel('Active/Reactive power')
% title('Inverter operation chart with Volt/VAR')
% grid on
% set(gcf,'color','w'); 
%% SQUARE: VOlt/VAR
figure(5);semilogy([Penetration_VARnew.store_VAR_Penet],Penetration_VARnew.store_VAR_PVGen)
hold on
semilogy(Penetration_VARnew.store_VAR_Penet,sum(PcQcTAble_VARnew.store_VAR_Pc,1))
legend({'PV_{av}','PV_{c}'}, 'Location', 'northwest')
ylim([1 10^4])
xlabel('Penetration')
ylabel('Active Power [kW]')
title('PV hosting capacity chart')
grid on
set(gcf,'color','w'); 
%% Pg 
figure(6);fig61 = plot(1:48, [PcQcTAble_OIDnew.store_Gug_Pg],'r');
hold on
fig62 = plot(1:48, [PcQcTAble_VARnew.store_VAR_Pg],'b--');
legend([fig61(1) fig62(1)],'P_g OID','P_g Volt/VAR', 'Location', 'southwest')
%ylim([1 10^4])
xlabel('Time')
ylabel('Active Power P_g [?]')
title('Grid Purchases')
grid on
set(gcf,'color','w'); 
% Qg 
figure(7); fig71 = plot(1:48, [PcQcTAble_OIDnew.store_Gug_Qg], 'r');
hold on
fig72 = plot(1:48, PcQcTAble_VARnew.store_VAR_Qg, 'b--');
legend([fig71(1) fig72(1)],'Q_g OID','Q_g Volt/VAR', 'Location', 'southwest')
%ylim([1 10^4])
xlabel('Time')
ylabel('Reactive Power [?]')
title('Reactive power from the grid')
grid on
set(gcf,'color','w');
%%
figure(8); fig81 = plot(1:48, [VoltageTable_OIDnew.store_Gug_V18], 'r');
hold on
fig82 = plot(1:48, VoltageTable_VARnew.store_VAR_V18, 'b--');
legend([fig81(1) fig82(1)],'V@H18 OID','V@H18 Volt/VAR', 'Location', 'southwest')
%ylim([1 10^4])
xlabel('Time')
ylabel('Voltage [p.u.]')
title('Voltage magnitude at different penetration levels @ Household 18')
grid on
set(gcf,'color','w');
%
figure(9); fig91 = plot(1:48, [VoltageTable_OIDnew.store_Gug_V9], 'r');
hold on
fig92 = plot(1:48, VoltageTable_VARnew.store_VAR_V9, 'b--');
legend([fig91(1) fig92(1)],'V@H9 OID','V@H9 Volt/VAR', 'Location', 'southwest')
%ylim([1 10^4])
xlabel('Time')
ylabel('Voltage [p.u.]')
title('Voltage magnitude at different penetration levels @ Household 9')
grid on
set(gcf,'color','w');
%
figure(10); fig101 = plot(1:48, [VoltageTable_OIDnew.store_Gug_V3], 'r');
hold on
fig102 = plot(1:48, VoltageTable_VARnew.store_VAR_V3, 'b--');
legend([fig101(1) fig102(1)],'V@H3 OID','V@H3 Volt/VAR', 'Location', 'southwest')
%ylim([1 10^4])
xlabel('Time')
ylabel('Voltage [p.u.]')
title('Voltage magnitude at different penetration levels @ Household 3')
grid on
set(gcf,'color','w');
%% ================================
solarCap = [1 2 2.5]; % capacity increase
testCase = IEEE_18BUS_PV; 
% Load Data
[solar, loadHH, loadTotal, timestamp, T, penetration, nBuses] = dataInput(testCase, solarCap(3));
% Reactive Power
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
% Total load versus total solar PV output
figure(37);plot([loadTotal sum(testCase.gen(2:end,9)).*solar])
legend({'load', 'solar'})
xlabel('Time')
ylabel('Active Power [kW]')
title('Load and solar PV output')
grid on
set(gcf,'color','w');  
%% Costs
tariff = 0.3 % $/kWh

