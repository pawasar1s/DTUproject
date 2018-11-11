load('PcQcTable_OID.mat');
load('VoltageTable_OID.mat');
load('PcQcTable_noCont.mat');
load('VoltageTable_noCont.mat');
%%
figure(21);bar(PcQcTAble_OID.store_Gug_Pc,'DisplayName','PcQcTAble_OID.store_Gug_Pg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([0 0.05])
title('PV output curtailment - OID')
set(gcf,'color','w');
grid on  
figure(22);bar(PcQcTAble_OID.store_Gug_Qc,'DisplayName','PcQcTAble_OID.store_Gug_Qg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([-0.5 0.05])
title('PV Reactive Power - OID')
set(gcf,'color','w');
grid on  
figure(23);bar(PcQcTAble_noCont.store_Gug_Pc,'DisplayName','PcQcTAble_noCont.store_Gug_Pc')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([0 0.05])
title('PV output curtailment - no control')
set(gcf,'color','w');
grid on  
figure(24);bar(PcQcTAble_noCont.store_Gug_Qc,'DisplayName','PcQcTAble_noCont.store_Gug_Pc')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([-0.5 0.05])
title('PV Reactive Power - no contol')
set(gcf,'color','w');
grid on  
%% grid reactive power
figure(25);bar(PcQcTAble_OID.store_Gug_Qg,'DisplayName','PcQcTAble_noCont.store_Gug_Qg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([-0.5 0.05])
title('Grid Reactive Power - OID')
set(gcf,'color','w');
grid on  
figure(26);bar(PcQcTAble_noCont.store_Gug_Qg,'DisplayName','PcQcTAble_noCont.store_Gug_Qg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
ylim([-0.5 0.05])
title('Grid Reactive Power - no contol')
set(gcf,'color','w');
grid on  
figure(27);bar(PcQcTAble_OID.store_Gug_Pg,'DisplayName','PcQcTAble_noCont.store_Gug_Pg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
%ylim([0 0.05])
title('Grid Active Power - OID')
set(gcf,'color','w');
grid on  
figure(28);bar(PcQcTAble_noCont.store_Gug_Pg,'DisplayName','PcQcTAble_noCont.store_Gug_Pg')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
%ylim([0 0.05])
title('Grid Active Power - no contol')
set(gcf,'color','w');
grid on  
%%
figure(29);plot(1:48, VoltageTable_OID.store_Gug_V18,'DisplayName','VoltageTable_OID.store_Gug_V18')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
%ylim([0 0.05])
title('Grid Active Power - OID')
set(gcf,'color','w');
ylim([0.95 1.1])
grid on  
figure(30);plot(1:48, VoltageTable_noCont.store_Gug_V18,'DisplayName','VoltageTable_noCont.store_Gug_V18')
legend({'SolarCap = 1','SolarCap = 1.5', 'SolarCap = 2'})
%ylim([0 0.05])
title('Grid Active Power - no contol')
ylim([0.95 1.1])
set(gcf,'color','w');
grid on  
%% Voltage duration curve 
ret = duration_plot([VoltageTable_OID.store_Gug_V18(:,3) VoltageTable_noCont.store_Gug_V18(:,3)], 'x_label','% of time', 'y_label','Voltage Mangnitude [p.u.]', 'Legend', {'OID','No Control'}, 'title', 'Voltage Duration Curve')