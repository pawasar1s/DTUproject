clear; clc; close all; 
addpath(genpath('matpower7.0b1')); addpath(genpath('cvx')); addpath(genpath('SC_Cases'));
addpath(genpath('IEEECases')); addpath(genpath('readMatPower')); addpath(genpath('Results'));
%addpath(genpath('C:\Users\plus0002\OneDrive\PhD\DTUproject\cvx')) 
rmpath('cvx/lib/narginchk_') % remove this function to avoid a potential name conflict
%checkcode Guggilam -cyc %% -> Applies McCabe Complexity (should keep below 10)
%% Settings
% No of time intervals
multiPer = 1; % 1 = multiperiod model, 2 = discretised model
per = 26; % for discretised model, which period to test? 26 = 1pm (30min intervals)
% Matpower settings
matpower = 0; % Matpower 1 - run Matpower for comparison 
% Plotting 
plotting = 0; % YES = 1, NO = 0
% Select test case
testCase = IEEE_18BUS_PV; % select case
% ==============
% to change the scenario with no control, change Vmax to 1.2 in
% IEEE_18BUS_PV.m  
% ===============
% testCase = matpower_LV_semiurban; % works fine 
% testCase = matpower_rural; % works fine 
%% network viz
%network2 = networkViz(testCase);
%% Results at varios penetration levels
solarScen = 1; % if various solar peentration scenarios are on 
solarCap = [1 1.5 2]; % capacity increase
%% 
T = 48;
store_Gug_V18 = complex(zeros(T, length(solarCap))); % Voltage vector 
store_Gug_V9 = complex(zeros(T, length(solarCap))); % Voltage vector 
store_Gug_V3 = complex(zeros(T, length(solarCap))); % Voltage vector 
store_Gug_Pc = complex(zeros(T, length(solarCap))); % PV curtailment
store_Gug_Qc = complex(zeros(T, length(solarCap))); % Q from PV
store_Gug_Qg = complex(zeros(T, length(solarCap))); % Q from grid
store_Gug_Pg = complex(zeros(T, length(solarCap))); % Q from grid
store_Gug_Penet = zeros(length(solarCap), 1); % 
for i = 1 : length(solarCap)
    testCase = IEEE_18BUS_PV; 
    %% Load Data
    [solar, loadHH, loadTotal, timestamp, T, penetration, nBuses] = dataInput(testCase, solarCap(i));
    %% MATPOWER power flow
    if matpower == 1
        ACOPFsolver = 1; % 1 = exact ACOPF, 2 = backward/forward sweep with I summation
        PF = 0.8;
        [ACOPF_struct, ACOPF_V, ACOPF_I2R, ACOPF_f, ACOPF_Qg] = ACOPF(testCase, T, solar, loadHH, ACOPFsolver, multiPer, per, PF);
    end
    %% Guggilam
    [V, Pc, Qc, Gug_V, Gug_I2R, Gug_PgTot, Gug_QgTot, Gug_PcTot, Gug_QcTot] = Guggilam(testCase, T, solar, loadHH, multiPer, per);
    %[V, Pc, Qc, Gug_V, Gug_I2R, Gug_PgTot, Gug_QgTot, Gug_PcTot, Gug_QcTot] = VoltVar(testCase, T, solar, loadHH, multiPer, per);    
    % 
    store_Gug_V18(:,i) = Gug_V(:,18); % Voltage vector Bus 18
    store_Gug_V9(:,i) = Gug_V(:,9); % Voltage vector Bus 18
    store_Gug_V3(:,i) = Gug_V(:,3); % Voltage vector Bus 18
    store_Gug_Pc(:,i) = Gug_PcTot; % P from PV (curtailment)
    store_Gug_Qc(:,i) = Gug_QcTot; % Q from PV (regulation)
    store_Gug_Qg(:,i) = Gug_QgTot; % Q from grid (total)
    store_Gug_Pg(:,i) = Gug_PgTot; % P from PV (curtailment)
    store_Gug_Penet(i) = penetration; % penetration level
end
    %% Save Results
% if exist('testCase.bus_name','var') == 0
%     %BusName = testCase.bus_name; 
%     VoltageTable_noCont = table( store_Gug_V18, store_Gug_V9, store_Gug_V3) % compare Guggilam and ACOPF
%     %VoltageTable.Properties.VariableNames = {'V_Guggilam' };
%     %writetable(VoltageTable_noCont, 'res_VoltageTable_noCont.mat');
%     save('VoltageTable_noCont.mat', 'VoltageTable_noCont')
%         
%     PcQcTAble_noCont = table(store_Gug_Pc, store_Gug_Qc, store_Gug_Pg, store_Gug_Qg) % compare Guggilam and ACOPF
%    % LineLossTable.Properties.VariableNames = {'Line_losses_Gug' };
%     %writetable(PcQcTAble_noCont, 'res_PcQcTable_noCont.mat');
%     save('PcQcTAble_noCont.mat', 'PcQcTAble_noCont')
% end
%% GRAPHS
if plotting == 1
    %nBuses = size(testCase.bus,1); % number of buses
    nLines = size(testCase.branch,1); % number of lines
    if multiPer == 0
        % ============  FEEDER VOLTAGE PROFILE
        f4 = figure(4);
        movegui(f4,'northwest');
        plot(1:nBuses,abs(Gug_V),'r--')
        hold on
        plot(1:nBuses,abs(ACOPF_V),'b--')
        %plot([3 6 9 12 15 18],Gug_V([3 6 9 12 15 18]),'k*')  % poles
        %plot([5 7 8 11 19],Gug_V([5 7 8 11 19]),'ko') % big solar
        %plot([3 6 9 12 15 18],ACOPF_V([3 6 9 12 15 18]),'k*') % poles
        %plot([5 7 8 11 19],ACOPF_V([5 7 8 11 19]),'ko') % big solar
        hold off
        xlabel('bus')
        ylabel('Voltage [p.u.]')
        xlim([1 nBuses]); ylim([0.90 1.07])
        xticks(0:1:nBuses);
        legend('V Gug','V MatPower','Poles','Large Solar','Location','SouthWest')
        title('Voltage Profile')
        set(gcf,'color','w');
        grid on
        %grid minor
        % ============  LINE LOSSES
        f5 = figure(5);
        movegui(f5,'southwest');
        %plot(1:size(ACOPF.branch,1),ACOPF_I2R(:,[25 40]),'+')
        %hold on
        plot([1 4 7 10 13 16], ACOPF_I2R([1 4 7 10 13 16]),'ko-')
        hold on
        plot([2 3 5 6 8 9 11 12 14 15 17 18], ACOPF_I2R([2 3 5 6 8 9 11 12 14 15 17 18]),'k*')
        %
        plot([1 4 7 10 13 16], Gug_I2R([1 4 7 10 13 16]),'ro-')
        hold on
        plot([2 3 5 6 8 9 11 12 14 15 17 18], Gug_I2R([2 3 5 6 8 9 11 12 14 15 17 18]),'r*')
        xlabel('Lines')
        ylabel('Line losses [kW]')
        legend('ACOPF Pole-to-pole','ACOPF Drop Lines',...
            'Gug Pole-to-pole', 'Gug Drop Lines','Location','NorthEast')
        title('Line Losses')
        xlim([0 nLines]);
        %xticks(1:1:size(ACOPF.branch,1));
        set(gcf,'color','w');
        grid on
        %grid minor
        % ============  ACTIVE/REACTIVE POWER PURCHASE
    elseif multiPer == 1
        UB = ones(size(solar,1),1)*1.05;
        timE = 0:minutes(30):hours(23.5);
        tt = 1:48;
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
        
        % ============  VOLTAGE OVER DAY FOR 3 BUSES
        f3 = figure(3) 
        movegui(f3,'southwest');
        plot(timE,ACOPF_V(:,3),'r')
        hold on
        plot(timE,ACOPF_V(:,9),'b')
        plot(timE,ACOPF_V(:,18),'k')
        plot(timE,abs(Gug_V(:,3)),'r--')
        plot(timE,abs(Gug_V(:,9)),'b--')
        plot(timE,abs(Gug_V(:,18)),'k--')
        %plot(timE,ones(T,1)*1.05,'Marker','*','Color','r')
        hold off
        %plot(t,ones(size(solarHH,1),1)*0.95,'Marker','*','Color','r')
        xlabel('Time (Hour)')
        ylabel('Voltage [p.u.]')
        ylim([0.94 1.08])
        %xticks(0:2:size(loadHH,1));
        title('Voltage profile over day')
        legend({'ACOPF: Pole 3','ACOPF: Pole 9','ACOPF: Pole 18','Gug: Pole 3','Gug: Pole 9','Gug: Pole 18'},'Location','northwest')
        set(gcf,'color','w');
        grid on
        
        % ============  FEEDER VOLTAGE PROFILE
        f4 = figure(4);
        movegui(f4,'northwest');
        %plot(1:nBuses,voltageVec)
        %hold on
        %plot([3 6 9 12 15 18],voltageVec([3 6 9 12 15 18]),'b*')  % poles
        %hold on
        plot(1:nBuses,ACOPF_V(26,:),'r')
        hold on
        plot(1:nBuses,ACOPF_V(40,:),'b')
        plot(1:nBuses,abs(Gug_V(26,:)),'r--')
        plot(1:nBuses,abs(Gug_V(40,:)),'b--')
        plot([3 6 9 12 15 18],ACOPF_V(26,[3 6 9 12 15 18]),'k*') % poles
        plot([5 7 8 11 19],ACOPF_V(26,[5 7 8 11 19]),'ro') % big solar
        %plot([3 6 9 12 15 18],ACOPF_V(40,[3 6 9 12 15 18]),'k*') % poles
        %plot([5 7 8 11 19],ACOPF_V(40,[5 7 8 11 19]),'ro') % big solar
        plot([3 6 9 12 15 18],abs(Gug_V(40,[3 6 9 12 15 18])),'k*') % poles
        plot([5 7 8 11 19],abs(Gug_V(40,[5 7 8 11 19])),'ro') % big solar
        hold off
        xlabel('bus')
        ylabel('Voltage [p.u.]')
        xlim([1 nBuses]); ylim([0.95 1.07])
        xticks(0:1:nBuses);
        legend('ACOPF: V @ 1.00pm','ACOPF: V @ 8.00pm','Gug: V @ 1.00pm','Gug: V @ 8.00pm','Location','SouthWest')
        title('Voltage Profile')
        set(gcf,'color','w');
        grid on
        grid minor
        
        % ============  LINE LOSSES
        f5 = figure(5);
        movegui(f5,'southeast');
        %plot(1:size(ACOPF.branch,1),ACOPF_I2R(:,[25 40]),'+')
        %hold on
        plot([1 4 7 10 13 16], ACOPF_I2R([26 40],[1 4 7 10 13 16]),'ro')
        hold on
        plot([1 4 7 10 13 16], abs(Gug_I2R([26 40],[1 4 7 10 13 16])),'bo')
        plot([2 3 5 6 8 9 11 12 14 15 17 18], ACOPF_I2R([26 40],[2 3 5 6 8 9 11 12 14 15 17 18]),'*')
        plot([2 3 5 6 8 9 11 12 14 15 17 18], abs(Gug_I2R([26 40],[2 3 5 6 8 9 11 12 14 15 17 18])),'*')
        xlabel('Lines')
        ylabel('Line losses [kW]')
        legend('ACOPF: Pole-to-pole @ 1.00pm','ACOPF: Pole-to-pole @ 8.00pm',...
            'Gug: Pole-to-pole @ 1.00pm','Gug: Pole-to-pole @ 8.00pm',...
            'Location','NorthEast')
        title('Line Losses')
        xlim([0 nLines]); %ylim([0.9 1.05])
        %xticks(1:1:size(ACOPF.branch,1));
        set(gcf,'color','w');
        grid on
        grid minor
        
        % ============  REACTIVE POWER PURCHASE
        f7 = figure(7);
        movegui(f7,'northeast');
        %plot(timE, ACOPF_PgQg)
        %hold on
        %plot(timE, loadTotal, 'g')
        %hold on
        %plot(timE, abs(sum(ACOPF_I2R,2)), 'k','LineWidth',2)
        plot(timE, abs(Gug_QgTot), 'r','LineWidth',2)
        % yyaxis right
        % plot(t, ACOPF_f,'*')
        hold off
        xlabel('Time [Hours]')
        ylabel('Grid Purchases [kW]')
        legend('Active Grid Power (P)','Reactive Grid Power (Q)','Total Load','Line Losses','Location','SouthWest')
        title('Reactive power from the grid')
        grid on
        set(gcf,'color','w');
        %% ============  Volt/Var
        f8 = figure(8);
        movegui(f8,'northeast');        
        scatter(Gug_V(:,18),abs(Gug_QcTot))
        xlim([0.93 1.1])
        %ylim([-0.01 0.01])
        xlabel('Voltage Magnitude [p.u.]')
        ylabel('Reactive Power [kVAR]')
        grid on
        set(gcf,'color','w');    
        % ============  VOlt/Watt
        f9 = figure(9);
        movegui(f9,'southeast');        
        scatter(Gug_V(:,18),abs(Gug_PcTot))
        xlim([0.93 1.1])
        %ylim([-0.01 0.01])
        xlabel('Voltage Magnitude [p.u.]')
        ylabel('Reactive Power [kVAR]')
        grid on
        set(gcf,'color','w');   
    end
% if solarScen == 1
%     f9 = figure(9) 
%         movegui(f9,'southwest');
%         %plot(timE,store_Gug_V(abs(Gug_V(:,3)),'r--')
%         %hold on
%         %plot(timE,abs(Gug_V(:,9)),'b--')
%         plot(timE,abs(store_Gug_V),'k--')
%         %plot(timE,ones(T,1)*1.05,'Marker','*','Color','r')
%         %hold off
%         %plot(t,ones(size(solarHH,1),1)*0.95,'Marker','*','Color','r')
%         xlabel('Time (Hour)')
%         ylabel('Voltage [p.u.]')
%         ylim([0.94 1.08])
%         %xticks(0:2:size(loadHH,1));
%         title('Voltage profile over day')
%         legend({'ACOPF: Pole 3','ACOPF: Pole 9','ACOPF: Pole 18','Gug: Pole 3','Gug: Pole 9','Gug: Pole 18'},'Location','northwest')
%         set(gcf,'color','w');
%         grid on    
% end
end
