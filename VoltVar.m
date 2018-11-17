function [V, Pc, Qc, Vmax, Gug_V, Gug_I2R, Gug_ITot, Gug_Vdrop, Gug_PgTot, Gug_QgTot, Gug_I2RTot, Gug_PcTot, Gug_QcTot] = VoltVar(testCase, T, T0, solar, loadHH, multiPer, per, plotting) 
[~, ZBus, Ysc, Aa, Ymn, Imax, nBuses, ~, nB] = readLinesMPC(testCase);
[~, Vnom, Vmin, Vmax, V0, Pd, Qd, Pav, Sinj, A, B, C, D, PF] = readGensMPC(testCase, nBuses);
%deadband = 0.02;
if multiPer == 0
    tic;
    for k = 1% : load_scenarios
        mpc = testCase;
        baseMVA = mpc.baseMVA;
        nPV = size(mpc.gen,1)-1; % No of PV systems
        idxPV = find(mpc.bus(:,2) == 2); % find index for buses with PV systems
        overSize = zeros(nBuses-1, 1);
        overSize(idxPV-1) = 0.1/baseMVA*mpc.gen(2:end,9); % inverter oversized by 10%
        % update solar data
        if nPV ~= 0
             Pav(idxPV-1) = mpc.gen(2:end,9)*solar(per)/baseMVA; % Pmax
             Sinj(idxPV-1) = mpc.gen(2:end,9)*solar(per)/baseMVA + oversize; % Smax
             Qmin = -tan(acos(PF))*(Pav) + oversize;
        end
        %slope = Qmin./(V0 - deadband/2 - Vmin(2));
        slope = 2*Qmin./(Vmax(2)-Vmin(2)); % droop control slope
        % update load data
        Pd(idxPV-1) = loadHH(per,idxPV-1)'/baseMVA; % Pg
        Qd(idxPV-1) = loadHH(per,idxPV-1)'*0.6/baseMVA; % Qg: 0.6*Pg - assumed
    % START ====================================================
        cvx_clear % clears any model being constructed
        cvx_begin
        variable V(nBuses-1,1) complex; % Voltage vector [p.u.]
        variable Pc(nBuses-1,1); % Curtailed PV vector [kW]
        variable Qc(nBuses-1,1); % Reative power absorbed/produced by inverter [kVA]
    % OBJECTIVE FUNCTION =======================================
        Obj1_Vec = [real(V0); real(V); imag(V0); imag(V)];
        Obj1 = 0.5 * quad_form(Obj1_Vec, Ymn); %Eq.1
        Obj2 = 0.5 * quad_form(Pc, A) + quad_form(Qc, C) + B' * Pc + D' * abs(Qc); %Eq.2
        Obj3_Vec = eye(nB) - (1/nB).*ones(nB,1)*ones(1,nB);
        Obj3 = norm(Obj3_Vec*diag(V),2); %Eq.3
        minimize(Obj1 + Obj2 + Obj3) % Eq.4
    % CONSTRAINTS =============================================
        subject to
        % Solar PV constraints
        0 <= Pc <= Pav; 
        Qc == slope.*(real(V)-Vnom(2)); % droop control constraint
        % Power balance eq.
        real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); 
        imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); 
        % Voltage magnitude limits
        Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; 
        % Line limit constraint
        Delta_V = Aa * [V0; V]; % Delta_V: votlage drop
        -Imax <= real(Ysc.*conj(Delta_V)) <= Imax; % Ysc: line admittance
        cvx_end
    % SAVE RESULTS ===========================================
        Gug_V = [V0; V];
        [Gug_I2R, Vdrop, Gug_I] = line_Losses2(mpc, Gug_V, nBuses); % losses, voltage drop, line current
        Gug_I2R = Gug_I2R'; Gug_V = Gug_V'; Gug_I = Gug_I'; 
        Gug_QgTot = sum(Qc-Qd); % Q from grid
        Gug_PgTot = sum(Pd-Pav+Pc); % P from grid
        Gug_PcTot = sum(Pc); % PV output curtailment
        Gug_QcTot = sum(Qc); % PV reactive power
        Gug_I2RTot = sum(Gug_I2R); % line loss
        Gug_ITot = Gug_I; % line current
        Gug_Vdrop = Vdrop; % voltage drop 
    end
    toc;
elseif multiPer == 1
    tic;
    mpc = testCase;
    baseMVA = mpc.baseMVA;
    nPV = size(mpc.gen,1)-1; % No of PV systems
    idxPV = find(mpc.bus(:,2) == 2); % find index for buses with PV systems
    overSize1 = 0.1/baseMVA*mpc.gen(2:end,9); % inverter oversized by 10%. vector (12,1)
    overSize2 = zeros(nBuses-1, 1); % inverter oversized by 10%. vector (18,1)
    overSize2(idxPV-1) = 0.1/baseMVA*mpc.gen(2:end,9); 
    % allocate empty matrices 
    Gug_V = complex(zeros(T-T0+1,size(mpc.bus,1))); % bus voltage
    Gug_QgTot = complex(zeros(T-T0+1,1)); % Q from grid
    Gug_PcTot = complex(zeros(T-T0+1,1)); % TOTAL PV output curtailment
    Gug_QcTot = complex(zeros(T-T0+1,1)); % TOTAL PV reactive power 
    Gug_PgTot = complex(zeros(T-T0+1,1)); % P from grid
    Gug_I2RTot = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % line loss
    Gug_ITot = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % line current 
    Gug_Vdrop = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % voltage drop
    Gug_QcInd = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % PV reactive power 
    Gug_QminInd = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % maximum PV reactive power (due to PF)
    Gug_Pinj = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % PV output injected in the grid
    Gug_check_Sinj = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % FOR GRAPHS: real inverter capacity considering Pc   
    Gug_actual_Sinj = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % FOR GRAPHS: inverter capacity based on Pav    
    Gug_check_PF = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % FOR GRAPHS: PF
    %
    Gug_check_slope = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % FOR GRAPHS: slope
    for t = T0 : T 
        % update solar data
        if nPV ~= 0
            Pav(idxPV-1) = mpc.gen(2:end,9)*solar(t)/baseMVA; % Pmax
            Sinj(idxPV-1) = mpc.gen(2:end,9)*solar(t)/baseMVA + overSize1; % Smax
            Qmin = -tan(acos(PF))*(Pav) + overSize2;
        end
        Gug_QminInd(t-T0+1,:) = Qmin; % minus to show maax absorbtion level
        Gug_actual_Sinj(t-T0+1,:) = Sinj; % FOR GRAPHS: inverter capacity based on Pav    
        % slope with deadband
        %slope = Qmin./(V0 + deadband/2 - Vmax(2));
        % slope without deadband
        slope = 2*Qmin./(Vmax(2)-Vmin(2));
        Gug_check_slope(t-T0+1,:) = slope; % FOR GRAPHS: slope
        % update load data
        Pd(idxPV-1) = loadHH(t,idxPV-1)'/baseMVA; % Pg
        Qd(idxPV-1) = loadHH(t,idxPV-1)'*0.6/baseMVA; % Qg: 0.6*Pg - assumed
    % START ====================================================
        cvx_clear % clears any model being constructed
        cvx_begin
        variable V(nBuses-1,1) complex; % Voltage vector [p.u.]
        variable Pc(nBuses-1,1); % Curtailed PV vector [kW]
        variable Qc(nBuses-1,1); % Reative power absorbed/produced by inverter [kVA]
    % OBJECTIVE FUNCTION =============================================
        Obj1_Vec = [real(V0); real(V); imag(V0); imag(V)];
        Obj1 = 0.5 * quad_form(Obj1_Vec, Ymn); %Eq.1
        Obj2 = 0.5 * quad_form(Pc, A) + quad_form(Qc, C) + B' * Pc + D' * abs(Qc); %Eq.2
        Obj3_Vec = eye(nB) - (1/nB).*ones(nB,1)*ones(1,nB);
        Obj3 = norm(Obj3_Vec*diag(V),2) %Eq.3
        minimize(Obj1 + Obj2 + Obj3) % Eq.4
    % CONSTRAINTS ====================================================
        subject to
        % Solar PV constraints
        0 <= Pc <= Pav; %Eq.9
        Qc == slope.*(real(V)-Vnom(2)); % droop control w/o deadband
        %Qc(idxPV-1) == slope(idxPV-1).*((V(idxPV-1))-Vnom(2)+deadband/2); % droop control with deadband
        % Power balance eq.
        real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); 
        imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); 
        % Voltage magnitude limits
        Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; 
        % Line limit constraint
        Delta_V = Aa * [V0; V]; % Delta_V: votlage drop
        -Imax <= real(Ysc.*conj(Delta_V)) <= Imax; % Ysc: line admittance
        cvx_end
    % SAVE RESULTS ===========================================    
        Vfull = [V0; V];
        Gug_V(t-T0+1,:) = Vfull; % voltage 
        [Gug_I2R, Gug_Vdrop(t-T0+1,:), Gug_ITot(t-T0+1,:)] = linesLoss(mpc, Vfull, nBuses); % losses, voltage drop, line current
        Gug_QgTot(t-T0+1,:) = sum(Qc-Qd); % Q from grid
        Gug_PgTot(t-T0+1,:) = sum(Pd-Pav+Pc); % P from grid
        Gug_PcTot(t-T0+1,:) = sum(Pc); % TOTAL PV output curtailment
        Gug_QcTot(t-T0+1,:) = sum(Qc); % TOTAL PV reactive power 
        Gug_QcInd(t-T0+1,:) = Qc; % PV reactive power 
        Gug_Pinj(t-T0+1,:) = Pav-Pc; % PV output injected in the grid
        Gug_I2RTot(t-T0+1,:) = Gug_I2R; % line loss
        % validate results
        Gug_check_Sinj(t-T0+1,:) = sqrt(Qc.^2+(Pav-Pc).^2); % FOR GRAPHS: real inverter capacity considering Pc
        Gug_check_PF(t-T0+1,:) = (Pav-Pc)./sqrt((Pav-Pc).^2+Qc.^2); % FOR GRAPHS: PF
        Gug_check_PF(Gug_check_PF<=0.75) = 1; % FOR GRAPHS: PF
    end
    toc;
    %% Graphs
    if plotting == 1
        %close all;
        select_house = 18;
        select_30minper = [18 21 24];
        % inverter available power vs injected power
        figure(200);
        plot(1:nB, (Gug_Pinj([26 30 34],:)),'s'); hold on;
        Pav(idxPV-1) = mpc.gen(2:end,9)*solar(26)/baseMVA;
        plot(1:nB, Pav, 'bx')
        Pav(idxPV-1) = mpc.gen(2:end,9)*solar(30)/baseMVA;
        plot(1:nB, Pav, 'rx')
        Pav(idxPV-1) = mpc.gen(2:end,9)*solar(34)/baseMVA;
        plot(1:nB, Pav, 'kx')
        xlim([0 nB+1])
        ylabel('Injected PV output [kW]'); xlabel('Bus')
        title('Inverter injected power')
        legend({'Pinj 1pm','Pinj 3pm','Pinj 5pm','Pav 1pm','Pav 3pm','Pav 5pm'},'Location','Southwest')
        grid on
        set(gcf,'color','w'); grid on
        % Inverter reactive power
        figure(201)
        plot(Gug_QminInd(:,select_house),'b*');hold on; plot(Gug_QcInd(:,select_house),'b-')
        plot(Gug_QminInd(:,select_house-6),'r*');hold on; plot(Gug_QcInd(:,select_house-6),'r--')
        ylabel('Reactive Power Q_{c} [kVAR]'); xlabel('Time')
        ylim([-0.2 0.05])
        legend({'HH18 Qc max','HH18 Qc actual', 'HH12 Qc max','HH12 Qc actual'},'Location','Southwest')
        title('Inverter reactive Power')
        set(gcf,'color','w'); grid on
        % voltage profile over feeder
        figure(202)
        plot(1:nBuses, Gug_V([26 30 34],:))%;hold on; plot(Gug_QcInd(26,:),'o')
        ylabel('Voltage Magnitude [p.u.]'); xlabel('Bus')
        ylim([0.99 1.05])
        legend({'V 1pm','V 3pm','V 5pm'},'Location','Southwest')
        title('Voltage profile')
        set(gcf,'color','w'); grid on
        % Inverter apparent power
        figure(203)
        plot(Gug_check_Sinj(:,select_house),'*');hold on; plot(Gug_actual_Sinj(:,select_house),'o')
        ylabel('Apparent power S_{inj} [kVA]'); xlabel('Bus')
        legend({'S = P_{av}*1.1', 'S = sqrt((P_{av}-P_c)^2 + Q_c^2)'},'Location','Northwest')
        title('Inverter apparent power')
        set(gcf,'color','w'); grid on
        % Inverter power curtailment
        figure(204)
        plot(T0:T, Gug_PcTot)%;hold on; plot(Gug_QcInd(26,:),'o')
        ylabel('Active Power P_{c} [kW]'); xlabel('Time')
        legend({'Pc'})
        title('Inverter power curtailment')
        ylim([0 0.45])
        set(gcf,'color','w'); grid on
        % Inverter Volt/Var ratio
        figure(205)
        plot(Gug_V(select_30minper,select_house), Gug_QcInd(select_30minper,select_house), 'r*')%;hold on; plot(Gug_QcInd(26,:),'o')
        xlim([0.95 1.05])
        hold on
        plot([Vmin(2) Vmax(2)], [-Gug_QminInd(select_30minper,select_house) Gug_QminInd(select_30minper,select_house)])
        xlim([Vmin(2) Vmax(2)])
        ylim([-0.3 0.01])
        ylabel('Reactive Power Q_{c} [kW]'); xlabel('Voltage [p.u.]')
        title('Inverter Volt/Var results')
        legend({'Q_c vs V'},'Location','Southwest')
        set(gcf,'color','w'); grid on
        % PF
        figure(206)
        plot(Gug_check_PF(:,select_house), '-')
        ylabel('Power Factor'); xlabel('Time')
        legend({'PF'},'Location','Southwest')
        title('Inverter power factor')
        ylim([0.5 1.05]); xlim([0 50])
        set(gcf,'color','w'); grid on
        % Line Current
        figure(207)
        plot(1:nB, Gug_ITot([26 36],:).*312.5, '*');hold on % 312.5 is Ibase
        plot([1 4 7 10 13 16], Gug_ITot([26 36],[1 4 7 10 13 16]).*312.5,'o')
        ylabel('Current [A]'); xlabel('Time')
        legend({'1pm Drop Line','6pm Drop Line','1pm Pole-to-pole Line', '6pm Pole-to-pole Line'},'Location','Southeast')
        title('Line current')
        set(gcf,'color','w'); grid on
    end
end
end