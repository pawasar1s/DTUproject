function [V, Pc, Qc, Vmax, Gug_V, Gug_I2R, Gug_ITot, Gug_PgTot, Gug_QgTot, Gug_I2RTot, Gug_PcTot, Gug_QcTot] = VoltVar(testCase, T, T0, solar, loadHH, multiPer, per) 
% inputs
[~, ZBus, Ysc, ~, Ymn, Imax, nBuses, ~, nB] = readLinesMPC(testCase);
[~, Vnom, Vmin, Vmax, deltaV, V0, Pd, Qd, Pav, Sinj, A, B, C, D, PF] = readGensMPC(testCase, nBuses);
%changePd = changePd; % between 2-10
inverterSize = 1.1;
%deadband = 0.02;
if multiPer == 0
    tic;
    for k = 1% : load_scenarios
        mpc = testCase;
        baseMVA = mpc.baseMVA;
        nPV = size(mpc.gen,1)-1; % No of PV systems
        idxPV = find(mpc.bus(:,2) == 2); % idx for buses with PV systems
        % update solar data
        if nPV ~= 0
             Pav(idxPV-1) = mpc.gen(2:end,9)*solar(per)/baseMVA; % Pmax
             Sinj(idxPV-1) = mpc.gen(2:end,9)*solar(per)*inverterSize/baseMVA; % Pmin
             Qmin = -tan(acos(PF))*(Pav)+0.1*Pav;
        end
        %slope = Qmin./(V0 - deadband/2 - Vmin(2));
        slope = 2*Qmin./(Vmax(2)-Vmin(2));
        % update load data
        Pd(idxPV-1) = loadHH(per,idxPV-1)'/baseMVA; % Pg
        Qd(idxPV-1) = loadHH(per,idxPV-1)'*0.6/baseMVA; % Qg, 0.6 assumed of Pd (we don't have Qd data)
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
        Obj3_Vec = eye(nB) - (1/nB).*ones(nB,1)*ones(1,nB);
        Obj3 = norm(Obj3_Vec*diag(V),2)
        minimize(Obj1 + Obj2 + Obj3) % Eq.4
        
        % CONSTRAINTS ===========================================
        subject to
        % Solar PV constraints
        0 <= Pc <= Pav; %Eq.9
        %(Qc).^2 <= (Sinj).^2-(Pav-Pc).^2; %Eq.10
        %abs(Qc) <= tan(acos(PF))*(Pav-Pc); %Eq.11
        Qc == slope.*(real(V)-Vnom(2));
        % Power balance eq.
        real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); %Eq.5
        imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); %Eq.5
        % Voltage magnitude limits
        Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; % Eq.7 & 8
        cvx_end
        
        % SAVE RESULTS ===========================================
        % Voltage vector Invalid constraint: {real affine} == {convex}
        Gug_V = [V0; V];
        % Line Losses
        [Gug_I2R, Vdrop, Gug_I] = linesLoss(mpc, Gug_V, nBuses);
        Gug_I2R = Gug_I2R'; Gug_V = Gug_V'; Gug_I = Gug_I'; 
        % Reactive Power 
        Gug_QgTot = sum(Qc-Qd);
        Gug_PgTot = sum(Pd-Pav+Pc);
        % Reactive Power 
        
        Gug_PcTot = sum(Pc);
        Gug_QcTot = sum(Qc);
        % Line losses
        Gug_I2RTot = sum(Gug_I2R);
        Gug_ITot = Gug_I;
    end
    toc;
elseif multiPer == 1
    tic;
    mpc = testCase;
    baseMVA = mpc.baseMVA;
    nPV = size(mpc.gen,1)-1; % No of PV systems
    idxPV = find(mpc.bus(:,2) == 2); % idx for buses with PV systems
    Gug_V = complex(zeros(T-T0+1,size(mpc.bus,1))); % Bus Voltage Vector
    %Gug_I2R = complex(zeros(T,size(mpc.bus,1)-1)); % Line Loss Vector
    Gug_QgTot = complex(zeros(T-T0+1,1)); % Line Loss Vector
    Gug_PgTot = complex(zeros(T-T0+1,1)); % Line Loss Vector
    Gug_PcTot = complex(zeros(T-T0+1,1)); % Line Loss Vector
    Gug_QcTot = complex(zeros(T-T0+1,1)); % Line Loss Vector
    Gug_I2RTot = complex(zeros(T-T0+1,1)); % Line Loss Vector
    Gug_ITot = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % Line current 
    %
    Gug_QcInd = complex(zeros(T-T0+1,size(mpc.bus,1)-1));  % Line Loss Vector
    Gug_QminInd = complex(zeros(T-T0+1,size(mpc.bus,1)-1));  % Line Loss Vector
    %
    Gug_check_Sinj = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); % check Sinj
    Gug_actual_Sinj = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); 
    Gug_check_PF = complex(zeros(T-T0+1,size(mpc.bus,1)-1));
    Gug_check_slope = complex(zeros(T-T0+1,size(mpc.bus,1)-1)); 
    for t = T0 : T 
        % INPUT PARAMETERS ====================================================
        % update solar data
        if nPV ~= 0
            Pav(idxPV-1) = mpc.gen(2:end,9)*solar(t)/baseMVA; % Pmax
            Sinj(idxPV-1) = mpc.gen(2:end,9)*solar(t)/baseMVA;%*inverterSize/baseMVA; % Pmin
            Qmin = -tan(acos(PF))*(Pav)+0.1*Pav;
        end
        Gug_QminInd(t-T0+1,:) = Qmin; % minus to show maax absorbtion level
        Gug_actual_Sinj(t-T0+1,:) = Sinj; 
        % with deadband
        %slope = Qmin./(V0 + deadband/2 - Vmax(2));
        % without deadband
        slope = 2*Qmin./(Vmax(2)-Vmin(2));
        Gug_check_slope(t-T0+1,:) = slope;
        % update load data
        Pd(idxPV-1) = loadHH(t,idxPV-1)'/baseMVA; % Pg
        Qd(idxPV-1) = loadHH(t,idxPV-1)'*0.6/baseMVA; % Qg, 0.6 assumed of Pd (we don't have Qd data)
         % BEGIN =========================================================
        cvx_clear % clears any model being constructed
        cvx_begin
        variable V(nBuses-1,1) complex; % Voltage vector [p.u.]
        variable Pc(nBuses-1,1); % Curtailed PV vector [MW]
        variable Qc(nBuses-1,1); % Reative power absorbed/produced by inverter [MW]
        % OBJECTIVE FUNCTION =============================================
        Obj1_Vec = [real(V0); real(V); imag(V0); imag(V)];
        Obj1 = 0.5 * quad_form(Obj1_Vec, Ymn); %Eq.1
        Obj2 = 0.5 * quad_form(Pc, A) + quad_form(Qc, C) + B' * Pc + D' * abs(Qc); %Eq.2
        %Obj3 = 0;
        Obj3_Vec = eye(nB) - (1/nB).*ones(nB,1)*ones(1,nB);
        Obj3 = norm(Obj3_Vec*diag(V),2)
        minimize(Obj1 + Obj2 + Obj3) % Eq.4
        
        % CONSTRAINTS ====================================================
        subject to
        % Solar PV constraints
        0 <= Pc <= Pav; %Eq.9
        %Qc(idxPV-1) == slope(idxPV-1).*((V(idxPV-1))-Vnom(2)+deadband/2)
        Qc == slope.*(real(V)-Vnom(2));
        % Power balance eq.
        real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); %Eq.5
        imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); %Eq.5
        % Voltage magnitude limits
        Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; % Eq.7 & 8
        cvx_end
        
        % Voltage vector
        Vfull = [V0; V];
        Gug_V(t-T0+1,:) = Vfull;
        % Line Losses
        [Gug_I2R, Gug_Vdrop, Gug_I] = linesLoss(mpc, Vfull, nBuses);
        % Reactive Power 
        Gug_QgTot(t-T0+1,:) = sum(Qc-Qd);
        Gug_PgTot(t-T0+1,:) = sum(Pd-Pav+Pc);
        % Reactive Power 
        Gug_PcTot(t-T0+1,:) = sum(Pc);
        Gug_QcTot(t-T0+1,:) = sum(Qc);
        Gug_QcInd(t-T0+1,:) = Qc;
        % Line losses
        Gug_I2RTot(t-T0+1,:) = sum(Gug_I2R);
        Gug_ITot(t-T0+1,:) = Gug_I;
        % validate results
        Gug_check_Sinj(t-T0+1,:) = sqrt(Qc.^2+(Pav-Pc).^2);
        Gug_check_PF(t-T0+1,:) = (Pav-Pc)./sqrt((Pav-Pc).^2+Qc.^2);
    end
    toc;
 %% plot Qmax vs Qc 
    % ================
    select_house = 18;
    % ================
    close all;
    figure(101)
    plot(Gug_QminInd(:,select_house),'*');hold on; plot(Gug_QcInd(:,select_house),'o')
    ylabel('Reactive Power [kVAR]')
    legend({'Qmin','Observed Qc'})
    % 
    figure(102)
    plot(1:nBuses, Gug_V(26,:))%;hold on; plot(Gug_QcInd(26,:),'o')
    ylabel('Reactive Power [kVAR]')
    legend({'Voltage'})   
    % 
    figure(103)
    plot(Gug_check_Sinj(:,select_house),'*');hold on; plot(Gug_actual_Sinj(:,select_house),'o')
    ylabel('Apparent power S_inj [kVA]')
    legend({'Calculated Apparent Power', 'Actual Apparent Power'})  
    % 
    figure(104)
    plot(T0:T, Gug_PcTot)%;hold on; plot(Gug_QcInd(26,:),'o')
    ylabel('Active Power [kW]')
    legend({'Active Power Curtailment'})     
    %
    % ===================
    select_30minper = [18 21 24];
    % ===================
    figure(105)
    plot(Gug_V(select_30minper,select_house), Gug_QcInd(select_30minper,select_house), 'r*')%;hold on; plot(Gug_QcInd(26,:),'o')
    xlim([0.95 1.05])
    hold on 
    plot([Vmin(2) Vmax(2)], [-Gug_QminInd(select_30minper,select_house) Gug_QminInd(select_30minper,select_house)])
    xlim([0.95 1.05])
    ylabel('Active Power [kW]')
    title('Volt/Var curve for different hours')
    legend({'Reactive power vs voltage'})    
    % 
    figure(106)
    plot(Gug_check_PF(:,select_house), '*')%;hold on; plot(1:48, 0.8*ones(1,48),'o')
    ylabel('Power Factor')
    legend({'Power Factor'})      
end
end