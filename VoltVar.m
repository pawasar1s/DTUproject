%function [V, Pc, Qc, Gug_V, Gug_I2R, Gug_PgTot, Gug_QgTot, Gug_I2RTot, Gug_PcTot, Gug_QcTot] = VoltVar(testCase, T, solar, loadHH, multiPer, per)
% inputs
[~, ZBus, Ymn, nBuses, ~, nB] = readLinesMPC(testCase);
[~, Vnom, Vmin, Vmax, ~, V0, Pd, Qd, Pav, Sinj, A, B, C, D, PF] = readGensMPC(testCase, nBuses);
%changePd = changePd; % between 2-10
inverterSize = 1.1;
deadband = 0.02;
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
             %Qmin = sqrt(Sinj.^2-Pav.^2);
             Qmin = -tan(acos(PF))*(Pav);%+0.1*Pav;
        end
        %slope = Qmin./(V0 - deadband/2 - Vmin(2));
        %sloppy = slope.*ones(nB,1)
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
        variable gg(nBuses-1,1); % Volt/VAR slope 
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
        (Qc).^2 <= (Sinj).^2-(Pav-Pc).^2; %Eq.10
        abs(Qc) <= tan(acos(PF))*(Pav-Pc); %Eq.11
        %Qc(idxPV-1)./Qmax(idxPV-1) == -slope(idxPV-1).*(V(idxPV-1))-0.98  
        %Qc(idxPV-1) == -slope(idxPV-1).*(real(V(idxPV-1))-Vnom(2)+deadband/2) 
        slope = 2*Qmin./(Vmax(2)-Vmin(2));
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
        [Gug_I2R, Vdrop] = linesLoss(mpc, Gug_V);
        Gug_I2R = Gug_I2R'; Gug_V = Gug_V';
        % Reactive Power 
        Gug_QgTot = sum(Qc-Qd);
        Gug_PgTot = sum(Pd-Pav+Pc);
        % Reactive Power 
        
        Gug_PcTot = sum(Pc);
        Gug_QcTot = sum(Qc);
        % Line losses
        Gug_I2RTot = sum(Gug_I2R);
    end
    toc;
elseif multiPer == 1
    tic;
    T = 19
    mpc = testCase;
    baseMVA = mpc.baseMVA;
    nPV = size(mpc.gen,1)-1; % No of PV systems
    idxPV = find(mpc.bus(:,2) == 2); % idx for buses with PV systems
    Gug_V = complex(zeros(T,size(mpc.bus,1))); % Bus Voltage Vector
    %Gug_I2R = complex(zeros(T,size(mpc.bus,1)-1)); % Line Loss Vector
    Gug_QgTot = complex(zeros(T,1)); % Line Loss Vector
    Gug_PgTot = complex(zeros(T,1)); % Line Loss Vector
    Gug_PcTot = complex(zeros(T,1)); % Line Loss Vector
    Gug_QcTot = complex(zeros(T,1)); % Line Loss Vector
    Gug_I2RTot = complex(zeros(T,1)); % Line Loss Vector
    %
    Gug_QcInd = complex(zeros(T,size(mpc.bus,1)-1));  % Line Loss Vector
    Gug_QminInd = complex(zeros(T,size(mpc.bus,1)-1));  % Line Loss Vector
    %
    Gug_check_Sinj = complex(zeros(T,size(mpc.bus,1)-1)); % check Sinj
    Gug_actual_Sinj = complex(zeros(T,size(mpc.bus,1)-1)); 
    Gug_check_PF = complex(zeros(T,size(mpc.bus,1)-1));
    Gug_check_slope = complex(zeros(T,size(mpc.bus,1)-1)); 
    for t = 14 : 26 
        % INPUT PARAMETERS ====================================================
        % update solar data
        if nPV ~= 0
            Pav(idxPV-1) = mpc.gen(2:end,9)*solar(t)/baseMVA; % Pmax
            Sinj(idxPV-1) = mpc.gen(2:end,9)*solar(t)*inverterSize/baseMVA; % Pmin
            %Qmin = sqrt(Sinj.^2-Pav.^2)+0.1*Pav;
            Qmin = -tan(acos(PF))*(Pav);%+0.1*Pav;
        end
        Gug_QminInd(t,:) = Qmin; % minus to show maax absorbtion level
        Gug_actual_Sinj(t,:) = Sinj; 
        %slope = Qmin./(V0 + deadband/2 - Vmax(2));
        slope = 2*Qmin./(Vmax(2)-Vmin(2));
        Gug_check_slope(t,:) = slope;
        %slope = -Qmin./(V0 - Vmax(2));
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
        (Qc).^2 <= (Sinj).^2-(Pav).^2;%+0.1*Pav; %Eq.10 * including inverter oversizing
        abs(Qc) <= tan(acos(PF))*(Pav-Pc); %Eq.11
        %Qc(idxPV-1) == slope(idxPV-1).*((V(idxPV-1))-Vnom(2)+deadband/2)
        %Qc(idxPV-1) == slope(idxPV-1).*((V(idxPV-1))-Vmin(2))
        Qc(idxPV-1) == slope(idxPV-1).*((V(idxPV-1))-Vnom(2))
        
        % Power balance eq.
        real(V) == Vnom + real(ZBus)*(Pav - Pc - Pd) + imag(ZBus)*(Qc - Qd); %Eq.5
        imag(V) == imag(ZBus)*(Pav - Pc - Pd) - real(ZBus)*(Qc - Qd); %Eq.5
        % Voltage magnitude limits
        Vmin <= Vnom + real(ZBus)*(Pav - Pc - Pd)+ imag(ZBus)*(Qc - Qd) <= Vmax; % Eq.7 & 8
        cvx_end
        
        % Voltage vector
        Vfull = [V0; V];
        Gug_V(t,:) = Vfull;
        % Line Losses
        [Gug_I2R, Vdrop] = linesLoss(mpc, Vfull);
        % Reactive Power 
        Gug_QgTot(t,:) = sum(Qc-Qd);
        Gug_PgTot(t,:) = sum(Pd-Pav+Pc);
        % Reactive Power 
        Gug_PcTot(t,:) = sum(Pc);
        Gug_QcTot(t,:) = sum(Qc);
        Gug_QcInd(t,:) = Qc;
        % Line losses
        Gug_I2RTot(t,:) = sum(Gug_I2R);
        % validate results
        Gug_check_Sinj(t,:) = sqrt(Qc.^2+(Pav-Pc).^2);
        Gug_check_PF(t,:) = (Pav-Pc)./sqrt((Pav-Pc).^2+Qc.^2);
    end
    toc;
 %% plot Qmax vs Qc 
 
    close all;
    figure(101)
    plot(Gug_QminInd(:,18),'*');hold on; plot(Gug_QcInd(:,18),'o')
    ylabel('Reactive Power [kVAR]')
    legend({'Qmin','Observed Qc'})
    % 
    figure(102)
    plot(1:19, Gug_V(26,:))%;hold on; plot(Gug_QcInd(26,:),'o')
    ylabel('Reactive Power [kVAR]')
    legend({'Voltage'})   
    % 
    figure(103)
    plot(Gug_check_Sinj(:,18),'*');hold on; plot(Gug_actual_Sinj(:,18),'o')
    ylabel('Apparent power S_inj [kVA]')
    legend({'Calculated Apparent Power', 'Actual Apparent Power'})  
    % 
    figure(104)
    plot(1:26, Gug_PcTot)%;hold on; plot(Gug_QcInd(26,:),'o')
    ylabel('Active Power [kW]')
    legend({'Active Power Curtailment'})     
    % 
    figure(105)
    plot(Gug_V(:,18), Gug_QcInd(:,18), '*')%;hold on; plot(Gug_QcInd(26,:),'o')
    ylabel('Active Power [kW]')
    legend({'Reactive power vs voltage'})    
    % 
    figure(106)
    plot(Gug_check_PF(:,18), '*')%;hold on; plot(1:48, 0.8*ones(1,48),'o')
    ylabel('Power Factor')
    legend({'Power Factor'})      
end
%end