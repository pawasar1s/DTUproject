function [ACOPF_struct, ACOPF_V, ACOPF_I2R, ACOPF_PVgen, ACOPF_f, ACOPF_Qg, ACOPF_I, ACOPF_Vdrop, ACOPF_struct_container] = ACOPF(testCase, T, T0, solar, loadHH, ACOPFsolver, multiPer, per, PF)
%% ACOPF
% choose solver
if ACOPFsolver == 1
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT', 'opf.flow_lim', 'I'); % ACOPF
else
    mpopt = mpoption('model','AC', 'pf.alg', 'ISUM', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT', 'opf.flow_lim', 'I'); % Backward/Forward Sweep with current summation
end
% run ACOPF
if multiPer == 0
    tic;
    mpc = testCase;
    nPV = size(mpc.gen,1)-1; % No of PV systems
    idxHH = find(mpc.bus(:,2) == 2); % find index for buses with PV systems
    % update solar PV data
    if nPV ~= 0
        mpc.gen(2:end,9) = mpc.gen(2:end,9)*solar(per); % Pmax
        mpc.gen(2:end,10) = mpc.gen(2:end,10)*solarHH(per); % Pmin
        mpc.gen(2:end,4) = mpc.gen(2:end,4)*solar(per)*tan(acos(PF)); % Qmax - ensures min PF
        mpc.gen(2:end,5) = mpc.gen(2:end,5)*solar(per)*(-tan(acos(PF))); % Qmin - ensures min PF
    end
    % update load data
    mpc.bus(idxHH,3) = loadHH(per,idxHH-1)'; % Pg 
    mpc.bus(idxHH,4) = loadHH(per,idxHH-1)'*0.6; % Qg: 0.6*Pg - assumed
    % save results
    ACOPF_struct = runopf(mpc,mpopt); % run ACOPF
    ACOPF_V = ACOPF_struct.bus(:,8); % save V
    %ACOPF_I2R= get_losses(ACOPF_struct);
    [ACOPF_I2R, ACOPF_Vdrop, ACOPF_I] = line_Losses2(ACOPF_struct); % save losses, voltage drop, line current
    ACOPF_f = ACOPF_struct.f; % save objective f value
    ACOPF_Qg = ACOPF_struct.var.val.Qg; % save Q from grid
    clear mpc;
    toc;
elseif multiPer == 1
    tic;
    % allocate empty matrices
    ACOPF_V = complex(zeros(T,size(testCase.bus,1))); % bus voltage 
    ACOPF_I2R = complex(zeros(T,size(testCase.bus,1)-1)); % line loss
    ACOPF_Vdrop = complex(zeros(T,size(testCase.bus,1)-1)); % voltage drop
    ACOPF_f = zeros(T,1); % objective f
    ACOPF_Qg = complex(zeros(T,size(testCase.gen,1))); % Q from grid
    ACOPF_I = complex(zeros(T,size(testCase.bus,1)-1)); % line current
    ACOPF_struct_container = cell(48,1); % save all ACOPF results
    ACOPF_PVgen = complex(zeros(T,size(testCase.gen,1))); % injected PV output 
    % solve Matpower
    for t = 1 : T
        mpc = testCase; % load data
        nPV = size(mpc.gen,1)-1; % No of PV systems
        idxHH = find(mpc.bus(:,2) == 2); % find index for buses with PV systems
        % update solar PV data
        if nPV ~= 0
            mpc.gen(2:end,10) = mpc.gen(2:end,10)*solar(t); % Pmax
            mpc.gen(2:end,9) = mpc.gen(2:end,9)*solar(t); % Pmin
            mpc.gen(2:end,4) = mpc.gen(2:end,4)*solar(t)*tan(acos(PF)); % Qmax - ensures min PF
            mpc.gen(2:end,5) = mpc.gen(2:end,5)*solar(t)*(-tan(acos(PF))); % Qmin - ensures min PF
        end
        % update load data
        mpc.bus(idxHH,3) = loadHH(t,idxHH-1)'; % Pg
        mpc.bus(idxHH,4) = loadHH(t,idxHH-1)'*0.6; % Qg: 0.6*Pg - assumed
        ACOPF_struct = runopf(mpc,mpopt); % run ACOPF
        ACOPF_V(t,:) = ACOPF_struct.bus(:,8); % save V
        ACOPF_PVgen(t,:) = ACOPF_struct.gen(:,2); % save P_inj from PV
        [ACOPF_I2R(t,:), ACOPF_Vdrop(t,:), ACOPF_I(t,:)] = line_Losses2(ACOPF_struct); % save losses, voltage drop, line current
        ACOPF_Qg(t,:) = ACOPF_struct.var.val.Qg; % Q from grid
        ACOPF_f(t) = ACOPF_struct.f; % save objective f
        ACOPF_struct_container{t} = ACOPF_struct; % save all ACOPF results
        clear mpc % clear input data for next iteration
    end
    toc;
    % plot P_inj (P_av-P_c) between 1 and 2pm to check results
    for i =  26 :28
        figure(201)
        plot(2:13, ACOPF_struct_container{i}.gen(2:end,2), '*');
        hold on
    end
    ylabel('Pav [kVA]'); xlabel('Time')
    legend({'Pav 1pm','1.30pm','2pm'})
    title('ACOPF PVinj')
    set(gcf,'color','w'); grid on
end
% ACOPF.om % shows the list of variables and constraints 
end