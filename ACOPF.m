function [ACOPF_struct, ACOPF_V, ACOPF_I2R, ACOPF_f, ACOPF_Qg] = ACOPF(testCase, T, solar, loadHH, ACOPFsolver, multiPer, per, PF);
%% ACOPF
% choose solver
if ACOPFsolver == 1
    mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); % ACOPF
else
    mpopt = mpoption('model','AC', 'pf.alg', 'ISUM', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT'); % Backward/Forward Sweep with current summation
end
% run ACOPF
if multiPer == 0
    tic;
    mpc = testCase;
    nPV = size(mpc.gen,1)-1; % No of PV systems
    idxHH = find(mpc.bus(:,2) == 2); % idx for buses with PV systems
    % update solar PV data
    if nPV ~= 0
        mpc.gen(2:end,9) = mpc.gen(2:end,9)*solar(per); % Pmax
        %mpc.gen(2:end,10) = mpc.gen(2:end,10)*solarHH(per); % Pmin
        mpc.gen(2:end,4) = mpc.gen(2:end,4)*solar(per)*tan(acos(PF)); % Qmax
        mpc.gen(2:end,5) = mpc.gen(2:end,5)*solar(per)*tan(acos(PF)); % Qmin
    end
    % update load data
    mpc.bus(idxHH,3) = loadHH(per,idxHH-1)'; % Pg 
    mpc.bus(idxHH,4) = loadHH(per,idxHH-1)'*0.6; % Qg - 0.6 as of Pg - assumed
    % save results
    ACOPF_struct = runopf(mpc,mpopt); %
    ACOPF_V = ACOPF_struct.bus(:,8);
    ACOPF_I2R= get_losses(ACOPF_struct);
    ACOPF_f = ACOPF_struct.f;
    ACOPF_Qg = ACOPF_struct.var.val.Qg;
    %ACOPF_PgQg = ACOPF_struct.bus(1,[2 3]);
    clear mpc;
    toc;
elseif multiPer == 1
    tic;
    % allocate empty matrices
    ACOPF_V = complex(zeros(T,size(testCase.bus,1))); % Bus Voltage Vector
    ACOPF_I2R = complex(zeros(T,size(testCase.bus,1)-1)); % Line Loss Vector
    ACOPF_f = zeros(T,1); % Objective Function Vector
    ACOPF_Qg = complex(zeros(T,size(testCase.gen,1))); % Bus Voltage Vector
    %ACOPF_I = complex(zeros(size(mpc.bus,1)-1,size(loadHH,1))); % Line Current
    % solve Matpower
    for t = 1 : T
        mpc = testCase; % AC OPF
        nPV = size(mpc.gen,1)-1; % No of PV systems
        idxHH = find(mpc.bus(:,2) == 2); % idx for buses with PV
        % update solar PV data
        if nPV ~= 0;
            mpc.gen(2:end,10) = mpc.gen(2:end,10)*solar(t);
            mpc.gen(2:end,9) = mpc.gen(2:end,9)*solar(t);
            mpc.gen(2:end,4) = mpc.gen(2:end,4)*solar(t);
            mpc.gen(2:end,5) = mpc.gen(2:end,5)*solar(t);
        end
        % update load data
        mpc.bus(idxHH,3) = loadHH(t,idxHH-1)'; % Pg
        mpc.bus(idxHH,4) = loadHH(t,idxHH-1)'*0.6; % Qg
        ACOPF_struct = runopf(mpc,mpopt); %
        ACOPF_V(t,:) = ACOPF_struct.bus(:,8);
        ACOPF_I2R(t,:) = get_losses(ACOPF_struct);
        ACOPF_Qg(t,:) = ACOPF_struct.var.val.Pg; ACOPF_struct.var.val.Qg;
        ACOPF_f(t) = ACOPF_struct.f;
        %ACOPF_I(:,l) = linesCurrent(mpc, ACOPF_V);
        clear mpc
    end
    toc;
end
% ACOPF.om % shows the list of variables and constraints 
end