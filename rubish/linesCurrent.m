function [ACOPF_I] = LineLoss(mpc, ACOPF_V)
%% Inputs
bus = mpc.bus;
baseMVA = mpc.baseMVA;
branch = mpc.branch; 
% if isstruct(baseMVA)
%     mpc = baseMVA;
%     [baseMVA, bus, branch] = deal(mpc.baseMVA, mpc.bus, mpc.branch);
% end
% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% create map of external bus numbers to bus indices
i2e = bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(bus, 1))';
out = find(branch(:, BR_STATUS) == 0);          %% out-of-service branches
% sizes of things
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches

% construct complex bus voltage vector
%Vfull = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
%V = ACOPF_V;
% parameters
Cf = sparse(1:nl, e2i(branch(:, F_BUS)), branch(:, BR_STATUS), nl, nb);
Ct = sparse(1:nl, e2i(branch(:, T_BUS)), branch(:, BR_STATUS), nl, nb);
tap = ones(nl, 1);                              %% default tap ratio = 1 for lines
xfmr = find(branch(:, TAP));                    %% indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  %% include transformer tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
Aa = spdiags(1 ./ tap, 0, nl, nl) * Cf - Ct;   
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X)); 
Vdrop = Aa * ACOPF_V;      %% vector of voltage drop across series impedance element
ACOPF_I = baseMVA * Ysc .* Vdrop;
%lineLoss = baseMVA * Ysc .* Vdrop .* conj(Vdrop);
end