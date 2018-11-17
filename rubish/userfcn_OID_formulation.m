function om = userfcn_OID_formulation(om, mpopt, args)
%% initialize some things
%define_constants;
mpc = om.get_mpc();
%r = mpc.reserves;
nGen = size(mpc.gen, 1); %% number of on-line gens
%% variable bounds
%Rmin = zeros(ng, 1); %% bound below by 0
%Rmax = r.qty; %% bound above by stated max reserve qty ...
%k = find(mpc.gen(:, RAMP_10) > 0 & mpc.gen(:, RAMP_10) < Rmax);
%Rmax(k) = mpc.gen(k, RAMP_10); %% ... and ramp rate
%Rmax = Rmax / mpc.baseMVA;
%% constraints
I = speye(nGen); %% identity matrix
Ar = [I I];
PF = tan(acos(0.8));
QinjMax = PF*mpc.gen(:,4);
QabsMax = -PF*mpc.gen(:,4);
%Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;

%lreq = r.req / mpc.baseMVA;
%% cost
%Cw = r.cost * mpc.baseMVA; %% per unit cost coefficients
%% add them to the model
%om.add_var('R', ng, [], Rmin, Rmax);
om.add_lin_constraint('QinjLim', Ar, QabsMax, QinjMax, {'Qg'});
%om.add_lin_constraint('Rreq', r.zones, lreq, [], {'R'});
%om.add_quad_cost('Rcost', [], Cw, 0, {'R'});
end