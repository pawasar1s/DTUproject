function om = userfcn_reserves_formulation(om, mpopt, args)
%% initialize some things
define_constants;
mpc = om.get_mpc();
r = mpc.reserves;
ng = size(mpc.gen, 1); %% number of on-line gens
%% variable bounds
Rmin = zeros(ng, 1); %% bound below by 0
Rmax = r.qty; %% bound above by stated max reserve qty ...
k = find(mpc.gen(:, RAMP_10) > 0 & mpc.gen(:, RAMP_10) < Rmax);
Rmax(k) = mpc.gen(k, RAMP_10); %% ... and ramp rate
Rmax = Rmax / mpc.baseMVA;
%% constraints
I = speye(ng); %% identity matrix
Ar = [I I];
Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;
lreq = r.req / mpc.baseMVA;
%% cost
Cw = r.cost * mpc.baseMVA; %% per unit cost coefficients
%% add them to the model
om.add_var(R, ng, [], Rmin, Rmax);
om.add_lin_constraint(Pg_plus_R, Ar, [], Pmax, {Pg, R});
om.add_lin_constraint(Rreq, r.zones, lreq, [], {R});
om.add_quad_cost(Rcost, [], Cw, 0, {R});
end