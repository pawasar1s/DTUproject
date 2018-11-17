function results=make_matpower_powerflow(mpc)
%Runs a power flow, and shows some figures

mpopt = mpoption('pf.alg','NR','pf.tol',1e-8,'pf.nr.max_it',100);
results=runpf(mpc,mpopt);

figure;
plot(sort(abs(results.bus(:,8))));
xlabel('Bus #');
ylabel('Voltage (p.u.)');

figure;
Sbranch_matpower=complex(results.branch(:,14),results.branch(:,15));
plot(sort(abs(Sbranch_matpower)));
xlabel('Branch #');
ylabel('Apparent power (MVA)');

figure;
plot(sort(abs(Sbranch_matpower./mpc.branch(:,6))));
xlabel('Branch #');
ylabel('Use factor: S/branch rating (p.u.)');


