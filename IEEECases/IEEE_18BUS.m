function mpc = IEEE_18BUS
%CASE9    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Joe H. Chow's book, p. 70.

%   This file is taken from MATPOWER available at http://www.pserc.cornell.edu/matpower/
%  R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas,
%  "MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education," 
%  Power Systems, IEEE Transactions on, vol. 26, no. 1, pp. 12-19, Feb. 2011.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system kVA base
mbase = 75; % [kVA]
mpc.baseMVA = mbase; % [kVA]
% based on 75kVA distribution transformer 
%% bus data
Pload = 4.9;
Qload = 3;
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
1	3	0       0       0	0	1	1	0	0.240   1	1.0	1.0; % grid
2	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
3	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
4	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
5	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
6	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
7	2	Pload   Qload	0	0	1	1	0	0.240	1	1.1	0.9;
8	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
9	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
10	2	Pload   Qload	0	0	1	1	0	0.240	1	1.1	0.9;
11	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
12	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
13	2	Pload   Qload	0	0	1	1	0	0.240	1	1.1	0.9;
14	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
15	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
16	2	Pload   Qload	0	0	1	1	0	0.240	1	1.1	0.9;
17	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
18	1	0       0       0	0	1	1	0	0.240	1	1.1	0.9;
19	2	Pload	Qload	0	0	1	1	0	0.240	1	1.1	0.9;
];
%     10 BASE_KV     baseKV, base voltage (kV)
%% generator data
PVeff = 0.77;
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	PmaxDC[kW]	PminDC[kW]	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	30	20   500	-500	1	1	mbase     1000        -1000       0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	5	-5      1	1	mbase     5.52*PVeff	5.52*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     4	0	0	5	-5      1	1	mbase     5.70*PVeff	5.70*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	5	-5      1	1	mbase     9.00*PVeff	9.00*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     7	0	0	5	-5      1	1	mbase     9.00*PVeff	9.00*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     8	0	0	5	-5      1	1	mbase     9.00*PVeff	9.00*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     10	0	0	5	-5      1	1	mbase     5.70*PVeff	5.70*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     11	0	0	5	-5      1	1	mbase     9.00*PVeff	9.00*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     13	0	0	5	-5      1	1	mbase     5.70*PVeff	5.70*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     14	0	0	5   -5      1	1	mbase    5.52*PVeff	5.52*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     16	0	0	5	-5      1	1	mbase     5.52*PVeff	5.52*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     17	0	0	5	-5      1	1	mbase     5.70*PVeff	5.70*PVeff	0	0	0	0	0	0	0	0	0	0	0;
%     19	0	0	5	-5      1	1	mbase    9.00*PVeff	9.00*PVeff	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
dropLineL = 20; % [m]
polePoleL = 50; % [m]
r_dropLine = 0.549; % [Ohms/km]
x_dropLine = 0.230 + 0.055; % [Ohms/km]
r_poleLine = 0.270; % [Ohms/km]
x_poleLine = 0.240 + 0.072; % [Ohms/km]
r_Dline = r_dropLine/1000*dropLineL; % [p.u.]
x_Dline = x_dropLine/1000*dropLineL; % [p.u.]
r_Pline = r_poleLine/1000*polePoleL; % [p.u.]
x_Pline = x_poleLine/1000*polePoleL; % [p.u.]
maxI = 320; % ACOPF works if set to 32000
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
1	3	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
2	3	r_Dline	x_Dline	    0	1*maxI	0	0	0	0	1	-360	360;
3	4	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
3	6	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
5	6	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
6	7	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
6	9	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
8	9	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
9	10	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
9	12	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
11	12	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
12	13	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
12	15	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
14	15	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
15	16	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
15	18	r_Pline	x_Pline 	0	1*maxI	0	0	0	0	1	-360	360;
17	18	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
18	19	r_Dline	x_Dline     0	1*maxI	0	0	0	0	1	-360	360;
];
%%
mpc.bus_name = {
	'Grid';
	'BUS NO2';
	'POLE NO3';
	'BUS NO4';
	'BUS NO5';
	'POLE  NO6';
	'BUS NO7';
	'BUS NO8';
	'POLE  NO9';
	'BUS NO10';
    'BUS NO11';
    'POLE  NO12';
    'BUS NO13';
    'BUS NO14';
    'POLE  NO15';
    'BUS NO16';
    'BUS NO17';
    'POLE  NO18';
    'BUS NO19';
};
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    2	0      0	3	0    0.3	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
%     2	0      0	3	0    0	 0;
]; % 80$/MWh
%% convert branch impedances from Ohms to p.u.
mpc.Vbase = mpc.bus(1,10) * 1000; % from kV to 240 Volts
mpc.Sbase = mpc.baseMVA * 1000; % from kVA to VA 
mpc.branch(:, [3, 4]) = mpc.branch(:, [3, 4]) /(mpc.Vbase^2 / mpc.Sbase);
mpc.Ibase = mpc.Sbase/mpc.Vbase; 
mpc.branch(:, 6) =  mpc.branch(:, 6) / mpc.Ibase; % in p.u. 
%% convert loads from kW to MW and kVAR to MVAR
%mpc.bus(:, [3, 4]) = mpc.bus(:, [3, 4]);
