function mpc = IEEE_18BUS_Radial_HHs
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
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
1	3	0   0	0	0	1	1	0	0.415	1	1	1; % grid
2	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
3	1	0   0	0	0	1	1	0	0.415	1	1.1	0.9;
4	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
5	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
6	1	0	0	0	0	1	1	0	0.415	1	1.1	0.9;
7	2	3   0.8	0	0	1	1	0	0.415	1	1.1	0.9;
8	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
9	1	0   0	0	0	1	1	0	0.415	1	1.1	0.9;
10	2	3   0.8	0	0	1	1	0	0.415	1	1.1	0.9;
11	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
12	1	0   0	0	0	1	1	0	0.415	1	1.1	0.9;
13	2	3   0.8	0	0	1	1	0	0.415	1	1.1	0.9;
14	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
15	1	0   0	0	0	1	1	0	0.415	1	1.1	0.9;
16	2	3   0.8	0	0	1	1	0	0.415	1	1.1	0.9;
17	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
18	1	0   0	0	0	1	1	0	0.415	1	1.1	0.9;
19	2	3	0.8	0	0	1	1	0	0.415	1	1.1	0.9;
];
%     10 BASE_KV     baseKV, base voltage (kV)
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0   20.000	-20.000	100	1	1   100     -100    0	0	0	0	0	0	0	0	0	0	0;
%     2	0	0	0.425	-0.425	100	1	1	4.25	4.25	0	0	0	0	0	0	0	0	0	0	0;
%     4	0	0	0.439	-0.439	100	1	1	4.39	4.39	0	0	0	0	0	0	0	0	0	0	0;
%     5	0	0	0.693	-0.693	100	1	1	6.93	6.93	0	0	0	0	0	0	0	0	0	0	0;
%     7	0	0	0.693	-0.693	100	1	1	6.93	6.93	0	0	0	0	0	0	0	0	0	0	0;
%     8	0	0	0.693	-0.693	100	1	1	6.93	6.93	0	0	0	0	0	0	0	0	0	0	0;
%     10	0	0	0.439	-0.439	100	1	1	4.39	4.39	0	0	0	0	0	0	0	0	0	0	0;
%     11	0	0	0.693	-0.693	100	1	1	6.93	6.93	0	0	0	0	0	0	0	0	0	0	0;
%     13	0	0	0.439	-0.439	100	1	1	4.39	4.39	0	0	0	0	0	0	0	0	0	0	0;
%     14	0	0	0.425   -0.425	100	1	1	4.25	4.25	0	0	0	0	0	0	0	0	0	0	0;
%     16	0	0	0.425	-0.425	100	1	1	4.25	4.25	0	0	0	0	0	0	0	0	0	0	0;
%     17	0	0	0.439	-0.439	100	1	1	4.39	4.39	0	0	0	0	0	0	0	0	0	0	0;
%     19	0	0	0.693	-0.693	100	1	1	6.93	6.93	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
1	3	13.5	1.56	0	0	0	0	0	0	1	-360	360;
3	2	10.98	5.7	    0	0	0	0	0	0	1	-360	360;
3	4	10.98	5.7     0	0	0	0	0	0	1	-360	360;
3	6	13.5	1.56	0	0	0	0	0	0	1	-360	360;
6	5	10.98	5.7     0	0	0	0	0	0	1	-360	360;
6	7	10.98	5.7     0	0	0	0	0	0	1	-360	360;
6	9	13.5	1.56	0	0	0	0	0	0	1	-360	360;
9	8	10.98	5.7     0	0	0	0	0	0	1	-360	360;
9	10	10.98	5.7     0	0	0	0	0	0	1	-360	360;
9	12	13.5	1.56	0	0	0	0	0	0	1	-360	360;
12	11	10.98	5.7     0	0	0	0	0	0	1	-360	360;
12	13	10.98	5.7     0	0	0	0	0	0	1	-360	360;
12	15	13.5	1.56	0	0	0	0	0	0	1	-360	360;
15	14	10.98	5.7     0	0	0	0	0	0	1	-360	360;
15	16	10.98	5.7     0	0	0	0	0	0	1	-360	360;
15	18	13.5	1.56	0	0	0	0	0	0	1	-360	360;
18	17	10.98	5.7     0	0	0	0	0	0	1	-360	360;
18	19	10.98	5.7     0	0	0	0	0	0	1	-360	360;
];
%%
mpc.bus_name = {
	'Grid';
	'BUS NO1';
	'BUS NO2';
	'BUS NO3';
	'BUS NO4';
	'BUS NO5';
	'BUS NO6';
	'BUS NO7';
	'BUS NO8';
	'BUS NO10';
    'BUS NO11';
    'BUS NO12';
    'BUS NO13';
    'BUS NO14';
    'BUS NO15';
    'BUS NO16';
    'BUS NO17';
    'BUS NO18';
};
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    2	0      0	3	0    250	 0;
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
mpc.Vbase = mpc.bus(1,10) * 1e6; % in Volts 415
mpc.Sbase = mpc.baseMVA * 1e6; % in VA 10e6
mpc.branch(:, [3, 4]) = mpc.branch(:, [3, 4]) /(mpc.Vbase^2 / mpc.Sbase) *10
%% convert loads from kW to MW and kVAR to MVAR
mpc.bus(:, [3, 4]) = mpc.bus(:, [3, 4]);
