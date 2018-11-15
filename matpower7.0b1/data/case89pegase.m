function mpc = case89pegase
%CASE89PEGASE  Power flow data for small part of European system.
%   Please see CASEFORMAT for details on the case file format.
%
%   This case accurately represents the size and complexity of part of the
%   European high voltage transmission network. The network contains 89
%   buses, 12 generators, and 210 branches and it operates at 380, 220, and
%   150 kV. Please note that the data are fictitious and do not correspond
%   to real world data. They can thus be used to validate methods and tools
%   but should not be used for operation and planning of the European grid.
%
%   The data stems from the Pan European Grid Advanced Simulation and State
%   Estimation (PEGASE) project, part of the 7th Framework Program of the
%   European Union (http://www.fp7-pegase.com/).
%
%   When publishing results based on this data, please cite:
%
%     C. Josz, S. Fliscounakis, J. Maeght, and P. Panciatici, "AC Power Flow
%     Data in MATPOWER and QCQP Format: iTesla, RTE Snapshots, and PEGASE"
%     http://arxiv.org/abs/1603.01533
%
%     S. Fliscounakis, P. Panciatici, F. Capitanescu, and L. Wehenkel,
%     "Contingency ranking with respect to overloads in very large power
%     systems taking into account uncertainty, preventive and corrective
%     actions", Power Systems, IEEE Trans. on, (28)4:4909-4917, 2013.
%     https://doi.org/10.1109/TPWRS.2013.2251015
%
%   Remarks:
%
%   1. Line flow limits are 20 MVA (at 1 p.u. voltage) greater than the
%   current flow limits found in PEGASE data.
%  
%   2. PEGASE data contains asymmetric shunt conductance and susceptance in
%   the PI transmission line model of branches. Thus total line charging
%   susceptance of branches is set to 0 p.u. and the nodal representation
%   of shunt condutance and susceptance is used. As a result, power flow
%   equations are left unchanged compared with original PEGASE data.
%   However, line flow constraints in the optimal flow problem are
%   modified.
%
%   3. Identical linear costs are used for all generators to form a loss
%   minimizing OPF objective function.
%
%   4. Since some parts of the network are aggregated, some generators
%   (e.g. with negative PMIN) represent aggregations of multiple loads
%   and generators.
%
%   Contacts:
%     Cédric Josz, Stéphane Fliscounakis, Jean Maeght, Patrick Panciatici
%     firstname.lastname@rte-france.com
%     Réseau de Transport d'Electricité (French Transmission System Operator)
%     Département Expertise Système, Immeuble "Le Colbert"
%     9 rue de la Porte de Buc, 78000 Versailles Cedex, France
%
%   March 18th, 2015

%   MATPOWER
%   Copyright (c) 2015, 2016 by Cédric Josz, Stéphane Fliscounakis, Jean Maeght,
%   and Patrick Panciatici
%   Licensed under the Creative Commons Attribution 4.0 International license,
%   http://creativecommons.org/licenses/by/4.0/

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	89	1	0	0	0.199111	4.530599	0	0.998636	-2.72828	380	2	1.1	0.9;
	228	1	-23.43	57.4	0	0	0	1.052188	-6.262599	150	2	1.1	0.9;
	271	1	96.7	26.8	0	0	0	1.039303	-7.776587	150	2	1.1	0.9;
	317	1	-0	-0	0.045946	-1.82203	0	1.052124	7.796992	380	2	1.1	0.9;
	659	1	-0	-0	0.174507	73.977052	0	1.051131	7.757217	380	2	1.1	0.9;
	792	1	295.3	41.9	0	0	0	1.034186	-0.893204	150	2	1.1	0.9;
	913	3	0	0	0	0	0	1.030951	0	380	2	1.1	0.9;
	955	1	244.2	152.7	0	0	0	1.002466	-9.134162	150	2	1.1	0.9;
	1037	1	-0	0	0	0.55	0	1.052997	6.607316	380	2	1.1	0.9;
	1163	1	-0	-0	0	16.04	0	1.009762	0.4517	380	2	1.1	0.9;
	1317	1	149.3	-2	0	20.77	0	1.038584	4.859365	380	2	1.1	0.9;
	1367	1	0	-0	0.175687	-1.509401	0	1.023381	1.554555	380	2	1.1	0.9;
	1445	1	238.4	210.9	0	0	0	1.034571	1.564741	220	2	1.1	0.9;
	1531	1	-0	0	0	0	0	1.018598	-9.561918	150	2	1.1	0.9;
	1579	1	-0	0	0	16.56	0	1.006544	-3.372817	380	2	1.1	0.9;
	1611	1	274.1	166.5	0	0.56	0	1.046304	-7.257075	150	2	1.1	0.9;
	1616	1	0	0	0	14.45	0	1.00541	-2.075628	380	2	1.1	0.9;
	1676	1	-0	0	0.171586	-1.4632	0	1.009368	0.392255	380	2	1.1	0.9;
	1815	1	63.66	9.7	0	0	0	1.057772	1.929189	220	2	1.1	0.9;
	1968	1	118.4	137	0	0	0	1.02308	-9.077175	150	2	1.1	0.9;
	2107	2	0	0	0	0	0	1.03888	1.910737	380	2	1.1	0.9;
	2154	1	-357.45	33.05	0.29	-5.6	0	1.038292	13.048352	380	2	1.1	0.9;
	2168	1	8.6	-1.48	0	0	0	1.016376	-10.987335	150	2	1.1	0.9;
	2267	2	0	0	0	0	0	1.038741	13.143978	380	2	1.1	0.9;
	2268	1	0	-0	0.183167	22.675025	0	0.99293	-3.22815	380	2	1.1	0.9;
	2299	1	402.1	110.2	0	0	0	1.041449	-9.18839	150	2	1.1	0.9;
	2441	1	0	-0	0.242534	1.143584	0	0.991005	-2.851769	380	2	1.1	0.9;
	2449	1	301.9	-103.2	0	0	0	1.086934	-3.975972	150	2	1.1	0.9;
	2520	1	-0	-0	0.057433	3.80797	0	0.999733	-2.864712	380	2	1.1	0.9;
	2870	1	-0	0	0	0.62	0	1.052997	6.607314	380	2	1.1	0.9;
	2908	1	203.8	61	0	0.6	0	1.022415	-8.627044	150	2	1.1	0.9;
	3097	1	361.91	3.66	0	9.39	0	1.048722	6.05208	380	2	1.1	0.9;
	3242	1	296.1	82.3	0	0	0	1.054504	-7.413077	150	2	1.1	0.9;
	3279	1	0	0	0	0	0	1.034258	-0.85281	150	2	1.1	0.9;
	3493	1	226.3	61.1	0	0	0	1.040629	-8.48147	150	2	1.1	0.9;
	3506	1	0	0	0	0.6	0	1.002708	-9.100433	150	2	1.1	0.9;
	3659	2	0	0	0	0	0	1.05217	3.046806	380	2	1.1	0.9;
	4014	1	0	0	0	0	0	1.033972	-11.212105	150	2	1.1	0.9;
	4423	1	-0	-0	0	0.33	0	1.041509	-9.18519	150	2	1.1	0.9;
	4427	1	296.8	77.7	0	0	0	1.036884	-8.989673	150	2	1.1	0.9;
	4495	1	336.3	104.4	0	0	0	1.038496	-8.944438	150	2	1.1	0.9;
	4586	2	0	0	0	0	0	1.047701	-8.937472	150	2	1.1	0.9;
	4665	1	390.5	62.5	0	0	0	1.038909	-9.573486	150	2	1.1	0.9;
	4929	1	37	8.7	0.070441	50.913893	0	1.052995	6.607328	380	2	1.1	0.9;
	5097	2	0	0	0	0	0	1.05444	-7.409995	150	2	1.1	0.9;
	5155	1	411.3	168	0	0	0	1.030998	-8.936324	150	2	1.1	0.9;
	5210	1	347.61	74.37	0	0	0	1.032893	1.841903	220	2	1.1	0.9;
	5416	1	0	0	0	25.13	0	1.038395	13.368566	380	2	1.1	0.9;
	5509	1	0	0	0.314057	17.296509	0	0.998445	-4.250016	380	2	1.1	0.9;
	5587	1	-0	-0	0	0.33	0	1.040763	-8.453527	150	2	1.1	0.9;
	5762	1	0	0	0	2.23	0	1.042316	2.397232	380	2	1.1	0.9;
	5776	1	222.7	80.5	0	0	0	1.074342	-5.370611	150	2	1.1	0.9;
	5848	1	0	-0	0.29	-0.42	0	1.009361	-3.019939	380	2	1.1	0.9;
	5996	1	-0	-0	0.164188	4.023827	0	1.042267	2.397462	380	2	1.1	0.9;
	6069	1	-456.66	-131.97	0.285753	35.371068	0	0.977261	-4.165509	380	2	1.1	0.9;
	6233	2	0	0	0	0	0	1.052521	7.854104	380	2	1.1	0.9;
	6293	1	-0	0	0	33.66	0	0.994291	-2.463618	380	2	1.1	0.9;
	6542	1	355.98	-147.98	0	0	0	1.070142	-0.548851	220	2	1.1	0.9;
	6704	1	-0	-0	0.211099	0.672832	0	1.003492	-1.525293	380	2	1.1	0.9;
	6798	2	0	0	0	0	0	1.053791	8.135392	380	2	1.1	0.9;
	6826	1	326.9	-82.3	0	0	0	1.076556	-5.786044	150	2	1.1	0.9;
	6833	1	0	0	0.246363	4.841277	0	0.968382	-4.882064	380	2	1.1	0.9;
	7051	1	0	0	0	21.19	0	1.014089	1.515243	380	2	1.1	0.9;
	7180	1	-0	0	0.234234	4.220723	0	0.970985	-4.970482	380	2	1.1	0.9;
	7279	2	0	0	0	0	0	1.033972	-11.212105	150	2	1.1	0.9;
	7526	1	-179.73	-63.08	0	0	0	1.015361	8.176427	380	2	1.1	0.9;
	7563	1	531.7	-12.8	0	0	0	1.040779	-11.087591	150	2	1.1	0.9;
	7637	1	0	0	0.29	6.18	0	1.035715	19.538796	380	2	1.1	0.9;
	7762	1	-0	0	0.304243	48.453399	0	1.010704	-1.708563	380	2	1.1	0.9;
	7829	1	-114.36	30.04	0	0	0	1.061047	-2.764068	220	2	1.1	0.9;
	7960	2	0	0	0	0	0	1.052291	7.889157	380	2	1.1	0.9;
	8103	1	39.34	-10.7	0	0	0	1.004552	0.141171	150	2	1.1	0.9;
	8179	1	-0	0	0.194464	19.1468	0	1.00098	-3.338785	380	2	1.1	0.9;
	8181	1	0	-0	0.210824	6.620599	0	1.002201	-2.957365	380	2	1.1	0.9;
	8229	1	0	0	0	0	0	1.074342	-5.370611	150	2	1.1	0.9;
	8329	1	0	0	0.149712	2.04987	0	1.044426	7.171947	380	2	1.1	0.9;
	8335	1	637.92	139.57	0	0	0	1.018566	-9.585336	150	2	1.1	0.9;
	8420	1	-0	0	0	0	0	1.046629	-7.203894	150	2	1.1	0.9;
	8574	1	-0	0	0.447839	20.994738	0	1.003438	-1.910009	380	2	1.1	0.9;
	8581	1	-1299.13	-140.85	0	0	0	1.039591	5.776595	380	2	1.1	0.9;
	8605	2	0	0	0	0	0	1.023822	1.617834	380	2	1.1	0.9;
	8847	1	-0	0	0.18076	-0.881251	0	1.012216	1.209962	380	2	1.1	0.9;
	8921	1	0	0	0	16.37	0	1.00666	-1.113996	380	2	1.1	0.9;
	8964	1	925.91	147.97	0	0	0	1.016233	-11.01836	150	2	1.1	0.9;
	9024	1	0	0	0.141019	25.246	0	1.057806	3.596649	380	2	1.1	0.9;
	9025	1	0	-0	0.205903	3.0768	0	1.004762	-2.614968	380	2	1.1	0.9;
	9064	1	0	-0	0	15.26	0	1.004902	-2.063811	380	2	1.1	0.9;
	9192	1	17.92	23.3	0	0	0	1.024174	-8.846156	150	2	1.1	0.9;
	9239	2	0	0	0	0	0	1.052319	7.893501	380	2	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	913	1249.4	684.16	772.11	-269.94	1.030951	100	1	1700	566.67	0	0	0	0	0	0	0	0	0	0	0;
	2107	1269.4	558.16	661.9	-257.62	1.03888	100	1	1500	500	0	0	0	0	0	0	0	0	0	0	0;
	2267	362	200.02	228.16	-78.83	1.038741	100	1	500	166.67	0	0	0	0	0	0	0	0	0	0	0;
	3659	1097.4	689.53	943.64	-282.01	1.05217	100	1	2000	666.67	0	0	0	0	0	0	0	0	0	0	0;
	4586	-545.7	41.95	92.96	-29.81	1.047701	100	1	100	-727.6	0	0	0	0	0	0	0	0	0	0	0;
	5097	21.1	2.52	8.91	-4.11	1.05444	100	1	21.23	7.08	0	0	0	0	0	0	0	0	0	0	0;
	6233	879.6	368.28	545.68	-189.97	1.052521	100	1	1200	400	0	0	0	0	0	0	0	0	0	0	0;
	6798	949.7	272.74	429.61	-183.48	1.053791	100	1	1000	333.33	0	0	0	0	0	0	0	0	0	0	0;
	7279	-681.7	5.38	14.78	-6.64	1.033972	100	1	100	-908.93	0	0	0	0	0	0	0	0	0	0	0;
	7960	419	153.39	274.2	-93.22	1.052291	100	1	600	200	0	0	0	0	0	0	0	0	0	0	0;
	8605	427	170.21	273.94	-93.65	1.023822	100	1	600	200	0	0	0	0	0	0	0	0	0	0	0;
	9239	419	153.17	274.2	-93.22	1.052319	100	1	600	200	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	3097	659	0.00073	0.00905	0	1205	0	0	0	0	1	-360	360;
	9024	4929	0.00192	0.02296	0	909	0	0	0	0	1	-360	360;
	1815	6542	0.01431	0.05899	0	0	0	0	0	0	1	-360	360;
	6542	1445	0.02607	0.07168	0	0	0	0	0	0	1	-360	360;
	6069	6833	0.00054	0.00465	0	1600	0	0	0	0	1	-360	360;
	6069	2268	0.00126	0.01124	0	1633	0	0	0	0	1	-360	360;
	6069	7180	0.00048	0.00527	0	1567	0	0	0	0	1	-360	360;
	6069	6293	0.00082	0.00938	0	1600	0	0	0	0	1	-360	360;
	1968	9192	0.00028	0.003502	0	0	0	0	0	0	1	-360	360;
	1968	3493	0.691467	4.414711	0	0	0	0	0	0	1	-360	360;
	1968	955	0.02004	0.0916	0	0	0	0	0	0	1	-360	360;
	1968	228	0.841911	4.768444	0	0	0	0	0	0	1	-360	360;
	1968	6826	0.026	0.08948	0	0	0	0	0	0	1	-360	360;
	1968	8964	0.002631	0.011889	0	319	0	0	0	0	1	-360	360;
	1968	5776	0.7232	1.944	0	0	0	0	0	0	1	-360	360;
	9192	955	0.156169	0.546756	0	0	0	0	0	0	1	-360	360;
	9192	6826	0.168809	0.580622	0	0	0	0	0	0	1	-360	360;
	9192	8964	0.019222	0.068209	0	0	0	0	0	0	1	-360	360;
	8574	1616	0.00039	0.00386	0	1435	0	0	0	0	1	-360	360;
	8574	1163	0.00096	0.01138	0	1468	0	0	0	0	1	-360	360;
	8574	9064	0.00037	0.00325	0	1435	0	0	0	0	1	-360	360;
	8574	8921	0.00029	0.00339	0	1600	0	0	0	0	1	-360	360;
	3493	2299	0.009458	0.158871	0	0	0	0	0	0	1	-360	360;
	2299	4423	0	0.000222	0	0	0	0	0	0	1	-360	360;
	6826	2299	0.312951	7.130667	0	0	0	0	0	0	1	-360	360;
	7563	2299	0.247649	1.075911	0	0	0	0	0	0	1	-360	360;
	4427	2299	0.06056	0.18612	0	0	0	0	0	0	1	-360	360;
	8335	2299	2.142489	4.347733	0	0	0	0	0	0	1	-360	360;
	5155	2299	0.07132	1.068889	0	0	0	0	0	0	1	-360	360;
	271	2299	0.02268	0.08652	0	0	0	0	0	0	1	-360	360;
	2449	2299	2.007778	4.494222	0	0	0	0	0	0	1	-360	360;
	4665	2299	0.006178	0.028938	0	0	0	0	0	0	1	-360	360;
	4495	2299	0.027711	0.260951	0	0	0	0	0	0	1	-360	360;
	3493	5587	0	0.000222	0	319	0	0	0	0	1	-360	360;
	4586	3493	0.00692	0.337378	0	0	0	0	0	0	1	-360	360;
	228	3493	0.025991	0.10048	0	0	0	0	0	0	1	-360	360;
	6826	3493	0.083711	0.8648	0	0	0	0	0	0	1	-360	360;
	7563	3493	0.447378	1.862089	0	0	0	0	0	0	1	-360	360;
	4427	3493	0.062169	0.153191	0	0	0	0	0	0	1	-360	360;
	8335	3493	2.492533	4.725778	0	0	0	0	0	0	1	-360	360;
	5155	3493	0.020271	0.08244	0	0	0	0	0	0	1	-360	360;
	271	3493	0.11416	1.777556	0	0	0	0	0	0	1	-360	360;
	5776	3493	0.574978	2.992578	0	0	0	0	0	0	1	-360	360;
	2449	3493	1.887067	4.964444	0	0	0	0	0	0	1	-360	360;
	4665	3493	0.250778	0.817422	0	0	0	0	0	0	1	-360	360;
	4495	3493	0.007049	0.029609	0	0	0	0	0	0	1	-360	360;
	5587	4586	0.003049	0.020218	0	0	0	0	0	0	1	-360	360;
	4423	4586	0.002911	0.020262	0	0	0	0	0	0	1	-360	360;
	4929	2870	4e-05	0.00041	0	777	0	0	0	0	1	-360	360;
	4929	659	0.00087	0.01087	0	1303	0	0	0	0	1	-360	360;
	4929	1037	4e-05	0.00046	0	777	0	0	0	0	1	-360	360;
	4929	659	0.00095	0.00924	0	1501	0	0	0	0	1	-360	360;
	1815	5210	0.528017	2.636157	0	0	0	0	0	0	1	-360	360;
	1815	1445	0.032031	0.128659	0	0	0	0	0	0	1	-360	360;
	1616	7762	0.00044	0.0032	0	1665	0	0	0	0	1	-360	360;
	1616	2520	0.00046	0.00489	0	1731	0	0	0	0	1	-360	360;
	9064	7762	0.0003	0.00352	0	1600	0	0	0	0	1	-360	360;
	9064	89	0.00051	0.00454	0	1665	0	0	0	0	1	-360	360;
	1317	659	0.00068	0.00746	0	1633	0	0	0	0	1	-360	360;
	1317	8605	0.00093	0.01022	0	1764	0	0	0	0	1	-360	360;
	8605	1367	4e-05	0.0004	0	1501	0	0	0	0	1	-360	360;
	8605	8921	0.00078	0.00698	0	1764	0	0	0	0	1	-360	360;
	1163	1676	4e-05	0.00043	0	1600	0	0	0	0	1	-360	360;
	1163	7051	0.00025	0.00309	0	1698	0	0	0	0	1	-360	360;
	913	7762	0.00047	0.00538	0	1633	0	0	0	0	1	-360	360;
	2107	7762	0.00046	0.00548	0	1764	0	0	0	0	1	-360	360;
	2107	5996	0.00016	0.00173	0	1863	0	0	0	0	1	-360	360;
	2107	6293	0.00166	0.01471	0	1764	0	0	0	0	1	-360	360;
	913	7762	0.00051	0.00512	0	1797	0	0	0	0	1	-360	360;
	6704	8921	0.00021	0.00248	0	1534	0	0	0	0	1	-360	360;
	228	4586	0.004391	0.021689	0	306	0	0	0	0	1	-360	360;
	2441	6293	0.00024	0.00289	0	1665	0	0	0	0	1	-360	360;
	955	3506	0	0.000222	0	358	0	0	0	0	1	-360	360;
	955	8964	0.01788	0.100471	0	0	0	0	0	0	1	-360	360;
	3506	2908	0.006458	0.039982	0	0	0	0	0	0	1	-360	360;
	228	6826	0.106151	0.993333	0	0	0	0	0	0	1	-360	360;
	228	5776	0.699778	3.232311	0	0	0	0	0	0	1	-360	360;
	8964	6826	0.257058	0.43592	0	0	0	0	0	0	1	-360	360;
	5776	6826	0.002769	0.011449	0	0	0	0	0	0	1	-360	360;
	659	8329	0.00034	0.00305	0	1600	0	0	0	0	1	-360	360;
	659	7051	0.00123	0.01255	0	1863	0	0	0	0	1	-360	360;
	659	9239	6e-05	0.00065	0	1435	0	0	0	0	1	-360	360;
	659	7960	6e-05	0.00063	0	1534	0	0	0	0	1	-360	360;
	659	6233	6e-05	0.00055	0	1764	0	0	0	0	1	-360	360;
	659	5416	0.00079	0.00968	0	1369	0	0	0	0	1	-360	360;
	659	6798	7e-05	0.00079	0	0	0	0	0	0	1	-360	360;
	792	3279	0	0.000222	0	397	0	0	0	0	1	-360	360;
	7563	792	2.064178	7.545333	0	0	0	0	0	0	1	-360	360;
	7279	792	0.153711	0.482489	0	0	0	0	0	0	1	-360	360;
	8335	792	0.015809	0.12368	0	0	0	0	0	0	1	-360	360;
	2449	792	0.057978	0.363422	0	0	0	0	0	0	1	-360	360;
	4665	792	0.6932	2.196889	0	0	0	0	0	0	1	-360	360;
	4495	792	0.560356	2.395867	0	0	0	0	0	0	1	-360	360;
	5416	2267	0.00024	0.00213	0	777	0	0	0	0	1	-360	360;
	5416	7637	0.00093	0.00884	0	1303	0	0	0	0	1	-360	360;
	8964	2168	0	0.000222	0	332	0	0	0	0	1	-360	360;
	3659	5996	0.00018	0.00204	0	1797	0	0	0	0	1	-360	360;
	5762	5996	0.00018	0.00208	0	777	0	0	0	0	1	-360	360;
	3242	5097	0	0.000222	0	0	0	0	0	0	1	-360	360;
	7563	3242	0.02068	0.093391	0	0	0	0	0	0	1	-360	360;
	4427	3242	0.777778	7.408	0	0	0	0	0	0	1	-360	360;
	3242	7279	0.00884	0.07292	0	0	0	0	0	0	1	-360	360;
	271	3242	0.131671	8.295556	0	0	0	0	0	0	1	-360	360;
	3242	1611	0.082462	0.568222	0	0	0	0	0	0	1	-360	360;
	5097	1611	0.0042	0.02624	0	0	0	0	0	0	1	-360	360;
	8181	8179	0.00031	0.00317	0	1271	0	0	0	0	1	-360	360;
	8181	7762	0.00043	0.00498	0	1501	0	0	0	0	1	-360	360;
	9025	7762	0.00042	0.00508	0	1731	0	0	0	0	1	-360	360;
	7563	4427	0.110671	0.773067	0	0	0	0	0	0	1	-360	360;
	7563	7279	0.00912	0.037609	0	0	0	0	0	0	1	-360	360;
	7563	8335	0.098049	0.598711	0	0	0	0	0	0	1	-360	360;
	7563	5155	0.019031	1.160356	0	0	0	0	0	0	1	-360	360;
	7563	271	0.024071	0.915733	0	0	0	0	0	0	1	-360	360;
	2449	7563	1.484267	4.479111	0	0	0	0	0	0	1	-360	360;
	4665	7563	0.011578	0.063271	0	0	0	0	0	0	1	-360	360;
	4495	7563	0.654889	1.839111	0	0	0	0	0	0	1	-360	360;
	7563	1611	0.022049	0.10076	0	0	0	0	0	0	1	-360	360;
	4427	7279	0.459867	2.604711	0	0	0	0	0	0	1	-360	360;
	4427	8335	2.320222	6.088	0	0	0	0	0	0	1	-360	360;
	4427	5155	0.006369	0.048089	0	0	0	0	0	0	1	-360	360;
	4427	271	0.006631	0.0432	0	0	0	0	0	0	1	-360	360;
	2449	4427	2.567244	6.767111	0	0	0	0	0	0	1	-360	360;
	4665	4427	0.309889	1.214178	0	0	0	0	0	0	1	-360	360;
	4495	4427	0.009191	0.043289	0	0	0	0	0	0	1	-360	360;
	4427	1611	0.691022	5.239556	0	0	0	0	0	0	1	-360	360;
	5210	1445	0.002669	0.01718	0	0	0	0	0	0	1	-360	360;
	8179	5509	0.00074	0.00694	0	1271	0	0	0	0	1	-360	360;
	8179	7762	0.00068	0.008	0	1567	0	0	0	0	1	-360	360;
	7279	4014	0	0.000222	0	0	0	0	0	0	1	-360	360;
	8335	7279	0.00524	0.030631	0	0	0	0	0	0	1	-360	360;
	5155	7279	0.33932	3.814667	0	0	0	0	0	0	1	-360	360;
	271	7279	0.22696	2.614933	0	0	0	0	0	0	1	-360	360;
	2449	7279	0.190671	0.719378	0	0	0	0	0	0	1	-360	360;
	4665	7279	1.050222	2.105556	0	0	0	0	0	0	1	-360	360;
	4495	7279	1.431644	4.294578	0	0	0	0	0	0	1	-360	360;
	1611	7279	0.022618	0.161089	0	0	0	0	0	0	1	-360	360;
	5509	1579	0.00088	0.00895	0	1830	0	0	0	0	1	-360	360;
	8335	1531	0	0.000222	0	319	0	0	0	0	1	-360	360;
	5155	8335	1.827778	7.768	0	0	0	0	0	0	1	-360	360;
	2449	8335	0.020849	0.153382	0	0	0	0	0	0	1	-360	360;
	4665	8335	0.156231	0.656267	0	0	0	0	0	0	1	-360	360;
	4495	8335	0.176498	0.891556	0	0	0	0	0	0	1	-360	360;
	7762	2268	0.00067	0.00777	0	1633	0	0	0	0	1	-360	360;
	5155	271	0.004578	0.061889	0	0	0	0	0	0	1	-360	360;
	5155	2908	0.01348	0.05932	0	0	0	0	0	0	1	-360	360;
	4665	5155	0.066391	2.802622	0	0	0	0	0	0	1	-360	360;
	4495	5155	0.194258	2.220444	0	0	0	0	0	0	1	-360	360;
	5155	1611	0.008151	0.072822	0	0	0	0	0	0	1	-360	360;
	271	2908	0.016502	0.055382	0	0	0	0	0	0	1	-360	360;
	4665	271	0.0806	4.633778	0	0	0	0	0	0	1	-360	360;
	4495	271	0.831556	3.729644	0	0	0	0	0	0	1	-360	360;
	271	1611	0.007951	0.074151	0	0	0	0	0	0	1	-360	360;
	1579	5848	0.0004	0.00362	0	1633	0	0	0	0	1	-360	360;
	7051	8847	0.00015	0.00177	0	1600	0	0	0	0	1	-360	360;
	5776	8229	0	0.000222	0	0	0	0	0	0	1	-360	360;
	317	6233	2e-05	0.00022	0	1698	0	0	0	0	1	-360	360;
	4665	2449	0.02864	0.172902	0	0	0	0	0	0	1	-360	360;
	2449	4495	0.017591	0.1404	0	0	0	0	0	0	1	-360	360;
	4665	4495	0.00752	0.031711	0	0	0	0	0	0	1	-360	360;
	1611	8420	0	0.000222	0	566	0	0	0	0	1	-360	360;
	9024	6542	0.000742	0.067219	0	0	0	0	0.994359	0	1	-360	360;
	9024	6542	0.000793	0.06388	0	0	0	0	0.998418	0	1	-360	360;
	6069	9192	0.000609	0.046809	0	711	0	0	0.933053	0	1	-360	360;
	6069	1968	0.00059	0.04521	0	711	0	0	0.91223	0	1	-360	360;
	7829	1968	0.519194	5.108058	0	0	0	0	0	0	1	-360	360;
	8574	2299	0.000559	0.039982	0	0	0	0	0.919523	0	1	-360	360;
	8574	5587	0.000537	0.045312	0	0	0	0	0.948262	0	1	-360	360;
	4929	1815	0.000889	0.061549	0	0	0	0	0.996867	0	1	-360	360;
	1815	792	1.467822	7.764888	0	0	0	0	0	0	1	-360	360;
	6704	4586	0.000586	0.049031	0	678	0	0	0.923407	0	1	-360	360;
	2441	3506	0.000417	0.047801	0	678	0	0	0.952413	0	1	-360	360;
	1367	6826	0.000755	0.05143	0	0	0	0	0.924009	0	1	-360	360;
	1676	228	0.000733	0.052399	0	678	0	0	0.934986	0	1	-360	360;
	7829	6826	0.00078	0.064311	0	0	0	0	0	0	1	-360	360;
	659	3279	0.00089	0.048632	0	0	0	0	0.987002	0	1	-360	360;
	5210	792	0.00855	0.08972	0	0	0	0	0	0	1	-360	360;
	1445	792	0.00223	0.04836	0	0	0	0	0	0	1	-360	360;
	7180	2168	0.000581	0.043111	0	678	0	0	0.924038	0	1	-360	360;
	6833	8964	0.000579	0.043272	0	711	0	0	0.901005	0	1	-360	360;
	3659	3242	0.000586	0.049031	0	711	0	0	0.923407	0	1	-360	360;
	8181	4427	0.000559	0.049192	0	0	0	0	0.924009	0	1	-360	360;
	9025	7563	0.000573	0.051129	0	678	0	0	0.934986	0	1	-360	360;
	2267	5210	0.000417	0.039229	0	744	0	0	0.967059	0	1	-360	360;
	5210	7279	0.230542	1.197733	0	0	0	0	0	0	1	-360	360;
	5210	8335	0.02284	0.2644	0	0	0	0	0	0	1	-360	360;
	5210	2449	0.08537	0.736798	0	0	0	0	0	0	1	-360	360;
	5210	4665	1.402768	3.640289	0	0	0	0	0	0	1	-360	360;
	5210	4495	1.050227	4.450207	0	0	0	0	0	0	1	-360	360;
	8179	7279	0.000716	0.04371	0	0	0	0	0.934986	0	1	-360	360;
	1445	7279	0.518622	2.317689	0	0	0	0	0	0	1	-360	360;
	5509	8335	0.00059	0.045662	0	0	0	0	0.955958	0	1	-360	360;
	5509	1531	0.000603	0.050959	0	0	0	0	0.968125	0	1	-360	360;
	1445	8335	0.047769	0.554178	0	0	0	0	0	0	1	-360	360;
	7762	5155	0.000591	0.045153	0	0	0	0	0.918868	0	1	-360	360;
	7762	271	0.0006	0.047778	0	711	0	0	0.948976	0	1	-360	360;
	7829	5776	0.02436	0.270151	0	0	0	0	0	0	1	-360	360;
	8329	1445	0.000802	0.05727	0	0	0	0	0.967021	0	1	-360	360;
	8329	1445	0.000823	0.056059	0	0	0	0	0.967021	0	1	-360	360;
	1445	2449	0.188099	1.291901	0	0	0	0	0	0	1	-360	360;
	1445	4665	2.466323	7.26343	0	0	0	0	0	0	1	-360	360;
	2268	2908	0.000523	0.050352	0	744	0	0	0.934622	0	1	-360	360;
	8847	5776	0.000682	0.049638	0	0	0	0	0.913034	0	1	-360	360;
	8847	8103	0.00099	0.04736	0	0	0	0	1.012146	0	1	-360	360;
	317	2449	0.000616	0.048759	0	678	0	0	0.933053	0	1	-360	360;
	7637	8581	9e-05	0.015499	0	1698	0	0	0	-0.428189	1	-360	360;
	5848	7526	9e-05	0.00963	0	1698	0	0	0	0.178581	1	-360	360;
	89	4495	0.000664	0.045609	0	0	0	0	0.924009	0	1	-360	360;
	2520	4665	0.000649	0.044862	0	711	0	0	0.933053	0	1	-360	360;
	5996	8420	0.00067	0.04326	0	744	0	0	0.923407	0	1	-360	360;
	2154	5996	9e-05	0.009901	0	1698	0	0	0	-0.153178	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
	2	0	0	3	0	1	0;
];