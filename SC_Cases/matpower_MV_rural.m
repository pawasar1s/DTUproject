function mpc = matpower_MV_rural
%Please see CASEFORMAT for details on the case file format.
% MATPOWER Case Format : Version 2
mpc.version = '2';
%-----  Power Flow Data  -----%%
% system MVA base
mpc.baseMVA = 1;
% bus data
%bus_i type Pd	Qd Gs Bs area Vm Va baseKV zone	Vmax Vmin
mpc.bus = [
1	3	0	0	0	0	1	1	0	132	1	1.1	0.9
2	1	0.2336	0.07008	0	0	1	1	0	20	1	1.1	0.9
3	1	0.168	0.0504	0	0	1	1	0	20	1	1.1	0.9
4	1	0.072	0.0216	0	0	1	1	0	20	1	1.1	0.9
5	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
6	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
7	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
8	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
9	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
10	1	0.0628	0.01884	0	0	1	1	0	20	1	1.1	0.9
11	1	0.1456	0.04368	0	0	1	1	0	20	1	1.1	0.9
12	1	0.282	0.0846	0	0	1	1	0	20	1	1.1	0.9
13	1	0.2524	0.07572	0	0	1	1	0	20	1	1.1	0.9
14	1	0.1648	0.04944	0	0	1	1	0	20	1	1.1	0.9
15	1	0.1948	0.05844	0	0	1	1	0	20	1	1.1	0.9
16	1	0.2852	0.08556	0	0	1	1	0	20	1	1.1	0.9
17	1	0.1112	0.03336	0	0	1	1	0	20	1	1.1	0.9
18	1	0.0544	0.01632	0	0	1	1	0	20	1	1.1	0.9
19	1	0.3184	0.09552	0	0	1	1	0	20	1	1.1	0.9
20	1	0.208	0.0624	0	0	1	1	0	20	1	1.1	0.9
21	1	0.2968	0.08904	0	0	1	1	0	20	1	1.1	0.9
22	1	0.322	0.0966	0	0	1	1	0	20	1	1.1	0.9
23	1	0.3248	0.09744	0	0	1	1	0	20	1	1.1	0.9
24	1	0.1608	0.04824	0	0	1	1	0	20	1	1.1	0.9
25	1	0.1048	0.03144	0	0	1	1	0	20	1	1.1	0.9
26	1	0.3252	0.09756	0	0	1	1	0	20	1	1.1	0.9
27	1	0.2088	0.06264	0	0	1	1	0	20	1	1.1	0.9
28	1	0.3388	0.10164	0	0	1	1	0	20	1	1.1	0.9
29	1	0.332	0.0996	0	0	1	1	0	20	1	1.1	0.9
30	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
31	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
32	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
33	1	0.2	0.06	0	0	1	1	0	20	1	1.1	0.9
34	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
35	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
36	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
37	1	0.0656	0.01968	0	0	1	1	0	20	1	1.1	0.9
38	1	0.0608	0.01824	0	0	1	1	0	20	1	1.1	0.9
39	1	0.1944	0.05832	0	0	1	1	0	20	1	1.1	0.9
40	1	0.1308	0.03924	0	0	1	1	0	20	1	1.1	0.9
41	1	0.07	0.021	0	0	1	1	0	20	1	1.1	0.9
42	1	0.0536	0.01608	0	0	1	1	0	20	1	1.1	0.9
43	1	0.126	0.0378	0	0	1	1	0	20	1	1.1	0.9
44	1	0.2436	0.07308	0	0	1	1	0	20	1	1.1	0.9
45	1	0.2868	0.08604	0	0	1	1	0	20	1	1.1	0.9
46	1	0.3372	0.10116	0	0	1	1	0	20	1	1.1	0.9
47	1	0.1088	0.03264	0	0	1	1	0	20	1	1.1	0.9
48	1	0.0716	0.02148	0	0	1	1	0	20	1	1.1	0.9
49	1	0.2596	0.07788	0	0	1	1	0	20	1	1.1	0.9
50	1	0.0848	0.02544	0	0	1	1	0	20	1	1.1	0.9
51	1	0.0704	0.02112	0	0	1	1	0	20	1	1.1	0.9
52	1	0.1728	0.05184	0	0	1	1	0	20	1	1.1	0.9
53	1	0.1748	0.05244	0	0	1	1	0	20	1	1.1	0.9
54	1	0.3064	0.09192	0	0	1	1	0	20	1	1.1	0.9
55	1	0.2388	0.07164	0	0	1	1	0	20	1	1.1	0.9
56	1	0.2392	0.07176	0	0	1	1	0	20	1	1.1	0.9
57	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
58	1	0.1	0.03	0	0	1	1	0	20	1	1.1	0.9
59	1	0.2912	0.08736	0	0	1	1	0	20	1	1.1	0.9
60	1	0.2004	0.06012	0	0	1	1	0	20	1	1.1	0.9
61	1	0.1556	0.04668	0	0	1	1	0	20	1	1.1	0.9
62	1	0.1136	0.03408	0	0	1	1	0	20	1	1.1	0.9
63	1	0.2532	0.07596	0	0	1	1	0	20	1	1.1	0.9
64	1	0.3352	0.10056	0	0	1	1	0	20	1	1.1	0.9
65	1	0.0376	0.01128	0	0	1	1	0	20	1	1.1	0.9
66	1	0.196	0.0588	0	0	1	1	0	20	1	1.1	0.9
67	1	0.2636	0.07908	0	0	1	1	0	20	1	1.1	0.9
68	1	0.1792	0.05376	0	0	1	1	0	20	1	1.1	0.9
69	1	0.0692	0.02076	0	0	1	1	0	20	1	1.1	0.9
70	1	0.2672	0.08016	0	0	1	1	0	20	1	1.1	0.9
71	1	0.2264	0.06792	0	0	1	1	0	20	1	1.1	0.9
72	1	0.328	0.0984	0	0	1	1	0	20	1	1.1	0.9
73	1	0.0764	0.02292	0	0	1	1	0	20	1	1.1	0.9
74	1	0.0392	0.01176	0	0	1	1	0	20	1	1.1	0.9
75	1	0.11	0.033	0	0	1	1	0	20	1	1.1	0.9
76	1	0.1196	0.03588	0	0	1	1	0	20	1	1.1	0.9
77	1	0.076	0.0228	0	0	1	1	0	20	1	1.1	0.9
78	1	0.0424	0.01272	0	0	1	1	0	20	1	1.1	0.9
79	1	0.114	0.0342	0	0	1	1	0	20	1	1.1	0.9
80	1	0.0516	0.01548	0	0	1	1	0	20	1	1.1	0.9
81	1	0.0316	0.00948	0	0	1	1	0	20	1	1.1	0.9
82	1	0.0876	0.02628	0	0	1	1	0	20	1	1.1	0.9
83	1	0.166	0.0498	0	0	1	1	0	20	1	1.1	0.9
84	1	0.1196	0.03588	0	0	1	1	0	20	1	1.1	0.9
85	1	0.2872	0.08616	0	0	1	1	0	20	1	1.1	0.9
86	1	0.078	0.0234	0	0	1	1	0	20	1	1.1	0.9
87	1	0.0324	0.00972	0	0	1	1	0	20	1	1.1	0.9
88	1	0.2952	0.08856	0	0	1	1	0	20	1	1.1	0.9
89	1	0.0884	0.02652	0	0	1	1	0	20	1	1.1	0.9
90	1	0.1296	0.03888	0	0	1	1	0	20	1	1.1	0.9
91	1	0.19	0.057	0	0	1	1	0	20	1	1.1	0.9
92	1	0.3112	0.09336	0	0	1	1	0	20	1	1.1	0.9
93	1	0.1704	0.05112	0	0	1	1	0	20	1	1.1	0.9
94	1	0.3352	0.10056	0	0	1	1	0	20	1	1.1	0.9
95	1	0.1244	0.03732	0	0	1	1	0	20	1	1.1	0.9
96	1	0.0392	0.01176	0	0	1	1	0	20	1	1.1	0.9
97	1	0.2324	0.06972	0	0	1	1	0	20	1	1.1	0.9
98	1	0.1444	0.04332	0	0	1	1	0	20	1	1.1	0.9
99	1	0.2992	0.08976	0	0	1	1	0	20	1	1.1	0.9
100	1	0.2228	0.06684	0	0	1	1	0	20	1	1.1	0.9
101	1	0.2852	0.08556	0	0	1	1	0	20	1	1.1	0.9
102	1	0.0628	0.01884	0	0	1	1	0	20	1	1.1	0.9
103	1	0.1056	0.03168	0	0	1	1	0	20	1	1.1	0.9
104	1	0.1116	0.03348	0	0	1	1	0	20	1	1.1	0.9
105	1	0.2028	0.06084	0	0	1	1	0	20	1	1.1	0.9
106	1	0.1176	0.03528	0	0	1	1	0	20	1	1.1	0.9
107	1	0.056	0.0168	0	0	1	1	0	20	1	1.1	0.9
108	1	0.1776	0.05328	0	0	1	1	0	20	1	1.1	0.9
109	1	0.158	0.0474	0	0	1	1	0	20	1	1.1	0.9
110	1	0.2296	0.06888	0	0	1	1	0	20	1	1.1	0.9
111	1	0.2392	0.07176	0	0	1	1	0	20	1	1.1	0.9
112	1	0.2072	0.06216	0	0	1	1	0	20	1	1.1	0.9
113	1	0.2416	0.07248	0	0	1	1	0	20	1	1.1	0.9
114	1	0.082	0.0246	0	0	1	1	0	20	1	1.1	0.9
115	1	0.306	0.0918	0	0	1	1	0	20	1	1.1	0.9
116	1	0.0228	0.00684	0	0	1	1	0	20	1	1.1	0.9
117	1	0	0	0	0	1	1	0	20	1	1.1	0.9
];
%% generator data
% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
mpc.gen = [
1 0 0 1000 -1000 1.0 1 1 1000 0 0 0 0 0 0 0 0 0 0 0 0];
%% branch data
%fbus tbus	r x b rA rB rC rat ang status angmin angmax
mpc.branch = [
6	47	7.42E-05	6.89E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
5	59	0.000151868	0.00014102	0	8.66025	8.66025	8.66025	0	0	1	-360	360
57	117	0.00657775	0.006108	0	8.66025	8.66025	8.66025	0	0	1	-360	360
58	109	0.00105785	0.0009823	0	8.66025	8.66025	8.66025	0	0	1	-360	360
34	95	7.42E-05	6.89E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
35	39	5.98E-05	5.56E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
36	80	0.00105883	0.0009832	0	8.66025	8.66025	8.66025	0	0	1	-360	360
30	110	0.000109525	0.000101703	0	8.66025	8.66025	8.66025	0	0	1	-360	360
31	63	0.000110603	0.000102703	0	8.66025	8.66025	8.66025	0	0	1	-360	360
32	97	0.00022227	0.000206393	0	8.66025	8.66025	8.66025	0	0	1	-360	360
33	57	0.000173067	0.000160705	0	8.66025	8.66025	8.66025	0	0	1	-360	360
7	108	0.000123815	0.000114973	0	8.66025	8.66025	8.66025	0	0	1	-360	360
8	117	0.00198905	0.00184697	0	8.66025	8.66025	8.66025	0	0	1	-360	360
9	83	0.0003606	0.00033485	0	8.66025	8.66025	8.66025	0	0	1	-360	360
2	39	0.000140848	0.000130785	0	8.66025	8.66025	8.66025	0	0	1	-360	360
3	35	0.00019955	0.000185298	0	8.66025	8.66025	8.66025	0	0	1	-360	360
4	48	0.00176148	0.00163565	0	8.66025	8.66025	8.66025	0	0	1	-360	360
37	117	0.00362325	0.0033645	0	8.66025	8.66025	8.66025	0	0	1	-360	360
38	51	0.000824975	0.00076605	0	8.66025	8.66025	8.66025	0	0	1	-360	360
39	41	9.80E-05	9.10E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
40	37	0.00146763	0.0013628	0	8.66025	8.66025	8.66025	0	0	1	-360	360
41	6	6.05E-05	5.62E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
42	40	0.00132213	0.0012277	0	8.66025	8.66025	8.66025	0	0	1	-360	360
43	4	0.00150473	0.00139725	0	8.66025	8.66025	8.66025	0	0	1	-360	360
44	49	0.0005546	0.000515	0	8.66025	8.66025	8.66025	0	0	1	-360	360
45	35	0.000104068	9.66E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
46	50	0.000129925	0.000120645	0	8.66025	8.66025	8.66025	0	0	1	-360	360
47	50	0.000124832	0.000115915	0	8.66025	8.66025	8.66025	0	0	1	-360	360
48	61	0.00167475	0.00155513	0	8.66025	8.66025	8.66025	0	0	1	-360	360
49	45	0.000190228	0.00017664	0	8.66025	8.66025	8.66025	0	0	1	-360	360
50	52	0.000162945	0.000151305	0	8.66025	8.66025	8.66025	0	0	1	-360	360
51	75	0.00140238	0.0013022	0	8.66025	8.66025	8.66025	0	0	1	-360	360
52	117	0.00496	0.00460575	0	8.66025	8.66025	8.66025	0	0	1	-360	360
53	10	0.000128925	0.000119718	0	8.66025	8.66025	8.66025	0	0	1	-360	360
54	47	0.000198502	0.000184325	0	8.66025	8.66025	8.66025	0	0	1	-360	360
55	52	0.00012576	0.000116778	0	8.66025	8.66025	8.66025	0	0	1	-360	360
56	17	0.000133675	0.000124128	0	8.66025	8.66025	8.66025	0	0	1	-360	360
10	56	0.000172013	0.000159725	0	8.66025	8.66025	8.66025	0	0	1	-360	360
11	16	0.000138643	0.00012874	0	8.66025	8.66025	8.66025	0	0	1	-360	360
12	11	0.000113712	0.00010559	0	8.66025	8.66025	8.66025	0	0	1	-360	360
13	10	0.000149187	0.00013853	0	8.66025	8.66025	8.66025	0	0	1	-360	360
14	13	0.000118453	0.000109993	0	8.66025	8.66025	8.66025	0	0	1	-360	360
15	22	0.000161757	0.000150203	0	8.66025	8.66025	8.66025	0	0	1	-360	360
16	15	0.000113185	0.0001051	0	8.66025	8.66025	8.66025	0	0	1	-360	360
17	19	0.000132618	0.000123145	0	8.66025	8.66025	8.66025	0	0	1	-360	360
18	43	0.0007365	0.0006839	0	8.66025	8.66025	8.66025	0	0	1	-360	360
19	20	0.000159652	0.00014825	0	8.66025	8.66025	8.66025	0	0	1	-360	360
20	16	9.58E-05	8.90E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
21	19	0.00013532	0.000125655	0	8.66025	8.66025	8.66025	0	0	1	-360	360
22	25	0.00015368	0.000142705	0	8.66025	8.66025	8.66025	0	0	1	-360	360
23	14	9.89E-05	9.18E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
24	26	0.000198055	0.000183908	0	8.66025	8.66025	8.66025	0	0	1	-360	360
25	62	0.000833475	0.00077395	0	8.66025	8.66025	8.66025	0	0	1	-360	360
26	27	0.000145278	0.0001349	0	8.66025	8.66025	8.66025	0	0	1	-360	360
27	60	0.00019337	0.000179557	0	8.66025	8.66025	8.66025	0	0	1	-360	360
28	26	0.000141605	0.00013149	0	8.66025	8.66025	8.66025	0	0	1	-360	360
29	31	0.000189768	0.000176213	0	8.66025	8.66025	8.66025	0	0	1	-360	360
59	27	0.000182353	0.000169327	0	8.66025	8.66025	8.66025	0	0	1	-360	360
60	18	0.000445425	0.0004136	0	8.66025	8.66025	8.66025	0	0	1	-360	360
61	8	0.000482875	0.000448375	0	8.66025	8.66025	8.66025	0	0	1	-360	360
62	117	0.003863	0.00358725	0	8.66025	8.66025	8.66025	0	0	1	-360	360
63	70	0.00013014	0.000120843	0	8.66025	8.66025	8.66025	0	0	1	-360	360
64	5	3.67E-05	3.41E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
65	117	0.00427	0.003965	0	8.66025	8.66025	8.66025	0	0	1	-360	360
66	63	0.000144945	0.00013459	0	8.66025	8.66025	8.66025	0	0	1	-360	360
67	64	0.00012897	0.000119758	0	8.66025	8.66025	8.66025	0	0	1	-360	360
68	57	0.00017983	0.000166985	0	8.66025	8.66025	8.66025	0	0	1	-360	360
69	73	5.77E-05	5.36E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
70	57	9.29E-05	8.63E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
71	33	5.77E-05	5.36E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
72	69	0.0001173	0.00010892	0	8.66025	8.66025	8.66025	0	0	1	-360	360
73	68	9.07E-05	8.42E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
74	68	0.000557375	0.00051755	0	8.66025	8.66025	8.66025	0	0	1	-360	360
75	72	0.0004038	0.000374975	0	8.66025	8.66025	8.66025	0	0	1	-360	360
76	117	0.0009739	0.00090435	0	8.66025	8.66025	8.66025	0	0	1	-360	360
77	74	0.0014885	0.00138218	0	8.66025	8.66025	8.66025	0	0	1	-360	360
78	80	0.0014083	0.0013077	0	8.66025	8.66025	8.66025	0	0	1	-360	360
79	81	0.00119422	0.00110893	0	8.66025	8.66025	8.66025	0	0	1	-360	360
80	76	0.00113108	0.0010503	0	8.66025	8.66025	8.66025	0	0	1	-360	360
81	103	0.00167903	0.0015591	0	8.66025	8.66025	8.66025	0	0	1	-360	360
82	117	0.00526125	0.0048855	0	8.66025	8.66025	8.66025	0	0	1	-360	360
83	36	0.000545625	0.00050665	0	8.66025	8.66025	8.66025	0	0	1	-360	360
84	97	0.00107868	0.00100162	0	8.66025	8.66025	8.66025	0	0	1	-360	360
85	88	0.000136708	0.000126943	0	8.66025	8.66025	8.66025	0	0	1	-360	360
86	84	0.00109015	0.00101228	0	8.66025	8.66025	8.66025	0	0	1	-360	360
87	82	0.00200335	0.00186025	0	8.66025	8.66025	8.66025	0	0	1	-360	360
88	92	0.00019118	0.000177525	0	8.66025	8.66025	8.66025	0	0	1	-360	360
89	82	0.0025825	0.00239813	0	8.66025	8.66025	8.66025	0	0	1	-360	360
90	85	0.000282975	0.00026275	0	8.66025	8.66025	8.66025	0	0	1	-360	360
91	9	0.000284825	0.000264475	0	8.66025	8.66025	8.66025	0	0	1	-360	360
92	93	0.00011109	0.000103155	0	8.66025	8.66025	8.66025	0	0	1	-360	360
93	91	8.93E-05	8.30E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
94	34	0.000132178	0.000122738	0	8.66025	8.66025	8.66025	0	0	1	-360	360
95	99	0.000131778	0.000122365	0	8.66025	8.66025	8.66025	0	0	1	-360	360
96	82	0.00178085	0.00165365	0	8.66025	8.66025	8.66025	0	0	1	-360	360
97	91	0.000115097	0.000106875	0	8.66025	8.66025	8.66025	0	0	1	-360	360
98	100	0.000387975	0.00036025	0	8.66025	8.66025	8.66025	0	0	1	-360	360
99	32	0.00015667	0.00014548	0	8.66025	8.66025	8.66025	0	0	1	-360	360
100	99	0.00019607	0.000182065	0	8.66025	8.66025	8.66025	0	0	1	-360	360
101	99	0.00017229	0.000159985	0	8.66025	8.66025	8.66025	0	0	1	-360	360
102	100	0.000121432	0.00011276	0	8.66025	8.66025	8.66025	0	0	1	-360	360
103	58	0.0005339	0.000495775	0	8.66025	8.66025	8.66025	0	0	1	-360	360
104	106	0.00111395	0.00103437	0	8.66025	8.66025	8.66025	0	0	1	-360	360
105	102	0.00093225	0.00086565	0	8.66025	8.66025	8.66025	0	0	1	-360	360
106	105	0.0006346	0.000589275	0	8.66025	8.66025	8.66025	0	0	1	-360	360
107	106	0.0014379	0.00133518	0	8.66025	8.66025	8.66025	0	0	1	-360	360
108	107	0.00089445	0.00083055	0	8.66025	8.66025	8.66025	0	0	1	-360	360
109	7	0.000442625	0.000411025	0	8.66025	8.66025	8.66025	0	0	1	-360	360
110	7	5.92E-05	5.50E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
111	112	0.000117647	0.000109245	0	8.66025	8.66025	8.66025	0	0	1	-360	360
112	30	0.000185383	0.000172143	0	8.66025	8.66025	8.66025	0	0	1	-360	360
113	30	6.13E-05	5.69E-05	0	8.66025	8.66025	8.66025	0	0	1	-360	360
114	115	0.00023271	0.000216087	0	8.66025	8.66025	8.66025	0	0	1	-360	360
115	113	0.000176723	0.0001641	0	8.66025	8.66025	8.66025	0	0	1	-360	360
116	89	0.00252475	0.0023443	0	8.66025	8.66025	8.66025	0	0	1	-360	360
31	59	0.000291175	0.000270375	0	8.66025	8.66025	8.66025	0	0	0	-360	360
42	44	0.0020407	0.00189493	0	8.66025	8.66025	8.66025	0	0	0	-360	360
55	56	0.000211033	0.000195957	0	8.66025	8.66025	8.66025	0	0	0	-360	360
1	117	3.75E-05	0.00125	0	80	80	80	0.963068	0	1	-360	360
];
