function model = model35
%%-----  Topo Data  -----%%																																																																																																																																																																																																														
% pipe	
model.pipe = [
    % from_node	to_node	length(m)	diameter(m)	roughness	conductivity m
    1	33	2	566.72	0.25	0.0004	0.321	66.04612219
    2	2	3	214.5	0.08	0.0004	0.21	8.812506012
    3	2	4	112.2	0.1	0.0004	0.21	11.92041343
    4	2	5	130.9	0.2	0.0004	0.327	45.31320276
    5	5	6	596.86	0.1	0.0004	0.189	8.856657287
    6	5	7	517.88	0.08	0.0004	0.236	9.079659871
    7	7	8	390.06	0.08	0.0004	0.21	8.86911535
    8	7	9	226.16	0.08	0.0004	0.21	8.849097137
    9	7	10	544.94	0.08	0.0004	0.21	8.888014674
    10	5	11	353.76	0.18	0.0004	0.327	27.3768856
    11	11	12	284.02	0.08	0.0004	0.21	8.838417868
    12	11	13	409.42	0.2	0.0004	0.327	37.53915068
    13	13	14	299.64	0.18	0.0004	0.278	37.53915068
    14	14	15	191.96	0.1	0.0004	0.219	13.36451188
    15	15	16	256.96	0.08	0.0004	0.189	6.679883146
    16	15	17	300.08	0.08	0.0004	0.189	6.68462873
    17	14	18	300.08	0.08	0.0004	0.189	6.678768025
    18	14	19	98.78	0.1	0.0004	0.278	10.8220877
    19	19	20	300.08	0.08	0.0004	0.189	6.686414994
    20	19	21	295.02	0.08	0.0004	0.189	6.685858112
    21	22	19	91.74	0.08	0.0004	0.236	2.550185403
    22	22	23	354.42	0.08	0.0004	0.189	8.826465819
    23	22	24	295.24	0.08	0.0004	0.189	8.819974643
    24	25	22	114.62	0.13	0.0004	0.236	20.19662586
    25	25	26	299.2	0.1	0.0004	0.189	8.813556575
    26	25	27	271.26	0.1	0.0004	0.189	8.810492745
    27	28	25	135.96	0.2	0.0004	0.21	37.82067518
    28	28	29	209.44	0.12	0.0004	0.189	8.799855999
    29	28	30	231.22	0.12	0.0004	0.189	8.802244518
    30	35	28	155.32	0.25	0.0004	0.321	55.4227757
    31	35	7	575.96	0.2	0.0004	0.321	26.37604637
    32	34	11	442.86	0.2	0.0004	0.321	30.9589869
    33	7	1	150	0.1	0.0004	0.321	8.849479077
    34	11	31	150	0.12	0.0004	0.321	11.95830395
    35	14	32	150	0.08	0.0004	0.321	6.673783076
    ];

% node																																																																																																																																																																																																														
model.node = [																																																																																																																																																																																																														
	% node load(MW) Ts Tr type																																																																																																																																																																																																						
    1	1.283999793	78.75	40	1	8.849479077
    2	0	78.75	40	1	0
    3	1.284000002	78.75	40	1	8.812506012
    4	1.740000003	78.75	40	1	11.92041343
    5	0	78.75	40	1	0
    6	1.283999971	78.75	40	1	8.856657287
    7	0	78.75	40	1	0
    8	1.28399972	78.75	40	1	8.86911535
    9	1.283999794	78.75	40	1	8.849097137
    10	1.283999606	78.75	40	1	8.888014674
    11	0	78.75	40	1	0
    12	1.283999997	78.75	40	1	8.838417868
    13	0	78.75	40	1	0
    14	0	78.75	40	1	0
    15	0	78.75	40	1	0
    16	0.965999975	78.75	40	1	6.679883146
    17	0.965999969	78.75	40	1	6.68462873
    18	0.965999977	78.75	40	1	6.678768025
    19	0	78.75	40	1	0
    20	0.965999761	78.75	40	1	6.686414994
    21	0.965999763	78.75	40	1	6.685858112
    22	0	78.75	40	1	0
    23	1.283999948	78.75	40	1	8.826465819
    24	1.283999955	78.75	40	1	8.819974643
    25	0	78.75	40	1	0
    26	1.283999981	78.75	40	1	8.813556575
    27	1.283999983	78.75	40	1	8.810492745
    28	0	78.75	40	1	0
    29	1.283999993	78.75	40	1	8.799855999
    30	1.283999992	78.75	40	1	8.802244518
    31	1.74	78.75	40	1	11.95830395
    32	0.965999982	78.75	40	1	6.673783076
    33	-9.45	78.75	40	2	-66.04612219
    34	-5.25	78.75	40	2	-30.9589869
    35	-11.267	78.75	40	0	-81.79882207
    ];


%%-----  System Data  -----%%																																																																																																																																																																																																														
																																																																																																																																																																																																														
model.water_c = 4182; % (J*kg^(-1)*K^(-1))																																																																																																																																																																																																														
model.water_dens = 958.4; % Density of water (kg/m^3) at 100 degrees C																																																																																																																																																																																																														
model.Pbase	=	1e6;	%	base	
model.g = 9.81; % Gravitational acceleration
model.viscosity=.294e-6; % temp is 100deg,unit:m2/s kinematic viscosity
