%%-----  Topo Data  -----%%																		
%% system MVA base																		
mpc.baseMVA = 100;																		
																		
%% bus data																		
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin	Pmax	Pmin			
mpc.bus = [  %% (Pd and Qd are for reference, and load power is shown below)																		
	1	3	0	0	0	0	1	1.05000000000000	0.999999999999999	230	1	1.1	0.9	2000	0			
	2	1	0	0	0	0	1	1.03281032665893	-0.583647141705036	230	1	1.1	0.9	500	0			
	3	1	0	0	0	0	1	1.02234275929398	-1.41712446831302	230	1	1.1	0.9	500	0			
	4	1	0	0	0	0	1	1.01587219394564	-2.29546702153801	230	1	1.1	0.9	500	0			
	5	1	0	0	0	0	1	1.00322687496341	-3.12715017662369	230	1	1.1	0.9	500	0			
	6	1	0	0	0	0	1	0.987963952891348	-4.21990737063525	230	1	1.1	0.9	500	0			
	7	1	0	0	0	0	1	1.00553061304068	-3.23704016691798	230	1	1.1	0.9	500	0			
	8	1	0	0	0	0	1	0.987504274773903	-4.25514836687740	230	1	1.1	0.9	500	0			
	9	2	0	0	0	0	1	1.05000000000000	-1.89699310043575	230	1	1.1	0.9	1500	0			
	10	2	0	0	0	0	1	1.05000000000000	-1.22569243814780	230	1	1.1	0.9	500	0			
	11	1	0	0	0	0	1	1.00412024382409	-2.58352109340641	230	1	1.1	0.9	500	0			
	12	1	0	0	0	0	1	0.994671498117181	-3.32777348982754	230	1	1.1	0.9	650	0			
	13	1	0	0	0	0	1	0.999818815711776	-3.02098781316980	230	1	1.1	0.9	500	0			
	14	1	0	0	0	0	1	1.03734880781526	-1.23098513900776	230	1	1.1	0.9	500	0			
	15	1	0	0	0	0	1	1.03126086915643	-0.977271050500901	230	1	1.1	0.9	500	0			
	16	1	0	0	0	0	1	1.03194474240816	-0.671032540253646	230	1	1.1	0.9	500	0			
	17	1	0	0	0	0	1	1.03165313162612	-0.690331228171748	230	1	1.1	0.9	500	0			
	18	1	0	0	0	0	1	1.03000759336129	-0.806211607379310	230	1	1.1	0.9	700	0			
	19	1	0	0	0	0	1	1.03155088733496	-0.710473706241255	230	1	1.1	0.9	500	0			
];																		
																		
%% generator data(here, status=1 means this gen is a CHP unit, status=2 means this gen is a thermal unit)																		
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10
mpc.gen = [																		
%	1	7.50314624400000	0	40000	-40000	1.05000000000000	100	1	40000	0	0	0	0	0	0	0	0	0
	9	50	0	40000	-40000	1.05000000000000	100	1	40000	0	0	0	0	0	0	0	0	0
	10	50	0	40000	-40000	1.05000000000000	100	1	40000	0	0	0	0	0	0	0	0	0
];																		
																		
%% branch data																		
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax					
mpc.branch = [  %% (r and x specified in ohms here)																		
	1	2	0.00370300000000000	0.0181010000000000	0	1000	1500	2000	0	0	1	-360	360					
	2	3	0.00688600000000000	0.0338820000000000	-0.0387930000000000	1000	1500	2000	0	0	1	-360	360					
	3	4	0.0179000000000000	0.0942000000000000	-0.0921500000000000	1000	1500	2000	0	0	1	-360	360					
	4	5	0.00730000000000000	0.0402570000000000	-0.0410000000000000	1000	1500	2000	0	0	1	-360	360					
	5	6	0.0143180000000000	0.0720270000000000	-0.108237000000000	1000	1500	2000	0	0	1	-360	360					
	5	7	0.0150290000000000	0.0601140000000000	-0.0365970000000000	1000	1500	2000	0	0	1	-360	360					
	7	8	0.0191000000000000	0.0859000000000000	-0.150900000000000	1000	1500	2000	0	0	1	-360	360					
	7	9	0.0186000000000000	0.0861000000000000	-0.164100000000000	1000	1500	2000	0	0	1	-360	360					
	4	10	0.0201960000000000	0.0807850000000000	-0.0612780000000000	1000	1500	2000	0	0	1	-360	360					
	1	11	0.0121000000000000	0.0918000000000000	-0.131500000000000	1000	1500	2000	0	0	1	-360	360					
	11	12	0.0146460000000000	0.0598950000000000	-0.0359570000000000	1000	1500	2000	0	0	1	-360	360					
	11	13	0.00415300000000000	0.0244830000000000	-0.0169640000000000	1000	1500	2000	0	0	1	-360	360					
	10	14	0.00898400000000000	0.0437190000000000	-0.0343100000000000	1000	1500	2000	0	0	1	-360	360					
	14	15	0.00600000000000000	0.0350000000000000	-0.0325000000000000	1000	1500	2000	0	0	1	-360	360					
	15	16	0.00450900000000000	0.0180340000000000	-0.00915000000000000	1000	1500	2000	0	0	1	-360	360					
	16	17	0.000250000000000000	0.00100000000000000	-0.0120000000000000	1000	1500	2000	0	0	1	-360	360					
	17	18	0.00250000000000000	0.0100000000000000	-0.0550000000000000	1000	1500	2000	0	0	1	-360	360					
	17	19	0.00384300000000000	0.0153730000000000	-0.0120420000000000	1000	1500	2000	0	0	1	-360	360					
	19	2	0.00275100000000000	0.0110050000000000	-0.00862000000000000	1000	1500	2000	0	0	1	-360	360					
	15	19	0.00860100000000000	0.0437190000000000	-0.0331840000000000	1000	1500	2000	0	0	1	-360	360					
	16	2	0.000475000000000000	0.00190000000000000	-0.00165000000000000	1000	1500	2000	0	0	1	-360	360					
];																		
																		
%%-----  OPF Data  -----%%																		
%% generator cost data (only thermal unit)																		
%	1	startup	shutdown	n	x1	y1	...	xn	yn									
%	2	startup	shutdown	n	c(n-1)	...	c0											
mpc.gencost = [																		
	9	0	0	0	0.102	30	0			%grid								
	10	0	0	0	0.115	26	0			%thermal gen								
];																		
																		
%% load active power (MW)																		
load = [																		
%bus(node)1	bus2	bus3	bus4	bus5	bus6	bus7	bus8	bus9	bus10	bus11	bus12	bus13	bus14	bus15	bus16	bus17	bus18	bus19
8.623923647	4.447394943	10.67010159	3.195690386	3.192718427	9.733798257	5.36294271	7.974634837	3.779672357	3.365291837	6.427920577	8.008691513	11.10311075	5.42592908	10.48401368	6.703578203	4.71600785	8.069986137	4.783296967
8.997052777	4.426408317	10.52321047	3.71370341	3.167834957	9.63623932	5.548969203	7.98115672	3.423540373	3.49872573	6.297139777	7.86222894	11.07774027	5.181514613	10.53012018	6.638866343	5.01385593	8.30090223	5.206715027
8.57740007	4.5612167	9.82378297	3.39667658	3.317028401	9.41104592	5.711906167	8.081447453	3.700230923	3.483253407	6.70720037	7.312133487	10.56469366	5.498683697	9.752618747	6.175174477	5.205719443	8.104611637	4.75473002
8.328824137	4.204183057	9.272914263	3.42828415	3.450989243	9.172443033	5.55507364	7.672886877	3.679809507	3.210122441	6.643644873	7.239031187	9.589737463	5.660657433	9.501160947	6.415806363	4.659459463	8.024806337	5.267837697
7.92125076	4.287103453	8.77408659	3.39552144	3.4632246	9.196787183	5.550751927	7.559197023	3.325404828	3.166697424	6.605330497	6.820674917	9.140761937	5.32548577	8.6864088	5.670709113	4.841666417	7.160300817	5.103791487
7.685470207	4.63283243	9.19023374	3.435834957	3.253086636	8.675471413	5.564703683	7.68711862	3.51281442	3.35775041	6.76856489	6.91727051	9.532641843	5.620799827	8.8451125	5.586533097	5.098546807	7.22518631	4.71514597
8.352977333	4.533982853	8.950294563	3.27679206	3.703942577	9.25856228	5.641050333	7.7494383	3.671828633	3.649717803	6.710225287	7.303605787	9.8364138	5.58005553	9.15491195	5.89002207	5.202730833	7.532470507	5.197918577
8.44463516	4.284170313	9.422006187	3.25326137	3.649364323	9.48981455	5.555188377	8.024939797	3.112869301	3.222208948	6.532408523	7.214976063	10.04380912	5.508746507	9.37172532	6.468821827	4.966665867	7.591278003	5.25386195
8.270651197	4.713133073	9.438993167	3.70197654	3.398651077	9.218094377	5.468364887	8.014586507	3.153500659	3.436144327	6.18062724	7.030885077	10.41746713	5.37306962	9.85839414	6.409078387	5.299449723	7.553516587	5.012399777
8.49654169	4.433442033	9.523704753	3.70280973	3.751719887	9.431747183	5.33248933	8.08343778	3.348154463	3.181221164	6.228399817	7.4448071	10.32767488	5.37388174	9.898458763	6.182689207	4.916903637	7.999928737	4.897025763
7.973518867	4.514249997	8.985875493	3.487803723	3.470886863	8.977055217	5.710968143	7.78097289	3.184837598	3.661794123	6.49706606	7.175159347	9.508335527	5.373725563	9.53923036	6.32524295	4.79786468	7.554577427	4.91403618
7.83430569	4.503188257	8.82938725	3.317229778	3.23002112	9.081319993	5.229358877	7.301377213	3.096470053	3.584835913	6.777500547	6.584305743	9.375758313	5.347057983	9.076284317	6.04779396	4.92342671	7.569173567	5.075912617
7.956511957	4.274742523	8.76475375	3.481720333	3.623053627	8.5453788	5.423814893	7.4402825	3.604159807	3.45889569	6.713728753	6.46550936	9.22671013	5.228960167	8.36901797	5.83166628	4.7244042	7.338677977	4.79053597
7.403664043	4.205924443	8.513680967	3.570984763	3.42617611	9.008100117	5.662016557	7.291296643	3.34349005	3.38503958	6.817950763	6.509327127	8.828608477	5.499647733	7.996023357	5.587264157	4.76589984	6.831036143	4.6789876
8.34318685	4.465292727	9.059255603	3.860352327	3.43062238	9.00723316	5.600533067	7.58661658	3.808960487	3.59073471	6.971950147	7.099133393	9.505918163	5.457667187	8.582464377	6.182411503	4.972386823	7.505304933	5.425469
8.217785817	4.180503193	9.44310259	3.528968637	4.10770331	9.16200171	6.391344893	7.881203733	3.58649099	4.062014353	6.765361037	7.87969109	10.28875676	6.204166053	9.467883667	6.408530767	5.498624243	7.59445724	5.657990173
8.9902426	4.45435914	10.44648236	4.00909319	3.744853187	9.328308627	6.527404293	8.003435323	4.10049811	4.253118477	6.897185217	8.462093537	11.39070066	6.20139665	10.32008838	6.281429933	5.57261821	7.96846877	5.686157703
9.3553866	4.54051604	11.15634415	4.085437847	4.14151348	10.26359547	7.058492267	8.521806483	4.38353821	4.13206386	7.342701623	9.506252777	11.86789886	7.033255673	10.92724956	6.765564433	5.954792477	8.488484647	5.969771393
9.888160257	4.55470142	11.04122923	4.27458233	4.51109676	10.26106212	6.882004273	8.848700983	4.36791698	4.215964633	7.055381143	9.855823157	12.60062435	7.22154984	11.43679534	7.12981989	5.90407051	8.74876396	5.713558637
10.26099943	4.08176383	11.30404659	4.19969725	4.25932708	10.53150101	7.330055477	8.79253554	3.980197143	4.43358048	7.57777925	9.685094963	12.68380863	6.93201085	11.80094491	7.242121923	5.764638607	8.945227477	5.468050663
10.4194198	4.462324743	11.47713362	3.915932147	4.17514569	10.37050964	7.331860793	9.08414459	4.00466077	4.15109732	7.557229157	9.997176233	12.90866559	7.362734797	11.49058865	7.057715897	5.65136847	8.82235378	5.66079765
10.22982028	4.197599743	11.64951443	3.94385743	4.187277993	10.09013064	7.122071507	8.665921837	4.11628714	4.391367813	7.18463029	9.731478493	13.17558286	7.36678794	11.67595327	7.310801833	5.598865573	8.56807636	5.902554517
10.55593707	4.74452585	11.98527283	4.15761485	4.53615021	10.83057602	6.83085178	9.417673097	4.329258487	4.59189283	7.190926403	10.2685538	13.57457786	7.42209102	11.95536378	7.673348957	6.057565007	9.341876657	6.004901693
10.81301963	5.136989203	12.98136162	4.590427137	5.15853375	10.89233691	7.06881975	9.893080737	4.775860107	4.864096427	7.46664453	10.48888947	14.14583903	6.929380907	13.17703958	8.25620552	5.84082963	9.573660337	6.360577193
11.36559396	5.66711869	13.99296724	5.206716123	5.242428067	11.51385124	7.074697813	10.1864989	5.04958831	5.13812629	7.529134453	11.1507423	15.01703299	6.892746357	14.3565041	9.600997883	6.092697043	10.54432297	6.52700681
11.21101464	5.57781392	15.23984795	5.918276927	5.55113577	11.44727017	7.132568623	10.91569857	5.69966569	5.800352413	7.630404877	11.26005697	15.39886117	7.120623597	14.80529384	9.888346477	6.73766332	10.95646153	6.659184393
11.87430456	6.51962943	15.63267453	6.01959659	6.078294113	12.25380531	7.172607967	11.32911517	6.099053553	5.91878353	7.571230983	11.85619817	16.08066908	7.25244228	15.92667982	10.88879493	6.722092117	11.42316486	6.668908917
12.0787534	7.65501754	16.47599089	7.455426663	7.224778403	12.8747158	7.260163687	12.06793229	7.06173793	7.174054067	7.680513573	11.81793536	16.40322362	7.41484231	16.56911798	11.72252195	7.716144643	12.23960618	7.553307807
12.39459408	9.353046687	17.38989525	9.188116123	8.82256954	13.26838529	8.089490247	13.35181271	8.950926893	9.100853743	8.427599843	11.89008573	16.62401398	7.998079457	17.48730426	13.50390928	8.489629977	13.02924649	8.84754205
13.01946098	10.86450751	17.90461061	10.32130547	10.30764504	13.63783785	8.043703617	14.47426308	9.92894087	9.973067847	8.52835943	11.75454293	16.41021178	7.8919904	18.36052864	14.57318329	9.670254527	14.49365246	9.451724817
12.70560683	12.68919955	18.8460585	11.72877277	11.72534972	13.78655286	8.38633126	15.02027536	11.48171963	11.79930906	9.455822327	12.0994002	16.35758232	8.29967942	18.80103843	16.14381881	10.4001525	15.16621097	10.21339274
13.16058535	13.07397792	18.86094048	11.89212278	11.99424772	14.42128014	8.413644903	15.54447212	12.19373944	12.02883175	9.78275114	12.22005984	16.87860461	8.376204257	19.26786907	16.97830135	11.01910922	15.3032777	11.11719245
13.21837096	13.3317791	18.98281919	12.31463101	12.40370959	14.71127041	9.005358527	15.51449629	12.19112106	12.2298129	10.00860974	12.41936506	16.91780345	8.573210287	18.96353004	17.11640476	11.37873901	15.50201917	11.00543923
13.68515308	13.41230339	19.50197587	12.37808656	11.98565749	15.0910249	9.13535809	15.79034353	11.78944889	12.27707994	10.02493883	12.25558833	17.36621641	9.081399467	19.68924899	16.72870408	10.88531707	15.81776111	10.94423758
14.20185831	13.71128133	19.46001877	12.32264558	12.01612698	15.30805048	9.27883235	15.97272446	11.94149727	11.9516268	10.65330145	12.12348669	17.30739328	8.644267667	19.26196786	17.00731176	11.27560345	16.46579483	10.96633046
13.84261489	13.83651336	19.58780071	12.1126532	12.05976169	15.90313279	9.23177956	16.08226424	11.96568705	11.74325782	10.94954364	12.47769931	17.21001333	8.87575127	19.42805402	16.87830742	11.15311273	16.20257795	11.52557524
14.54122963	14.15746045	19.45816786	12.14865312	11.75336817	15.79593223	9.39155105	16.38090366	11.80981371	12.26915671	10.92083545	12.2529875	17.18655866	9.166154037	19.83299283	17.31139944	11.68520122	16.51876906	11.448214
14.14843039	14.23914084	20.05728398	12.31970001	11.86393549	15.96673066	9.260976393	16.50491043	11.75937975	11.95589613	10.95999123	12.43956444	17.55777601	9.598695363	19.7259123	16.79151297	11.93195431	16.95764021	11.57465681
14.25788503	13.93090583	19.72313484	11.9925485	11.78275044	16.50197427	9.484962657	16.65320257	12.0109217	12.14465265	11.69243885	12.70886195	17.43252123	9.278625013	19.63091889	17.22602179	11.47450366	16.54336107	11.71112234
14.61072671	13.91303802	19.36714795	11.88877384	11.76145336	16.63282037	9.413875633	16.54565337	11.94252959	12.29274661	11.83394552	12.3770067	17.50004178	9.855128987	19.84586724	16.75995784	11.85463998	16.75474997	11.85208862
14.21590097	14.49696685	19.65492319	12.13040076	12.30266286	16.67244363	9.79655622	17.07932426	12.2588509	11.8918786	11.83217149	12.57279044	17.27479575	9.964170453	19.25277893	17.00195591	11.96919715	16.82546126	11.80085289
14.05191064	14.16432233	18.8318987	11.91963967	11.74495116	16.70989066	9.910714227	16.62995676	12.08104888	12.38795743	12.1125139	11.99724626	16.62046968	9.642208227	19.0490759	16.483024	12.13871444	16.68375154	12.00082962
14.28592106	14.42704549	18.60601518	12.1523195	12.29699031	16.78693797	9.52326957	16.40029011	12.28264437	12.29386923	12.31884825	11.98976171	16.27974933	10.12580198	18.50474478	16.2529349	11.96428617	16.58564406	11.7686809
13.72428232	14.76063351	18.20117742	11.7898488	12.09633715	16.72508137	9.705305437	16.16766763	12.06507988	12.28572179	12.16762277	11.70884196	15.60191204	9.638975807	17.95500058	16.05067848	12.25391986	16.17955258	11.96695864
14.39852083	14.28951125	17.84302492	12.45049902	11.86323672	16.785963	9.766959473	16.2841488	11.79678163	12.22516756	12.68215638	11.36692525	16.03944805	10.13829146	18.26462116	16.54928727	12.06256863	16.32355863	12.13260431
14.25494913	14.46666466	18.58165119	11.96126448	12.15038548	16.7695781	10.09947828	16.46716105	11.82674385	12.18587039	12.53932746	11.69887976	15.97780394	10.0696532	18.322011	16.2899708	12.42697788	16.73400522	12.48544931
14.31291654	14.64495398	18.62979696	12.27259049	12.267309	17.17374059	10.2358348	16.98517242	11.97849904	12.12546874	12.43687839	12.18759126	16.06661206	10.16646908	18.23339838	16.60253566	12.52188141	16.98918046	12.44357226
14.691561	14.81079246	18.7688182	12.14369424	11.97208478	16.81207149	10.37756017	16.85288281	12.15589991	12.46009541	12.5129302	11.74373515	16.25854823	10.1479202	18.58549277	16.34207041	12.42993408	16.83571359	12.59095194
14.96397624	15.09848767	19.12785283	12.44432097	12.61797894	17.09254682	9.901012443	17.04368894	12.44463657	12.36271926	13.0148448	12.0554232	16.87863024	10.03538242	18.6624568	16.56449498	12.29871725	16.84176652	12.73005408
14.59772841	14.68805168	18.97728201	12.06273924	12.47942289	16.8658547	10.03649824	16.66648988	12.15035947	12.56517126	12.6443115	12.0048127	16.74512353	9.997842653	18.97945531	17.06842765	12.50335109	16.95619879	12.7901507
14.75591823	15.22993975	18.5609022	12.1668857	12.42974594	17.11581972	10.22674211	17.02141754	12.30817217	12.78077423	12.47695103	11.98292429	16.68293278	10.33746394	18.68136379	16.80905194	12.64700222	16.83417173	12.32811328
14.61395541	14.92610576	18.72954913	12.81705585	12.77485183	17.14183216	10.18426854	17.08327269	12.62897205	12.72451562	13.00066189	12.01288304	16.49332016	10.46694888	18.73600797	16.48162492	12.49928263	16.74830639	12.72257074
14.61941401	15.21153283	18.75512053	12.33183377	12.44853865	16.63838915	10.1721611	16.91224062	12.4647887	12.84193319	12.79412077	12.3036453	16.36423214	10.39166078	18.70270842	16.74062371	12.87059412	16.56661676	12.81874256
14.45844767	15.1147485	18.51478456	12.53154455	12.85925414	16.80270306	10.06344788	16.78795742	12.52456611	12.6336874	13.13033722	11.94424334	16.03923526	10.33528239	18.6607759	16.90233788	12.47452652	16.54677927	12.86082237
14.29937682	15.06146235	18.23662674	12.67061404	12.53914365	16.63481992	10.60196482	16.97220751	12.42343797	12.90233083	12.74441433	12.0632977	15.81747544	10.72405961	18.36333604	16.85526657	12.63496312	16.75099033	12.87092162
14.58870202	15.20435873	17.68028922	12.48833469	12.9514403	16.81647329	10.42477311	17.04337402	12.5823487	12.54647133	13.30154981	11.66239422	15.87997124	10.5506523	17.72767808	16.3898938	13.16155236	16.45947385	12.56199176
14.13626103	15.47124921	17.95200338	12.69125466	12.4920604	17.1691791	10.67045293	16.58078938	12.92647711	12.46619954	13.2134543	11.62823279	15.34925507	10.37209163	18.09864722	16.25150294	13.04359277	16.98952435	12.73688736
14.44897607	15.62859332	17.61691145	12.86214783	12.97964271	16.95186016	10.60119977	16.84104522	12.93341385	12.67114065	13.30388436	11.87551899	15.08706609	10.68971179	17.96694295	16.57980932	13.12774121	16.40863934	12.8462551
14.5723032	15.60360828	17.8467681	12.67203002	12.66465841	17.10657801	10.18589707	16.52716507	12.58651732	12.34491441	12.86766719	11.3102236	15.16544267	10.36535419	17.81176387	16.37206708	13.12149664	16.444034	13.00947989
14.09631822	15.08692202	18.01597984	12.93602963	12.69111625	17.01712841	10.28492069	16.37730488	12.76503484	12.86644595	13.13861713	11.27326087	15.65842433	10.57920243	18.15217092	16.48164216	12.57080517	16.84585736	12.68345664
14.42311828	15.05776193	18.22783745	12.41393334	12.97833611	17.20692453	10.66674038	16.48965225	12.48400149	12.49203031	13.1074348	11.75829053	15.68244015	10.68853934	17.91177847	16.36592884	12.77693807	16.87681539	12.55173852
14.48585591	15.57298319	17.92964362	12.60454302	12.84068852	16.73387048	10.62987739	16.74613703	12.96002436	12.77836544	12.99160949	12.01324736	15.99468324	10.18871619	18.51101666	16.46156552	12.97678047	16.76600137	13.14684734
14.44460948	15.42451876	18.61420322	12.50240834	12.52474771	16.95881545	10.29385922	16.84540684	12.73200122	12.84478932	12.75517245	11.5670711	15.87313734	10.66308782	18.09456879	17.04629457	12.60818091	16.70200122	12.65332421
14.52329521	15.03038772	18.45279437	12.70021306	12.34006126	17.17072601	10.48336782	16.66910402	12.43944481	12.75604209	13.11708597	11.92352368	16.03060113	10.53584572	18.39674738	17.08284037	12.95926652	16.66258796	13.1930639
14.82768042	15.34742943	18.08602788	12.30923005	12.31429601	17.41416494	10.46436721	16.76387359	12.35801692	12.71875297	12.85084688	11.56626785	15.59551777	10.15285113	18.13359996	16.65055456	12.71470417	16.90425468	12.95706254
14.83420263	15.61467088	18.16375412	12.31080511	12.81031009	17.54407249	10.66080933	16.95037434	12.64558052	12.42657526	13.41954459	11.91136901	16.15372874	10.27101697	18.20600096	16.91120047	13.12271689	17.17991946	12.72805441
14.71512851	15.23951222	18.16461663	12.46033134	12.51979185	16.98103208	10.72598274	17.18511449	12.42796526	12.56334532	13.22124901	11.90630514	15.94139003	10.71062195	18.18989725	17.02011084	13.09271896	16.97212517	12.8802437
14.4991477	15.12298496	18.58898083	12.26718485	12.45770795	17.0610876	10.05712455	17.15761312	12.24904267	12.82587868	13.19415121	11.54515875	16.09851887	10.32664359	18.60793617	17.0334008	12.76958386	16.76985896	12.94275628
14.26657124	15.01120556	18.71033698	12.27148508	12.50687302	17.21963356	10.49904791	17.19028201	12.26375487	12.31885767	13.03557534	12.07473946	16.38648276	10.49538347	18.78808813	17.05161061	12.6442789	17.25114495	12.55435243
14.54518862	14.82094042	18.37878825	12.26721951	12.61454296	17.33347529	10.14601913	16.92104711	12.05279051	12.15918623	12.85387869	12.14684217	16.09902473	10.29058581	18.56548242	16.86081172	12.87007809	17.05416318	12.30547189
14.24853636	14.09073649	18.25994916	11.23752228	11.48162829	16.84787064	9.753264927	16.591659	11.66493703	11.39266255	12.88578497	11.5370399	15.91156652	9.6358841	17.63679088	15.5739118	11.73589918	16.28038946	12.09689331
13.63052678	13.22606698	16.56780529	10.71065878	10.64075741	16.48843468	9.62105504	15.73109417	10.75900202	10.75504309	12.18309125	11.03186219	14.94028282	9.987261007	16.43185369	14.31651754	11.41710583	15.42737346	11.29741454
13.4756166	11.23594976	15.01937088	8.926073407	8.56171487	15.5605283	9.76289281	14.39287097	9.036755963	8.719453443	12.16750085	10.36284019	14.13413929	9.749583047	14.8462941	12.38446788	10.49988646	14.12224671	10.382822
12.79883401	9.71399976	13.39400185	6.906438717	7.014385663	15.33960512	9.53935858	12.76892301	6.96710054	7.192447947	12.21175988	10.02397601	13.54552263	9.71996233	13.24123855	10.04806923	9.632368987	12.88102815	9.8690332
12.82689021	8.568070097	12.89589315	5.626857773	5.68198437	15.26958375	9.40651056	12.52218398	6.18268941	5.70809972	11.79322423	10.35195305	13.64897219	9.402652957	12.87831749	9.159754207	8.57034253	12.34994549	8.729796193
12.89453704	7.424870217	13.70524989	5.25247959	5.144626997	15.07612289	9.542350633	12.18444462	5.38669737	5.03967921	11.77914042	10.72411026	14.19247263	9.314091053	14.00243503	9.435501107	8.221247237	12.05469734	8.756953623
13.22179082	7.32030085	14.44673067	4.849668617	4.793607773	15.95866853	9.182394927	12.59380257	4.55056056	4.746903587	11.20766479	11.04716916	15.41369576	8.868664107	14.76591359	9.12685986	8.2505044	12.72627958	8.088061223
13.97090404	6.818815153	15.79297696	4.338734133	4.383831007	16.61906074	8.759935907	12.99515568	4.432970433	4.35753687	11.33728746	12.02852411	17.3218492	8.835865193	16.02209404	9.66648336	7.685458407	12.752657	7.54680108
14.59353519	6.355332293	17.60281454	4.562078817	4.126111757	16.57901322	8.68831887	13.73175794	4.33815538	4.58850359	11.24338388	12.42058873	18.70624828	8.601496833	16.99334852	10.24994603	7.387494387	13.44206053	7.514785763
14.68916652	6.236355297	18.10268517	4.195296947	4.21915744	16.51670063	8.345908153	13.31672237	4.215178417	4.421259597	10.20402776	13.26782847	19.26118056	8.574208907	17.56823527	10.57841805	6.974887823	13.28398131	7.04960098
14.57061451	5.28349682	17.48182725	3.786688123	4.30448948	15.58019171	7.944601837	13.13321523	4.069170633	3.792598327	9.580356447	12.61878605	18.93234773	8.06482437	17.57520031	10.19322288	6.64352452	12.87247319	7.195874093
13.61808497	5.256131213	17.64559446	4.39713897	3.936530697	15.04664574	7.874739333	12.62752938	4.33825844	4.19328604	8.51413599	12.67462208	18.54933062	7.618558977	17.84653657	10.34351167	6.482615313	12.23180293	6.442037513
13.21695655	4.714584517	17.40126949	4.12507041	4.42061046	13.96606002	7.496965367	12.29359466	3.919658327	3.852318013	7.95666246	12.80997835	18.79757737	7.40136457	17.80954214	10.00074727	6.392814323	12.47480658	6.178114503
13.30583988	4.43372117	17.74266007	3.890427317	4.348882847	14.18505861	7.125626313	11.94077289	4.50162945	4.135011763	7.881862	12.74224847	18.51992891	6.998472807	17.92583804	9.983042747	6.15559297	12.08096316	5.81964033
13.56575893	4.72496012	18.18465006	3.92447124	4.211215527	13.65965331	6.892536567	12.16187048	4.409465033	4.20134584	7.849613467	13.03274607	19.23489482	7.406482503	18.06680684	10.64815259	5.651138217	12.25019168	6.149531757
13.5467285	4.870197433	18.64935366	4.52207084	4.50309744	13.91793952	6.99627316	12.35066443	4.083188087	4.41668026	7.418690873	13.23960165	19.81649469	7.36338148	18.67032618	11.01767026	5.56594343	12.57421936	5.59137036
14.14312177	4.674845547	19.40134534	4.43845037	4.237487547	14.69838803	7.2502868	12.57991373	4.27172481	4.27758522	7.17781892	13.78940838	20.57885298	7.11127781	19.1974574	11.06127199	5.958571623	12.73733286	6.059316783
14.38274245	4.86385866	19.84823671	4.27365871	4.00864811	14.83372541	6.982598823	12.80509152	4.066025623	4.033657963	7.760006623	13.56804109	20.99618274	6.716608767	19.77311167	11.34860845	5.83147446	12.83563485	5.794786997
14.06828907	4.417471517	19.22021697	3.870970367	4.42484957	14.35210724	7.21598748	12.56745804	3.912486557	4.456134867	7.454445543	13.55179547	20.37177262	7.314579053	19.21595859	10.89802393	5.82698432	12.98976581	6.110268833
13.25463831	4.67190659	17.95464062	4.049821457	4.40947258	14.00393283	6.890513007	12.21745051	4.182490237	4.3463401	7.72361238	13.27011966	19.24968532	7.148968027	17.94783024	10.75902638	5.684710477	12.17608731	5.827729473
12.82044929	4.351581703	17.01132594	4.12578951	4.025473487	13.37595558	7.102439097	12.05110494	4.36679971	4.19171569	7.313384093	12.41213548	18.32126447	6.7250025	17.19576056	9.9951519	5.671335447	11.88845691	6.03622456
12.24829909	4.408638767	16.23802066	3.84533199	3.997309973	12.95306341	7.21195466	10.94182193	4.288604913	4.377422857	7.518903973	12.23216974	17.2463278	7.286132663	16.0214324	9.74947259	6.02481193	11.36967576	5.653951597
12.11645088	4.224769223	16.18573015	4.048520473	4.227870647	12.48986294	6.53111777	11.2016959	4.29962125	3.913247983	6.947559763	11.79865354	16.81497218	6.600993533	16.04786446	9.51692549	5.60211916	11.25906922	5.56866302
11.93036894	4.684037027	16.19324771	3.728268853	3.83937398	12.75356482	6.381077337	11.39378776	3.733383637	3.736172793	7.107758223	11.09011088	17.13982631	6.60031957	16.25743409	9.805815047	5.637696757	10.97204536	5.537221483
11.77740142	4.85503917	16.28539502	3.63833029	3.372603723	12.50710983	5.77330123	11.42924451	3.597827653	3.750572303	6.80149483	10.81347273	16.66106112	6.005609927	16.43786338	9.364747547	5.08264611	11.04130168	5.566822507
12.16783946	4.72359	16.08266912	3.62774683	3.761196547	12.5657935	5.57842367	10.99314634	3.342792887	3.300447927	6.69598502	10.84159623	17.23329179	5.88160489	16.51202217	9.69778907	4.913828523	11.39187062	4.945871857
];																		
