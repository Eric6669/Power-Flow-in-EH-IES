% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-06 09:15:22
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-06 16:58:39
%  */
model.pipe = [																																																																																																																																																																																																														
	69	2	5328	1.2	0.0004	2.08776
	2	3	5650	0.2	0.0004	2.08917
	2	4	7664	1.1	0.0004	2.08788
	4	5	1461	0.2	0.0004	2.08884
	4	6	1099	0.3	0.0004	2.0884
	4	7	910	1.1	0.0004	2.08797
	7	8	705	1.1	0.0004	2.08799
	8	9	973	1.1	0.0004	2.08801
	9	10	618	0.5	0.0004	2.08839
	10	11	1739	0.5	0.0004	2.08982
	9	12	1793	0.8	0.0004	2.08804
	12	13	405	0.15	0.0004	2.08817
	12	14	312	1.1	0.0004	2.08809
	14	15	350	1.1	0.0004	2.08808
	15	16	922	1.1	0.0004	2.08812
	16	17	737.5	0.8	0.0004	2.08817
	17	18	423.26	1	0.0004	2.08821
	18	19	315	1	0.0004	2.08824
	19	20	600	0.8	0.0004	2.08828
	20	21	231	0.8	0.0004	2.08832
	21	22	169	0.8	0.0004	2.08835
	22	23	64	0.8	0.0004	2.08836
	23	24	100	0.8	0.0004	2.08839
	24	25	125	0.8	0.0004	2.08841
	25	26	166	0.8	0.0004	2.08852
	26	27	532	0.2	0.0004	2.08895
	26	28	967	0.8	0.0004	2.08814
	55	49	354	0.8	0.0004	2.08703
	53	54	134.5	0.8	0.0004	2.10099
	51	53	281	1.2	0.0004	2.10097
	51	52	491	0.4	0.0004	2.10112
	50	51	296.8	1.2	0.0004	2.10098
	49	50	520	1.2	0.0004	2.101
	45	49	334	1.2	0.0004	2.10102
	44	45	109	1.2	0.0004	2.10103
	45	47	145	0.1	0.0004	2.10116
	43	44	373	1	0.0004	2.10105
	45	46	225.14	0.2	0.0004	2.10123
	45	48	813.2	0.1	0.0004	2.10251
	41	43	46.2	1	0.0004	2.10106
	41	42	204	0.2	0.0004	2.10115
	38	41	228.45	1	0.0004	2.10107
	38	39	115	0.2	0.0004	2.10113
	38	40	92	0.2	0.0004	2.1011
	36	38	388	1	0.0004	2.10109
	68	36	175	1	0.0004	2.10112
	67	66	281.2	0.6	0.0004	2.1012
	66	65	314.33	0.6	0.0004	2.10121
	65	64	72	0.6	0.0004	2.10121
	64	63	575	0.6	0.0004	2.10123
	63	62	244	0.6	0.0004	2.10125
	62	60	94	0.6	0.0004	2.10126
	60	61	250	0.2	0.0004	2.10128
	60	59	76	0.6	0.0004	2.10127
	37	62	909.38	0.2	0.0004	2.1013
	1	37	536	0.2	0.0004	2.10137
	28	1	51	0.2	0.0004	2.1014
	28	29	85	0.2	0.0004	2.10141
	28	30	97	0.2	0.0004	2.10161
	31	32	244.1	0.6	0.0004	2.10144
	32	33	144	0.2	0.0004	2.10158
	32	34	332	0.6	0.0004	2.10155
	34	35	413.5	0.6	0.0004	2.10174
	59	58	253	0.8	0.0004	2.10131
	58	57	90	0.8	0.0004	2.10137
	57	56	143.3	0.8	0.0004	2.10146
	56	55	197.5	0.8	0.0004	2.10175
	28	31	36	0.6	0.0004	2.10254
	35	36	36.5	0.7	0.0004	2.0881																																																																																																																																																																																																	
];

% node 类型1是负荷节点、类型2是非平衡热源节点、类型0是平衡热源节点																																																																																																																																																																																																													
model.node = [																																																																																																																																																																																																														
	% node	load(kJ)	Tsmin	Tsmax	Trmin	Trmax	Prmin	Prmax type																																																																																																																																																																																																					
	1	1 	0	100	0	100	-Inf	Inf	1																																																																																																																																																																																																					
	2	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	3	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	4	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																					
	5	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	6	1	0	100	0	100	-Inf	Inf 1 																																																																																																																																																																																																						
	7	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	8	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	9	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	10	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	11	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	12	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	13	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	14	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	15	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	16	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	17	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	18	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	19	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	20	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	21	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	22	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	23	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	24	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	25	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	26	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	27	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	28	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	29	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	30	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	31	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	32	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	33	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	34	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	35	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	36	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	37	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	38	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	39	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	40	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	41	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	42	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	43	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	44	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	45	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	46	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	47	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	48	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																					
	49	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	50	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	51	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	52	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	53	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	54	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	55	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	56	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	57	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	58	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	59	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	60	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	61	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	62	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	63	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	64	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	65	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	66	1	0	100	0	100	-Inf	Inf 1																																																																																																																																																																																																						
	67	1	0	100	0	100	-Inf	Inf 2																																																																																																																																																																																																						
	68	1	0	100	0	100	-Inf	Inf 2																																																																																																																																																																																																						
	69	1	0	100	0	100	-Inf	Inf 0																																																																																																																																																																																																						
];			


model.water_c = 4182; % (J*kg^(-1)*K^(-1))																																																																																																																																																																																																														
model.water_dens = 1e3; %(kg*m^(-3))																																																																																																																																																																																																														
model.Pbase	=	1e6;	%	base	
model.g = 9.81; % Gravitational acceleration
model.viscosity=.294e-6; % temp is 100deg,unit:m2/s kinematic viscosity

[pipe,node,Cp,rho,Pbase,g,viscosity] = deal(model.pipe,model.node,model.water_c,model.water_dens,model.Pbase,model.g,model.viscosity);
nnodes = length(node(:,1));
npipes = length(pipe(:,1));
nloads = 0;
for i = 1:nnodes
    if node(i,9) == 1
        nloads = nloads + 1;
    end
end

D = pipe(:,4);
lamda = pipe(:,6); 
R = 1./lamda; % heat resistance
rough = pipe(:,5);
s = pi*D.*D/4; % area of cross-section of pipe   
len = pipe(:,3);

% network incidence matrix Ah
A = [];
for i = 1:npipes
    A = [A; pipe(i,1) i -1];
end
for i = 1:npipes
    A = [A; pipe(i,2) i 1];
end
A=A';  
Ah = sparse(A(1,:),A(2,:),A(3,:),nnodes,npipes);

Aa = full(Ah);
index_refNode = find(node(:,9) == 0); % reference node
Aa(index_refNode,:) = []; % basic node-branch incidence matrix

nloops = 2 - nnodes + npipes -1; % number of basic loops; Euler's Planar Formula

if nloops > 1e-06 
	isAnnular = 1; % Annular network
	isDendritic = 0;
	Arref = rref(Aa); % Reduced row echelon form
	At = Arref; % block matrix A ~ [At Al]
	Al = Arref;
	index_t = find(sum(abs(Arref)) == 1); % the column index of At and Bt
	index_l = setdiff(1:npipes,index_t); % the column index of Al and Al_b_l
	At(:,index_l) = []; % tree matrix
	Al(:,index_t) = []; % lian zhi matrix
	Bt = (-At\Al)'; % Bt' = -At\Al
	% change the column order
	Bh = zeros(npipes-rank(Aa),npipes);
	Bh(:,index_t) = Bt;
	Bh(:,index_l) = eye(numel(index_l));
else 
	isAnnular = 0;
	isDendritic = 1; % Dendritic network
	Bh = 0;
end
% test = 0?
% verify = Aa*Bh'; 