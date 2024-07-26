function mpc = case5
    %CASE5  Power flow data for modified 5 bus, 5 gen case based on PJM 5-bus system
    %   Please see CASEFORMAT for details on the case file format.
    %
    %   Based on data from ...
    %     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
    %     Proceedings of the 2010 IEEE Power & Energy Society General Meeting
    
    %   Created by Rui Bo in 2006, modified in 2010, 2014.
    %   Distributed with permission.
    
    %   MATPOWER
    
    %% MATPOWER Case Format : Version 2
    mpc.version = '2';
    
    %%-----  Power Flow Data  -----%%
    %% system MVA base
    mpc.baseMVA = 100;
    
    %% bus data
    %	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
    mpc.bus = [
        1	3	0	0	0	0	1	1.06	0	230	1	1.1	0.9;
        2	1	-20	-20	0	0	1	1	0	230	1	1.1	0.9;
        3	1	45	15	0	0	1	1	0	230	1	1.1	0.9;
        4	1	40	5	0	0	1	1	0	230	1	1.1	0.9;
        5	1	60	10	0	0	1	1	0	230	1	1.1	0.9;
    ];
    
    %% generator data
    %	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    mpc.gen = [
        1	40	0	200	-200	1.06	100	1	300	0	0	0	0	0	0	0	0	0	0	0	0;
    ];
    
    %% branch data
    %	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
    mpc.branch = [
        1	2	0.02	0.06	0	400	400	400	0	0	1	-360	360;
        1	3	0.08	0.24	0	0	0	0	0	0	1	-360	360;
        2	3	0.06	0.18	0	0	0	0	0	0	1	-360	360;
        2	4	0.06	0.18	0	0	0	0	0	0	1	-360	360;
        2	5	0.04	0.12	0	0	0	0	0	0	1	-360	360;
        3	4	0.01	0.03	0	240	240	240	0	0	1	-360	360;
        4	5	0.08	0.24	0	0	0	0	0	0	1	-360	360;
    ];
    
    