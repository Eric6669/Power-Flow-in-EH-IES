% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-06 09:10:17
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-08 19:59:05
%  */
clear;clc;
HEAT_DATA;
ELECTRICITY_DATA;
INIT_DATA;

% Phi load
Phi = -[10.0294448	9.957627727	40.273569375	13.10859375	11.393594125...	
10.923405292	2.343467	5.127948	10.021023333	13.09284013...	
20.405475779	6.8384421	3.187027754	44.06466183	0.094882667...	
2.3901983	2.06179995	2.482797158	3.764163375	0.695273367	...
5.7730113	2.376980375	2.930970992	5.905707488	7.852104925...
8.6056397	7.272569967	1.201022375	0.7732263	5.063259234...
50.8837266	6.75106659	28.33800991	5.64171148	8.965308196	...
7.086144086	2.059494	25.58327568	7.388423982	3.10164321	4.07140726 13.15310534	6.009421514...
8.74262627	10.9273357	8.317739864	8.028660155	7.694059215	12.673768	11.01184383	5.346804418...
8.556597763	27.45703723	9.469866555	2.347831324	23.47174068	3.731597295	4.793815673	2.123330851...
5.818911986	2.782336939	28.27255014	29.89465381	6.255552299	11.91517176	13.20652273]';
% P load
pf = 1*ones(n,1); % power factor at each bus except slack bus
qf = sqrt(1-pf.*pf);
Pd = [0 8.57740007	4.5612167	9.82378297	3.39667658	3.317028401...	
9.41104592	5.711906167	0 0 8.081447453	3.700230923...	
3.483253407	6.70720037	7.312133487	10.56469366	5.498683697	9.752618747	6.175174477]';
Qd = Pd./pf(1:n).*qf(1:n); % Qd the reactive power of PQ bus(nPQ variables)


% CHP1 electricity slack node
Phi_chp1_eslack0 = 100;

% CHP3
P_chp3 = 10;
Cm3 = 1.0/.79;
Phi_chp3 = P_chp3*Cm3;

% Error
Total_Herr = inf*ones(20,1);
Total_Perr = Total_Herr;
t = 1;

Phichp1_totalerr = 1;
while Phichp1_totalerr > 1e-3
    tic
    %% -----------------------------------HS power flow------------------------------------- %%
    hydraulic_err = 1;
    thermal_err = 1;
    while hydraulic_err>1e-3||thermal_err>1e-3
        % update dT
        dT = [Ts(1:nloads)-To(1:nloads);Ts(nloads+1:nnodes-1)-Tr(nloads+1:nnodes-1)];      
        mfeedback = sign(m);
        % update Ah,Bh,m
        Ah = Ah.*repmat(mfeedback',nnodes,1);
        Bh = Bh.*mfeedback';
        m = abs(m);
        % update m_node
        m_node = Ah*m;
        
        % nodal flow  mismatches 
        dPhi = Cp/Pbase*m_node(1:nnodes-1).*dT-[-Phi;-Phi_chp1_eslack0;-Phi_chp3];
        % nodal flow  mismatches Jacobian 
        JPhi = Cp/Pbase*Ah(1:nnodes-1,:).*repmat(dT,1,npipes);

        % Kf
        v = m./s/rho;
        Re = abs(v).*D./viscosity;
        for ire = 1:npipes
            if Re(ire) < 2300
                factor(ire) = 64/Re(ire); 
            else
                factor(ire) = colebrook(Re(ire),rough(ire)./D(ire));
            end
        end
        Kf = factor'.*len./D./s./s/2/g/rho/rho;
        
        % loop pressure equation
        dhf = Bh*(Kf.*abs(m).*m);
        % loop pressure equation Jacobian
        Jhf = 2*Bh.*(Kf'.*abs(m'));

        % Jacobian update
        d_hydraulic = [dPhi;dhf];
        J_hydraulic = [JPhi;Jhf];
        delta_m = -J_hydraulic\d_hydraulic;
        hydraulic_err = max(abs(d_hydraulic));
        m = m + delta_m;
        % pipe temperature fall coefficient
        Kt = exp(-lamda.*len./Cp./m);
        % network topology description nd
        nd = Net_Topo(npipes,nnodes,Ah);
        [newTs,newTo,newTr]=STEADY_TSOL(Kt,Ah,m,Ts,To,nd,nloads);
        thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
        Tr=newTr;
        Ts=newTs;
        To=newTo;
        Herr = max(hydraulic_err,thermal_err);
        Total_Herr(t,:) = Herr;
        t = t+1;
    end

    m_node = Ah*m;
    Phi_chp2_hslack = Cp/Pbase*(-m_node(nnodes))*(Ts(nnodes)-Tr(nnodes));

    % CHP2
    z2 = 8.1; 
    Pcon2 = 80;
    P_chp2_hslack = Pcon2-Phi_chp2_hslack/z2;

    %% -----------------------------------ES power flow------------------------------------- %%
    % Active power and reactive power in each bus
    P = zeros(n,1);
    Q = zeros(n,1);
    P_CHP = [P_chp2_hslack; 0; P_chp3];
    for k = 1:length(gen(:,1))
        Q(gen(k,1)) = gen(k,3);
        P(gen(k,1)) = P_CHP(k);
    end
    P = P-Pd;
    Q = Q-Qd;

    % delete slackbus P&Q and PVbus Q
    P_n_1 = P;
    Q_nPQ = Q;
    P_n_1(slackbus) = [];
    Q_nPQ(slackbus) = -inf;
    for k = 1:nPV
        Q_nPQ(PVbus(k)) = -inf;
    end
    Q_nPQ([find(Q_nPQ(:) == -inf)]) = [];
    [P_n_1,Q_nPQ]=deal(P_n_1'/baseMVA,Q_nPQ'/baseMVA);
    % disp('P_n_1,Q_nPQ matrix:');disp(P_n_1);disp(Q_nPQ);

    % delete slackbus U&cita and PVbus U
    U_nPQ = U;
    cita_n_1 = cita;
    U_nPQ(slackbus) = -inf;
    cita_n_1(slackbus) = [];
    for k = 1:nPV
        U_nPQ(PVbus(k)) = -inf;
    end
    U_nPQ([find(U_nPQ(:) == -inf)]) = [];
    % disp('U_nPQ,cita_n_1 matrix:');disp(U_nPQ);disp(cita_n_1);

    % Ubalanced dP&dQ caculation
    Pi = zeros(1,n);
    Qi = zeros(1,n);
    for i = 1:n
        for j = 1:n
            Pi(i) = Pi(i) + U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
            Qi(i) = Qi(i) + U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
        end
    end

    % delete slackbus Pi&Qi and PVbus Qi
    Pi_n_1 = Pi;
    Qi_nPQ = Qi;
    Pi_n_1(slackbus) = [];
    Qi_nPQ(slackbus) = -inf;
    for k = 1:nPV
        Qi_nPQ(PVbus(k)) = -inf;
    end
    Qi_nPQ([find(Qi_nPQ(:) == -inf)]) = [];
    dP_n_1 = P_n_1-Pi_n_1; % dP n-1  PQ + PV
    dQ_nPQ = Q_nPQ-Qi_nPQ; % dQ nPQ  PQ
    % disp('Ubalanced active power:dP_n_1');disp(dP_n_1);
    % disp('Ubalanced reactive power:dQ_nPQ');disp(dQ_nPQ);

    Perr = 1;
    while Perr > 1e-3
        % fprintf('%d iteration\n',iter_time);
        J=Jacobi_matrix(n,nPV,U,cita,B,G,Pi,Qi,slackbus,PVbus);
        dPQ = [dP_n_1 dQ_nPQ]';
        dUcita = (-inv(J)*dPQ)';      % inv(J)*dPQ  J\dPQ

        % cita correction
        dcita = dUcita(1:n-1);        % dcita1-dcita2-...-dcita(n-1)
        cita_n_1 = cita_n_1 + dcita; 

        % voltage correction
        U_amplitude = zeros(nPQ,nPQ);
        for i = 1:nPQ
            U_amplitude(i,i) = U_nPQ(i);
        end
        dU = (U_amplitude*dUcita(n:n+nPQ-1)')';  % 后nPQ对应电压的修正量
        U_nPQ = U_nPQ + dU;
       
        % Append slack bus and PV bus
        % citaslack = cita + slackbus cita 
        citaslack = zeros(1,n);
        index = (1:n)';index(slackbus) = [];
        for i = 1:length(index)
            citaslack(index(i)) = cita_n_1(i);
        end
        citaslack(slackbus) = 0;
        cita = citaslack;

        % UPVslack = U + PVbus U + slackbus U
        UPVslack = zeros(n,1);
        index = (1:n)';
        index(slackbus) = -inf;
        for k = 1:nPV
            index(PVbus(k)) = -inf;
        end
        index(index(:) == -inf) = [];

        for i = 1:length(index)
            UPVslack(index(i)) = U_nPQ(i);
        end
        UPVslack(slackbus) = 1.05;
        for k = 1:nPV
            UPVslack(PVbus(k)) = 1.05;
        end
        U = UPVslack';

        % Error calculation
        Pi = zeros(1,n);
        Qi = zeros(1,n);
        for i = 1:n
            for j = 1:n
                Pi(i) = Pi(i) + U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
                Qi(i) = Qi(i) + U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
            end
        end
        % delete slackbus Pi&Qi and PVbus Qi
        Pi_n_1 = Pi;
        Qi_nPQ = Qi;
        Pi_n_1(slackbus) = [];
        Qi_nPQ(slackbus) = -inf;
        for k = 1:nPV
            Qi_nPQ(PVbus(k)) = -inf;
        end
        Qi_nPQ([find(Qi_nPQ(:) == -inf)]) = [];
        dP_n_1 = P_n_1-Pi_n_1; % dP n-1  PQ + PV
        dQ_nPQ = Q_nPQ-Qi_nPQ; % dQ nPQ  PQ

        % if (max(abs(dP_n_1))<1e-3 && max(abs(dQ_nPQ))<1e-3 )
        %     disp('power flow converge!!!');
        %     break
        % end
        Perr = max([abs(dP_n_1)';abs(dQ_nPQ)']);
        Total_Perr(t,:) = Perr;
        t = t+1;
    end
    disp('power flow converge!!!');
    % % electricity slack node 
    % U_cita = U.*exp(1j*cita);
    % S = U_cita.*conj(Y*U_cita');
    % S_slack = S(slackbus);
    % P_chp1_eslack = real(S_slack);
    P_chp1_eslack = Pi(slackbus)*baseMVA;
    
    % CHP1 electricity slack node
    Cm1=1.3;
    Phi_chp1_eslack = P_chp1_eslack*Cm1;

    % thermal power of the electrical slack node
    Phichp1_totalerr = abs(Phi_chp1_eslack0 - Phi_chp1_eslack);
    Phi_chp1_eslack0 = Phi_chp1_eslack;
end
toc
Ts=Ts+Ta;Tr=Tr+Ta;To=To+Ta;
Ploss = P_chp2_hslack + P_chp1_eslack + P_chp3 + sum(-Pd); 
Philoss = Phi_chp1_eslack0 + Phi_chp2_hslack + Phi_chp3 + sum(Phi);


disp('mass flow rates within each pipe');disp(m);
disp('flow injection at the node');disp(m_node);
disp('supply temperature');disp(Ts);
disp('return temperature');disp(Tr);  
disp('U:');disp(U');
disp('cita:');disp(rad2deg(cita)');
fprintf('CHP1(Electricity slack node) Phi:%f\t,P:%f\n',Phi_chp1_eslack,P_chp1_eslack);
fprintf('CHP2(Heat slack node)        Phi:%f\t,P:%f\n',Phi_chp2_hslack,P_chp2_hslack);
fprintf('CHP3(None slack node)        Phi:%f\t,P:%f\n',Phi_chp3,P_chp3);
fprintf('electricity losses (MWe)     %f\n',Ploss);
fprintf('heat losses (MWth)           %f\n',Philoss);
Vis_EH_plot(m,m_node,Ts,Tr,U,cita,Total_Herr,Total_Perr);