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

% CHP1 electricity slack node 
% Phi_chp1_eslack = P_chp1_eslack*Cm1;
Phi_chp1_eslack0 = 150;
Cm1 = 1.3;

% CHP3
P_chp3 = 83.74;
Cm3 = 1.0/.79;
Phi_chp3 = P_chp3*Cm3;

% Error
Total_Herr = inf*ones(20,1);
Total_Perr = Total_Herr;
t = 1;

Ah0=Ah;
Phichp1_totalerr = 1;
while Phichp1_totalerr > 1e-3
    tic
    %% -----------------------------------HS power flow------------------------------------- %%
    hydraulic_err = 1;
    thermal_err = 1;
    while hydraulic_err>1e-3||thermal_err>1e-3
        % update dT
        dT = [Ts(1:nloads)-To(1:nloads);Ts(nloads+1:nnodes-1)-Tr(nloads+1:nnodes-1)];      
        % mfeedback = sign(m);
        % % update Ah,Bh,m
        % Ah = Ah.*repmat(mfeedback',nnodes,1);
        % Bh = Bh.*mfeedback';
        % m = abs(m);
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
    Pcon2 = 200;
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
    Phi_chp1_eslack = P_chp1_eslack*Cm1;

    % thermal power of the electrical slack node
    Phichp1_totalerr = abs(Phi_chp1_eslack0 - Phi_chp1_eslack);
    Phi_chp1_eslack0 = Phi_chp1_eslack;
end
toc
Ts=Ts+Ta;Tr=Tr+Ta;To=To+Ta;
[row,sign_flow]=find(Ah-Ah0); % return the correct flow direction
m(sign_flow)=-m(sign_flow);
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
EH_STEADY_VIS(m,m_node,Ts,Tr,U,cita,Total_Herr,Total_Perr);