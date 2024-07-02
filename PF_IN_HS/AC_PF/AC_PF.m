% /*
%  * @Descripttion: This function calculates the AC power flow for a given power system.
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 20:12:59
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 20:13:03
%  */

function [U,cita,Pi,Qi,Sij] = AC_PF(mpc)

    %% Data init
    [baseMVA,bus,gen,branch]=deal(mpc.baseMVA,mpc.bus,mpc.gen,mpc.branch);
    n = length(bus(:,1));
    slackbus = find(mpc.bus(:,2)==3); 
    nslack = length(slackbus);
    PVbus = find(mpc.bus(:,2)==2);
    nPV = length(PVbus);
    nPQ = n-nPV-nslack;

    %% Active power and reactive power in each bus
    P = zeros(n,1);
    Q = zeros(n,1);
    for k = 1:length(gen(:,1))
        P(gen(k,1)) = gen(k,2);
        Q(gen(k,1)) = gen(k,3);
    end
    P = P-bus(:,3);
    Q = Q-bus(:,4);

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
    disp('P_n_1,Q_nPQ matrix:');disp(P_n_1);disp(Q_nPQ);

    %% voltage and its angle
    U = bus(:,8)'; 
    cita = bus(:,9)';
    cita = (deg2rad(cita));  % degree2rad

    % delete slackbus U&cita and PVbus U
    U_n_1 = U;
    cita_nPQ = cita;
    U_n_1(slackbus) = -inf;
    cita_nPQ(slackbus) = [];
    for k = 1:nPV
        U_n_1(PVbus(k)) = -inf;
    end
    U_n_1([find(U_n_1(:) == -inf)]) = [];
    disp('U_n_1,cita_nPQ matrix:');disp(U_n_1);disp(cita_nPQ);

    %% node admittance matrix 
    Y = Y_matrix(n,branch); 
    G = real(Y);
    B = imag(Y);
    disp('node admittance matrix:');disp(Y);

    %% Ubalanced dP&dQ caculation
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
    disp('Ubalanced active power:dP_n_1');disp(dP_n_1);
    disp('Ubalanced reactive power:dQ_nPQ');disp(dQ_nPQ);

    

    %% Newton-Rapshon method 
    iter_time = 1;
    limit = 50;
    epr = 1e-5; % 迭代收敛精度
    while iter_time < limit
        fprintf('%d iteration\n',iter_time);



        J=Jacobi_matrix(n,nPQ,nPV,U,cita,B,G,Pi,Qi,slackbus,PVbus);
        dPQ = [dP_n_1 dQ_nPQ]';
        dUcita = (-inv(J)*dPQ)';      % inv(J)*dPQ  J\dPQ



        % cita correction
        dcita = dUcita(1:n-1);        % dcita1-dcita2-...-dcita(n-1)
        cita_nPQ = cita_nPQ + dcita; 


        % voltage correction
        U_amplitude = zeros(nPQ,nPQ);
        for i = 1:nPQ
            U_amplitude(i,i) = U_n_1(i);
        end
        dU = (Us*dUcita(n:n+nPQ-1)')';  % 后m对应电压的修正量
        U_n_1 = U_n_1 + dU;
       
        % Append slack bus and PV bus
        % citaslack = cita + slackbus cita 
        citaslack = zeros(n,1);
        index = (1:n)';index(slackbus) = [];
        for i = 1:length(index)
            citaslack(index(i)) = cita_nPQ(i);
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
            UPVslack(index(i)) = U(i);
        end
        UPVslack(slackbus) = 1;
        for k = 1:nPV
            UPVslack(PVbus(k)) = bus(PVbus(k),8);
        end
        U = UPVslack;

        % Error calculation
        Pi = zeros(1,n);
        Qi = zeros(1,n);
        for i = 1:n
            for j = 1:n
                Pi(i) = Pi(i) + U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
                Qi(i) = Qi(i) + U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
            end
        end
        Pi(slackbus) = [];
        Qi(slackbus) = -inf;
        for k = 1:nPV
            Qi(PVbus(k),:) = -inf;
        end
        Qi([find(Qi(:) == -inf)], :) = [];
        dP = P-Pi; % dP n-1  PQ + PV
        dQ = Q-Qi; % dQ nPQ  PQ

        if (max(abs(dP))<epr && max(abs(dQ))<epr )
            disp('power flow converge!!!');
            break
        end
        iter_time = iter_time+1;
    end

    

    y=zeros(n,n);
    for i=1:n
        for j=1:n
            if i~=j
                y(i,j)=-Y(i,j);
            else
                y(i,j)=sum(Y(i,:));
            end
        end
    end
    Sij = line_power(n,y,U,cita);

    fprintf('迭代总次数：%d\n', iter_time);
    disp('节点电压幅值：');
    disp(U);
    disp('节点电压相角：');
    disp(rad2deg(cita));
    % disp('节点注入有功计算结果：');
    % disp(Pi);
    % disp('节点注入无功计算结果：');
    % disp(Qi);
    disp('支路功率计算结果：');
    disp(sparse(Sij))
end


