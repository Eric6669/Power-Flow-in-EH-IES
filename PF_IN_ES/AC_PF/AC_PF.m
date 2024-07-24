% /*
%  * @Descripttion: This function calculates the AC power flow for a given power system.
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 20:12:59
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 20:13:03
%  */

function [U,cita,Sij] = AC_PF(mpc)

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
    % disp('P_n_1,Q_nPQ matrix:');disp(P_n_1);disp(Q_nPQ);

    %% voltage and its angle
    U = bus(:,8)'; 
    cita = bus(:,9)';
    cita = (deg2rad(cita));  % degree2rad

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

    %% node admittance matrix 
    Y = Y_matrix(n,branch); 
    G = real(Y);
    B = imag(Y);
    % disp('node admittance matrix:');disp(Y);

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
    % disp('Ubalanced active power:dP_n_1');disp(dP_n_1);
    % disp('Ubalanced reactive power:dQ_nPQ');disp(dQ_nPQ);

    %% Newton-Rapshon method 
    iter_time = 1;
    limit = 50;
    epr = 1e-5; % ������������
    while iter_time < limit
        fprintf('%d iteration\n',iter_time);
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
        dU = (U_amplitude*dUcita(n:n+nPQ-1)')';  % ��nPQ��Ӧ��ѹ��������
        U_nPQ = U_nPQ + dU;
       
        % Append slack bus and PV bus
        % citaslack = cita + slackbus cita 
        citaslack = zeros(n,1);
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
        UPVslack(slackbus) = 1.04;
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

        if (max(abs(dP_n_1))<epr && max(abs(dQ_nPQ))<epr )
            disp('power flow converge!!!');
            break
        end
        iter_time = iter_time+1;
    end

    % admittance matrix y
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
    Sij = Line_power(n,y,U,cita);

    disp('U:');disp(U);
    disp('cita:');disp(rad2deg(cita));
    disp('power flow:');disp(sparse(Sij))
end


