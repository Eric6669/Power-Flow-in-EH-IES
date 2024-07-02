% /*
%  * @Descripttion: This function calculates the DC power flow for a given power system.
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 18:44:03
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 19:31:26
%  */

function [cita,Sij] = DC_PF(mpc)

    %% Data init
    [baseMVA,bus,gen,branch]=deal(mpc.baseMVA,mpc.bus,mpc.gen,mpc.branch);
    n = length(bus(:,1));
    slackbus = find(mpc.bus(:,2)==3);
    
    %% Each node's active power
    P = zeros(n,1);
    for k=1:length(gen(:,1))
        P(gen(k,1)) = gen(k,2);
    end
    P = P-bus(:,3);
    P = deal(P/baseMVA);
    P(slackbus,:) = [];
    disp('Node power matrix:');disp(P);

    %% Construct B matrix (Susceptance matrix)
    B=zeros(n);
    for k=1:length(branch(:,1))
        i=branch(k,1);j=branch(k,2);
        B(i,j) = -1/branch(k,4);
        B(j,i) = B(i,j); 
    end
    for i = 1:n
        B(i,i) = -sum(B(i,:));
    end
    B(slackbus,:) = [];
    B(:,slackbus) = [];
    disp('Susceptanc(电纳) matrix:');disp(B);

    %% P = B*cita Solve for voltage angles
    cita = B\P;
    
    %% Append slack bus
    citaslack = zeros(n,1); % cita include slackbus
    index = (1:n)';index(slackbus) = [];
    for i = 1:length(index)
        citaslack(index(i)) = cita(i);
    end
    citaslack(slackbus) = 0;
    disp('Node voltage angles:');disp(rad2deg(citaslack));

    %% Sij = [cita(i) - cita(j)] / xij
    Sij = zeros(n, n);
    for k=1:length(branch(:,1))
        i=branch(k,1);j=branch(k,2);
        Sij(i,j) = (citaslack(i)-citaslack(j))/branch(k,4);
    end
    disp('Branch power flow results:');disp(sparse(Sij))
    
end

