clear;clc;
HEAT_DATA;ELECTRICITY_DATA;INIT_DATA;
% CHP1 electricity slack node
Phi_chp1_eslack0 = 200;

% CHP3
P_chp3 = 20;
Cm3 = 1.0/.79;
Phi_chp3 = P_chp3*Cm3;

% Error
Total_Herr = inf*ones(20,1);
Total_Perr = Total_Herr;
t = 1;

Phichp1_totalerr = 1;
while Phichp1_totalerr > 1e-3
    hydraulic_err = 1;
    thermal_err = 1;
    tic
    % -----------------------------------------HS-----------------------------------------% 
    while hydraulic_err>1e-3||thermal_err>1e-3
        %% hydraulic model
        % update dT
        dT = [Ts(1:nloads)-To(1:nloads);Ts(nloads+1:nnodes-1)-Tr(nloads+1:nnodes-1)];      
        mfeedback = sign(m);
        % update Ah,Bh,m
        Ah = Ah.*repmat(mfeedback',nnodes,1);
        Bh = Bh.*mfeedback';
        m = abs(m);
        update m_node
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
        
        %% thermal model
        % pipe temperature fall coefficient
        Kt = exp(-lamda.*len./Cp./m);
        %% network topology description nd
        nd = Net_Topo(npipes,nnodes,Ah);
        [newTs,newTo,newTr]=STEADY_TSOL(Kt,Ah,m,Ts,To,nd,nloads);
        thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
        Tr=newTr;
        Ts=newTs;
        To=newTo;
        Total_Herr(t,:) = max(hydraulic_err,thermal_err);
        t = t+1;
    end
    
    % CHP2
    m_node = Ah*m;
    Phi_chp2_hslack = Cp/Pbase*(-m_node(nnodes))*(Ts(nnodes)-Tr(nnodes));
    z2 = 8.1; 
    Pcon2 = 100;
    P_chp2_hslack = Pcon2-Phi_chp2_hslack/z2;
    
    Perr = 1;
    while Perr > 1e-3
        U_cita = U.*exp(1j*cita);
        S = U_cita.*conj(Y*U_cita);
        dP = real(S(1:n-1))-[P;P_chp2_hslack;P_chp3];
        dQ = imag(S(1:nPQ))-[Q;]; % CHP power factor P_coupling*qf(n-1)/pf(n-1)
        % compute H N K L
        HK = -1j*diag(U_cita)*conj(Y*diag(U_cita))+diag(1j*S);
        EJ1 = real(HK(1:n-1,1:n-1));
        EJ3 = imag(HK(1:nPQ,1:n-1));
        e = exp(1j*cita);
        NL = diag(U_cita)*conj(Y*diag(e))+diag(e.*conj(Y*U_cita));
        EJ2 = real(NL(1:n-1,1:nPQ));
        EJ4 = imag(NL(1:nPQ,1:nPQ));

        J = [EJ1 EJ2;EJ3 EJ4];
        x = [cita(1:n-1); U(1:nPQ)];
        dx = -J\[dP;dQ;];
        Perr = max(abs([dP;dQ;]));
        x=x+dx;
        cita(1:n-1) = x(1:n-1);
        U(1:nPQ) = x(n:n-1+nPQ);
        
        Total_Herr(t,:)=Perr;
        Total_Perr(t,:)=Perr;
        t = t+1;
    end
    %calculate P_slack and Phi_slack
    U_cita = U.*exp(1j*cita);
    S = U_cita.*conj(Y*U_cita);
    S_slack = S(n); % power of slack bus S_slack=V(n)*conj(Y(n,:)*V)
    P_chp1_eslack = real(S_slack);
    
    % CHP1 electricity slack node
    Cm1=1.3;
    Phi_chp1_eslack = P_chp1_eslack*Cm1;

    % thermal power of the electrical slack node
    Phichp1_totalerr = abs(Phi_chp1_eslack0 - Phi_chp1_eslack);
    Phi_chp1_eslack0 = Phi_chp1_eslack;
end
toc
Ts=Ts+Ta;Tr=Tr+Ta;To=To+Ta;
Ploss = P_chp2_hslack + P_chp1_eslack + P_chp3 + sum(P); 
Philoss = Phi_chp1_eslack0 + Phi_chp2_hslack + Phi_chp3 + sum(Phid);

    
disp('mass flow rates within each pipe');disp(m);
disp('flow injection at the node');disp(m_node);
disp('supply temperature');disp(Ts);
disp('return temperature');disp(Tr);  
disp('U:');disp(U);
disp('cita:');disp(rad2deg(cita));
disp('power flow:');disp(sparse(Sij));
fprintf('CHP1(Electricity slack node) Phi:%f\t,P:%f\n',Phi_chp1_eslack,P_chp1_eslack);
fprintf('CHP2(Heat slack node)        Phi:%f\t,P:%f\n',Phi_chp2_hslack,P_chp2_hslack);
fprintf('CHP3(None slack node)        Phi:%f\t,P:%f\n',Phi_chp3,P_chp3);
fprintf('electricity losses (MWe)\n%f\n heat losses (MWth)\n%f\n',Ploss,Philoss);
Vis_EH_plot(m,m_node,Ts,Tr,U,cita,Herr,Perr);