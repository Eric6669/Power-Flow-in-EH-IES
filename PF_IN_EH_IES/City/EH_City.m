% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-03 15:12:31
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-04 14:31:04
%  */

% CHP1 electricity slack node
Phi_chp1_eslack0 = 200;

% CHP3
P_chp3 = 20;
Cm3 = 1.0/.79;
Phi_chp3 = P_chp3*Cm3;

Phichp1_totalerr = 1;
while Phichp1_totalerr > 1e-3
    % volume adjustment in HS
    [m,m_node,Ts,Tr,To,Herr,Phi_chp2_hslack] = EH_STEADY_PF(Data_HS_City,1,Phi_chp1_eslack0,Phi_chp3,m,Ts,Tr,To,Ta);
    % [m,m_node,Ts,Tr,Herr,Phi_chp2_hslack] = EH_DYNAMIC_PF(Data_HS_City,1,Phi_chp1_eslack0,Phi_chp3);
    
    % CHP2 heat slack node
    z2 = 8.1; 
    Pcon2 = 100;
    P_chp2_hslack = Pcon2-Phi_chp2_hslack/z2;

    % power flow in ES
    [U,cita,Sij,P_chp1_eslack,Perr] = EH_AC_PF(Data_ES_City,P_chp2_hslack,P_chp3,U,cita);

    % CHP1 electricity slack node
    Cm1=1.3;
    Phi_chp1_eslack = P_chp1_eslack*Cm1;

    % thermal power of the electrical slack node
    Phichp1_totalerr = abs(Phi_chp1_eslack0 - Phi_chp1_eslack);
    Phi_chp1_eslack0 = Phi_chp1_eslack;
end 

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
Vis_EH_plot(m,m_node,Ts,Tr,U,cita,Herr,Perr);
% save('./EH_City_result.mat','m','m_node','Ts','Tr','Herr','U','cita','Perr');