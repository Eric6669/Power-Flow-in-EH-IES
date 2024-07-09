% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-04 13:36:22
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-04 13:36:37
%  */

function EH_STEADY_VIS(m,m_node,Ts,Tr,U,cita,Total_Herr,Total_Perr)
    
    figure(1);
    subplot(2,1,1);
    plot(m,'-ob');
    xlabel('Pipe number');
    ylabel('mass flow rates');
    grid on;
    subplot(2,1,2);
    plot(m_node,'-db','MarkerFaceColor','b');
    xlabel('Node number');
    ylabel('flow injection');
    grid on;


    figure(2);
    subplot(2,1,1);
    plot(Ts,'-^b');
    xlabel('Node number');
    ylabel('Supply temperature');
    grid on;
    subplot(2,1,2);
    plot(Tr,'-^b','MarkerFaceColor','b');
    xlabel('Node number');
    ylabel('Return temperature');
    grid on;

    figure(3);
    subplot(2,1,1);
    plot(rad2deg(cita),'-^m');
    xlabel('Node number');
    ylabel('voltage angle (deg)');
    grid on;
    subplot(2,1,2);
    plot(U,'-vm','MarkerFaceColor','m');
    xlabel('Node number');
    ylabel('voltage magnitude');
    grid on;

    figure(4);
    plot(Total_Herr,'-*r');
    hold on;
    plot(Total_Perr,'-db');
    hold on;
    xlabel('Iteration number');
    ylabel('error');
    legend('Heat error','Electricity error');
    grid on;





end