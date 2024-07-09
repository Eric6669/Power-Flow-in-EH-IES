% /*
%  * @Descripttion: Vis for dynamic
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-09 10:47:32
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-09 10:49:55
%  */
function EH_DYNAMIC_VIS(m,m_node,Ts,Tr,U,cita,Pipe_Ts,Pipe_Tr,npipes,M)
    
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

    % figure(4);
    % plot(Total_Herr,'-*r');
    % hold on;
    % plot(Total_Perr,'-db');
    % hold on;
    % xlabel('Iteration number');
    % ylabel('error');
    % legend('Heat error','Electricity error');
    % grid on;

    figure(4);
    for k = 1:npipes
        plot(Pipe_Ts{k}(M(k)+1,:),':','LineWidth',1.8); 
        hold on;
    end
    figure(5);
    for k = 1:npipes
        plot(Pipe_Tr{k}(M(k)+1,:),'-','LineWidth',1.8); 
        hold on;
    end

end