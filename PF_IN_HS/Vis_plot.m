% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 14:47:11
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-02 21:17:43
%  */
function Vis_plot(m,m_node,Ts,Tr,Herr)
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
    plot(Herr,'-*r');
    xlabel('Iteration number');
    ylabel('Heat error');
    grid on;
end