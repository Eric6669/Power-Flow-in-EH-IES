% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-06 16:58:25
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-06 20:08:29
%  */
m = ones(npipes,1);
% node temperature
Ta = -10*ones(nnodes,1);
% assume load supply temperature 
Ts = 80*ones(nnodes,1)-Ta; 
% assume source return temperature
To = 40*ones(nnodes,1)-Ta;
Tr = To; 

U = 1.05*ones(1,n);
cita = zeros(1,n);
cita = (deg2rad(cita));  % degree2rad