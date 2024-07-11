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

% Phi load
Phi = -[
    0
    0
    6.306307429
    0
    4.036036755
    12.36036256
    0
    0
    0
    0
    30.52252796
    0
    3.3456
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    3.3456
    0
    6.306307429
    6.306307429
    0
    0
    3.3456
    0
    0
    0
    0
    0
    3.3456
    3.3456
    0
    3.3456
    0
    0
    0
    3.3456
    3.3456
    3.3456
    0
    0
    0
    4.036036755
    0
    228.2694063
    0
    0
    0
    0
    0
    0
    6.306307429
    0
    0
    0
    0
    0    
];

% P load
pf = 1*ones(n,1); % power factor at each bus except slack bus
qf = sqrt(1-pf.*pf);
Pd = 3*[0 8.57740007	4.5612167	9.82378297	3.39667658	3.317028401...	
9.41104592	5.711906167	0 0 8.081447453	3.700230923...	
3.483253407	6.70720037	7.312133487	10.56469366	5.498683697	9.752618747	6.175174477]';
Qd = Pd./pf(1:n).*qf(1:n); % Qd the reactive power of PQ bus(nPQ variables)