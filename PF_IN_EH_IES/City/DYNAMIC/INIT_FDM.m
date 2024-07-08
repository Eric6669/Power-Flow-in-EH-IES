% /*
%  * @Descripttion: FDM parameters for EH-IES CITY power flow
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-08 15:48:36
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-08 19:59:49
%  */

v = [
    2.9700
    1.2605
    3.2113
    0.5575
    0.1737
    2.9844
    2.9696
    2.9372
    0.2147
    0.4423
    2.6591
    0.2040
    2.5951
    2.3161
    2.3155
    4.3492
    7.6879
    7.6350
   13.4307
    7.5399
    7.4169
    7.3662
    7.3037
    7.1777
    7.0103
    0.9057
    6.6695
    1.8672
    0.3631
    0.2365
    0.2675
    0.3259
    0.3329
    0.0766
    0.2650
    1.5529
    0.4490
    1.0314
    1.5125
    0.4953
    1.6246
    0.6282
    0.4654
    0.1959
    1.1185
    0.0962
    0.8816
    0.7810
    0.6902
    0.6426
    0.4144
    2.2601
    0.0869
    2.9888
    0.7609
    0.7805
    1.9701
    0.0168
    0.1084
    3.4624
    2.4222
    2.7128
    2.5923
    2.9236
    2.7762
    2.6614
    1.9394
    4.5484
    1.7637];
tau = 60; % dt = 1min
h = 50; % dx = 100m
M = round(len./h); % 1<=x<=M


tsim = 3600*1; % total simulation time 1h
N = tsim/tau; % 1<=t<=N  N=60

alpha = v.*tau./h;
belta = v.*tau./(2.*m.*Cp.*R);
omega1 = (ones(npipes,1)+alpha-belta)./(ones(npipes,1)+alpha+belta);
omega2 = (alpha-ones(npipes,1)-belta)./(ones(npipes,1)+alpha+belta);
omega3 = (ones(npipes,1)-alpha-belta)./(ones(npipes,1)+alpha+belta);

% omega2^M for each pipe
W = zeros(npipes,1);
for k = 1:npipes
    W(k) = omega2(k)^M(k);
end

% A for each pipe
Ak = cell(npipes,1);
for k = 1:npipes
    Ak{k} = zeros(1,M(k));
    for x = 1:M(k)
        Ak{k}(x) = omega2(k)^(M(k)-x);
    end
end

%% pipe temperature
Ta = 0;
% cell Pipe_Ts Pipe_Tr Pipe_To
Pipe_Ts = cell(npipes,1);
Pipe_Tr = cell(npipes,1);
Pipe_To = cell(npipes,1);
for k = 1:npipes
    % supply temperature
    Pipe_Ts{k} = zeros(M(k)+1,N+1);
    % return temperature
    Pipe_Tr{k} = zeros(M(k)+1,N+1);
    Pipe_To{k} = zeros(M(k)+1,N+1);
end

% Initial conditions(t=0)
% differ in different system 
for k = 1:npipes
    Tsi0 = 80*ones(M(k)+1,1)-Ta;
    Pipe_Ts{k}(:,1) = Tsi0;
    Tri0 = 40*ones(M(k)+1,1)-Ta;
    Pipe_Tr{k}(:,1) = Tri0;
end

% Boundary conditions
% Ts,source temperature
for k = 1:npipes
    % supply source temperature
    Pipe_Ts{k}(1,2:N+1) = 80*ones(1,N)-Ta;
end

% Tr,load temperature
for k = 1:npipes
    % return loads temperature
    Pipe_Tr{k}(1,2:N+1) = 40*ones(1,N)-Ta;
end

% To temperature
% Pipe_To{k}(:,1:N+1) = Pipe_Tr{k}(:,1:N+1);

Bs = cell(npipes,1);
Br = cell(npipes,1);
b_supply = cell(npipes,1);
b_return = cell(npipes,1);
for k = 1:npipes
    Bs{k} = zeros(M(k),1);
    Br{k} = zeros(M(k),1);
    b_supply{k} = zeros(M(k),1);
    b_return{k} = zeros(M(k),1);
end