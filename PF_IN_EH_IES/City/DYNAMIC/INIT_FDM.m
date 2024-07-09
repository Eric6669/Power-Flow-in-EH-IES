% /*
%  * @Descripttion: FDM parameters for EH-IES CITY power flow
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-08 15:48:36
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-08 19:59:49
%  */

tau = 60; % dt = 1min
h = 50; % dx = 50m
M = round(len./h); % 1<=x<=M


tsim = 3600*2; % total simulation time 1h
N = tsim/tau; % 1<=t<=N  N=120


%% pipe temperature
Pipe_Ta = -10;
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
    Tsi0 = 80*ones(M(k)+1,1)-Pipe_Ta;
    Pipe_Ts{k}(:,1) = Tsi0;
    Tri0 = 40*ones(M(k)+1,1)-Pipe_Ta;
    Pipe_Tr{k}(:,1) = Tri0;
end

% Boundary conditions
% Ts,source temperature
for k = 1:npipes
    % supply source temperature
    Pipe_Ts{k}(1,2:N+1) = 80*ones(1,N)-Pipe_Ta;
end

% Tr,load temperature
for k = 1:npipes
    % return loads temperature
    Pipe_Tr{k}(1,2:N+1) = 40*ones(1,N)-Pipe_Ta;
end

% To temperature
% Pipe_To{k}(:,1:N+1) = Pipe_Tr{k}(:,1:N+1);
