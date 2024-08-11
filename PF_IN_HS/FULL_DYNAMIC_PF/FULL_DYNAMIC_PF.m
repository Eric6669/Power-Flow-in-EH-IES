% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-08-11 15:20:40
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-08-11 15:21:18
%  */
clc;clear;
% -----------------------------------------------------Data init----------------------------------------------------------------
nnodes = 3;
npipes = 3;
nloads = 2;
D=ones(3,1)*150e-3;%pipe diameter
rough=ones(3,1)*1.25e-3;%pipe roughness
s=pi*D.*D/4;%area of cross-section of pipe
len=[400 400 600]';%length of pipe
rho=958.4; % Density of water (kg/m^3) at 100 degrees C
g=9.81;%Gravitational acceleration
viscosity=.294e-6;%temp is 100deg,unit:m2/s kinematic viscosity
Ah=[1 -1 0; 0 1 1; -1 0 -1];
lamda=[.2 .2 .2]';%thermal conductivity
Pbase=1e6;%base value of power
Cp=4182; %J/kg/K specific capacity of water


Node_m_load = [2; 3];
msteady = [2.712; 0.712; 2.289];
Node_m_all = [2; 3; -5];
Hssteady = [1e6/rho/g; 1e6/rho/g; 1*1e6/rho/g];

a = 1200;
tau = 60; % dt = 1min
tsim = 3600*2; % total simulation time
N = tsim/tau; % 1<=t<=N  N=120
h = 50; % dx = 50m
M = round(len./h); % 1<=x<=M
N_piece = 0;
for k = 1:npipes
    N_piece = N_piece + M(k);
end
%% ----------------------------------- Center-implicit difference method init ----------------------------------- %%
% pipe mass flow
Pipe_m = cell(npipes,1);
% node mass flow
Node_m = cell(nnodes,1);
% pipe pressure H
Pipe_Hs = cell(npipes,1);
Pipe_Hr = cell(npipes,1);
% node pressure H
Node_Hs = cell(nnodes,1);
Node_Hr = cell(nnodes,1);

for k = 1:nnodes
    Node_m{k} = zeros(N+1,1);
    Node_Hs{k} = zeros(N+1,1);
    Node_Hr{k} = zeros(N+1,1);
end

for k = 1:npipes
    Pipe_m{k} = zeros(M(k)+1,N+1);
    Pipe_Hs{k} = zeros(M(k)+1,N+1);
    Pipe_Hr{k} = zeros(M(k)+1,N+1);
end
% ----------------------------------- Initial conditions(t=1) ----------------------------------- %
for k = 1:npipes
    Pipe_m{k}(:,1) = msteady(k);
    Pipe_Hs{k}(:,1) = Hssteady(k);
end
% ----------------------------------- Boundary conditions ----------------------------------- %
for k = 1:npipes
    Pipe_m{k}(1,2:N+1) = msteady(k);
    Pipe_Hs{k}(1,2:N+1) = Hssteady(k);
end
% ----------------------------------- Iteration initial value ----------------------------------- %
for k = 1:npipes
    Pipe_m{k}(2:M(k)+1,2:N+1) = msteady(k);
    Pipe_Hs{k}(2:M(k)+1,2:N+1) = Hssteady(k);
end
for k = 1:nnodes
    Node_m{k}(:,1) = Node_m_all(k);
    Node_Hs{k}(:,1) = Hssteady(k);
end


%% ------------------------------------------------- FDM init ---------------------------------------------- %%
% pipe temperature
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

% ----------------------------------- Initial conditions(t=0) ----------------------------------- %
% differ in different system 
for k = 1:npipes
    Tsi0 = 80*ones(M(k)+1,1)-Pipe_Ta;
    Pipe_Ts{k}(:,1) = Tsi0;
    Tri0 = 40*ones(M(k)+1,1)-Pipe_Ta;
    Pipe_Tr{k}(:,1) = Tri0;
end

% ----------------------------------- Boundary conditions ----------------------------------- %
% Ts,source temperature
% Tr,load temperature
for k = 1:npipes
    % supply source temperature
    Pipe_Ts{k}(1,2:N+1) = 80*ones(1,N)-Pipe_Ta;
    % return loads temperature
    Pipe_Tr{k}(1,2:N+1) = 40*ones(1,N)-Pipe_Ta;
end

% To temperature
% Pipe_To{k}(:,1:N+1) = Pipe_Tr{k}(:,1:N+1);

% node temperature
Node_Ta = -10;
Node_Ts = cell(nnodes,1);
Node_To = cell(nnodes,1);
Node_Tr = cell(nnodes,1);

for k = 1:nnodes
    Node_Ts{k} = 80*ones(N,1)-Node_Ta;
    Node_To{k} = 40*ones(N,1)-Node_Ta;
    Node_Tr{k} = Node_To{k};
end


for t = 1:N
    hydraulic_err = 1; thermal_err = 1;
    while hydraulic_err>1e-3 || thermal_err>1e-3
        % Kf
        for k = 1:npipes
            m(k) = sum(Pipe_m{k}(:,t))/M(k); 
        end
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

        % ----------------------------------- Center-implicit para ----------------------------------- %
        % f1 f2 
        df1 = cell(npipes,1); df2 = cell(npipes,1);
        for k = 1:npipes
            df1{k} = zeros(M(k),1);
            df2{k} = zeros(M(k),1);
        end
        for k = 1:npipes
            for i = 1:M(k)
                df1{k}(i) = (Pipe_Hs{k}(i,t+1)-Pipe_Hs{k}(i,t)+Pipe_Hs{k}(i+1,t+1)-Pipe_Hs{k}(i+1,t))/tau + ...
                            (Pipe_m{k}(i+1,t)-Pipe_m{k}(i,t)+Pipe_m{k}(i+1,t+1)-Pipe_m{k}(i,t+1))/h*a*a/rho/g/s(k);
                df2{k}(i) = (Pipe_Hs{k}(i+1,t)-Pipe_Hs{k}(i,t)+Pipe_Hs{k}(i+1,t+1)-Pipe_Hs{k}(i,t+1))*g/h + ...
                            (Pipe_m{k}(i+1,t+1)-Pipe_m{k}(i+1,t)+Pipe_m{k}(i,t+1)-Pipe_m{k}(i,t))/tau/rho/s(k) + ...
                            (Pipe_m{k}(i,t)+Pipe_m{k}(i,t+1)+Pipe_m{k}(i+1,t)+Pipe_m{k}(i+1,t+1))^2*Kf(k)/16/D(k)/rho/rho/s(k)/s(k);  
            end
        end      
        f1 = vertcat(df1{:});
        f2 = vertcat(df2{:});

        % Jf1H Jf1m Jf2H Jf2m
        Pipe_Jf1H = cell(npipes,1); Pipe_Jf2H = cell(npipes,1);
        Pipe_Jf1m = cell(npipes,1); Pipe_Jf2m = cell(npipes,1);
        for k = 1:npipes
            for i = 1:M(k)
                Pipe_Jf1H{k}(i,i) = 1/tau;
                Pipe_Jf1H{k}(i,i+1) = 1/tau; 
                Pipe_Jf2H{k}(i,i) = -g/h;
                Pipe_Jf2H{k}(i,i+1) = g/h;
            end
        end 
        Jf1H = blkdiag(Pipe_Jf1H{:});
        Jf2H = blkdiag(Pipe_Jf2H{:});

        for k = 1:npipes
            for i = 1:M(k)
                Pipe_Jf1m{k}(i,i) = -a*a/rho/g/s(k)/h;
                Pipe_Jf1m{k}(i,i+1) = a*a/rho/g/s(k)/h;
                % Pipe_Jf2m{k}(i,i) = 1/rho/s(k)/tau+Kf(k)/D(k)/rho/rho/s(k)/s(k)/2*Pipe_m{k}(i,t+1);
                % Pipe_Jf2m{k}(i,i+1) = 1/rho/s(k)/tau+Kf(k)/D(k)/rho/rho/s(k)/s(k)/2*Pipe_m{k}(i+1,t+1);
                Pipe_Jf2m{k}(i,i) = 1/rho/s(k)/tau+Kf(k)/D(k)/rho/rho/s(k)/s(k)/8*(Pipe_m{k}(i,t)+Pipe_m{k}(i,t+1)+Pipe_m{k}(i+1,t)+Pipe_m{k}(i+1,t+1));
                Pipe_Jf2m{k}(i,i+1) = 1/rho/s(k)/tau+Kf(k)/D(k)/rho/rho/s(k)/s(k)/8*(Pipe_m{k}(i,t)+Pipe_m{k}(i,t+1)+Pipe_m{k}(i+1,t)+Pipe_m{k}(i+1,t+1));
            end
        end
        Jf1m = blkdiag(Pipe_Jf1m{:});
        Jf2m = blkdiag(Pipe_Jf2m{:});

        % f3
        % H continuation equation
        Pipe_df31 = cell(npipes,1); Pipe_df32 = cell(npipes,2);
        for k = 1:npipes
            Pipe_df31{k} = zeros(1,M(k)+1);
            Pipe_df32{k} = zeros(1,M(k)+1);
        end
        for k = 1:npipes
            pipe2node = find(Ah(:,k)==1); % Pipe_end - Node = 0
            node2pipe = find(Ah(:,k)==-1); % Pipe_start - Node = 0
            Pipe_df31{k} = Pipe_Hs{k}(end,t+1)-Node_Hs{pipe2node}(t+1);
            Pipe_df32{k} = Pipe_Hs{k}(end,t+1)-Node_Hs{node2pipe}(t+1);
        end
        f31 = vertcat(Pipe_df31{:}); f32 = vertcat(Pipe_df32{:});
        f3 = [f31;f32];

        % Jf3H Jf3Hn
        Pipe_Jf31H = cell(npipes,1); Pipe_Jf32H = cell(npipes,1);
        Jf31Hn = zeros(npipes,nnodes); Jf32Hn = zeros(npipes,nnodes);
        for k = 1:npipes
            Pipe_Jf31H{k} = zeros(1,M(k)+1);
            Pipe_Jf32H{k} = zeros(1,M(k)+1);
        end
        for k = 1:npipes
            pipe2node = find(Ah(:,k)==1); % Pipe_end - Node = 0
            node2pipe = find(Ah(:,k)==-1); % Pipe_start - Node = 0
            Pipe_Jf31H{k}(end) = 1;
            Pipe_Jf32H{k}(1) = 1;
            Jf31Hn(k,pipe2node) = -1;
            Jf32Hn(k,node2pipe) = -1;
        end
        Jf31H = blkdiag(Pipe_Jf31H{:});
        Jf32H = blkdiag(Pipe_Jf32H{:});
        Jf3H = [Jf31H;Jf32H];
        Jf3Hn = [Jf31Hn;Jf32Hn];
        
        % f4
        % node conservation of mass
        df4 = zeros(nnodes,1); 
        for k = 1:nnodes
            m_branch2node = 0; m_node2branch = 0;
            branch2node = find(Ah(k,:)==1);
            node2branch = find(Ah(k,:)==-1);
            for i = 1:length(branch2node)
                m_branch2node = m_branch2node + Pipe_m{branch2node(i)}(end,t+1);
            end
            for i = 1:length(node2branch)
                m_node2branch = m_node2branch + Pipe_m{node2branch(i)}(1,t+1);
            end
            df4(k) = m_branch2node - m_node2branch - Node_m{k}(t+1);
        end
        f4 = df4;

        % Jf4m Jf4mn
        Pipe_Jf4m = cell(1,npipes); Jf4mn = zeros(nnodes,nnodes);
        for k = 1:npipes
            Pipe_Jf4m{k} = zeros(nnodes,M(k)+1);
        end
        for k = 1:nnodes
            branch2node = find(Ah(k,:)==1);
            node2branch = find(Ah(k,:)==-1);
            for i = 1:length(branch2node)
                Pipe_Jf4m{branch2node(i)}(k,end) = 1;
            end
            for i = 1:length(node2branch)
                Pipe_Jf4m{node2branch(i)}(k,1) = -1;
            end
            Jf4mn(k,k) = -1;
        end 
        Jf4m = horzcat(Pipe_Jf4m{:});

        
        % f5
        % m for load node
        df5 = zeros(nnodes,1);
        for k = 1:nloads
            df5(k) =  Node_m{k}(t+1) - Node_m_load(k);
        end
        % H for source node
        for k = (nloads+1):nnodes
            df5(k) = Node_Hs{k}(t+1) - 1*1e6/rho/g; % supply H for source node --highest(MPa) P = rho*g*h(m)
        end
        f5 = df5;

        % Jf5Hn Jf5mn
        Jf5mn = zeros(nnodes,nnodes); Jf5Hn = zeros(nnodes,nnodes);
        for k = 1:nloads
            Jf5mn(k,k) = 1; 
        end
        for k = (nloads+1):nnodes
            Jf5Hn(k,k) = 1; 
        end 

        d_hydraulic = [f1;f2;f3;f4;f5];
        J_hydraulic = [Jf1H Jf1m zeros(N_piece,nnodes) zeros(N_piece,nnodes);
                        Jf2H Jf2m zeros(N_piece,nnodes) zeros(N_piece,nnodes);
                        Jf3H zeros(2*nnodes,N_piece+npipes) Jf3Hn zeros(2*nnodes,nnodes);
                        zeros(nnodes,N_piece+npipes) Jf4m zeros(nnodes,nnodes) Jf4mn;
                        zeros(nnodes,N_piece+npipes) zeros(nnodes,N_piece+npipes) Jf5Hn Jf5mn];
        delta_X = -J_hydraulic\d_hydraulic;
        hydraulic_err = max(abs(d_hydraulic));

        delta_Pipe_Hs = delta_X(1:N_piece+npipes);
        delta_Pipe_m = delta_X(N_piece+npipes+1:2*N_piece+2*npipes);
        delta_Node_Hs = delta_X(2*N_piece+2*npipes+1:2*N_piece+2*npipes+nnodes);
        delta_Node_m = delta_X(2*N_piece+2*npipes+nnodes+1:2*N_piece+2*npipes+2*nnodes);
        
        for k = 1:npipes
            for i = 1:(M(k)+1)
                delta_id = 0 ;
                for i1 = 1:(k-1)
                    delta_id = delta_id + (M(i1)+1);
                end
                delta_id = delta_id + i;
                Pipe_Hs{k}(i,t+1) = Pipe_Hs{k}(i,t+1) + delta_Pipe_Hs(delta_id);
                Pipe_m{k}(i,t+1) = Pipe_m{k}(i,t+1) + delta_Pipe_m(delta_id);
            end
        end
        for k = 1:nnodes
            Node_Hs{k}(t+1) = Node_Hs{k}(t+1) + delta_Node_Hs(k);
            Node_m{k}(t+1) = Node_m{k}(t+1) + delta_Node_m(k);
        end
    end
    % ----------------------------------- FDM para ----------------------------------- %
    for k = 1:npipes
        m(k) = sum(Pipe_m{k}(:,t))/M(k); 
    end
    v = m./s/rho;
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

    % Bs\Br for each pipe  bs = Ak*Bs br = Ak*Br
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
    for k = 1:npipes
        for x = 1:M(k)
            Bs{k}(x) = omega1(k)*Pipe_Ts{k}(x,t) + omega3(k)*Pipe_Ts{k}(x+1,t);
            Br{k}(x) = omega1(k)*Pipe_Tr{k}(x,t) + omega3(k)*Pipe_Tr{k}(x+1,t);
        end
        b_supply{k} = Ak{k}*Bs{k};
        b_return{k} = Ak{k}*Br{k};
    end
    % network topology
    nd = Net_Topo(npipes,nnodes,Ah);
    [newTs,newTo,newTr] = FULL_DYNAMIC_TSOL(W,b_return,b_supply,Node_Ts,Node_To,Pipe_m,Node_m,Ah,nd,nloads);
    thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
    Tr=newTr;
    Ts=newTs;
    To=newTo;
end

% nodeT = PipeT
for k = 1:npipes
    Pipe_Ts{nd(k).k}(1,t+1) = Ts(nd(k).j);
    Pipe_Tr{nd(k).k}(1,t+1) = Tr(nd(k).i);
end

% all T for t+1 
for k = 1:npipes
    Aki = cell(npipes,M(k));
    Bksi = cell(npipes,M(k));
    Bkri = cell(npipes,M(k));
end

for k = 1:npipes
    for x = 2:M(k)+1
        Aki{k,x} = zeros(1,x-1);
        Bksi{k,x} = zeros(x-1,1);
        Bkri{k,x} = zeros(x-1,1);
        for i = 1:x-1
            Aki{k,x}(i) = omega2(k)^(i-1);
            Bksi{k,x}(i) = omega1(k)*Pipe_Ts{k}(i,t) + omega3(k)*Pipe_Ts{k}(i+1,t);
            Bkri{k,x}(i) = omega1(k)*Pipe_Tr{k}(i,t) + omega3(k)*Pipe_Tr{k}(i+1,t);
        end
    end
end

for k = 1:npipes
    for x = 2:M(k)+1
        Pipe_Ts{k}(x,t+1) = omega2(k)^(x-1)*Pipe_Ts{k}(1,t+1) + Aki{k,x}*Bksi{k,x};
        Pipe_Tr{k}(x,t+1) = omega2(k)^(x-1)*Pipe_Tr{k}(1,t+1) + Aki{k,x}*Bkri{k,x};
    end
end


