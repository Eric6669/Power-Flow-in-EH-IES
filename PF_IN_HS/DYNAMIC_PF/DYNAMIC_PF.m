% /*
%  * @Descripttion: model in DATA , adjustment 1 or 2
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 23:33:20
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-03 09:09:50
%  */

function [m,m_node,Ts,Tr] = DYNAMIC_PF(model,adjustmethod)
     
    % -----------------------------------------------------Data init----------------------------------------------------------------
    [pipe,node,Cp,rho,Pbase,g,viscosity] = deal(model.pipe,model.node,model.water_c,model.water_dens,model.Pbase,model.g,model.viscosity);
    nnodes = length(node(:,1));
    npipes = length(pipe(:,1));
    nloads = 0;
    for i = 1:nnodes
        if node(i,5) == 1
            nloads = nloads + 1;
        end
    end

    D = pipe(:,5);
    % lamda = 0.35*ones(npipes,1);
    lamda = pipe(:,7); 
    R = 1./lamda; % heat resistance
    rough = pipe(:,6);
    s = pi*D.*D/4; % area of cross-section of pipe   
    % node heat power (p.u.) base value 1e6
    Phi = node(1:nnodes-1,2);
    m_node = node(:,6); 


    %% network topology
    % network incidence matrix Ah
    A = [];
    for i = 1:npipes
        A = [A; pipe(i,2) pipe(i,1) -1];
    end
    for i = 1:npipes
        A = [A; pipe(i,3) pipe(i,1) 1];
    end
    A=A';  
    Ah = sparse(A(1,:),A(2,:),A(3,:),nnodes,npipes);

    %% isAnnular or isDendritic
    % type judgement plot digraph
    fnode = pipe(:,2);
    tnode = pipe(:,3);
    pipe_number = pipe(:,1);
    DIG = digraph(fnode, tnode, pipe_number);
    plot(DIG, 'EdgeLabel', DIG.Edges.Weight, 'linewidth', 2, 'EdgeColor','r');

    Aa = full(Ah);
    index_refNode = find(node(:,5) == 0); % reference node
    Aa(index_refNode,:) = []; % basic node-branch incidence matrix

    nloops = 2 - nnodes + npipes -1; % number of basic loops; Euler's Planar Formula

    if nloops > 1e-06 
        isAnnular = 1; % Annular network
        isDendritic = 0;
        Arref = rref(Aa); % Reduced row echelon form
        At = Arref; % block matrix A ~ [At Al]
        Al = Arref;
        index_t = find(sum(abs(Arref)) == 1); % the column index of At and Bt
        index_l = setdiff(1:npipes,index_t); % the column index of Al and Al_b_l
        At(:,index_l) = []; % tree matrix
        Al(:,index_t) = []; % lian zhi matrix
        Bt = (-At\Al)'; % Bt' = -At\Al
        % change the column order
        Bh = zeros(npipes-rank(Aa),npipes);
        Bh(:,index_t) = Bt;
        Bh(:,index_l) = eye(numel(index_l));
    else 
        isAnnular = 0;
        isDendritic = 1; % Dendritic network
        Bh = 0;
    end
    % test = 0?
    % verify = Aa*Bh'; 

    switch adjustmethod 
        case 1 % volume adjustment
            m = pipe(:,8);
            % len = 100*floor(pipe(:,4));
            % len = 10*pipe(:,4);
            len = pipe(:,4);

        case 2 % quality adjustment
            % m = filelocation (m result of case 1)
            % load model35_result.mat m;len = 100*floor(pipe(:,4));
            % load model51_result.mat m;len = 10*pipe(:,4);
            load model225_result.mat m;len = pipe(:,4);
    end

    %% node temperature
    Ta = 0;
    % assume load supply temperature is 80
    Ts = node(:,3)-Ta;
    % assume source return temperature is 40
    To = node(:,4)-Ta;
    Tr = To; 

    % -----------------------------------------------FDM init-------------------------------------------------------------------------
    % h = 50;
    h = 100; % dx = 100m
    M = len./h; % 1<=x<=M

    tau = 60*6; % dt = 6min
    % tau = 60; % dt = 1min
    tsim = 3600*3; % total simulation time 3h
    N = tsim/tau; % 1<=t<=N  N=200

    v = m./s./rho;
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
    A = cell(npipes,1);
    for k = 1:npipes
        A{k} = zeros(1,M(k));
        for x = 1:M(k)
            A{k}(x) = omega2(k)^(M(k)-x);
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

    % 初始条件（管道初始温度分布）-all pipes 
    % differ in different system 
    for k = 1:npipes
        % Tsi0 = 78.75*ones(M(k)+1,1)-Ta;
        % Tsi0 = 80*ones(M(k)+1,1)-Ta;
        Tsi0 = 84.05*ones(M(k)+1,1)-Ta;
        Pipe_Ts{k}(:,1) = Tsi0;
        Tri0 = 40*ones(M(k)+1,1)-Ta;
        Pipe_Tr{k}(:,1) = Tri0;
    end

    % 边值条件（入口管段热媒温度随时间变化情况）
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

    
   % --------------------------------------------Dynamic power flow-------------------------------------------------------------------

    % 考虑时间断面 
    for t = 1:N
        for k = 1:npipes
            for x = 1:M(k)
                Bs{k}(x) = omega1(k)*Pipe_Ts{k}(x,t) + omega3(k)*Pipe_Ts{k}(x+1,t);
                Br{k}(x) = omega1(k)*Pipe_Tr{k}(x,t) + omega3(k)*Pipe_Tr{k}(x+1,t);
            end
            b_supply{k} = A{k}*Bs{k};
            b_return{k} = A{k}*Br{k};
        end

        % 每个时间断面求解
        switch adjustmethod
            % ----------------------------------------- volume adjustment -----------------------------------------------------------------
            case 1 
                thermal_err = 1;
                hydraulic_err = 1;
                % -----------------------------------------Annular network-----------------------------------------
                if isAnnular == 1    
                    while hydraulic_err>1e-3||thermal_err>1e-3
                        %% hydraulic model
                        % update dT
                        dT = [Ts(1:nloads)-To(1:nloads);Ts(nloads+1:nnodes-1)-Tr(nloads+1:nnodes-1)];      
                        mfeedback = sign(m);
                        % update Ah,Bh,m
                        Ah = Ah.*repmat(mfeedback',nnodes,1);
                        Bh = Bh.*mfeedback';
                        m = abs(m);
                        % update m_node
                        m_node = Ah*m;
                        
                        % nodal flow  mismatches 
                        dPhi = Cp/Pbase*m_node(1:nnodes-1).*dT-Phi;
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
                        % d_hydraulic = [dm;dhf];
                        d_hydraulic = [dPhi;dhf];
                        J_hydraulic = [JPhi;Jhf];
                        delta_m = -J_hydraulic\d_hydraulic;
                        hydraulic_err = max(abs(d_hydraulic));
                        m = m + delta_m;
                        
                        %% thermal model
                        % network topology
                        nd = Net_Topo(npipes,nnodes,Ah);
                        [newTs,newTo,newTr] = DYNAMIC_TSOL(W,b_return,b_supply,Ts,To,m,Ah,nd,nloads);
                        thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                        Tr=newTr;
                        Ts=newTs;
                        To=newTo;
                    end
                % -----------------------------------------Dendritic network-----------------------------------------
                elseif isDendritic ==1
                    while thermal_err>1e-3
                        %% hydraulic model
                        % update dT
                        dT = [Ts(1:nloads)-To(1:nloads);Ts(nloads+1:nnodes-1)-Tr(nloads+1:nnodes-1)];     
                        % heat power equation->m_node
                        m_node(1:nnodes-1) = Phi./Cp./dT.*Pbase;
                        m = linsolve(full(Ah(1:nnodes-1,:)),m_node(1:nnodes-1));
                        
                        %% thermal model
                        nd = Net_Topo(npipes,nnodes,Ah);
                        [newTs,newTo,newTr] = DYNAMIC_TSOL(W,b_return,b_supply,Ts,To,m,Ah,nd,nloads);
                        thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                        Tr=newTr;
                        Ts=newTs;
                        To=newTo;
                    end
                end
            % ----------------------------------------- quality adjustment -----------------------------------------------------------------
            case 2
                thermal_err = 1;
                while thermal_err > 1e-3
                    % thermal model
                    mq = Ah*m;
                    nd = Net_Topo(npipes,nnodes,Ah);
                    [newTs,newTo,newTr] = DYNAMIC_TSOL(W,b_return,b_supply,Ts,To,m,Ah,nd,nloads);
                    newTs(nloads+1:nnodes-1) = Phi(nloads+1:nnodes-1)./mq(nloads+1:nnodes-1).*Pbase./Cp + Tr(nloads+1:nnodes-1);
                    thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                    Tr=newTr;
                    Ts=newTs;
                    To=newTo;
                end
        end

        %% nodeT = PipeT
        for k = 1:npipes
            Pipe_Ts{nd(k).k}(1,t+1) = Ts(nd(k).j);
            Pipe_Tr{nd(k).k}(1,t+1) = Tr(nd(k).i);
        end

        %% all T for t+1 
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
    end

    figure(1);
    for k = 1:npipes
        plot(Pipe_Ts{k}(M(k)+1,:),':','LineWidth',1.8); hold on;
    end
    figure(2);
    for k = 1:npipes
        plot(Pipe_Tr{k}(M(k)+1,:),'-','LineWidth',1.8); hold on;
    end

end


