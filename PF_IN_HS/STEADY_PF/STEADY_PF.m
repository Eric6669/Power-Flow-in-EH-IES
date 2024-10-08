% /*
%  * @Descripttion: Steady power flow in HS with method 1 volume_adjustment and 2 quality_adjustment
%  * @version: 2.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 12:16:34
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-02 12:18:04
%  */

function [m,m_node,Ts,Tr] = STEADY_PF(model,adjustmethod)

    %% -----------------------------------Data init-------------------------------------------------------------------------------------
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
    % R = 1./lamda; % heat resistance
    rough = pipe(:,6);
    s = pi*D.*D/4; % area of cross-section of pipe   
    len = pipe(:,4);
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
            % m = ones(npipes,1);
        case 2 % quality adjustment
            % m = filelocation (m result of case 1)
            % load model35_result.mat m;
            % load model51_result.mat m;
            load model225_result.mat m;
    end
    
    %% node temperature
    Ta = 0;
    % assume load supply temperature 
    Ts = node(:,3)-Ta; 
    % assume source return temperature
    To = node(:,4)-Ta;
    Tr = To; 
    % ----------------------------------------------------------------------------------------------------------------------------------
    
    % --------------------------------------------Steady power flow-------------------------------------------------------------------
    % iteration times
    t = 0;
    Herr = inf*ones(20,1);

    switch adjustmethod
        % ----------------------------------------- volume adjustment -----------------------------------------------------------------
        case 1
            hydraulic_err = 1;
            thermal_err = 1;
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
                    d_hydraulic = [dPhi;dhf];
                    J_hydraulic = [JPhi;Jhf];
                    delta_m = -J_hydraulic\d_hydraulic;
                    hydraulic_err = max(abs(d_hydraulic));
                    m = m + delta_m;
                    
                    %% thermal model
                    % pipe temperature fall coefficient
                    Kt = exp(-lamda.*len./Cp./m);
                    %% network topology description nd
                    nd = Net_Topo(npipes,nnodes,Ah);
                    [newTs,newTo,newTr]=STEADY_TSOL(Kt,Ah,m,Ts,To,nd,nloads);
                    thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                    Tr=newTr;
                    Ts=newTs;
                    To=newTo;
                    t = t+1;
                    Herr(t,:) = max(hydraulic_err,thermal_err);
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
                    % pipe temperature fall coefficient
                    Kt = exp(-lamda.*len./Cp./m);
                    %% network topology description nd
                    nd = Net_Topo(npipes,nnodes,Ah);
                    [newTs,newTo,newTr]=STEADY_TSOL(Kt,Ah,m,Ts,To,nd,nloads);
                    thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                    Tr=newTr;
                    Ts=newTs;
                    To=newTo;
                    t = t+1;
                    Herr(t,:) = thermal_err;
                end
            end
        % ----------------------------------------- quality adjustment -----------------------------------------------------------------
        case 2
            thermal_err = 1;
            mq = Ah*m;
            Kt = exp(-lamda.*len./Cp./m);
            while thermal_err > 1e-3
                % thermal model
                % network topology description nd
                nd = Net_Topo(npipes,nnodes,Ah);
                [newTs,newTo,newTr] = STEADY_TSOL(Kt,Ah,m,Ts,To,nd,nloads);
                newTs(nloads+1:nnodes-1) = Phi(nloads+1:nnodes-1)./mq(nloads+1:nnodes-1).*Pbase./Cp + Tr(nloads+1:nnodes-1);
                thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
                Tr=newTr;
                Ts=newTs;
                To=newTo;
            end
    end
    % --------------------------------------------------------------------------------------------------------------------------------

    Ts=Ts+Ta;
    Tr=Tr+Ta;
    % To=To+Ta;

    disp('mass flow rates within each pipe');disp(m);
    disp('flow injection at the node');disp(m_node);
    disp('supply temperature');fprintf('%f\n',Ts);
    disp('return temperature');fprintf('%f\n',Tr);  
    Vis_plot(m,m_node,Ts,Tr,Herr);
    
    Phi_hslack = Cp/Pbase*(-m_node(nnodes))*(Ts(nnodes)-Tr(nnodes));
    disp('heat power of hslack Node (MWth)');fprintf('%f\n',Phi_hslack);
    
    % volume adjustment, save m result for quaility adjustment
    % save('./RESULT/model35_result.mat','m','m_node','Ts','Tr','Herr');
    % save('./RESULT/model51_result.mat','m','m_node','Ts','Tr','Herr');
    % save('./RESULT/model225_result.mat','m','m_node','Ts','Tr','Herr');

end

