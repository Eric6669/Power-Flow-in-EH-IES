% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-03 15:37:38
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-03 17:28:07
%  */
function [m,m_node,Ts,Tr,Herr,Phi_chp2_hslack] = STEADY_PF(model,adjustmethod,Phi_chp1_eslack0,Phi_chp3)

    %% -----------------------------------Data init-------------------------------------------------------------------------------------
    [pipe,node,Cp,rho,Pbase,g,viscosity] = deal(model.pipe,model.node,model.water_c,model.water_dens,model.Pbase,model.g,model.viscosity);
    nnodes = length(node(:,1));
    npipes = length(pipe(:,1));
    nloads = 0;
    for i = 1:nnodes
        if node(i,9) == 1
            nloads = nloads + 1;
        end
    end

    D = pipe(:,4);
    lamda = pipe(:,6); 
    % R = 1./lamda; % heat resistance
    rough = pipe(:,5);
    s = pi*D.*D/4; % area of cross-section of pipe   
    len = pipe(:,3);
    
    % node heat power (p.u.) base value 1e6
    % load96 = deal(model.load);
    % load_chosen = load96(48,:)';
    % Phi = load_chosen(1:nloads,1);
    Phi = [10.0294448	9.957627727	40.273569375	13.10859375	11.393594125	10.923405292	2.343467	5.127948	10.021023333	13.09284013	20.405475779	6.8384421	3.187027754	44.06466183	0.094882667	2.3901983	2.06179995	2.482797158	3.764163375	0.695273367	5.7730113	2.376980375	2.930970992	5.905707488	7.852104925	8.6056397	7.272569967	1.201022375	0.7732263	5.063259234	50.8837266	6.75106659	28.33800991	5.64171148	8.965308196	7.086144086	2.059494	25.58327568	7.388423982	3.10164321	4.07140726	13.15310534	6.009421514	8.74262627	10.9273357	8.317739864	8.028660155	7.694059215	12.673768	11.01184383	5.346804418	8.556597763	27.45703723	9.469866555	2.347831324	23.47174068	3.731597295	4.793815673	2.123330851	5.818911986	2.782336939	28.27255014	29.89465381	6.255552299	11.91517176	13.20652273]';
    
    % m_node = node(:,6); 
    m_node = ones(nnodes,1);

    %% network topology
    % network incidence matrix Ah
    A = [];
    for i = 1:npipes
        A = [A; pipe(i,1) i -1];
    end
    for i = 1:npipes
        A = [A; pipe(i,2) i 1];
    end
    A=A';  
    Ah = sparse(A(1,:),A(2,:),A(3,:),nnodes,npipes);

    %% isAnnular or isDendritic

    % type judgement plot digraph
    % fnode = pipe(:,1);
    % tnode = pipe(:,2);
    % pipe_number = zeros(npipes,1);
    % for k = 1:npipes
    %     pipe_number(k)=k;
    % end
    % DIG = digraph(fnode, tnode, pipe_number);
    % plot(DIG, 'EdgeLabel', DIG.Edges.Weight, 'linewidth', 2, 'EdgeColor','r');

    Aa = full(Ah);
    index_refNode = find(node(:,9) == 0); % reference node
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
            % m = pipe(:,8);
            m = 50*ones(npipes,1);
        case 2 % quality adjustment
            % m = filelocation (m result of case 1)
            % load model35_result.mat m;
            % load model51_result.mat m;
            % load model225_result.mat m;
    end
    
    %% node temperature
    Ta = -10*ones(nnodes,1);
    % Ta = ones(nnodes,1);
    % assume load supply temperature 
    Ts = 90*ones(nnodes,1)-Ta; 
    % assume source return temperature
    To = 70*ones(nnodes,1)-Ta;
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
                    dPhi = Cp/Pbase*m_node(1:nnodes-1).*dT-[Phi;-Phi_chp1_eslack0;-Phi_chp3];
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
    m_node = Ah*m;
    Phi_chp2_hslack = Cp/Pbase*(-m_node(nnodes))*(Ts(nnodes)-Tr(nnodes));
end

