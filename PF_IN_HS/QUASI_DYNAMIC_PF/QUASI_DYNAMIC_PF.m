% /*
%  * @Descripttion: Quasi dynamic power flow in quailty adjustment (only in 2)
%  * @version: 2.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 23:32:08
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-02 23:34:21
%  */

function [m,m_node,Ts,Tr] = QUASI_DYNAMIC_PF(model,adjustmethod)

   %% -----------------------------------------Data init--------------------------------------------------------------------------------
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
           % m = pipe(:,8);
       case 2 % quality adjustment
           % m = filelocation (m result of volume adjustment)
           load model35_result.mat m;
           % load model51_result.mat m;
           % load model225_result.mat m;
   end
   
   %% node temperature
   Ta = 0;
   % assume load supply temperature 
   Ts = node(:,3)-Ta; 
   % assume source return temperature
   To = node(:,4)-Ta;
   Tr = To; 
   % --------------------------------------------------------------------------------------------------------------------------------------------

    % time
    delta_t = 0.1; % dt = 12min
    
    % node method
    gama = rho.*s.*len./m;
    phi = gama + 1;
    Jk = exp(-delta_t./(s.*Cp.*rho.*lamda).*(gama+0.5));
    
    %% network topology
    nd = Net_Topo(npipes,nnodes,Ah);

    % --------------------------------------------Quasi dynamic power flow-------------------------------------------------------------------
    thermal_err = 1;
    while thermal_err > 1e-3
        % thermal model
        mq = Ah*m;
        [newTs,newTo,newTr] = QUASI_DYNAMIC_TSOL(Jk,Ah,m,Ts,To,nd,nloads);
        newTs(nloads+1:nnodes-1) = Phi(nloads+1:nnodes-1)./mq(nloads+1:nnodes-1).*Pbase./Cp + Tr(nloads+1:nnodes-1);
        thermal_err = max([abs(newTs-Ts);abs(newTr-Tr)]);
        Tr=newTr;
        Ts=newTs;
        To=newTo;
    end
    % --------------------------------------------------------------------------------------------------------------------------------------

    Ts=Ts+Ta;
    Tr=Tr+Ta;
    To=To+Ta;

    disp('supply temperature');fprintf('%f\n',Ts);
    disp('return temperature');fprintf('%f\n',Tr);
    Vis_plot(m,m_node,Ts,Tr,0)