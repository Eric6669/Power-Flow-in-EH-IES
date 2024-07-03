% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 23:32:26
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-03 09:09:25
%  */
function [newTs,newTo,newTr] = QUASI_DYNAMIC_TSOL(Jk,Ah,m,Ts,To,nd,nloads)
    
    mq = Ah*m;
    npipes = size(Ah,2);
    nnodes = size(Ah,1);
    %% initialize Cs,bs
    Cs = zeros(nloads);
    bs = zeros(nloads,1);
    Cr = zeros(nloads);
    br = zeros(nloads,1);

    %% supplynetwork j——k——>i
    for index = 1:npipes
        if nd(index).i <= nloads
            if nd(index).mix_supply == 1
                % Cs(i,i)
                Cs(nd(index).i,nd(index).i) = Cs(nd(index).i,nd(index).i) + m(nd(index).k);
                % judge j load node?
                if nd(index).j <= nloads
                    % is load node
                    Cs(nd(index).i,nd(index).j) = -m(nd(index).k)*Jk(nd(index).k);
                else
                    bs(nd(index).i) = bs(nd(index).i) + m(nd(index).k)*Ts(nd(index).j)*Jk(nd(index).k);
                end
            else
                Cs(nd(index).i,nd(index).i) = 1;
                bs(nd(index).i) = Ts(nd(index).j)*Jk(nd(index).k);
            end
        end
    end

    newTs = [Cs\bs;Ts(nloads+1:nnodes)];    

    %% returnnetwork i——k——>j
    % load nodes
    runonce = zeros(nnodes,1);
    for index1 = 1:npipes
        for index2 = index1+1:npipes
            if nd(index1).j == nd(index2).j
                runonce(index2) = 1;
            end
        end
    end

    for index = 1:npipes
        if nd(index).j <= nloads
            % mix load nodes
            if nd(index).mix_return == 1
                Cr(nd(index).j,nd(index).i) = -m(nd(index).k)*Jk(nd(index).k);
                br(nd(index).j) = mq(nd(index).j)*To(nd(index).j);
                % 从j流出的所有m
                if runonce(index) == 0
                    for index2 = 1:npipes
                        if nd(index2).i == nd(index).j
                            Cr(nd(index).j,nd(index).j) = Cr(nd(index).j,nd(index).j) + m(nd(index2).k);
                        end
                    end
                end
            else
                Cr(nd(index).j,nd(index).j) = 1;
                br(nd(index).j) = To(nd(index).j);
            end
        end

    end

    for index = 1:npipes
        J(index) = nd(index).j;
    end
    for j =1:nloads
        if ~ismember(j,J)
            Cr(j,j) = 1;
            br(j) = To(j);
        end
    end
    newTr = zeros(1,nnodes);
    newTr(1:nloads) = Cr\br;

    % new Tr,sources
    for index = 1:npipes
        if nd(index).j > nloads
            if nd(index).mix_return == 1
                newTr(nd(index).j) = newTr(nd(index).j) + (m(nd(index).k)*newTr(nd(index).i)*Jk(nd(index).k))/(-mq(nd(index).j));
            else
                newTr(nd(index).j) = newTr(nd(index).i)*Jk(nd(index).k);
            end
        end
    end
  
    newTr = [newTr(1:nloads)';newTr(nloads+1:nnodes)'];
    newTo = [To(1:nloads);newTr(nloads+1:nnodes)];
    
end





