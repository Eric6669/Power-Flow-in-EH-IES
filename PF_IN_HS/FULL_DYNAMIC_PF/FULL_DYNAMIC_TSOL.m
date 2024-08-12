% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 23:33:33
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-03 09:11:02
%  */
function [newTs,newTo,newTr] = FULL_DYNAMIC_TSOL(W,b_return,b_supply,Node_Ts,Node_To,Pipe_m,Node_m,Ah,nd,nloads,t)
    npipes = size(Ah,2);
    nnodes = size(Ah,1);

    Cs = zeros(nloads);
    bs = zeros(nloads,1);
    Cr = zeros(nloads);
    br = zeros(nloads,1);

    %% supply 
    % Ts loads
    for k =1:npipes
        if nd(k).i <= nloads
            if nd(k).mix_supply == 1
                % Cs(i,i)
                Cs(nd(k).i,nd(k).i) = Cs(nd(k).i,nd(k).i) + Pipe_m{nd(k).k}(1,t+1);
                % judge j load node?
                if nd(k).j <= nloads
                    % is loadnode
                    Cs(nd(k).i,nd(k).j) = -Pipe_m{nd(k).k}(end,t+1)*W(nd(k).k);
                    bs(nd(k).i) = bs(nd(k).i) + Pipe_m{nd(k).k}(end,t+1)*b_supply{nd(k).k};
                else
                    bs(nd(k).i) = bs(nd(k).i) + Pipe_m{nd(k).k}(end,t+1)*(W(nd(k).k)*Node_Ts{nd(k).j}(t)+b_supply{nd(k).k});
                end
            else
                Cs(nd(k).i,nd(k).i) = 1;
                bs(nd(k).i) = W(nd(k).k)*Node_Ts{nd(k).j}(t)+b_supply{nd(k).k};
            end
        end
    end

    newTs = Cs\bs;
    for k = nloads+1:nnodes
        newTs = [newTs;Node_Ts{k}(t)];
    end

    %% return
    % Tr loads
    runonce = zeros(npipes,1);
    for k1 = 1:npipes
        for k2 = k1+1:npipes
            if nd(k1).j == nd(k2).j
                runonce(k2) = 1;
            end
        end
    end

    for k = 1:npipes
        if nd(k).j <= nloads
            % mix load nodes
            if nd(k).mix_return == 1
                Cr(nd(k).j,nd(k).i) = -Pipe_m{nd(k).k}(end,t+1)*W(nd(k).k);
                br(nd(k).j) = br(nd(k).j) + Pipe_m{nd(k).k}(end,t+1)*b_return{nd(k).k};
                if runonce(k) == 0
                    br(nd(k).j) = br(nd(k).j) + Node_m{nd(k).j}(t+1)*Node_To{nd(k).j}(t);
                end
                % 从j流出的所有m
                if runonce(k) == 0
                    for k2 = 1:npipes
                        if nd(k2).i == nd(k).j
                            Cr(nd(k).j,nd(k).j) = Cr(nd(k).j,nd(k).j) + Pipe_m{nd(k).k}(1,t+1);
                        end
                    end
                end
            else
                Cr(nd(k).j,nd(k).j) = 1;
                br(nd(k).j) = Node_To{nd(k).j}(t);
            end
        end
    end

    for index = 1:npipes
        J(index) = nd(index).j;
    end
    for j =1:nloads
        if ~ismember(j,J)
            Cr(j,j) = 1;
            br(j) = Node_To{j}(t);
        end
    end
    newTr = zeros(nnodes,1);
    newTr(1:nloads) = Cr\br;

    % Tr sources
    for k = 1:npipes
        if nd(k).j > nloads
            if nd(k).mix_return == 1
                newTr(nd(k).j) = newTr(nd(k).j) + (Pipe_m{nd(k).k}(end,t+1)*(newTr(nd(k).i)*W(nd(k).k)+b_return{nd(k).k}))/(-Node_m{nd(k).j}(t+1));
            else
                newTr(nd(k).j) = W(nd(k).k)*newTr(nd(k).i)+b_return{nd(k).k};
            end
        end
    end
    newTr = [newTr(1:nloads);newTr(nloads+1:nnodes)];
    newTo = newTr(nloads+1:nnodes);
    for k = 1:nloads
        newTo = [Node_To{k}(t);newTo];
    end

end