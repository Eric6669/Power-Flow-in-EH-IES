% /*
%  * @Descripttion: 
%  * @version: 
%  * @Author: Ke Wang
%  * @Date: 2024-07-02 11:55:55
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-02 12:17:31
%  */
function nd = Net_Topo(npipes,nnodes,Ah)
    %% network topology description nd
    % supply j--k-->i  
    % return i--k-->j
    nd = struct('i', {}, 'j', {}, 'k', {},'mix_supply',{},'mix_return',{});
    for k = 1:npipes
        for i = 1:nnodes
            if Ah(i,k) == 1
                for j = 1:nnodes
                    if Ah(j,k) == -1
                        nd(end+1).i = i;
                        nd(end).j = j;
                        nd(end).k = k;
                        break;
                    end
                end
            end
        end
    end
    for k1 = 1:npipes-1
        for k2 = k1+1:npipes
            if nd(k1).i == nd(k2).i
                nd(k1).mix_supply = 1;
                nd(k2).mix_supply = 1;
            end
            if nd(k1).j == nd(k2).j
                nd(k1).mix_return = 1;
                nd(k2).mix_return = 1;
            end
        end
    end

end