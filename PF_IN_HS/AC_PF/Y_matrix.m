% /*
%  * @Descripttion: node admittance matrix
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 20:47:52
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 20:47:59
%  */

function [Y] = Y_matrix(nodenumber,branch)
    Y = zeros(nodenumber);
    for k=1:length(branch(:,1))
        start=branch(k,1);
        final=branch(k,2);
        r=branch(k,3);
        x=branch(k,4);
        b=branch(k,5)/2;
        % K=branch(k,6);  % Remove this line
        Y(start,start)=Y(start,start)+1/(r+1i*x)+1i*b;
        Y(start,final)=Y(start,final)-1/(r+1i*x);
        Y(final,start)=Y(final,start)-1/(r+1i*x);
        Y(final,final)=Y(final,final)+1/(r+1i*x)+1i*b;
    end
end