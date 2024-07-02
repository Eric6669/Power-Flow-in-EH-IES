% /*
%  * @Descripttion: Jacobi_matrix solution
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 20:04:36
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 20:06:33
%  */
function J = Jacobi_matrix(n,nPV,U,cita,B,G,Pi,Qi,slackbus,PVbus)
        % H N M L
        % H = zeros(n-1,n-1);
        % N = zeros(n-1,nPQ);
        % M = zeros(nPQ,n-1);
        % L = zeros(nPQ,nPQ);
        H = zeros(n,n);
        N = zeros(n,n);
        M = zeros(n,n);
        L = zeros(n,n);
        % H
        for i=1:n
            for j=1:n
                if i~=j
                    H(i,j)=-U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
                else
                    H(i,i)=U(i).^2*B(i,i)+Qi(i);
                end
            end
        end
        H(slackbus,:) = [];
        H(:,slackbus) = [];

        % N
        for i=1:n
            for j=1:n
                if i~=j
                    N(i,j)=-U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
                else
                    N(i,i)=-U(i).^2*G(i,i)-Pi(i);
                end
            end
        end
        N(slackbus,:) = [];
        N(:,slackbus) = -inf;
        for k = 1:nPV
            N(:,PVbus(k)) = -inf;
        end
        N(:,[find(N(1,:) == -inf)]) = [];

        % Minf
        for i=1:n
            for j=1:n
                if i~=j
                    M(i,j)=U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
                else
                    M(i,i)=U(i).^2*G(i,i)-Pi(i);
                end
            end
        end
        M(:,slackbus) = [];
        M(slackbus,:) = -inf;
        for k = 1:nPV
            M(PVbus(k),:) = -inf;
        end
        M([find(M(:,1) == -inf)],:) = [];


        % L
        for i=1:n
            for j=1:n
                if i~=j
                    L(i,j)=-U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
                else
                    L(i,i)=U(i).^2*B(i,i)-Qi(i);
                end
            end
        end
        L(slackbus,:) = -inf;
        for k = 1:nPV
            L(PVbus(k),:) = -inf;
        end
        L([find(L(:,1) == -inf)],:) = [];

        L(:,slackbus) = -inf;
        for k = 1:nPV
            L(:,PVbus(k)) = -inf;
        end
        L(:,[find(L(1,:) == -inf)]) = [];
    
        J=[H N;M L]; 
end
    