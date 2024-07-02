% /*
%  * @Descripttion: power on line solution
%  * @version: 1.0
%  * @Author: Ke Wang
%  * @Date: 2024-07-01 20:05:07
%  * @LastEditors: Ke Wang
%  * @LastEditTime: 2024-07-01 20:11:03
%  */
function Sij = Line_power(n,y,U,cita)
    for i=1:n
        U1(i)=complex(U(i)*cos(cita(i)),U(i)*sin(cita(i)));
        for j=1:n
            U1(j)=complex(U(j)*cos(cita(j)),U(j)*sin(cita(j)));
            Sij(i,j)=U1(i)*conj((U1(i)-U1(j))*y(i,j));
        end
    end
end