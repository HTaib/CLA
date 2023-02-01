%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of research works on Combination Lock Algorithm %
% Paper co-author and programmer: Hasnanizan Taib                %
%                                                                                                %
% Email: hasnanizan.taib@gmail.com                                         %
%           hasnanizan.taib@utb.edu.bn                                        %
%                                                                                                %
% Main paper:                                                                             %
% The novel combination lock algorithm for improving the            %
% performance of metaheuristic optimizers,                                  % 
% by A. Bahreininejad & H. Taib                                                  %
% Advances in Engineering Software                                           %
% DOI: 10.1016/j.advengsoft.2022.103177                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g,geq]=constraint(x)
% Inequality constraints

%==========================================================================
%%%%% Three-bar truss problem
p=2;
sigma=2;
g(1)=((sqrt(2)*x(1)+x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*p-sigma;
g(2)=((x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*p-sigma;
g(3)= (1/(sqrt(2)*x(2)+x(1)))*p-sigma;
%==========================================================================
%%%%% Speed reducer design problem
% g(1)=27/(x(1)*x(2)^2*x(3))-1;
% g(2)=397.5/(x(1)*x(2)^2*x(3)^2)-1;
% g(3)=(1.93*x(4)^3)/(x(2)*x(3)*x(6)^4)-1;
% g(4)=(1.93*x(5)^3)/(x(2)*x(3)*x(7)^4)-1;
% g(5)=((sqrt(((745*x(4))/(x(2)*x(3)))^2+16.9e6))/(110*x(6)^3))-1;
% g(6)= ((sqrt(((745*x(5))/(x(2)*x(3)))^2+157.5e6))/(85*x(7)^3))-1;
% g(7)=((x(2)*x(3))/40)-1;
% g(8)=(5*x(2)/x(1))-1;
% g(9)= (x(1)/12*x(2))-1;
% g(10)=((1.5*x(6)+1.9)/x(4))-1;
% g(11)=((1.1*x(7)+1.9)/x(5))-1;
%==========================================================================
% %%%% Welded beam problem
% R=sqrt((x(2)^2)/4+((x(1)+x(3))/2)^2);
% P=6000;
% L=14;
% E=30e6;
% G=12e6;
% Px=((4.013*E*sqrt((x(3)^2*x(4)^6)/36))/(L^2))*(1-(x(3)/(2*L))*sqrt(E/(4*G)));
% tap=P/(sqrt(2)*x(1)*x(2));
% Q=P*(L+(x(2)/2));
% J=2*(sqrt(2)*x(1)*x(2)*((x(2)^2/12)+((x(1)+x(3))/2)^2));
% tapp=(Q*R)/J;
% ta=sqrt(tap^2+((2*tap*tapp*x(2))/(2*R))+tapp^2);
% sigma=(6*P*L)/(x(4)*(x(3)^2));
% delta=(4*P*(L^3))/(E*(x(3)^3)*x(4));
% g(1)=ta-13600;
% g(2)=sigma-30000;
% g(3)=x(1)-x(4);
% g(4)=0.10471*x(1)^2+0.04811*x(3)*x(4)*(14+x(2))-5;
% g(5)=0.125-x(1);
% g(6)=delta-0.25;
% g(7)=P-Px;
%==========================================================================

% If no equality constraint at all, put geq=[] as follows
geq=[];