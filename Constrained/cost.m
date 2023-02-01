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

function f=cost(x)

%==========================================================================
% Three-bar truss problem
f=(2*sqrt(2)*x(1)+x(2))*100;
%==========================================================================
% Speed reducer problem  
% f=0.7854*x(1)*x(2)^2*(3.3333*x(3)^2+14.9334*x(3)-43.0934)-1.508*x(1)*(x(6)^2+x(7)^2)+7.4777*(x(6)^3+x(7)^3)+0.7854*(x(4)*x(6)^2+x(5)*x(7)^2);
%==========================================================================
% Welded Beam
% f=1.10471*x(1)^2*x(2)+0.04811*x(3)*x(4)*(14+x(2));
%==========================================================================
end