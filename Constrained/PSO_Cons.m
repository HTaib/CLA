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
   
function [ConvergenceCurve, bestfun, best_variables]=PSO_Cons(nPop, iterMax, fhandle, fnonlin, nVar, wmax, wmin, c1, c2, LB, UB, x, v, pbest, gbest, f0, fmin0);
ffmin=zeros(iterMax,1);
f=zeros(nPop,1);

ite=1;
tolerance=1;
while ite<=iterMax && tolerance > 10^-6	%%%%% Start WHILE loop

	w=wmax-(wmax-wmin)*ite/iterMax;  %%%%% update inertial weight

%%%%% PSO velocity updates---------------
	for i=1:nPop
		for j=1:nVar
			v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...
			+c2*rand()*(gbest(1,j)-x(i,j));
		end
    end
    
%%%%% PSO position update
	for i=1:nPop
		for j=1:nVar
			x(i,j)=x(i,j)+v(i,j);
		end
    end
    
%%%%% Handling boundary violations
	for i=1:nPop
		for j=1:nVar
			if x(i,j)<LB(j)
				x(i,j)=LB(j);
			elseif x(i,j)>UB(j)
				x(i,j)=UB(j);
			end
		end
    end

%%%%% Evaluating fitness
	for i=1:nPop
% 		f(i,1)=fun(x(i,:));
        f(i,1)=Fun(fhandle,fnonlin,x(i,:));  
	end

%%%%% Updating pbest and fitness
	for i=1:nPop
		if f(i,1)<f0(i,1)
			pbest(i,:)=x(i,:);
			f0(i,1)=f(i,1);
		end
	end
	[fmin,index]=min(f0);		%%%%% Finding the best particle
	ffmin(ite,1)=fmin;		       %%%%% Storing best fitness
	ffite(1)=ite;			            %%%%% Storing iteration count

%%%%% Updating gbest and best fitness
	if fmin<fmin0
		gbest=pbest(index,:);
		fmin0=fmin;
	end

%%%%% Calculating tolerance
	if ite>1000
        tolerance=abs(fmin0-ffmin(ite-300,1));
	end
%%%%% Displaying iterative results
	if ite==1
		fprintf('\nIteration\t Fitness Value\n');
	end
	fprintf('%8g\t %12.4f\n',ite,fmin0);	
		ite=ite+1;
	end		% End WHILE loop
%%%%% End of PSO algorithm

	gbest;
%%%%%---------You must use your objective (fitness) function---------
%   fvalue=0.6224*gbest(1)*gbest(3)*gbest(4)+1.7781*gbest(2)*gbest(3)^2+3.1661*gbest(1)^2*gbest(4)+19.84*gbest(1)^2*gbest(3);
    
fvalue=Fun(fhandle,fnonlin,gbest);
%  f(i,1)=Fun(fun,fnonlin,gbest);  
%%%%%----------------------------------------------------------------
	fff(1)=fvalue;
	rgbest(1,:)=gbest;

%%%%% End of PSO main program---------------

%toc	%%%%% Time out for calculations

[bestfun,bestrun]=min(fff);
best_variables=rgbest(bestrun,:);
ConvergenceCurve=ffmin;
end