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

clc
clear all
close all

rng('default');  

%%%%% List of Benchmark Engineering Problems %%%%%
%==========================================================================
% Three-bar truss problem
nVar=2;
lb=[0 0];
ub=[1 1];
%==========================================================================
% Speed reducer problem
% nVar=7;
% lb=[2.6 0.7 17 7.3 7.3 2.9 5];
% ub=[3.6 0.8 28 8.3 8.3 3.9 5.5];
%==========================================================================
% Welded beam problem
% nVar=4;
% lb=[0.1 0.1 0.1 0.1];
% ub=[2 10 10 2];
%==========================================================================

%%%%% Experimental Settings
nPop=50; % Number of search agents for PSO/ACOR/GWO                    - CHANGE HERE
iterMax=1000; % Maximum number of iterations for PSO/ACOR/GWO      - CHANGE HERE 
fhandle=@cost;
fnonlin=@constraint;
%%%%% Metaheuristics' Parameters Settings
%PSO's Values of Parameters
wmax=0.9;	% maximum inertia weight
wmin=0.4;	 % minimum inertia weight
c1=2;		    % personal acceleration factor
c2=2;		    % global acceleration factor
%%%%%----------------------------------------------------------------

%%%%% Start PSO main program
n=50; % No. of Simulation Runs
for run=1:n
    rng(run);  
    fprintf('\n\nRun  ------>');
	fprintf('\t%8g\n',run);
    
    %%%%% Start PSO initialization
    x0=zeros(nPop,nVar);     %%%%% Initialize x0 vector
    f0=zeros(nPop,1);

    x0=initializationPSO(nPop,ub,lb); %%%%% Random generation of initial starting individuals

    x=x0;		%%%%% initial population
    v=0.1*x0;	%%%%% initial velocity

    for i=1:nPop
        f0(i,1)=Fun(fhandle,fnonlin,x0(i,:)); 
    end

    [fmin0,index0]=min(f0);

    pbest=x0;	%%%%% initial pbest
    gbest=x0(index0,:); %%%%% initial gbest

    %%%%% End of PSO initialization

    %%%%% Start PSO algorithm
    %tic	% Begin timing
    [ConvergenceCurve, bestfun, best_variables]=PSO_Cons(nPop, iterMax, fhandle, fnonlin, nVar, wmax, wmin, c1, c2, lb, ub, x, v, pbest, gbest, f0, fmin0);
    %BestTime(run)=toc;% End timing
    BestConvergenceCurve(run,:)=ConvergenceCurve;
    BestSolution(run,:)=best_variables; %Storing Best Solution in every run
    BestFunVal(run)=bestfun;              %Storing Best Function Values (from the Best Solution) in every run
end
%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:)
%%%%% Plot PSO convergence characteristic
ConvergenceCurveOfTheOverallBestFunction=BestConvergenceCurve(index,:);
plot(ConvergenceCurveOfTheOverallBestFunction,'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('Convergence characteristic of PSO')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)