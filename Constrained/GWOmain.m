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
%------------------------------------------------------------------------------------------------------------------------

n=50;
for run=1:n	%%%%% Start main FOR loop

    rng(run);	

    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    %Initialize the positions of search agents
    Positions=initializationGWO(nPop,ub,lb);

    Convergence_curve=zeros(1,iterMax);

    %%%%% End of GWO initialization (ultimate npop)

    %%%%% Start GWO algorithm
    %tic	% Begin timing
    [Best_score,Best_pos,GWO_cg_curve]=GWO(nPop,iterMax,lb,ub,nVar,fhandle,fnonlin, Positions);
    %BestTime(run)=toc;% End timing
    BestConvergenceCurve(run,:)=GWO_cg_curve;
    BestSolution(run,:)=Best_pos; %Storing Best Solution in every run
    BestFunVal(run)=Best_score;               %Storing Best Function Values (from the Best Solution) in every run

end
%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:)
%%%%% Plot GWO convergence characteristic
ConvergenceCurveOfTheOverallBestFunction=BestConvergenceCurve(index,:);
plot(ConvergenceCurveOfTheOverallBestFunction,'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('Convergence characteristic of GWO')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)