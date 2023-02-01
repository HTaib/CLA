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
VarSize = [1 nVar];   % Variables Matrix Size
%%%%% Metaheuristics' Parameters Settings
%ACOR's Values of Parameters
nSample = 40;         % Sample Size
q = 0.5;              % Intensification Factor (Selection Pressure)
zeta = 1;             % Deviation-Distance Ratio

%%%%%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

n=50;
for run=1:n	%%%%% Start main FOR loop

    rng(run);	

    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    % Initialization

    % Create Empty Individual Structure
    empty_individual.Position = [];
    empty_individual.Cost = [];

    % Create Population Matrix
    pop = repmat(empty_individual, nPop, 1);

    % Initialize Population Members
    for i = 1:nPop
    % Create Random Solution
    pop(i).Position=initializationACOR(VarSize,ub,lb);
    % Evaluation
    pop(i).Cost = Fun(fhandle, fnonlin, pop(i).Position);
    end

    % Sort Population
    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);

    % Update Best Solution Ever Found
    BestSol = pop(1);

    % Array to Hold Best Cost Values
    BestCost = zeros(iterMax, 1);

    % Solution Weights
    w = 1/(sqrt(2*pi)*q*nPop)*exp(-0.5*(((1:nPop)-1)/(q*nPop)).^2);

    % Selection Probabilities
    p = w/sum(w);

    %%%%% Start ACOR algorithm
    %tic	% Begin timing
    [BestSol, BestCost, ConvergenceCurve]=ACOR(fhandle, fnonlin, nVar,VarSize,lb, ub,iterMax,nPop,nSample,zeta,empty_individual,pop,BestSol,BestCost,p);
    %BestTime(run)=toc;% End timing
    BestConvergenceCurve(run,:)=BestCost;
    BestSolution(run,:)=BestSol.Position; %Storing Best Solution in every run
    BestFunVal(run)=BestCost(end);             %Storing Best Function Values (from the Best Solution) in every run
end
%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:)
%%%%% Plot ACOR convergence characteristic
ConvergenceCurveOfTheOverallBestFunction=BestConvergenceCurve(index,:);
plot(ConvergenceCurveOfTheOverallBestFunction,'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('Convergence characteristic of ACOR')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)