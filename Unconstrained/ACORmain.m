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

%-----------------------------------------------------------------------------------------------------------------------
%%%%% List of Test Functions/Objective Functions %%%%%
%Unimodal Functions
% F1 - Sphere
% F2 - Powell Sum
% F3 - Schwefel 2.20
% F4 - Schwefel 2.21
% F5 - Step
% F6 - Schwefel 2.22
% F7 - Schwefel 2.23
% F8 - Rosenbrock
% F9 - Brown
% F10 - Dixon & Price
% F11 - Powell Singular
% F12 - Perm 0,D,Beta
% F13 - Sum Squares
%Multimodal Functions
% F14 - Schwefel 2.26
% F15 - Rastrigin
% F16 - Periodic
% F17 - Qing
% F18 - Alpine N. 1
% F19 - Xin-She Yang
% F20 - Ackley
% F21 - Trignometric 2
% F22 - Salomon
% F23 - Styblinski-Tang
% F24 - Griewank
% F25 - Xin-She Yang N. 4
% F26 - Xin-She Yang N. 2
% F27 - Gen. Pendlized
% F28 - Pendlized
% F29 - Michalewics
% F30 - Quartic Noise
% F31 - Shifted & Rotated Weierstrass
% F32 - Shifted & Rotated Happy Cat
%-----------------------------------------------------------------------------------------------------------------------
%%%%% Experimental Settings
Function_name='F1'; % Name of the test function                                    - CHANGE HERE
nPop=20; % Number of search agents for PSO/ACOR/GWO                    - CHANGE HERE
iterMax=1000; % Maximum number of iterations for PSO/ACOR/GWO      - CHANGE HERE
% Load details of the selected test / objective function
[lb,ub,nVar,fun]=Test_Function(Function_name);
LB=repmat(lb,1,nVar);
UB=repmat(ub,1,nVar);
VarSize = [1 nVar];   % Variables Matrix Size
%------------------------------------------------------------------------------------------------------------------------
%%%%% Metaheuristics' Parameters Settings
% ACOR's Values of Parameters
nSample = 40;         % Sample Size
q = 0.5;                    % Intensification Factor (Selection Pressure)
zeta = 1;                   % Deviation-Distance Ratio
%-----------------------------------------------------------------------------------------------------------------------

%%%%% Start Main program
n=50; % No. of Simulation Runs                                                                - CHANGE HERE
for run=1:n
    rng(run); % This will ensure that both version, i.e. with and without CLA have same seeds at each run
    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    %%%%% Start of ACOR initialization
    % Create Empty Individual Structure
    empty_individual.Position = [];
    empty_individual.Cost = [];
    % Create Population Matrix
    pop = repmat(empty_individual, nPop, 1);
    % Initialize Population Members
    for i = 1:nPop
        % Create Random Solution
        pop(i).Position = unifrnd(lb, ub, VarSize);
        % Evaluation
        pop(i).Cost = fun(pop(i).Position);
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
    %%%%% End of ACOR initialization

    %%%% Start ACOR main program
    %tic	%%%%% Begin timing
    [BestSol, BestCost]=ACOR(fun, nVar,VarSize,lb, ub,iterMax,nPop,nSample,zeta,empty_individual,pop,BestSol,BestCost,p);
    %BestTime(run)=toc;% End timing
    BestConvergenceCurve(run,:)=BestCost;
    BestSolution(run,:)=BestSol; %Storing Best Solution in every run
    BestFunVal(run)=BestCost(end);               %Storing Best Function Values (from the Best Solution) in every run
end
%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:).Position
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
