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
nInitPop=1000; % CLA Initial Population  
%------------------------------------------------------------------------------------------------------------------------
%%%%% Metaheuristics' Parameters Settings
%ACOR's Values of Parameters
nSample = 40;         % Sample Size
q = 0.5;              % Intensification Factor (Selection Pressure)
zeta = 1;             % Deviation-Distance Ratio

%%%%%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%% Start Main program
n=50; % No. of Simulation Runs                                                 - CHANGE HERE
for run=1:n	
    rng('default');
    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    %%%% Start ACOR initialization (preprocessing with initpop)

    % Create Empty Individual Structure for Efficiency
    empty_individual.Position = [];
    empty_individual.Cost = [];
    % Create Population Matrix
    pop = repmat(empty_individual, nPop, 1);
    pop0= repmat(empty_individual, nPop, 1);
    % Initialize Population Members
    for i = 1:nInitPop
        pop(i).Position=initializationACOR(VarSize,ub,lb);
        pop(i).Cost = Fun(fhandle, fnonlin, pop(i).Position);
    end

    % --------------------------------------------------------------------------------------------------------------
    %%%% The finding of the elite individual begin
    for VarNo=1:nVar
        % Randomizing the VarNo's column only while maintaining the values in other column
        for i=1:nInitPop
            pop0(i).Position=initializationACOR(VarSize,ub,lb);
            pop0(i).Position(VarNo)=unifrnd(lb(VarNo), ub(VarNo));
            % Re-evaluating the objective function values of all the revised initial individuals
            pop0(i).Cost =Fun(fhandle, fnonlin, pop0(i).Position(1,:));
        end
        % Updating the best initial individual based on evaluation of objective functions of each individuals 
        [~, SortOrder] = sort([pop0.Cost]);
        pop0 = pop0(SortOrder);
        % Replacing the values of the first 'VarNo' of all initial individuals with the the first 'VarNo'
        % values of the newly obtained best initial

        for i=1:nInitPop
            pop0(i).Position(1:VarNo)=pop0(1).Position(1:VarNo);
        end
    end % END of the column/decision variable. Here the best initial has been found.

    % Update Best Solution Ever Found
    BestSol = pop0(1);
    pop(1) =BestSol;

    %%%%% Start of ACOR initialization (ultimate nPop)
    % Obtaining (n-1) individuals to make up initial population for ACOR
    %(not within the neighbourhood of the previously found best initial elite)

    rng(run) % This will ensure that both version, i.e. with and without CLA have same (n-1) seeds at each run 

    for i=2:nPop
        pop(i).Position=initializationACOR(VarSize,ub,lb);
    pop(i).Cost = Fun(fhandle, fnonlin, pop(i).Position); %%%%% calculating initial population's fitness function for ACOR
    end

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
title('Convergence characteristic of CLA+ACOR')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)