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

clear all 
clc
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
nInitPop=1000; % CLA Initial Population  
%------------------------------------------------------------------------------------------------------------------------

n=50;
for run=1:n	%%%%% Start main FOR loop

    rng('default');	
    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    Positions=initializationGWO(nPop, ub,lb);
    Positions0=initializationGWO(nInitPop, ub,lb);
    % --------------------------------------------------------------------------------------------------------------
    
    %%%% The finding of the elite individual begin
    for VarNo=1:nVar
        % Randomizing the VarNo's column only while maintaining the values in other column
        for i=1:nInitPop
            Positions0(i,VarNo)=rand().*(ub(VarNo)-lb(VarNo))+lb(VarNo);
        % Re-evaluating the objective function values of all the revised initial individuals
        f0(i)=Fun(fhandle, fnonlin, Positions0(i,:));
        end

        % Updating the best initial individual based on evaluation of objective functions of each individuals 
        [initialbestfun,index0]=min(f0);
        initialbest=Positions0(index0,:);
        % Replacing the values of the first 'VarNo' of all initial individuals with the the first 'VarNo' values of the newly obtained best initial
        Positions0(:,1:VarNo)=repmat(initialbest(1:VarNo),nInitPop,1);
    end
    % END of the column/decision variable. Here the best initial (the elite) individual has been found.

    % Preparing to store the final best individual before combining it with (n-1) other individuals to make up initial population for GWO
    Positions=zeros(nPop,nVar); 
    Positions(1,:)=Positions0(1,:);	

    %%%%% Start of GWO initialization (ultimate nPop)
    % Obtaining (n-1) individuals to make up initial population for GWO
    %(not within the neighbourhood of the previously found best initial elite)

    rng(run);	   

    for i=2:nPop
            Positions(i,:)=rand().*(ub-lb)+lb;
              f(i)=Fun(fhandle, fnonlin, Positions(i,:));
    end

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
title('Convergence characteristic of CLA+GWO')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)