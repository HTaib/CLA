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
%%%%% CLA setting
nInitPop=1000; % CLA Initial Population                                                                                      - CHANGE HERE                                                        
%------------------------------------------------------------------------------------------------------------------------

%%%%% Start Main program
n=50; % No. of Simulation Runs
for run=1:n	%%%%% Start main FOR loop
    rng('default');
    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);

    %%%%% Start of GWO initialization
    % Initialize alpha, beta, and delta_pos
    Alpha_pos=zeros(1,nVar);
    Alpha_score=inf; 

    Beta_pos=zeros(1,nVar);
    Beta_score=inf; 
    
    Delta_pos=zeros(1,nVar);
    Delta_score=inf; 

    Convergence_curve=zeros(1,iterMax);

    %%%%% Random generation of initial starting individuals for CLA algorithm
    for i=1:nInitPop
        Positions0(i,:)=rand(1,nVar).*(ub-lb)+lb;
    end   
% Evaluating the fitness function values of each individuals      
    for i=1:nInitPop
        f0(i)=fun(Positions0(i,:));
    end   

    % --------------------------------------------------------------------------------------------------------------
    %%%% The finding of the elite individual begin
    for VarNo=1:nVar
        % Randomizing the VarNo's column only while maintaining the values in other column
        for i=1:nInitPop
            Positions0(i,VarNo)=rand().*(UB(VarNo)-LB(VarNo))+LB(VarNo);
            % Re-evaluating the objective function values of all the revised initial individuals
            f0(i)=fun(Positions0(i,:));
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
    Positions(1,:)=Positions0(1,:);	%This elite individual

    %%%%% Start of GWO initialization (ultimate nPop)
    % Obtaining (n-1) individuals to make up initial population for GWO
    %(not within the neighbourhood of the previously found best initial elite)

    rng(run);	% This will ensure that both version, i.e. with and without CLA have same (n-1) seeds at each run 

    for i=2:nPop
        Positions(i,:)=rand().*(ub-lb)+lb; % Creating initial population for GWO
         f(i)=fun(Positions(i));
    end

    %%%%% Start GWO algorithm
    %tic	% Begin timing
    [Best_score,Best_pos,GWO_cg_curve]=GWO(nPop,iterMax,LB,UB,nVar,fun,Positions,Alpha_pos,Alpha_score,Beta_pos,Beta_score,Delta_pos,Delta_score,Convergence_curve);
    %BestTime(run)=toc;  %%%%End timing
    BestGWO_cg_curve(run,:)=GWO_cg_curve;
    BestSolution(run,:)=Best_pos;
    BestFunVal(run)=Best_score;
end

%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:)
%%%%% Plot GWO convergence characteristic
ConvergenceCurveOfTheOverallBestFunction=BestGWO_cg_curve(index,:);
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
