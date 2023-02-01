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
%%%%% Metaheuristics' Parameters Settings
%PSO's Values of Parameters
wmax=0.9;	% maximum inertia weight
wmin=0.4;	 % minimum inertia weight
c1=2;		    % personal acceleration factor
c2=2;		    % global acceleration factor
%-----------------------------------------------------------------------------------------------------------------------

%%%%% Start Main program
n=50; % No. of Simulation Runs
for run=1:n	
    rng('default');
    fprintf('\n\nRun  ------>');
    fprintf('\t%8g\n',run);
    %%%% Start PSO preprocessing with initpop

    %%%% Initialize x0 and f0 vectors for efficiency purposes
    x0=zeros(nInitPop,nVar);     
    f0=zeros(nInitPop,1);

    %%%%% Random generation of initial starting individuals for CLA algorithm
    for i=1:nInitPop
        for j=1:nVar
            x0(i,j)=(LB(j)+rand()*(UB(j)-LB(j))); 
        end
    end   

    % --------------------------------------------------------------------------------------------------------------
    %%%% The finding of the elite individual begin
    for VarNo=1:nVar
        % Randomizing the VarNo's column only while maintaining the values in other column
        for i=1:nInitPop
            x0(i,VarNo)=(LB(VarNo)+rand()*(UB(VarNo)-LB(VarNo)));
        end
        % Re-evaluating the objective function values of all the revised initial individuals
        for i=1:nInitPop
            f0(i,1)=fun(x0(i,:));
        end
        % Updating the best initial individual based on evaluation of objective functions of each individuals 
        [initialbestfun,index0]=min(f0);
        initialbest=x0(index0,:);
        % Replacing the values of the first 'VarNo' of all initial individuals with the the first 'VarNo' values of the newly obtained best initial
        x0(:,1:VarNo)=repmat(initialbest(1:VarNo),nInitPop,1);
    end
    % END of the column/decision variable. Here the best initial (the elite) individual has been found.

    % Preparing to store the final best individual before combining it with (n-1) other individuals to make up initial population for PSO
    x=zeros(nPop,nVar); 
    x(1,:)=x0(1,:);	% The elite individual

    fmin0=initialbestfun;

    %%%%% Start of PSO initialization (ultimate nPop)
    % Obtaining (n-1) individuals to make up initial population for PSO
    %(not within the neighbourhood of the previously found best initial elite)

    rng(run) % This will ensure that both version, i.e. with and without CLA have same (n-1) seeds at each run 

    for i=2:nPop
        for j=1:nVar
            x(i,j)=(lb+rand()*(ub-lb)); % Creating initial population for PSO
        end
    end

    v=0.1*x;	%%%%% initial velocity

    f0=zeros(nPop,1); 
    for i=1:nPop
        f0(i,1)=fun(x(i,:)); %%%%% calculating initial population's fitness function for PSO
    end

    [fmin0,index0]=min(f0); %%%%% initial best function value for PSO
    pbest=x;	%%%%% initial pbest
    gbest=x(index0,:); %%%%% initial gbest

    ffmin=zeros(iterMax,1);
    f=zeros(nPop,1);

    %%%%% End of PSO initialization (ultimate npop)

    %%%%% Start PSO algorithm
    %tic	% Begin timing
    [ConvergenceCurve, bestfun, best_variables]=myPSO(nPop, iterMax, fun, nVar, wmax, wmin, c1, c2, LB, UB, x, v, pbest, gbest, ffmin, f, f0, fmin0);
    %BestTime(run)=toc;% End timing
    BestConvergenceCurve(run,:)=ConvergenceCurve;
    BestSolution(run,:)=best_variables; %Storing Best Solution in every run
    BestFunVal(run)=bestfun;               %Storing Best Function Values (from the Best Solution) in every run

end  
%%%%% Overall Optimum Values
[bestfunction, index]=min(BestFunVal);
TheOverallBestPosition=BestSolution(index,:)
%%%%% Plot PSO convergence characteristic
ConvergenceCurveOfTheOverallBestFunction=BestConvergenceCurve(index,:);
plot(ConvergenceCurveOfTheOverallBestFunction,'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('Convergence characteristic of CLA+PSO')

%%%%% Analysing the overall results from the simulation runs
AverageValue=mean(BestFunVal)
StandardValue=std(BestFunVal)
MinimumValue=min(BestFunVal)
MaximumValue=max(BestFunVal)
MedianValue=median(BestFunVal)
