%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA104
% Project Title: Ant Colony Optimization for Continuous Domains (ACOR)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, ACO for Continuous Domains in MATLAB (URL: https://yarpiz.com/67/ypea104-acor), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

%% ACOR Main Loop
function [BestSol, BestCost, ConvergenceCurve]=ACOR(CostFunction, fnonlin, nVar,VarSize,VarMin, VarMax,MaxIt,nPop,nSample,zeta,empty_individual,pop,BestSol,BestCost,p)
% pop(1).Position
% pop(2).Position
for it = 1:MaxIt
    
    % Means
    s = zeros(nPop, nVar);
    for l = 1:nPop
        s(l, :) = pop(l).Position;
    end
    
    % Standard Deviations
    sigma = zeros(nPop, nVar);
    for l = 1:nPop
        D = 0;
        for r = 1:nPop
            D = D+abs(s(l, :)-s(r, :));
        end
        sigma(l, :) = zeta*D/(nPop-1);
    end
    
    % Create New Population Array
    newpop = repmat(empty_individual, nSample, 1);
    for t = 1:nSample
        
        % Initialize Position Matrix
        newpop(t).Position = zeros(VarSize);
        
        % Solution Construction
        for i = 1:nVar
            
            % Select Gaussian Kernel
            l = RouletteWheelSelection(p);
            
            % Generate Gaussian Random Variable
            newpop(t).Position(i) = s(l, i)+sigma(l, i)*randn;
            
        end
        
        % Apply Variable Bounds
        newpop(t).Position = max(newpop(t).Position, VarMin);
        newpop(t).Position = min(newpop(t).Position, VarMax);
        
        % Evaluation
%         newpop(t).Cost = CostFunction(newpop(t).Position);
        newpop(t).Cost = Fun(CostFunction,fnonlin, newpop(t).Position);
        
    end
   
    % Merge Main Population (Archive) and New Population (Samples)
    pop = [pop
         newpop]; %#ok
     
    % Sort Population
    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);
    
    % Delete Extra Members
    pop = pop(1:nPop);
    
    % Update Best Solution Ever Found
    BestSol = pop(1);
    
    % Store Best Cost
    BestCost(it) = BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    ConvergenceCurve=BestCost;
end