function [ PlMin ] = localLMI( physics, vref )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = physics.m;
g = physics.g;
kv = physics.kv;
k = physics.k;
muS = physics.muS;
muC = physics.muC;
vs = physics.vs;

[vref1, vref2] = A0Hurwitz(physics);
if vref >= vref1 && vref <= vref2
    PlMin = nan;
    return
end

% Define the decision variables
Pl = sdpvar(2);
tau = sdpvar(1);
eta = sdpvar(1);
gamma = sdpvar(1);
tau1 = sdpvar(1);

% Creation of Theta
A = [-kv/m -k/m; 1 0];
B = [-1/m; 0];
C = [1 0];
Gamma = -2*m*g*(muS-muC)*vref*exp(-vref^2/vs^2)/vs^2;
A0 = A + B*Gamma*C;

PlMin = nan;
etaMin = Inf;
for rl = linspace(vref, 1e-5)
    
    lambda1 = lambdaMin(physics,vref,rl);
    
    %ThetaLocal = [A0'*Pl + Pl*A0 - 2*tau*Gamma*(Gamma + lambda)*(C'*C), Pl*B - tau*(2*Gamma + lambda)*C';...
    %    B'*Pl - tau*(2*Gamma + lambda)*C, -2*tau];
    
    ThetaLocal = [A0'*Pl + Pl*A0 - 2*Gamma*(tau*Gamma + gamma)*(C'*C), Pl*B - (2*tau*Gamma + gamma)*C';...
        B'*Pl - (2*tau*Gamma + gamma)*C, -2*tau] + tau1*[zeros(2), C'; C, 0];
    
    %ThetaLocalA = [A'*Pl + Pl*A, Pl*B - tau*lambda*C';...
    %    B'*Pl - tau*lambda*C, -2*tau];
    
    ThetaIncl = [Pl, C'; C, rl^2];
    
    % Soving the optimization problem
    Constraints = [(Pl >= 0):'Positivity',...
        (tau >= 0):'S-variable',...
        (tau1 >= 0):'S-variable',...
        (ThetaLocal <= -1e-5):'Negativity',...
        (ThetaIncl >= 0):'Inclusion',...
        (Pl <= eta*eye(2)):'Objective function',...
        (gamma >= tau*lambda1):'lambdaMin'];
    
    warning('off', 'YALMIP:strict');
    option = sdpsettings('verbose', 0, 'solver', 'mosek');
    diagnostics = optimize(Constraints, eta, option);
    
    if diagnostics.problem == 0
        if value(eta) >= etaMin
            break
        end
        PlMin = value(Pl);
        etaMin = value(eta); 
    end
end

end

