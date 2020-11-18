function [ PgMax ] = globalLMI( physics, vref )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    m = physics.m;
    g = physics.g;
    kv = physics.kv;
    k = physics.k;
    FS = physics.muS*m*g;
    FC = physics.muC*m*g;

    % Define the decision variables
    Pg = sdpvar(2);
    tau1 = sdpvar(1);
    tau1Prime = sdpvar(1);
    tau2 = sdpvar(1);
    tau3 = sdpvar(1);
    tau4 = sdpvar(1);
    
    % Creation of Theta
    A = [-kv/m -k/m; 1 0];
    B = [-1/m; 0];
    Fnl = @(thetaDot) friction(thetaDot, physics) - kv*thetaDot;
    D = [A B -Fnl(vref)*B];
    F = eye(2, 4);
    
    Theta = D'*Pg*F + F'*Pg*D;
    
    % Creation of Pi matrices
    pi1 = [1 0 0 0];
    pi3 = [0 0 1 0];
    pi4 = [0 0 0 1];
    
    Pi0 = (pi4'*pi4) - F'*Pg*F;
    Pi1 = pi3'*pi3 - FS^2*(pi4'*pi4);
    Pi2 = -pi3'*pi3 + FC^2*(pi4'*pi4);
    Pi3 = -(pi1 + vref*pi4)'*pi3;
    Pi3 = Pi3' + Pi3;
    Pi4 = (pi + vref*pi4)'*(pi + vref*pi4);
    % Creation of ThetaBar
    
    tau0List = linspace(0, -2*max(real(eig(A))), 30);
    
    geomeanMax = 0;
    PgMax = nan;
    for tau0 = tau0List

        ThetaBar = Theta - tau0*Pi0 - tau1*Pi1 - tau2*Pi2 - tau3*Pi3;
        ThetaBar2 = Theta - tau0*Pi0 - tau1Prime*Pi1 - tau4*Pi4;
        
        % Soving the optimization problem    
        Constraints = [(Pg >= 1e-5):'Positivity',...
            (tau1 >= 0):'S-variable F_{max}',...
            (tau1Prime >= 0):'S-variable F_{max}',...
            (tau2 >= 0):'S-variable F_{min}'...
            (tau3 >= 0):'S-variable sign'...
            (ThetaBar <= -1e-5):'Negativity',...
            (ThetaBar2 <= -1e-5):'Negativity'];

        warning('off', 'YALMIP:strict');
        option = sdpsettings('verbose', 0, 'solver', 'mosek');
        diagnostics = optimize(Constraints, -geomean(Pg), option);

        if diagnostics.problem == 0
            if geomean(value(Pg),'all') >= geomeanMax
                geomeanMax = geomean(value(Pg),'all');
                PgMax = value(Pg);
            end
        end
    end

end

