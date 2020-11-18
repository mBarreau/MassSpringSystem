function [globalAsympStab] = globalAsympLMI(physics, vref)
    globalAsympStab = 0;
    
    m = physics.m;
    g = physics.g;
    kv = physics.kv;
    k = physics.k;
    muS = physics.muS;
    muC = physics.muC;
    vs = physics.vs;
    
    FS = muS*m*g;
    FC = muC*m*g;

    % Define the decision variables
    Pg = sdpvar(2);
    Pl = sdpvar(2);
    tau1 = sdpvar(1);
    tau1Prime = sdpvar(1);
    tau2 = sdpvar(1);
    tau3 = sdpvar(1);
    tau4 = sdpvar(1);
    tau = sdpvar(1);
    gamma = sdpvar(1);
    
    % Creation of Theta
    A = [-kv/m -k/m; 1 0];
    B = [-1/m; 0];
    C = [1 0];
    
    Fnl = @(thetaDot) friction(thetaDot, physics) - kv*thetaDot;
    D = [A B -Fnl(vref)*B];
    F = eye(2, 4);
    
    Gamma = -2*m*g*(muS-muC)*vref*exp(-vref^2/vs^2)/vs^2;
    A0 = A + B*Gamma*C;
    
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
    
    for rl = linspace(vref, 1e-5, 30)
        
        lambda = lambdaMin3(physics,vref,rl);
        ThetaLocal = [A0'*Pl + Pl*A0 - 2*Gamma*(tau*Gamma + gamma)*(C'*C), Pl*B - (2*tau*Gamma + gamma)*C';...
                B'*Pl - (2*tau*Gamma + gamma)*C, -2*tau];
        ThetaIncl = [Pl, C'; C, rl^2];
        
        for tau0 = tau0List
            
            ThetaBar = Theta - tau0*Pi0 - tau1*Pi1 - tau2*Pi2 - tau3*Pi3;
            ThetaBar2 = Theta - tau0*Pi0 - tau1Prime*Pi1 - tau4*Pi4;
        
            % Soving the optimization problem    
            Constraints = [(Pl >= 0):'Positivity 1',...
                (Pg - Pl >= 0):'Attractor in Basin',...
                (tau1 >= 0):'S-variable F_{max}',...
                (tau1Prime >= 0):'S-variable F_{max}',...
                (tau2 >= 0):'S-variable F_{min}'...
                (tau3 >= 0):'S-variable sign'...
                (ThetaBar <= -1e-5):'Negativity',...
                (ThetaBar2 <= -1e-5):'Negativity',...
                (tau >= 0):'S-variable',...
                (ThetaLocal <= -1e-5):'Negativity',...
                (gamma >= tau*lambda):'lambdaMin',...
                (ThetaIncl >= 0):'Inclusion'];

            warning('off', 'YALMIP:strict');
            option = sdpsettings('verbose', 0, 'solver', 'mosek');
            diagnostics = optimize(Constraints, {}, option);

            if diagnostics.problem == 0
                globalAsympStab = 1;
                return
            end
        end
    end

end

