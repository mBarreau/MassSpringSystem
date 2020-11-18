function [lambdamin] = lambdaMin(physics, vref, vDesired)

    m = physics.m;
    g = physics.g;
    muS = physics.muS;
    muC = physics.muC;
    vs = physics.vs;
    kv = physics.kv;

    f = @(x) 1 - exp((x.^2 + 2*x*vref)/vs^2) + 2*x.*(x+vref)/vs^2;
    
    epsilon = 3; % error in e-epsilon
    i = max(0, ceil(-log10(vref/10)));
    epsilon = max(epsilon, i); % error in e-epsilon
    v = -vref;
    while i <= epsilon
        if f(v)*f(v+10^(-i)) <= 0
            if v*(v+10^(-i)) >= 1e-5
                i = i + 1;
                continue
            end
        end
        v = v+10^(-i);
    end
    
    vStarMinus = v;
    vStarPlus = v+10^(-epsilon);
    
    if vref < vs/sqrt(2)
        
        if vDesired >= vStarMinus
            lambdamin = 2*m*g*(muS - muC)*(vStarMinus + vref)*exp(-(vStarMinus + vref)^2/vs^2)/vs^2;
        else
            phi = friction(vDesired+vref,physics) - kv*vDesired - friction(vref,physics);
            lambdamin = -phi/vDesired;
        end
        
    else
        
        if -vDesired <= vStarPlus
            lambdamin = 2*m*g*(muS - muC)*(vStarPlus + vref)*exp(-(vStarPlus + vref)^2/vs^2)/vs^2;
        else
            phi = friction(-vDesired+vref,physics) + kv*vDesired - friction(vref,physics);
            lambdamin = phi/vDesired;
        end
        
    end
    
%    lambdamax = min()
    
    if lambdamin <= 0
        lambdamin = 2*m*g*(muS - muC)*(vs/sqrt(2))*exp(-1/2)/vs^2;
    end

end

