function [lambdaMin] = lambdaMin(physics, vref, rl)

    m = physics.m;
    g = physics.g;
    muS = physics.muS;
    muC = physics.muC;
    vs = physics.vs;

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
    
    partialPhi = @(eps1) -2*m*g*(muS - muC)*(eps1+vref)/vs^2*exp(-(eps1+vref)^2/vs^2);
    phi = @(eps1) m*g*(muS - muC)*(exp(-((eps1+vref)/vs)^2) - exp(-(vref/vs)^2));

    if max([abs(vStarMinus), abs(vStarPlus)]) <= rl
        lambdaMin = max([-partialPhi(vStarMinus), -partialPhi(vStarPlus)]);
    else
        lambdaMin = -min([-phi(-rl)/rl, phi(rl)/rl]);
    end

end

