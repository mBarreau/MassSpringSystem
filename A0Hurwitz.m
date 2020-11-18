function [vref1,vref2] = A0Hurwitz(physics)

    m = physics.m;
    g = physics.g;
    kv = physics.kv;
    vs = physics.vs;
    muS = physics.muS;
    muC = physics.muC;

    theta = @(v) v.*exp(-v.^2/vs^2)-kv*vs^2/(2*m*g*(muS-muC));
    [vref1, ~, flag] = fsolve(theta, 0, optimoptions('fsolve','Display','off'));
    
    if flag <= 0
       vref1 = nan;
       vref2 = nan;
       return
    end
    
    i = 0;
    v0 = vs*sqrt(5);
    while i <= 5
        [vref2, ~, flag] = fsolve(theta, v0, optimoptions('fsolve','Display','off'));
        if vref2 <= 1.001*vref1
            v0 = v0*1.1;
        else
            if flag <= 0
                v0 = (v0 + vref1)/2;
                i = i+1;
            else
               return 
            end
        end
    end
    vref2 = vref1;
end

