function [ t, x, v, z, Force ] = simulation( physics, simu, vref )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    m = physics.m;
    g = physics.g;
    k = physics.k;
    l0 = physics.l0;
    muS = physics.muS;

    deltaT = simu.deltaT;
    tmax = simu.tmax;
    kmax = ceil(tmax/deltaT);
    
    t = (0:(kmax-1))*deltaT;
    x = zeros(1, kmax);
    v = zeros(1, kmax);
    z = zeros(1, kmax);
    Force = zeros(1, kmax);
    
    x(1) = simu.x0;
    v(1) = simu.v0;
    z(1) = simu.z0;
    Force(1) = friction(v(1), physics);
   
    for i=2:kmax
        % Model 1
        %v(i) = v(i-1) + deltaT*(-friction(v(i-1), physics) + m*g*sin(alpha));
        %x(i) = x(i-1) + deltaT*v(i);
        
        % Model 2
        z(i) = z(i-1) + deltaT*(v(i-1) - vref - l0);
        Force(i) = friction(v(i-1), physics);
        if Force(i) == 0
            Force(i) = Force(i-1);
        end
        v(i) = v(i-1) + deltaT*(-Force(i) - k*z(i))/m;
        
        if v(i)*v(i-1) <= 0
            F = v(i-1)*m/deltaT - k*z(i);
            if F <= muS*m*g && F >= -muS*m*g
                v(i) = 0;
                Force(i) = F;
            end
        end
        
        x(i) = x(i-1) + deltaT*v(i);
        
    end

end

