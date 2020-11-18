function [ frictionForce ] = friction( thetaDot, physics )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    muS = physics.muS;
    muC = physics.muC;
    kv = physics.kv;
    vs = physics.vs;
    RN = physics.m*physics.g;

    frictionForce = RN*(muC + (muS - muC)*exp(-(abs(thetaDot)/vs).^2)).*sign(thetaDot)...
        + kv*thetaDot;

end

