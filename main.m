% Travail sur les forces de friction 
clear
close all
clc

%% Parameters
vref = 1.45;

% Physical parameters
m = 1;
muC = 0.2997;
muS = 0.5994;
kv = 1;
vs = 0.8;
alpha = 0;
g = 9.81;
k = 2;
l0 = 0;
kp = 0;

% Simulation parameters
tmax = 20;
deltaT = 0.003*min(1, vref);

x0 = 2*l0;
v0 = 1.6;
z0 = -2.1;

v02 = 1.7;
z02 = 5;

% Generates structure
physics = struct('m', m, 'muS', muS, 'muC', muC, 'kv', kv, 'vs', vs,...
    'g', g, 'k', k, 'l0', l0);
simu = struct('deltaT', deltaT, 'tmax', tmax, 'x0', x0, 'v0', v0, 'z0', z0);
simu2 = struct('deltaT', deltaT, 'tmax', tmax, 'x0', x0, 'v0', v02, 'z0', z02);

%% Plot of the friction force
figure
thetaDot1 = -5:0.01:-0.01;
thetaDot2 = 0.01:0.01:5;
thetaDot = [thetaDot1 0 thetaDot2];
plot(thetaDot1, friction(thetaDot1, physics), '-b', 'LineWidth', 2);
hold on
plot(thetaDot, friction(thetaDot, physics)-kv*thetaDot, '--r', 'LineWidth', 2);
plot(thetaDot, friction(thetaDot, physics), '-b', 'LineWidth', 2);
grid on
xlabel('$\dot{x}$ [m.s${}^{-1}$]','interpreter','latex','FontSize',16);
ylabel('Force value [N]', 'interpreter','latex','FontSize',16);
h = legend('$F$', '$F_{nl}$');
set(h, 'interpreter','latex','FontSize',11,'location','northwest');
set(gcf, 'Position', [20, 20, 560, 320]);

figure
Fnl = @(theta) friction(theta, physics) - kv*theta;
phi = @(theta) Fnl(theta+vref) - Fnl(vref);
thetaT = linspace((-vref*1.05),0,100);
thetaT = [thetaT, linspace(0,(vref*1.05))];
lambda = lambdaMin3(physics, vref, vref);
lambda2 = lambdaMin3(physics, vref, vref*0.8);
plot(thetaT, phi(thetaT), 'LineWidth', 2);
hold on
plot(thetaT, -lambda*thetaT, '--', 'LineWidth', 2);
plot(thetaT, -lambda2*thetaT, '-.', 'LineWidth', 2);
xlabel('$\varepsilon_1$ [m.s${}^{-1}$]','interpreter','latex','FontSize',16);
ylabel('Force value [N]', 'interpreter','latex','FontSize',16);
h = legend('$\phi_{v_{ref}}$',...
    strcat(['$\theta \mapsto -\lambda_{min} \theta, \ r_l = ', num2str(vref), '$']),...
    strcat(['$\theta \mapsto -\lambda_{min} \theta, \ r_l = ', num2str(vref*0.8), '$']));
set(h, 'interpreter','latex','FontSize',11);
set(gcf, 'Position', [20, 20, 560, 320]);
grid on
axis([-vref*1.05, vref*1.05, min(-lambda*thetaT), max(-lambda*thetaT)])


% Global assessment
% Pg = globalLMI(physics, vref);
Pg = globalLMI(physics, vref);
geomean(Pg,'all')

% Local assessment
[vref1, vref2] = A0Hurwitz(physics)
Pl = localLMI(physics, vref);
eig(Pl)

% Simulation
[ t, x, v, z, F ] = simulation( physics, simu, vref );
[ ~, x2, v2, z2, F2 ] = simulation( physics, simu2, vref );

figure
plot(t, v, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'LineStyle', '-')
hold on
plot(t, v2, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2, 'LineStyle', '-.')
grid on
xlabel('$t$ [s]','interpreter','latex','FontSize',16);
ylabel('$v$ [m.s${}^{-1}$]', 'interpreter','latex','FontSize',16);
set(gcf, 'Position', [20, 20, 560, 320]);
h = legend('$v_{ref} = 1$ m.s${}^{-1}$', '$v_{ref} = 1.6$ m.s${}^{-1}$');
set(h,'interpreter','latex','FontSize',11,'Location','northeast');

figure
plot(t, F, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
hold on
plot(t, F2, 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2)
grid on
xlabel('$t$ [s]','interpreter','latex','FontSize',16);
ylabel('Force [N]', 'interpreter','latex','FontSize',16);
set(gcf, 'Position', [20, 20, 560, 320]);


figure
hold on
zinf = -friction(vref, physics)/k; 
axisLimitPlus = sqrt(1/min(eig(Pg)))+max(vref, zinf);
fimplicit(@(v,z) [v-vref;z-zinf]'*Pg*[v-vref;z-zinf] - 1, [-axisLimitPlus axisLimitPlus], 'LineWidth', 1, 'LineStyle', '--');
try
    axisLimitPlus = sqrt(1/min(eig(Pl)))+max(vref, zinf); %[-0.01 0.15 -3.1 -2.8]
    fimplicit(@(v,z) [v-vref;z-zinf]'*Pl*[v-vref;z-zinf] - 1, [-2 2 -3 -3], 'LineWidth',1, 'LineStyle', '-.');
catch
    plot(0,0,'LineWidth',1, 'LineStyle', '-.');
    fprintf('No local stability analysis.')
end
plot(v, z, 'Color', [0.4660 0.6740 0.1880]	, 'LineWidth',1)
plot(v2, z2, 'Color', [0.3010 0.7450 0.9330], 'LineWidth',1)

xlabel('$v$ [m.s${}^{-1}$]','interpreter','latex','FontSize',16);
ylabel('z [m]', 'interpreter','latex','FontSize',16);
set(gcf, 'Position', [20, 20, 560, 320]);
grid on

h = legend('Attractor', 'Basin of attraction', 'Trajectory 1', 'Trajectory 2');
set(h,'interpreter','latex','FontSize',11,'AutoUpdate','off','Location','southwest');

plot(v0,z0,'Color',[0.4660 0.6740 0.1880]	,'Marker','d') % Initial point 1
plot(v02,z02,'Color',[0.3010 0.7450 0.9330],'Marker','d') % Initial point 2
plot(vref,zinf,'r*') % Eq. point
space = 0;
space2 = 0;
maxSpace = 5;
for i = 1:(length(t)-1)
    if space >= maxSpace
        arrowh([v(i) v(i+1)],[z(i), z(i+1)], [0.4660 0.6740 0.1880]);
        space = 0;
    end
    space = space + sqrt((v(i+1)-v(i))^2 + (z(i+1)-z(i))^2);
    
    if space2 >= maxSpace
        arrowh([v2(i) v2(i+1)],[z2(i), z2(i+1)], [0.3010 0.7450 0.9330]);
        space2 = 0;
    end
    space2 = space2 + sqrt((v2(i+1)-v2(i))^2 + (z2(i+1)-z2(i))^2);
end
