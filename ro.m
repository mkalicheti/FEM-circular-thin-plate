%{
FEA code to analyse the deflections and stresses in a circular thin plate
-static
-linear
-isotropic
%}

clc;
clear;
close all;

%% Properties and constants
E = 90e9; %(in Pa) Young's Modulus
v = 0.23; % Poisson's ratio
g = 9.8; % (in m/s^2) acceleration due to gravity
rho = 2.7*10^3; % (in kg/m^3) density

%% Dimensions (in m)
Ri = 0.01; %inner radius
Ro = 1; %outer radius
h = 0.05; %thickness

nelm = 20; %number of elements
nnod = nelm+1;%number of nodes
le = (Ro-Ri)/nelm; %length of element

D = E*h^3/(12*(1 - v^2)); % (in Nm) Bending Stress
q0 = -rho*h*g; %distributed load
GDof = 2*nnod; %total degrees of freedom

K = zeros(GDof, GDof); %Global stiffness matrix initialisation
f = zeros(GDof, 1); %Global force vector initialisation

%% Finite Element Assembly
for i = 1:nelm
    a = Ri + (i-1)*le; b = Ri + i*le; %inner and outer radii of element
    [Ke, fe] = element_calc(a, b, D, v, q0);
    dof = 2*i-1 : 2*(i+1);
    K(dof, dof) = K(dof, dof) + Ke;
    f(dof, 1) = f(dof, 1) + fe;
end

%% Boundary conditions
% simply supported
fixedDof = [2*nnod-1];
%clamped
%fixedDof = [2*nnod-1 2*nnod];

%% FEM Solution
activeDof = setdiff(1:GDof, fixedDof);
U = zeros(GDof, 1);
U(activeDof) = K(activeDof, activeDof)\f(activeDof);
w = U(1:2:end); %deflection
theta = U(2:2:end); %slope
F = K*U;
R = F(fixedDof) - f(fixedDof); %Reaction force

%% Stress calculation
sigmarr = zeros(nelm, 1); %Radial stress
sigmatheta = zeros(nelm, 1); %Hoope's stress
for i = 1:nelm
    a = Ri + (i-1)*le; b = Ri + i*le;
    ue = U(2*i-1:2*(i+1), 1);
    [sigmarr(i), sigmatheta(i)] = Stress(a, b, v, ue, E, h);
end

%% Analytical solution
% simply supported
syms x;
beta = Ri/Ro;
kappa = (beta^2/(1-beta^2))*log(beta);
alpha1 = (3+v)*(1-beta^2)-4*(1+v)*beta^2*kappa;
alpha2 = (3+v)+4*(1+v)*kappa;
w_exact = (q0*Ro^4/(64*D))*(-(1-(x/Ro)^4)+(2*alpha1/(1+v))*(1-(x/Ro)^2)-(4*alpha2*beta^2/(1-v))*log(x/Ro));
theta_exact = (q0*Ro^4/(64*D))*(4*x^3/Ro^4 - ((4*alpha1)/(1+v))*(x/Ro^2) - ((4*alpha2*beta^2)/(1-v))*(1/x));
theta_exact_ = (q0*Ro^4/(64*D))*(12*x^2/Ro^4 - ((4*alpha1)/(1+v))*(1/Ro^2) - ((4*alpha2*beta^2)/(1-v))*(-1/x^2));
sigmarr_exact = (E*h/(2*(1-v^2))) * (theta_exact_ + (v/x)*theta_exact);
sigmatheta_exact = (E*h/(2*(1-v^2))) * (v*theta_exact_ + (1/x)*theta_exact);

%% Plots
%Deflection
figure;
fplot(w_exact,[Ri Ro], 'b'); %exact
hold on; grid on;
r = linspace(Ri, Ro, nnod);
plot(r,w , 'ro-'); %FEM
xlabel("Radial distance($r$) in metres", 'Interpreter', 'Latex');
ylabel("Deflection($w$) in metres", 'Interpreter','Latex');
title("$w$ vs $r$", 'Interpreter', 'Latex');
legend("Exact", "FEM", 'location', 'northwest');

%Theta
figure;
fplot(theta_exact, [Ri Ro], 'b'); %exact
hold on; grid on;
plot(r,theta , 'ro-'); %FEM
xlabel("Radial distance($r$) in metres", 'Interpreter', 'Latex');
ylabel("$\theta$(in radians)", 'Interpreter', 'Latex');
title("$\theta$ vs $r$", 'Interpreter', 'Latex');
legend("Exact", "FEM", 'location', 'northwest');

%Radial stress
figure;
r = linspace(Ri + le, Ro-le, nelm);
fplot(sigmarr_exact, [Ri Ro], 'b'); %exact
hold on; grid on;
plot(r, sigmarr, 'ro'); %FEM
xlabel("Radial distance($r$) in metres", 'Interpreter', 'Latex');
ylabel("Radial stress($\sigma_{rr}$) in $N/m^2$", 'Interpreter', 'Latex');
title("$\sigma_{rr}$ vs $r$", 'Interpreter', 'Latex');
legend("Exact", "FEM");

%Hoope's stress
figure;
fplot(sigmatheta_exact, [Ri Ro], 'b');
hold on; grid on;
plot(r, sigmatheta, 'ro');
xlabel("Radial distance($r$)", 'Interpreter', 'Latex');
ylabel("Hoope's stress($\sigma_{\theta\theta}$) in $N/m^2$", 'Interpreter', 'Latex');
title("$\sigma_{\theta\theta}$ vs $r$", 'Interpreter', 'Latex');
legend("Exact", "FEM");

%% Outputs
Dof_w = 1:2:GDof;
Dof_theta = 2:2:GDof;
format
deflections = [Dof_w' w]
slopes = [Dof_theta' theta]
fprintf("reactions = \n\n     %d       %0.0f\n", fixedDof, R);

%% Helper functions

function [sigmarr, sigmatheta] = Stress(a, b, v, we, E, h)
    syms r;
    N1(r) = (3*a - b - 2*r)*(b - r)^2/(a-b)^3;
    N2(r) = -(a - r)*(b - r)^2/(a-b)^2;
    N3(r) = (a - r)^2*(a - 3*b + 2*r)/(a - b)^3;
    N4(r) = -(a - r)^2*(b - r)/(a - b)^2;
    Ne(r) = [N1 N2 N3 N4];
    N_e(r) = diff(Ne(r), r);
    Be(r) = diff(N_e(r), r);
    p(r) = E*h*(Be*we + (v/r)*((N_e)*(we)))/(2*(1-v^2));
    q(r) = (E*h/(2*(1-v^2))) * (v*Be*we + (1/r)*((N_e)*(we)));
    sigmarr = p((a+b)/2);
    sigmatheta = q((a+b)/2);
end

function [Ke, fe] = element_calc(a, b, D, v, q)     
     syms r real;
     N1(r) = (3*a - b - 2*r)*(b - r)^2/(a-b)^3;
     N2(r) = -(a - r)*(b - r)^2/(a-b)^2;
     N3(r) = (a - r)^2*(a - 3*b + 2*r)/(a - b)^3;
     N4(r) = -(a - r)^2*(b - r)/(a - b)^2;
     N(r) = [N1 N2 N3 N4];
     N_(r) = diff(N(r), r);
     B(r) = diff(N_(r), r);
     P(r) = D*(B'*B*r + v*((B')*(N_) + (N_)'*(B)) + (1/r)*(N_)'*(N_));
     Q(r) = N'*q*r;
     Ke = int(P, a, b);
     fe = int(Q, a, b);
end

