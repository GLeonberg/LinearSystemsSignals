% Gregory Leonberg

%% Problem 1

% Clear variables and command window
clear all;
clc;

% Declare constants and unit step function
t = linspace(-1,1,1001);
u = @(t) (t>=0);
e = [0.1, 0.05, 0.001, 0.002];

%% Rectangular Pulse

% Declare equations 
eq1 = (1/e(2))*(u(t+0.5*e(2))-u(t-0.5*e(2)));
eq2 = (1/e(1))*(u(t+0.5*e(1))-u(t-0.5*e(1)));

% Plot equations
figure % create new figure
plot(t, eq1, 'b', t, eq2, 'r:') % plot equations
axis([-1 1 0 25]) % set axis limits
title('Rectangular Pulse') % add title
legend('\epsilon = 0.05', '\epsilon = 0.10') % add legend

%% Gaussian

% Declare equations 
eq1 = (1/sqrt(2*pi*e(3)))*exp(-(t.^2)/(2*e(3)));
eq2 = (1/sqrt(2*pi*e(4)))*exp(-(t.^2)/(2*e(4)));

% Plot equations
figure % create new figure
plot(t, eq1, 'b', t, eq2, 'r:') % plot equations
axis([-1 1 0 15]) % set axis limits
title('Gaussian') % add title
legend('\epsilon = 0.001', '\epsilon = 0.002') % add legend

%% Lorentzian

% Declare equations 
eq1 = (1./pi).*(e(2)./(e(2).^2 + t.^2));
eq2 = (1./pi).*(e(1)./(e(1).^2 + t.^2));

% Plot equations
figure % create new figure
plot(t, eq1, 'b', t, eq2, 'r:') % plot equations
title('Lorentzian') % add title
legend('\epsilon = 0.05', '\epsilon = 0.1') % add legend

%% Sinc

% Declare equations 
eq1 = sinc(t ./ e(1)) ./ (pi.*e(1));
eq2 = sinc(t ./ e(2)) ./ (pi.*e(2));

% Plot equations
figure % create new figure
plot(t, eq1, 'b', t, eq2, 'r:') % plot equations
title('Sinc') % add title
legend('\epsilon = 0.1', '\epsilon = 0.05') % add legend

%% Problem II

clear all;
clc;

%% Part A - Verifying Equality of Equations
% y(t) = -4v(t) + 3x(t)
% y'(t) = -4v'(t) + 3x'(t)
% v'(t) = -2v(t) + x(t)

% y'(t) = -4(-2v(t)+x(t)) + 3x'(t)
% y'(t) = 8v(t) - 4x(t) + 3x'(t)
% y'(t) = 8v(t) - 6x(t) + 2x(t) + 3x'(t)
% y'(t) = -2(-4v(t) + 3x(t)) + 2x(t) + 3x'(t)
% y'(t) = -2y(t) + 2x(t) + 3x'(t)

% y'(t) + 2y(t) = 2x(t) + 3x'(t)

%% Part B - Transfer Function Inverse Laplace

% y'(t) + 2y(t) = 3x'(t) + 2x(t)
% Using form 
% y'(t) + ay(t) = B0*x'(t) + B1*x(t) => H(s) = (B0*s + B1) / (s + a)
% a = 2, B0 = 3, B1 = 2

syms s; % declare s symbolic
H = (3*s + 2)/(s + 2);
h = ilaplace(H)

%% Part C - Initial Values

% y(t) = -4v(t) + 3x(t)
% 4v(t) = 3x(t) - y(t)
% 4v(t) = 3e^-t*u(t) - y(t)

% How is v(0-) = v(0+) related to y(0-) where y(0-) = -1?

% 4v(0-) = 3e^-(0-)*u(0-) - y(0-)
% 4v(0-) = -y(0-)
% 4v(0-) = 1

% v(0-) = 1/4
% v(0+) = 1/4

% 4v(0+) = 3e^-(0+)*u(0+) - y(0+)
% 1 = 3e^-(0+)*u(0+) - y(0+)
% 1 = 3 * 1 * 1 - y(0+)

% 2 = y(0+)

%% Part D - Laplace, Partial Fractions, Inverse Laplace
syms Y y0 t s;
x = exp(-1*t);       % input signal, exponential, b~=a 
X = laplace(x);      % X(s)

Y1 = solve(s*Y+1 + 2*Y == 3*s*X + 2*X, Y);  % y(t) - Laplace method
y1 = ilaplace(Y1)

%% Part E - Alternative Method using DSolve

syms t y(t);
x = exp(-t);
dy = diff(y,t); 
dx = diff(x,t);

y2 = dsolve(dy + 2*y(t) == 3*dx + 2*x, y(0) == 2)
[y2 - y1] % show that dsolve and laplace method get same answer

%% Part F - State Equation 2

syms t v(t) V;

x = exp(-1*t);       % input signal, exponential, b~=a 
x0 = subs(x,t,0);    % x(0+), needed for dsolve
X = laplace(x);      % X(s)

V1 = solve(s*V-0.25 + 2*V == X, V);  % y(t) - Laplace method
v1 = ilaplace(V1)

x = exp(-t);
dv = diff(v,t); 
dx = diff(x,t);

v2 = dsolve(dv == -2*v + x, v(0) == 0.25)
[v2 - v1] % show that dsolve and laplace method get same answer

%% Part G - Plotting y(t) and v(t)
clear t v2 y2;
t = linspace(0,5,1001);
y2 = 3*exp(-2*t) - exp(-t);
v2 = exp(-t) - (3*exp(-2*t))/4;
figure;
plot(t, y2, 'b', t, v2, 'r:') % plot equations
title('Exact Output, y(0-) = -1') % add title
legend('y(t)', 'v(t)') % add legend

%% Part H - Using lsim
% Lots of inspiration from week 2 source code

t = linspace(0,5,201)';
x = exp(-1*t); 

num = [3,2]; den = [1,2];    % numerator and denominator coeffs of H(s)
H = tf(num,den);             % class(H) is tf, can be passed to lsim for yzs
[a,b,c,d] = tf2ss(num,den);  % convert to state-space controller form
                             %   dv/dt = a*v + b*x 
                             %       y = c*v + d*x
                             % in our example, a=-1, b=1, c=-2, d=3
S = ss(a,b,c,d);             % class(S) is ss, and can be passed to lsim

s = tf('s');                 % class(s) is tf
Hv = 1/(s+1);                % class(Hv) is tf
                             % transfer function of the state equation
                             % convert Hv to state-space form
Sv = ss(Hv);                 % class(Sv) is ss and can be passed to lsim

a = 2; b0 = 3; b1 = 2; 
y0 = -1; v0 = 0.25;  % initial conditions, y(0-) and v(0)

y   = lsim(S,x,t,v0);        % full solution, must use v0 as initial state
v   = lsim(Sv,x,t,v0);       % state v(t)

figure; 
plot(t,y,'b-', t,v,'r--')
xlabel('\itt'); 
title('lsim outputs,  y(0-) = -1'), 
legend(' {\ity}({\itt})', ' {\itv}({\itt})')

%% Part I - Discrete Time Implementation
% Lots of inspiration from week 2 source code

clear
T = 0.1;                    % discretization time step
a = 2; B0 = 3; B1 = 2;       % continuous-time parameters
y0 = -1; v0 = 0.25;  % initial conditions, y(0-), v(0-)

a1 = -exp(-a*T);         % discrete-time parameters, ZOH-method
b0 = B0;
b1 = B1*(1-exp(-a*T))/a - B0;

% [a1, b0, b1]

tmax = 5;   
tn = 0:T:tmax;           % time range, sampled in steps of T
N = length(tn);

x = exp(-1*tn);          % sampled input, x(tn)

w = y0; v = 0;           % initialize w,v
                         % calculate full solution y(t)
for n=0:N-1,             
    y(n+1) = -a1*w + b0*x(n+1) + b1*v;       % (n+1) is MATLAB index
    w = y(n+1);
    v = x(n+1);
end

w = 0; v = 0;            % initialize w,v
                         % calculate zero-state solution yzs(t)
for n=0:N-1,             
    yzs(n+1) = -a1*w + b0*x(n+1) + b1*v;     % (n+1) is MATLAB index
    w = yzs(n+1);
    v = x(n+1);
end

clear v
a1 = -exp(-a*T);          % discretization parameters for state v(t)
b1 = (1 - exp(-a*T))/a;
w = v0; u = 0;            % initialize w,u (using letter u instead of v)
                          % calculate state variable v(t), with v(0-)=v0
for n=0:N-1,             
    v(n+1) = -a1*w + b1*u;     % (n+1) is MATLAB index
    w = v(n+1);
    u = x(n+1);
end

figure; plot(tn,y,'b-', tn,v,'r--', tn,yzs, 'k:')
xlabel('\itt_n'); 
title('discrete-time outputs, y(0) = -1, T = 0.1'), 
legend(' {\ity}({\itt_n})', ' {\itv}({\itt_n})', 'exact')

%% Part J - Repeat C-I with y(0-) = 0

%% Part C - Initial Values

% y(t) = -4v(t) + 3x(t)
% 4v(t) = 3x(t) - y(t)
% 4v(t) = 3e^-t*u(t) - y(t)

% How is v(0-) = v(0+) related to y(0-) where y(0-) = 0?

% 4v(0-) = 3e^-(0-)*u(0-) - y(0-)
% 4v(0-) = -y(0-)
% 4v(0-) = 0

% v(0-) = 0
% v(0+) = 0

% 4v(0+) = 3e^-(0+)*u(0+) - y(0+)
% 0 = 3e^-(0+)*u(0+) - y(0+)
% 0 = 3 * 1 * 1 - y(0+)

% 3 = y(0+)

%% Part D - Laplace, Partial Fractions, Inverse Laplace
syms Y y0 t s;
x = exp(-1*t);       % input signal, exponential, b~=a 
X = laplace(x);      % X(s)

Y1 = solve(s*Y+0 + 2*Y == 3*s*X + 2*X, Y);  % y(t) - Laplace method
y1 = ilaplace(Y1)

%% Part E - Alternative Method using DSolve

syms t y(t);
x = exp(-t);
dy = diff(y,t); 
dx = diff(x,t);

y2 = dsolve(dy + 2*y(t) == 3*dx + 2*x, y(0) == 3)
[y2 - y1] % show that dsolve and laplace method get same answer

%% Part F - State Equation 2

syms t v(t) V s;

x = exp(-1*t);       % input signal, exponential, b~=a 
x0 = subs(x,t,0);    % x(0+), needed for dsolve
X = laplace(x);      % X(s)

V1 = solve(s*V-0 + 2*V == X, V);  % y(t) - Laplace method
v1 = ilaplace(V1)

x = exp(-t);
dv = diff(v,t); 
dx = diff(x,t);

v2 = dsolve(dv == -2*v + x, v(0) == 0)
[v2 - v1] % show that dsolve and laplace method get same answer

%% Part G - Plotting y(t) and v(t)

clear t v2 y2;
t = linspace(0,5,1001);
y2 = 4*exp(-2*t) - exp(-t);
v2 = exp(-t) - exp(-2*t);
figure;
plot(t, y2, 'b', t, v2, 'r:') % plot equations
title('Exact Output, y(0-) = 0') % add title
legend('y(t)', 'v(t)') % add legend

%% Part H - Using lsim
% Lots of inspiration from week 2 source code

t = linspace(0,5,201)';
x = exp(-1*t); 

num = [3,2]; den = [1,2];    % numerator and denominator coeffs of H(s)
H = tf(num,den);             % class(H) is tf, can be passed to lsim for yzs
[a,b,c,d] = tf2ss(num,den);  % convert to state-space controller form
                             %   dv/dt = a*v + b*x 
                             %       y = c*v + d*x
                             % in our example, a=-1, b=1, c=-2, d=3
S = ss(a,b,c,d);             % class(S) is ss, and can be passed to lsim

s = tf('s');                 % class(s) is tf
Hv = 1/(s+1);                % class(Hv) is tf
                             % transfer function of the state equation
                             % convert Hv to state-space form
Sv = ss(Hv);                 % class(Sv) is ss and can be passed to lsim

a = 2; b0 = 3; b1 = 2; 
y0 = 0; v0 = 0;  % initial conditions, y(0-) and v(0)

y   = lsim(S,x,t,v0);        % full solution, must use v0 as initial state
v   = lsim(Sv,x,t,v0);       % state v(t)

figure; 
plot(t,y,'b-', t,v,'r--')
xlabel('\itt'); 
title('lsim outputs,  y(0-) = 0'), 
legend(' {\ity}({\itt})', ' {\itv}({\itt})')

%% Part I - Discrete Time Implementation
% Lots of inspiration from week 2 source code

clear
T = 0.1;                    % discretization time step
a = 2; B0 = 3; B1 = 2;       % continuous-time parameters
y0 = 0; v0 = 0;  % initial conditions, y(0-), v(0-)

a1 = -exp(-a*T);         % discrete-time parameters, ZOH-method
b0 = B0;
b1 = B1*(1-exp(-a*T))/a - B0;

% [a1, b0, b1]

tmax = 5;   
tn = 0:T:tmax;           % time range, sampled in steps of T
N = length(tn);

x = exp(-1*tn);          % sampled input, x(tn)

w = y0; v = 0;           % initialize w,v
                         % calculate full solution y(t)
for n=0:N-1,             
    y(n+1) = -a1*w + b0*x(n+1) + b1*v;       % (n+1) is MATLAB index
    w = y(n+1);
    v = x(n+1);
end

w = 0; v = 0;            % initialize w,v
                         % calculate zero-state solution yzs(t)
for n=0:N-1,             
    yzs(n+1) = -a1*w + b0*x(n+1) + b1*v;     % (n+1) is MATLAB index
    w = yzs(n+1);
    v = x(n+1);
end

clear v
a1 = -exp(-a*T);          % discretization parameters for state v(t)
b1 = (1 - exp(-a*T))/a;
w = v0; u = 0;            % initialize w,u (using letter u instead of v)
                          % calculate state variable v(t), with v(0-)=v0
for n=0:N-1,             
    v(n+1) = -a1*w + b1*u;     % (n+1) is MATLAB index
    w = v(n+1);
    u = x(n+1);
end

figure; plot(tn,y,'b-', tn,v,'r--', tn,yzs, 'k:')
xlabel('\itt_n'); 
title('discrete-time outputs, y(0) = 0, T = 0.1'), 
legend(' {\ity}({\itt_n})', ' {\itv}({\itt_n})', 'exact')