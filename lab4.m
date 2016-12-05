%% Lab 4

% Gregory Leonberg

%% Part A -- PID Response to Step Input

%% original case

% sample code given to calculate transfer function
a = 2; kp = 10; ki = 5; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 10, ki = 5, kd = 3')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% ki = 10

% sample code given to calculate transfer function
a = 2; kp = 10; ki = 10; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 10, ki = 10, kd = 3')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% kp = 20

% sample code given to calculate transfer function
a = 2; kp = 20; ki = 5; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 20, ki = 5, kd = 3')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% kd = 6

% sample code given to calculate transfer function
a = 2; kp = 10; ki = 5; kd = 6;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 10, ki = 5, kd = 6')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% kp = 20, ki = 10

% sample code given to calculate transfer function
a = 2; kp = 20; ki = 10; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 20, ki = 10, kd = 3')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% kp = 20, kd = 6

% sample code given to calculate transfer function
a = 2; kp = 20; ki = 5; kd = 6;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 20, ki = 5, kd = 6')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% ki = 10, kd = 6

% sample code given to calculate transfer function
a = 2; kp = 10; ki = 10; kd = 6;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 10, ki = 10, kd = 6')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% kp = 20, ki = 10, kd = 6

% sample code given to calculate transfer function
a = 2; kp = 20; ki = 10; kd = 6;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% determine poles
p = pzmap(H);

% determine 40dB time constant
tau = -log(100) / real(p(3))

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% calculate step response
vals = lsim(H, u(t), t);

% plot step response
figure;
plot(t, vals)
title('step response, kp = 20, ki = 10, kd = 6')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 1.4])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.2:1.4)
grid on;

%% Part B - Tracking

clear all;

a = 2; kp = 10; ki = 5; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% creating reference input signals
r1 = @(t) u(t) + u(t-10);
r2 = @(t) 0.1 .* t .* u(t);
r3 = @(t) atan(0.1 .* t) .* u(t);
r4 = @(t) 0.04.*t.*u(t) ... 
          - 0.04.*t.*u(t-10) ...
          + (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-10) ...
          - (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-14) ...
          + (0.8 + 0.2.*(t-14)).*u(t-14) ...
          - (0.8 + 0.2.*(t-14)).*u(t-20);

% generating output with lsim for each input signal
y1 = lsim(H, r1(t), t);
y2 = lsim(H, r2(t), t);     
y3 = lsim(H, r3(t), t);      
y4 = lsim(H, r4(t), t);    

% generate errors for each input signal
ye1 = lsim(Herr, r1(t), t);
ye2 = lsim(Herr, r2(t), t);
ye3 = lsim(Herr, r3(t), t);
ye4 = lsim(Herr, r4(t), t);

% plot outputs wrt inputs, and errors for each input signal

% input signal r1
figure;
plot(t, y1, 'b-', t, r1(t), 'r:');
title('tracking step changes')
legend('y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)
grid on;

figure;
plot(t, ye1);
title('tracking error')
xlabel('t')
ylabel('e(t)')
xlim([0, 20])
ylim([-0.2, 1.2])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -0.2:0.2:1.2)
grid on; 
      
% input signal r2
figure;
plot(t, y2, 'b-', t, r2(t), 'r:');
title('ramp tracking')
legend('y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)
grid on;

figure;
plot(t, ye2);
title('tracking error')
xlabel('t')
ylabel('e(t)')
xlim([0, 20])
ylim([-0.1, 0.1])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -0.1:0.05:0.1)
grid on;   
  
% input signal r3
figure;
plot(t, y3, 'b-', t, r3(t), 'r:');
title('ramp tracking with correct angle')
legend('y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)
grid on;

figure;
plot(t, ye3);
title('tracking error')
xlabel('t')
ylabel('e(t)')
xlim([0, 20])
ylim([-0.1, 0.1])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -0.1:0.05:0.1)
grid on;       
      
% input signal r4
figure;
plot(t, y4, 'b-', t, r4(t), 'r:');
title('accelerating case')
legend('y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)
grid on;

figure;
plot(t, ye4);
title('tracking error')
xlabel('t')
ylabel('e(t)')
xlim([0, 20])
ylim([-0.1, 0.1])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -0.1:0.05:0.1)
grid on;       

% ki = 0 plot

% recalculate transfer functions for ki = 0
a = 2; kp = 10; ki = 0; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));

% reapply lsim for new transfer functions
y2 = lsim(H, r2(t), t);  
ye2 = lsim(Herr, r2(t), t);

% plot new transfer functions for input r2
figure;
plot(t, y2, 'b-', t, r2(t), 'r:');
title('ramp tracking, ki = 0')
legend('y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)
grid on;

figure;
plot(t, ye2);
title('tracking error, ki = 0')
xlabel('t')
ylabel('e(t)')
xlim([0, 20])
ylim([-0.1, 0.1])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -0.1:0.05:0.1)
grid on;   

%% Part C - Torque

clear all;

s = tf('s');

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% creating reference input signals
r1 = @(t) u(t) + u(t-10);
r2 = @(t) 0.1 .* t .* u(t);
r3 = @(t) atan(0.1 .* t) .* u(t);
r4 = @(t) 0.04.*t.*u(t) ... 
          - 0.04.*t.*u(t-10) ...
          + (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-10) ...
          - (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-14) ...
          + (0.8 + 0.2.*(t-14)).*u(t-14) ...
          - (0.8 + 0.2.*(t-14)).*u(t-20);
      
% sample code given to calculate transfer function
a = 2; kp = 10; ki = 5; kd = 3;

% declare given tau
tau = 0.05;

% recalculate Gc with new tau
Gc = kp + ki/s + (kd*s)/(tau*s+1);

% recalculate H with new tau
H = ((kd + tau*kp)*s*s + (kp + tau*ki)*s + ki) / ...
    (tau*s*s*s*s + (tau*a+1)*s*s*s + (a+kd+tau*kp)*s*s ...
    + (kp + tau*ki)*s + ki);

% f(t) = [r(t) - y(t)] * gc(t)

% generate each y(t) vector
yt1 = lsim(H, r1(t), t);
yt2 = lsim(H, r2(t), t);
yt3 = lsim(H, r3(t), t);
yt4 = lsim(H, r4(t), t);

% generate f(t) vectors
ft1 = lsim(Gc, (r1(t)'-yt1), t);
ft2 = lsim(Gc, (r2(t)'-yt2), t);
ft3 = lsim(Gc, (r3(t)'-yt3), t);
ft4 = lsim(Gc, (r4(t)'-yt4), t);

% plot torque step
figure;
plot(t, ft1)
title('torque f(t) -- step changes')
xlabel('t')
ylabel('f(t)')
xlim([0,20])
ylim([-20, 80])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', -20:20:80)

% plot torque ramp tracking
figure;
plot(t, ft2)
title('torque f(t) -- ramp tracking')
xlabel('t')
ylabel('f(t)')
xlim([0,20])
ylim([0, 0.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.1:0.5)
grid on;

% plot torque ramp angle
figure;
plot(t, ft3)
title('torque f(t) -- ramp with correct angle')
xlabel('t')
ylabel('f(t)')
xlim([0,20])
ylim([0, 0.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.1:0.5)
grid on;

% plot torque accel
figure;
plot(t, ft4)
title('torque f(t) -- accelerating case')
xlabel('t')
ylabel('f(t)')
xlim([0,20])
ylim([0,0.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.1:0.5)
grid on;

% set tau to zero
tau = 0;


%% Part D - Disturbance

clear all;

% declare x region for plotting
t = linspace(0, 20, 1001);

% unit step function
u = @(t) double(t>=0);

% creating reference input signals
r1 = @(t) u(t) + u(t-10);
r2 = @(t) 0.1 .* t .* u(t);
r3 = @(t) atan(0.1 .* t) .* u(t);
r4 = @(t) 0.04.*t.*u(t) ... 
          - 0.04.*t.*u(t-10) ...
          + (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-10) ...
          - (-2 + 0.69.*t - 0.07.*t.*t + 0.0025.*t.*t.*t).*u(t-14) ...
          + (0.8 + 0.2.*(t-14)).*u(t-14) ...
          - (0.8 + 0.2.*(t-14)).*u(t-20);
      
% sample code given to calculate transfer function
a = 2; kp = 10; ki = 5; kd = 3;
s = tf('s');
G = 1/(s*(s+a));
Gc = kp + ki/s + kd*s;
H = minreal(Gc*G/(1+Gc*G));
Herr = minreal(1/(1+Gc*G));
Hdist = minreal(G/(1+Gc*G));     

% generate noise
fdist1 = 2*(u(t-4)-u(t-6)); % wind gust
seed=2016; rng(seed); % initialize random number generator
fdist2 = randn(size(t)); % zero-mean, unit-variance noise

% case 1 subcase 1 calculations
ydist11 = lsim(Hdist,fdist1,t);
y11 = lsim(H,r1(t),t);
ytot11 = y11 + ydist11;

% case 1 subcase 2 calculations
ydist12 = lsim(Hdist,fdist1,t);
y12 = lsim(H,r2(t),t);
ytot12 = y12 + ydist12;

% case 1 subcase 3 calculations
ydist13 = lsim(Hdist,fdist1,t);
y13 = lsim(H,r3(t),t);
ytot13 = y13 + ydist13;

% case 1 subcase 4 calculations
ydist14 = lsim(Hdist,fdist1,t);
y14 = lsim(H,r4(t),t);
ytot14 = y14 + ydist14;

% case 2 subcase 1 calculations
ydist21 = lsim(Hdist,fdist2,t);
y21 = lsim(H,r1(t),t);
ytot21 = y21 + ydist21;

% case 2 subcase 2 calculations
ydist22 = lsim(Hdist,fdist2,t);
y22 = lsim(H,r2(t),t);
ytot22 = y22 + ydist22;

% case 2 subcase 3 calculations
ydist23 = lsim(Hdist,fdist2,t);
y23 = lsim(H,r3(t),t);
ytot23 = y23 + ydist23;

% case 2 subcase 4 calculations
ydist24 = lsim(Hdist,fdist2,t);
y24 = lsim(H,r4(t),t);
ytot24 = y24 + ydist24;

% plot case 1 subcase 1
figure;
plot(t, ytot11, 'b-', t, r1(t), 'r--', t, fdist1, 'r:');
title('wind gust -- step changes')
legend('location', 'se', 'y(t)', 'r(t)', 'fdist')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 1 subcase 2
figure;
plot(t, ytot12, 'b-', t, r2(t), 'r--', t, fdist1, 'r:');
title('wind gust -- ramp')
legend('location', 'se', 'y(t)', 'r(t)', 'fdist')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 1 subcase 3
figure;
plot(t, ytot13, 'b-', t, r3(t), 'r--', t, fdist1, 'r:');
title('wind gust -- ramp with correct angle')
legend('location', 'ne', 'y(t)', 'r(t)', 'fdist')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 1 subcase 4
figure;
plot(t, ytot14, 'b-', t, r4(t), 'r--', t, fdist1, 'r:');
title('wind gust -- accelerating')
legend('location', 'se', 'y(t)', 'r(t)', 'fdist')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 2 subcase 1
figure;
plot(t, ytot21, 'b-', t, r1(t), 'r--');
title('wind noise -- step changes')
legend('location', 'nw', 'y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 2 subcase 2
figure;
plot(t, ytot22, 'b-', t, r2(t), 'r--');
title('wind noise -- ramp')
legend('location', 'nw', 'y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 2 subcase 3
figure;
plot(t, ytot23, 'b-', t, r3(t), 'r--');
title('wind noise -- ramp with correct angle')
legend('location', 'nw', 'y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)

% plot case 2 subcase 4
figure;
plot(t, ytot24, 'b-', t, r4(t), 'r--');
title('wind noise -- accelerating')
legend('location', 'nw', 'y(t)', 'r(t)')
xlabel('t')
ylabel('y(t)')
xlim([0, 20])
ylim([0, 2.5])
set(gca, 'xtick', 0:2:20)
set(gca, 'ytick', 0:0.5:2.5)