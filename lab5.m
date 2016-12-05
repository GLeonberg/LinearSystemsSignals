%% LSS Lab 5
% Gregory Leonberg

%% Problem 1 -- AC/DC half-wave rectifier/converter

close all; clear; clc;

wo = 2*pi;
u = @(t) (t >= 0);
x = @(t) sin(wo.*t);
f = @(t) x(t).*u(t);
T = 2*pi/wo;

%% Part A -- M-term Approximations

P = @(s) wo.*(1+exp(-1.*s.*T./2)) ./ (s.^2 + wo.^2); 

t = linspace(0, 3, 1001);

for M = [10, 30]
    for tau = [T, 5*T, 10*T]

        H = @(s) 1./(1 + s.*tau);

        c0 = 1/pi;
        c1 = 1/(4*j);
        ck = @(k) (1/T).*P(j*k*wo);

        fm = c0 + 2*real(c1*exp(j*1*wo*t));
        ymsteady = c0.*H(0) + 2*real(c1*H(j*1*wo)*exp(j*1*wo*t));

        for k = 2:1:M
            fm = fm + 2*real(ck(k)*exp(j*k*wo*t));
            ymsteady = ymsteady + 2*real(ck(k)*H(j*k*wo)*exp(j*k*wo*t));
        end;

        figure;
        plot(t, fm, t, ymsteady);
        xlim([0, 3]);
        ylim([-0.2, 1.2]);
        set(gca, 'xtick', 0:0.5:3);
        set(gca, 'ytick', -0.2:0.2:1.2);
        xlabel('t/T');
        legend('half-wave f(t)', 'steady output');
        title( horzcat('half-wave rectifier, M = ', num2str(M), ', \tau = ', num2str(tau), 'T'));
        grid on;
        
    end;
end;

%% Part B -- M-term Approximations Continued

t = linspace(0, 24, 1001);
A = @(ta) ((-1.*ta.^-1).*P(-1.*ta.^-1)) ./ (exp(T./ta) - 1);

for M = 30
    for tau = [5*T, 10*T]
        
        fm = c0 + 2*real(c1*exp(j*1*wo*t));
        ymsteady = c0.*H(0) + 2*real(c1*H(j*1*wo)*exp(j*1*wo*t));   
        ym = A(tau)*exp(-t/tau) + c0*H(0) + 2*real(c1*H(j*1*wo)*exp(j*1*wo*t));
        
        for k = 2:1:M
            fm = fm + 2*real(ck(k)*exp(j*k*wo*t));
            ymsteady = ymsteady + 2*real(ck(k)*H(j*k*wo)*exp(j*k*wo*t));
            ym = ym + 2*real(ck(k)*H(j*k*wo)*exp(j*k*wo*t));
        end;
        
        figure;
        plot(t, fm, t, ymsteady, t, ym);
        xlim([0, 24]);
        ylim([-0.2, 1.2]);
        set(gca, 'xtick', 0:2:24);
        set(gca, 'ytick', -0.2:0.2:1.2);
        xlabel('t/T');
        legend('half-wave f(t)', 'steady output', 'exact output');
        title( horzcat('half-wave rectifier, M = ', num2str(M), ', \tau = ', num2str(tau), 'T'));
        grid on;
        
    end;
end;

%% Part C -- Error norms using lsim

for M = 30
    for tau = [5*T, 10*T]
        tran = tf([0, 0, 1], [0, tau, 1]);
        norm(lsim(tran, f(t), t)' - ym)    
    end;
end;

%% Problem 2

close all;
clear;
clc;

omega0 = 0.1*pi;
omega1 = 0.2*pi;
omega2 = 0.3*pi;

s = @(n) sin(omega0*n);
v = @(n) sin(omega1*n) + sin(omega2*n);
x = @(n) s(n) + v(n);

%% Part A - Plotting x(n) and s(n)

n = 0:1:199;
xn = x(n);
sn = s(n);

figure;
plot(n, xn, n, sn);
title('x(n) and s(n)');
xlabel('time samples, n');
set(gca, 'xtick', 0:50:200);
set(gca, 'ytick', -3:1:3);
ylim([-3, 3])
grid on

%% Part B - Using Filter

omegaA = 0.15*pi;
omegaB = 0.25*pi;

filt = fir1(50, [omegaA, omegaB]);
freqz(filt, x(n));
















