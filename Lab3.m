% Gregory Leonberg

%% Problem 1

syms s t;

wo = 4;
Hnum = (s + 3);
Hden = (s.^2 + s + 1.25);
H = Hnum ./ Hden;

xt = sin(wo*t);

%% Part A -- Pole, 40dB, h(t) calculation

time = linspace(0, 10, 1001);

% find poles and roots
poles = roots([1, 1, 1.25])
ht = feval(symengine, 'partfrac', H);
ht = ilaplace(ht)

% 40dB = 1 percent of max
fortyval = 0.01*max(subs(ht, time));
fortylocs = find(subs(ht, time) >= fortyval);
fortyloc = fortylocs(1, end);
fortytime = fortyloc/100
% 40dB tau is 8.58

% graph filter response
figure;
plot(time, subs(ht, time));
xlim([0, 10]);
ylim([-2, 2]);
set(gca, 'xtick', 0:1:10);
set(gca, 'ytick', -2:1:2);
title('impulse response, h(t)');
xlabel('t (sec)');
grid on;

%% Part B -- Find H(jwo), |H(jwo)|, arg H(jwo)

freqResp = subs(H, 's', 4)
freMag = 20*log(freqResp)
phase = angle(subs(H, 's', 4*j))

yst = freqResp .* sin(wo.*time + phase);
plot(time, subs(yst, 't', time))

%% Part C -- Exact response for complex sinusoidal input

% define sinusoid input
sinusoid = exp(j.*wo.*t);

% laplace sinusoid and multiply by input
yc = H.*laplace(sinusoid)

% partial fraction expand
yc = feval(symengine, 'partfrac', yc)

% inverse laplace to get yc
yc = ilaplace(yc);

% use reside and imag to find imaginary parts
[r,p,k] = residue([1, 3], vertcat([4*i], poles))
yt = imag(subs(yc, 't', time));

% plot yt, yst, x
figure;
plot(time, yt, time, subs(yst, 't', time), ...
    time, subs(xt, 't', time));

xlim([0, 10]);
ylim([-1.3, 1.3]);
set(gca, 'xtick', 0:1:10);
set(gca, 'ytick', -1:0.5:1);
grid on;
title('x(t), y(t), yst(t)');
legend('location', 'sw', 'exact', 'steady', 'input');
hold on

%% Part D -- Phase Delay and Max

tph = phase./wo

% declare range to check for max
temp = linspace(8, 9, 101);
% create vector of values in that time period
xtemp = subs(xt, 't', time);
xtemp = xtemp(1, 800:950);
ytemp = yst(1, 800:950);
% find maxes for each function in that time period
yt1 = max(ytemp);
xt2 = max(xtemp);

% find index of each of those maxes
t1 = time(find(max(ytemp) == ytemp) + 800);
t2 = time(find(max(xtemp) == xtemp) + 800);

t_est = t2 - t1

plot(t1, yt1, 'g*', t2, xt2, 'g*');
legend('location', 'sw', 'exact', 'steady', 'input', 'max');

%% Part E -- Transient Plot

ytr = yt-yst;
figure;
plot(time, ytr);
xlim([0, 10]);
ylim([-1.3, 1.3]);
set(gca, 'xtick', 0:1:10);
set(gca, 'ytick', -1:0.5:1);
grid on;
title('transients, t40 = 9.21 sec');
xlabel('t (sec)');

% The decrease IS consistent with our calculated value of 8.58 sec

%% Problem 2









