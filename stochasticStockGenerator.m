clc
clear variables
close all

figure;
hold on
for i = 1:5
[P, t] = generateStock(10);
plot(t,P);
end

function [P, t] = generateStock(T)
N = 10;
P(1) = 100;
t(1) = 0;

mu = 0.03;
sigma = 0.15;
deltat = 1/260;
dt = T/N;
for i = 1:N
    phi = normrnd(0,1)
    P(i + 1) = P(i) + mu*P(i)*deltat + sigma*P(i)*sqrt(deltat)*phi;
    t(i + 1) = i*deltat;
end

end