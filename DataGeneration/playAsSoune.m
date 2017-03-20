clear;

n = csvread('NoisySignal.csv');

f = n(:,1);
z = n(:,2) + n(:,3) .* sqrt(-1);

nt = ifft(z, 'symmetric');

h = floor(length(f)/2) + 2;

sz2 = floor(length(f));
df = f(h,1);
%df = 0.01;
T = 1 / df;
dt = T / (sz2);

t = (0:dt:T-dt);
nt = nt ./ dt;
soundsc(nt);
plot(t,nt);