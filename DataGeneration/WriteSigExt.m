clear;

n = csvread('NoisySignal.csv');

f = n(:,1);
z = n(:,2) + n(:,3) .* sqrt(-1);

nt = ifft(z,'symmetric');

h = floor(length(f)/2) + 2;

sz2 = floor(length(f));
df = f(h,1);
T = 1 / df;
dt = T / (sz2);

t = (0:dt:(T-dt));

tdub = repelem(t,2);
nt = nt .* (sqrt(sz2));
nq = transpose(nt);
a1 = (nq(1,:));
a2 = zeros(1,length(nq));

a3=[a1; a2];
a3=a3(:)';

signal(1,:) = tdub(1,:);
signal(2,:) = a3;
csvwrite('NoisySignals/30_30_whit_01.csv',signal);
plot(t,nt);
soundsc(nt);
