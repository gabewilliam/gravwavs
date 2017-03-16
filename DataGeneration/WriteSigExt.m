n = csvread('NoisySignal.csv');

f = n(:,1);
z = n(:,2) + n(:,3) .* sqrt(-1);

nt = ifft(z, 'symmetric');

h = ceil(length(f+1)/2) + 1;

sz2 = floor(length(f));
df = f(h,1);
T = 1 / df;
dt = T / (sz2);

t = (0:dt:T-dt);

tdub = repelem(t,2);

nq = transpose(nt);
a1 = (nq(1,:));
a2 = zeros(1,length(nq));

a3=[a1; a2];
a3=a3(:)';

signal(1,:) = tdub(1,:);
signal(2,:) = a3;

csvwrite('NoisySignals/noggit.csv',signal);