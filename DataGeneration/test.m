clear;

n = csvread('NoisySignal.csv');
L =length(n(:,1));
H =(length(n(:,1))+1)/2;
n=transpose(n);

f = [n(1,H:L) n(1,1:H)];
f=transpose(f);
re = [n(2,H:L) n(2,1:H)];
re=transpose(re);
im = [n(3,H:L) n(3,1:H)];
im=transpose(im);
z = re + im .* sqrt(-1);

nt = ifft(z,'symmetric');

sz2 = floor(length(f));
df = f(2,1);
T = 1 / df;
dt = T / (sz2);
nt = nt .* (sqrt(sz2));
t = (0:dt:(T-dt));
nt = flipud(nt);
plot(t,nt);
soundsc(nt);

tdub = repelem(t,2);
nt = nt .* (sqrt(sz2));
nq = transpose(nt);
a1 = (nq(1,:));
a2 = zeros(1,length(nq));

a3=[a1; a2];
a3=a3(:)';

signal(1,:) = tdub(1,:);
signal(2,:) = a3;
%csvwrite('NoisySignals/43_24_zdhp_01.csv',signal);