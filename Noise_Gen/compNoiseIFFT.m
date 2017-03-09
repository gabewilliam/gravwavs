noise = csvread('noise2.csv');
%noise(:,2) = noise(:,2) .* sqrt(3E8);
Z = noise(:,2) + noise(:,3) * sqrt(-1);
F = noise(:,1);

figure;


loglog(F,(imag(Z)).^2);
hold on;
loglog(F,real(Z).^2);
xlim([0 8000]);
ylim([-2E-21 2E-21]);
nt = ifft(Z, 'symmetric');

h = ((length(F) - 1) / 2) + 2 ;

sz2 = floor(length(F(:,1)));
df = F(h,1);
T = 1 / df;
dt = T / (sz2);

t = (0:dt:T-dt);

figure;

plot(t, nt);

