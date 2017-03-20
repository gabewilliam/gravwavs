n = csvread('Noise.csv');

f = n(:,1);
z = n(:,2) + n(:,3) .* sqrt(-1);

nt = ifft(z, 'symmetric');

h = floor(length(f)/2) + 2;

sz2 = floor(length(f));
df = f(h,1);
T = 1 / df;
dt = T / (sz2);
nt = nt .* 1/dt;
t = (0:dt:(T-dt));
figure;
subplot(1,2,1);
plot(t,nt);
xlabel('$t$/s','Interpreter','latex');
ylabel('$n(t)$','Interpreter','latex');
grid on;
subplot(1,2,2);
loglog(f,abs(z),'r');
xlim([0 f(length(f))]);
xlabel('$f$/Hz','Interpreter','latex');
ylabel('$\left|\tilde{n}(f)\right|$/Hz$^{-1/2}$','Interpreter','latex');
grid on;

fig = gcf;
set(gcf, 'PaperPosition', [0 0 20 8]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 8]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'succ', 'pdf') %Save figure