set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

clear;

n = csvread('Noise.csv');
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

subplot(1,2,1);

plot(t,nt);
xlabel('$t$/s','interpreter','latex');
ylabel('$n(t)$','interpreter','latex');
xlim([0.0 1.0]);
grid on;
subplot(1,2,2);

loglog(f,abs(z),'r');
xlabel('$f$/Hz','interpreter','latex');
ylabel('$\tilde{n}(f)$','interpreter','latex');
grid on;
fig = gcf;
set(gcf, 'PaperPosition', [0 0 20 8]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [20 8]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'whit', 'pdf') %Save figure
