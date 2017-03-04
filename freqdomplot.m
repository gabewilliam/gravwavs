load('30_26.csv');

freq = X30_26(1, :);
dat = X30_26(2, :);
sz = floor(length(dat)/2);

df = 0.05;
T = 1/df;
dt = T/sz;
% 
t = (0:dt:T-dt);

% figure;
% plot(abs(freq), dat);
% ylim([-1E-13 1E-13]);

amp = zeros(1, sz);
f = zeros(1, sz);
im = zeros(1, sz);
re = zeros(1, sz);
z = zeros(1, sz);

for d=1:sz
    f(d) = freq(d*2);
    re(d) = dat(d*2-1);
    im(d) = dat(d*2);
    amp(d) = sqrt(im(d)^2 + re(d)^2);
    
    sum(d) = im(d) + re(d);
    z(d) = re(d) + im(d)*1i;
end

figure;
subplot(1, 2, 1);
semilogx(f, re);
xlabel('Frequency, Hz');
ylabel('Re(h)');
grid on;
%ylim([-1E-9 1E-9]);
%xlim([20 5E2]);
subplot(1, 2, 2);
loglog(f, amp, '.');
xlabel('Frequency, Hz');
ylabel('|h|');
xlim([20 5E2]);
ylim([0 2E-14]);
grid on;

% filename='FreqDomPlots';
% fig = gcf; 
% u = fig.Units;
% fig.Units = 'inches';
% set(fig, 'PaperPositionMode', 'manual');
% set(fig, 'PaperPosition', [0 0 25 10]);
% fig_pos = fig.PaperPosition; 
% fig.PaperSize = [25 10];
% print(fig,[filename,'.pdf'],'-dpdf')

% IFFTAMP = ifft(amp);
% IFFTRE = ifft(re);
% IFFTIM = ifft(im);
% IFFTALL = ifft(dat);
% 
IFFTZ = ifft(z, 'symmetric');
figure;

h = (-6.103552878131E-6 - IFFTZ);

plot(t, h);
grid on;
xlabel('Time');
ylabel('Amplitude');
%ylim([-0.6 0.6]*10^(-16));
%xlim([19.75 19.98]);

% filename='TimeDomPlots';
% fig = gcf; 
% u = fig.Units;
% fig.Units = 'inches';
% set(fig, 'PaperPositionMode', 'manual');
% set(fig, 'PaperPosition', [0 0 25 10]);
% fig_pos = fig.PaperPosition; 
% fig.PaperSize = [25 10];
% print(fig,[filename,'.pdf'],'-dpdf')

% figure;
% plot(IFFTAMP, '.');
% figure;
% plot(IFFTRE);
% figure;
% plot(IFFTIM);
% figure;
% plot(IFFTALL, '.');

% figure
% subplot(1, 2, 1);
% semilogy(freq, dat, '.');
% grid on;
% xlabel('Frequency, MHz');
% ylabel('Amplitude');
% title('Semi-log plot');
% subplot(1, 2, 2);
% loglog(freq, dat, '.');
% grid on;
% xlabel('Frequency, MHz');
% ylabel('Amplitude');
% title('Log-log plot');
% % 
% filename='FreqDomPlots';
% fig = gcf; 
% u = fig.Units;
% fig.Units = 'inches';
% set(fig, 'PaperPositionMode', 'manual');
% set(fig, 'PaperPosition', [0 0 25 10]);
% fig_pos = fig.PaperPosition; 
% fig.PaperSize = [25 10];
% print(fig,[filename,'.pdf'],'-dpdf')