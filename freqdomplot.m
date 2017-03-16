load('36.2_29.1.csv');
load('23_13.csv');
load('14.2_7.5.csv');
load('sig_30_39.csv');

freq = sig_30_39(1, :);
dat = sig_30_39(2, :);

freq1 = X36_2_29_1(1, :);
dat1 = X36_2_29_1(2, :);

freq2 = X23_13(1, :);
dat2 = X23_13(2, :);

freq3 = X14_2_7_5(1, :);
dat3 = X14_2_7_5(2, :);

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

amp1 = zeros(1, sz);
f1 = zeros(1, sz);
im1 = zeros(1, sz);
re1 = zeros(1, sz);
z1 = zeros(1, sz);
amp2 = zeros(1, sz);
f2 = zeros(1, sz);
im2 = zeros(1, sz);
re2 = zeros(1, sz);
z2 = zeros(1, sz);
amp3 = zeros(1, sz);
f3 = zeros(1, sz);
im3 = zeros(1, sz);
re3 = zeros(1, sz);
z3 = zeros(1, sz);

for d=1:sz
    
    f(d) = freq(d*2);
    re(d) = dat(d*2-1);
    im(d) = dat(d*2);
    amp(d) = sqrt(im(d)^2 + re(d)^2);
    
    sum(d) = im(d) + re(d);
    z(d) = re(d) + im(d)*1i;
    
%     f1(d) = freq1(d*2);
%     re1(d) = dat1(d*2-1);
%     im1(d) = dat1(d*2);
%     amp1(d) = sqrt(im1(d)^2 + re1(d)^2);
%     
%     sum1(d) = im1(d) + re1(d);
%     z1(d) = re1(d) + im1(d)*1i;
%     
%     f2(d) = freq2(d*2);
%     re2(d) = dat2(d*2-1);
%     im2(d) = dat2(d*2);
%     amp2(d) = sqrt(im2(d)^2 + re2(d)^2);
%     
%     sum2(d) = im2(d) + re2(d);
%     z2(d) = re2(d) + im2(d)*1i;
%     
%     f3(d) = freq3(d*2);
%     re3(d) = dat3(d*2-1);
%     im3(d) = dat3(d*2);
%     amp3(d) = sqrt(im3(d)^2 + re3(d)^2);
%     
%     sum3(d) = im3(d) + re3(d);
%     z3(d) = re3(d) + im3(d)*1i;
end

% figure;
% subplot(1, 2, 1);
% semilogx(f, re);
% xlabel('Frequency, Hz');
% ylabel('Re(h)');
% grid on;
% %ylim([-1E-9 1E-9]);
% %xlim([20 5E2]);
% subplot(1, 2, 2);
% loglog(f, amp, '.');
% xlabel('Frequency, Hz');
% ylabel('|h|');
% %xlim([20 5E2]);
% %ylim([0 2E-14]);
% grid on;

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

% filename='MultiTimeDomPlots';
% fig = gcf; 
% u = fig.Units;
c = 3*10^8;
b = 1*10^21;
% 
IFFTZ = b*ifft(z, 'symmetric')/sqrt(c);
% IFFTZ2 = b*ifft(z2, 'symmetric')/sqrt(c);
% IFFTZ3 = b*ifft(z3, 'symmetric')/sqrt(c);
figure;
% subplot(3,1,3);
plot(t, IFFTZ);
% legend('M = 23 + 13, d = 1000', 'location', 'northwest');
% grid on;
% xlabel('Time (s)');
% ylabel('Strain amplitude (10^{-21})');
% %ylim([-5 5]*10^(-21));
% xlim([9.5 11]);
% subplot(3,1,2);
% plot(t, IFFTZ2);
% legend('M = 14.2 + 7.5, d = 440', 'location', 'northwest');
% grid on;
% xlabel('Time (s)');
% ylabel('Strain amplitude (10^{-21})');
% %ylim([-5 5]*10^(-21));
% xlim([9.5 11]);
% subplot(3,1,1);
% plot(t, IFFTZ3);
% legend('M = 36.2 + 29.1, d = 420', 'location', 'northwest');
% grid on;
% xlabel('Time (s)');
% ylabel('Strain amplitude (10^{-21})');
% %ylim([-5 5]*10^(-21));
% xlim([9.5 11]);

% fig.Units = 'inches';
% set(fig, 'PaperPositionMode', 'manual');
% set(fig, 'PaperPosition', [0 0 25 30]);
% fig_pos = fig.PaperPosition; 
% fig.PaperSize = [25 30];
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