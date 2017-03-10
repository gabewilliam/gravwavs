load('30_26.csv');
load('Noisy_Sig.csv');

freq = X30_26(1, :);
dat = X30_26(2, :);

datn = Noisy_Sig(2, :);

sz = floor(length(dat)/2);

im = zeros(1, sz);
re = zeros(1, sz);
imn = zeros(1, sz);
ren = zeros(1, sz);
z = zeros(1, sz);
zn = zeros(1, sz);

for d=1:sz

    re(d) = dat(d*2-1);
    im(d) = dat(d*2);

    z(d) = re(d) + im(d)*1i;
    
    ren(d) = datn(d*2-1);
    imn(d) = datn(d*2);

    zn(d) = ren(d) + imn(d)*1i;
end

df = 0.05;
T = 1/df;
dt = T/sz;
t = (0:dt:T-dt);

b = 1*10^21;

cleanData = b*ifft(z, 'symmetric');

duttyData = b*ifft(zn, 'symmetric');

figure;
subplot(2,2,1);
plot(t, cleanData);
grid on;
xlabel('Time (s)');
ylabel('Strain amplitue, h(t) (10^{-21})');
subplot(2,2,2);
plot(t, duttyData);
grid on;
xlabel('Time (s)');
ylabel('Strain amplitue, h(t) (10^{-21})');
subplot(2,2,3);
plot(t, cleanData)
xlim([10.0 10.5]);
ylim([-3.5 3.5]);
grid on;
xlabel('Time (s)');
ylabel('Strain amplitue, h(t) (10^{-21})');
subplot(2,2,4);
plot(t, duttyData);
xlim([10.0 10.5]);
ylim([-3.5 3.5]);
grid on;
xlabel('Time (s)');
ylabel('Strain amplitue, h(t) (10^{-21})');