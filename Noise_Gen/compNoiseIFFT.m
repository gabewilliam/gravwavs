noise = csvread('noise.csv');

z = noise(:,2) + sqrt(-1) * noise(:,3);

f = transpose(linspace(-5000,5000,200000));

figure;

subplot(2,1,1);
plot(f, real(z));
subplot(2,1,2);
plot(f, imag(z));