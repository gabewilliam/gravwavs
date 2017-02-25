noise = csvread('noise.csv');

z = noise(:,2) + sqrt(-1) * noise(:,3);

%f = transpose(linspace(-10000,10000,200000));

t = ifft(z);

plot(t);