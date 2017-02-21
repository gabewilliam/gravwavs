re = csvread('real.csv');
im = csvread('imag.csv');

z = re + sqrt(-1) * im;

f = transpose(linspace(-10000,10000,200000));

t = ifft(z);

plot(t);