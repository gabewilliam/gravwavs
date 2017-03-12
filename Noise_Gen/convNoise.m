clear;
noise = csvread('noise.csv');

z = noise(:,2) + sqrt(-1) * noise(:,3);

T = ifft(z);

f1=noise(:,1);

f2=transpose(repelem(f1,2));

a1 = transpose(noise(:,2));
a2 = transpose(noise(:,3));

a3=[a1; a2];
a3=a3(:)';

n2(1,:) = f2;
n2(2,:) = a3;

i = length(f2) / 2 + 3;

H1 = n2(:, i:length(f2));
H2 = n2(:, 1:(i-1));

n3 = [H1, H2];

plot(imag(z));

%sig = csvread('sig.csv');