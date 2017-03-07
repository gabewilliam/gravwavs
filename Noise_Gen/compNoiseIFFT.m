sig = load('30_26.csv');

noise = csvread('noise.csv');
%noise(:,2) = noise(:,2) .* sqrt(3E8);
n = transpose(noise);

nsig = sig + n;

plot(nsig(1,:),nsig(2,:));

csvwrite('Noisy_Sig.csv',nsig);

%nt = IFFT(nsig);

%z = noise(:,2) + sqrt(-1) * noise(:,3);

%t = size(z);
%f = transpose(linspace(-5000,5000,t(1)));