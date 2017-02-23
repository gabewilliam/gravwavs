noise = csvread('noise.csv');

Z = noise(:,1).*exp(noise(:,2)*sqrt(-1));

plot(Z);

nt = 1/0.001 .* ifft(Z, 340446);
figure;
plot(real(nt));
figure;
plot(imag(nt));
%ns(:,1) = s(:,1);

%ns(:,2) = s(:,2) + real(nt(:,1));
%figure;
%plot(ns(:,1), ns(:,2));

%csvwrite('bigassnoisyasscoolsignal.csv', ns);