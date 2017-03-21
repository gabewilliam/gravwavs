clear;

n = csvread('NoisySignal.csv');

f = n(:,1);
re = n(:,2);
im = n(:,3);

fdub = repelem(f,2)';

a1 = re;
a2 = im;

z = a1 + a2*1i;

IFFT = ifft(z, 'symmetric');
%plot(IFFT);

a3=[a1; a2];
a3=a3(:)';

freqsig(1,:) = fdub(1,:);
freqsig(2,:) = a3;

plot(f,re);

fig(1,:) = [freqsig(1,(length(freqsig)/2:length(freqsig))) freqsig(1,1:(length(freqsig)/2)+1)];
fig(2,:) = [freqsig(2,(length(freqsig)/2:length(freqsig))) freqsig(2,1:(length(freqsig)/2)+1)];

csvwrite('NoisySignals/31_46_zdhp_freq_01.csv',fig);