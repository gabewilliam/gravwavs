nsig = cleanData + nt;

plot(t, nsig);

tdub = repelem(t,2);


a1 = (nsig(1,:));
a2 = zeros(1,length(nsig));

a3=[a1; a2];
a3=a3(:)';

signal(1,:) = tdub(1,:);
signal(2,:) = a3;

csvwrite('SpecialNoiseSignalForEverybody.csv',signal);