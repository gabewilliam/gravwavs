nsig = cleanData + nt;

plot(t, nsig);

signal(1,:) = t;
signal(2,:) = nsig;

csvwrite('SpecialNoiseSignalForEverybody.csv',signal);