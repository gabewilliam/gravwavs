noise = csvread('noise2.csv');
%noise(:,2) = noise(:,2) .* sqrt(3E8);
Z = noise(:,2) + noise(:,3) * sqrt(-1);
F = noise(:,1);

nt = ifft(Z, 'symmetric');

h = ((length(F) - 1) / 2) + 2 ;

sz2 = floor(length(F(:,1)));
df = F(h,1);
T = 1 / df;
dt = T / (sz2);

t = (0:dt:T-dt);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 900, 650]);

subplot(2,1,1);
loglog(F,abs(Z));
grid minor;
% Create xlabel
xlabel({'frequency f / Hz'});

% Create title
title({'Loglog Plot of Noise in Frequency Domain'});

% Create ylabel
ylabel({'|ASD| |N^~(f)| / Hz^-^0^.^5'});

xlim([8 8000]);

%ylim([0 3E-22]);


subplot(2,1,2);
plot(t,nt);
grid minor;
% Create xlabel
xlabel({'time t / s'});

% Create title
title({'Corresponding Time Domain Noise'});

% Create ylabel
ylabel({'Noise N(t)'});


