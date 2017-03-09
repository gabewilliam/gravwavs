load('30_26.csv');

freq = X30_26(1, :);
dat = X30_26(2, :);

dat = dat(:, 1:end-2);
freq = freq(:, 1:end-2);

sz = floor(length(dat)/2);

im = zeros(1, sz);
re = zeros(1, sz);
z = zeros(1, sz);

for d=1:(sz)

    re(d) = dat(d*2-1);
    im(d) = dat(d*2);

    z(d) = re(d) + im(d)*1i;
end

df = 0.01;
T = 1/df;
dt = T/(sz);

T = 200;

t = (0:dt:(T)-dt);

b = 1;

cleanData = b*ifft(z, 'symmetric');
cleanData = [cleanData, zeros(1 , 819201)];
plot(t,cleanData);