Xnoise = load('noise.csv');
Xnoise = transpose(Xnoise);
freq2 = Xnoise(1, :);
dat2 = Xnoise(2, :);

sz2 = floor(length(dat2)/2);

im2 = zeros(1, sz2);
re2 = zeros(1, sz2);
z2 = zeros(1, sz2);

for d=1:(sz2)

    re2(d) = dat2(d*2-1);
    im2(d) = dat2(d*2);

    z2(d) = re2(d) + im2(d)*1i;
end


df2 = 0.005;
T2 = 1 / df2;
dt2 = T2 / (sz2);

t2 = (0:dt2:T2-dt2);

figure;

nt = ifft(z2, 'symmetric');

plot(t2, nt);