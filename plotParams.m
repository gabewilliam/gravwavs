function[best1, best2] = plotParams(param1Val, param1Name, param1Unit, param2Val, param2Name, param2Unit)

figure

%Makes a figure featuring a (param1, param2) scatter and a histogram.
%Also now features a nice colourmap.
%It only really looks good when you maximise the plot window.
%This first block makes the scatter.
subplot(2,2,1)
samp = scatter(param1Val,param2Val,'.');
%xlim([min(param1Val) max(param1Val)]);
%ylim([min(param2Val) max(param2Val)]);
xlabel(param1Unit)
ylabel(param2Unit)
title('Data Point Plot')
hold on
hold off;

%This makes the histogram.
subplot(2,2,2)
hist3([param1Val, param2Val],[20 20])
xlabel(param1Unit)
ylabel(param2Unit)
zlabel('Sample Density')
title('Posterior Histogram')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

%This extracts some data from the histogram.
% n is the number of counts per bin and
% C is the value of the bin centers.
% M is the number of counts of the highest bin in each column of n,
% and I is the index of each of said bins within its column of n.
% N is then the number of counts in the absolute highest bin,
% and I2 is its index (ie. which column of n it lies in).
% I1 is then the row in which it lies.
subplot(2,2,3)
[n,C] = hist3([param1Val, param2Val],[20 20]);
[M,I] = max(n);
[N,I2] = max(M);
I1 = I(I2);

%Finds the (param1Val,param2Val) combination of the highest histogram bar.
best1 = C{1}(I1);
best2 = C{2}(I2);
DispmaMax = sprintf('%e is the best value of %s', best1 , param1Unit);
DispmbMax = sprintf('%e is the best value of %s', best2 , param2Unit);
disp(DispmaMax);
disp(DispmbMax);

%This plots the density of the histogram as a flat colourmap.
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(param1Val),max(param1Val),size(n,1)+1);
yb = linspace(min(param2Val),max(param2Val),size(n,1)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(jet)
hold on
xlabel(param1Unit)
ylabel(param2Unit)
title('Histogram Density Plot')
hold off

%Makes a kernel density plot of the data
subplot(2,2,4)
ksdensity([param1Val, param2Val])
xlabel(param1Unit)
ylabel(param2Unit)
zlabel('p')

end