%Loads the data from the text file into the workspace
fid=fopen('2DMonte.txt');

%Extracts the data from the file.
MaMb = fscanf(fid,'%g,%g\n',[2 Inf]);
MaMb = MaMb';
Ma = MaMb(:,1);
Mb = MaMb(:,2);

%Makes a figure featuring an (ma,mb) scatter and a histogram.
%Also now features a nice colourmap.
%It only really looks good when you maximise the plot window.
%This first block makes the scatter.
subplot(2,2,1)
samp = scatter(Ma,Mb,'.');
xlim([min(Ma) max(Ma)]);
ylim([min(Mb) max(Mb)]);
xlabel('m_a/kg')
ylabel('m_b/kg')
title('Data Point Plot')
hold on
legend([samp],{'(m_a,m_b) Samples'})
%daspect([1 1 1])
hold off;

%This makes the histogram.
subplot(2,2,[2,4])
hist3(MaMb,[20 20])
xlabel('m_a/kg')
ylabel('m_b/kg')
zlabel('Sample Density')
title('Posterior Histogram')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

%This extracts some data from the histogram. n is the number of counts per
%bin and C is the mass value of the bin centers. M is the number of counts
%of the highest bin in each column of n, and I is the index of each of said
%bins within its column of n. N is then the number of counts in the 
%absolute highest bin, and Imb is its index (ie. which column of n it lies 
%in). Ima is then the row in which it lies.
subplot(2,2,3)
[n,C] = hist3(MaMb,[20 20]);
[M,I] = max(n);
[N,Imb] = max(M);
Ima = I(Imb);

%Finds the (m1,m2) combination of the highest histogram bar.
maMax = C{1}(Ima);
mbMax = C{2}(Imb);
DispmaMax = sprintf('The best value of m_a is %e kg',maMax);
DispmbMax = sprintf('The best value of m_b is %e kg',mbMax);
disp(DispmaMax);
disp(DispmbMax);

%This plots the density of the histogram as a flat colourmap.
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(Ma),max(Ma),size(n,1)+1);
yb = linspace(min(Mb),max(Mb),size(n,1)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hsv)
hold on
xlabel('m_a/kg')
ylabel('m_b/kg')
title('Histogram Density Plot')

fclose(fid);

%Some of the code is lifted from the MATLAB documentation on hist3
