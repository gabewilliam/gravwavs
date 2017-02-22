%Loads the data from the text file into the workspace
fid=fopen('2DMonte.txt');

%Extracts the data from the file. The m1 values are stored in X. The m2
%values are stored in Y. This convention is followed throughout.
constants = fscanf(fid,'%g',1);
XY = fscanf(fid,'%g,%g\n',[2 Inf]);
XY = XY';
X = XY(:,1);
Y = XY(:,2);

%Makes a figure featuring an (m1,m2) scatter and a histogram.
%Also now features a nice colourmap.
%It only really looks good when you maximise the plot window.
%This first block makes the scatter.
subplot(2,2,1)
samp = scatter(X,Y,'.');
xlabel('m_1/kg')
ylabel('m_2/kg')
title('Data Point Plot')
hold on
legend([samp],{'(m_1,m_2) Samples'})
%daspect([1 1 1])
hold off;

%This makes the histogram.
subplot(2,2,[2,4])
hist3(XY,[20 20])
xlabel('m_1/kg')
ylabel('m_2/kg')
zlabel('Sample Density')
title('Posterior Histogram')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

%This extracts some data from the histogram.
subplot(2,2,3)
[n,C] = hist3(XY,[20 20]);
[M,I] = max(n);
[N,Iy] = max(M);
Ix = I(Iy);

%Finds the highest histogram bar and the coresponding (m1,m2) combination.
xMax = C{1}(Ix);
yMax = C{2}(Iy);
m1Max = sprintf('The best value of m_1 is %e kg',xMax);
m2Max = sprintf('The best value of m_2 is %e kg',yMax);
disp(m1Max);
disp(m2Max);

%This plots the density of the histogram as a flat colourmap.
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(X),max(X),size(n,1)+1);
yb = linspace(min(Y),max(Y),size(n,1)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hsv)
hold on
xlabel('m_1/kg')
ylabel('m_2/kg')
title('Histogram Density Plot')

fclose(fid);

%Some of the code is lifted from the MATLAB documentation on hist3
