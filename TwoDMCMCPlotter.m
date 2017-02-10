%Loads the data from the text file into the workspace
fid=fopen('2DMonte.txt');

test=fscanf(fid,'%g,%g,%g',3)
XY=fscanf(fid,'%g,%g\n',[2 Inf])

s=test(1);
a=test(2);
b=test(3);

%Creates circles of radius 1sigma and 2sigma, for ilustrative purposes
x12 = linspace(-s+a,s+a,10000);
y1 = b + sqrt(s^2 - (x12-a).^2);
y2 = b - sqrt(s^2 - (x12-a).^2);
x34 = linspace(-2*s+a,2*s+a,10000);
y3 = b + sqrt(4*s^2 - (x34-a).^2);
y4 = b - sqrt(4*s^2 - (x34-a).^2);

%Makes a figure featuring an (x,y) scatter and a histogram.
%Also now features a nice colourmap, and rings showing 1 and 2 sigma.
%It only really looks good when you maximise the plot window.
subplot(2,2,1)
samp = scatter(XY(1,:),XY(2,:),'.')
xlabel('x')
ylabel('y')
hold on
c1 = plot(x12,y1,'r')
plot(x12,y2,'r')
c2 = plot(x34,y3,'g')
plot(x34,y4,'g')
legend([samp,c1,c2],{'(x,y) Samples','1\sigma Contour','2\sigma Contour'})
%daspect([1 1 1])
hold off;

subplot(2,2,[2,4])
hist3(XY',[20 20])
xlabel('x')
ylabel('y')
zlabel('Sample Density')
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')

subplot(2,2,3)
n = hist3(XY',[20 20]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(XY(1,:)),max(XY(1,:)),size(n,1)+1);
yb = linspace(min(XY(2,:)),max(XY(2,:)),size(n,1)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hsv)
hold on
plot(x12,y1,'k')
plot(x12,y2,'k')
plot(x34,y3,'k')
plot(x34,y4,'k')
xlabel('x')
ylabel('y')

fclose(fid)
%Some of the code is lifted from the MATLAB documentation on hist3
