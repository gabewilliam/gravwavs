%Reads in the full and mass limited stellar populations, pre evolution
fullPop = csvread('pop.csv');


%Extracts the columns of m1 and m2 values from the files
m1 = fullPop(:,1);
m2 = fullPop(:,2);
a = fullPop(:,3);


%Plots the systems m2 vs m1
%scatter(m1,m2,'.b');
%hold on
colormap(parula(3));
scatter(m1,m2,5, log(a));
xlabel('M_1/M_{Solar}')
ylabel('M_2/M_{Solar}')
c = colorbar;
c.Label.String = 'Log_{10}(Binary separation/AU)';
axis([0 350 0 350])
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MySavedFile','-dpdf')

%ksdensity([m1 m2])