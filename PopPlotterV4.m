%Reads in the full and mass limited stellar populations, pre evolution
fullPop = csvread('pop.csv');


%Extracts the columns of m1 and m2 values from the files
m1 = fullPop(:,1);
m2 = fullPop(:,2);


%Plots the systems m2 vs m1
%scatter(m1f,m2f,'.b');
%hold on
scatter(m1,m2,'.r');
xlabel('m1/M_{Solar}')
ylabel('m2/M_{Solar}')

ksdensity([m1 m2])