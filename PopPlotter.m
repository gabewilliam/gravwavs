%Reads in the full and mass limited stellar populations, pre evolution
fullPop = csvread('fullpop.csv');
massLimitPop = csvread('masslimitpop.csv');

%Extracts the columns of m1 and m2 values from the files
m1f = fullPop(:,1);
m2f = fullPop(:,2);
m1ml = massLimitPop(:,1);
m2ml = massLimitPop(:,2);

%Plots the systems m2 vs m1
scatter(m1f,m2f,'.b');
hold on
scatter(m1ml,m2ml,'.r');
xlabel('m1/M_{Solar}')
ylabel('m2/M_{Solar}')
legend('Full Population','20M_{Solar} \leq m \leq 300M_{Solar}')
