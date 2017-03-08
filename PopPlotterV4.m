%Reads in the full and mass limited stellar populations, pre evolution
fullPop = csvread('pop.csv');


%Extracts the columns of m1 and m2 values from the files
m1 = fullPop(:,1);
m2 = fullPop(:,2);
tm = fullPop(:,3);


%Plots the systems m2 vs m1
%hold on
%scatter(m1,m2,'.r');
xlabel('m1/M_{Solar}')
ylabel('m2/M_{Solar}')

%ksdensity([m1 m2])

%scatter((m1+m2),(m2./m1),'.b')
%ksdensity([(m1+m2) (m2./m1)])

%hist(tm,20)

%scatter(m1.^(0.6).*m2.^(0.6).*(m1+m2).^(-0.2) tm)
