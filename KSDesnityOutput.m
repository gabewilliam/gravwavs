set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%Reads in the full and mass limited stellar populations, pre evolution
fullPop = csvread('pop.csv');


%Extracts the columns of m1 and m2 values from the files
m1 = fullPop(:,1);
m2 = fullPop(:,2);

mChirp = ((m1.*m2).^0.6)./((m1+m2).^0.2);
mRatio = m2./m1;


[P12, m12] = ksdensity([m1 m2]);
[PChirpRat, mChirpRat] = ksdensity([mChirp mRatio]);

fid1 = fopen('m12KernelGrid.csv','w');
fid11 = fopen('m12KernelProb.txt','w');
fid2 = fopen('ChirpRatKernelGrid.csv','w');
fid22 = fopen('ChirpRatKernelProb.txt','w');
fid3 = fopen('m12KernelEven.csv','w');

dat3 = [m12(:,1)'; m12(:,2)'; P12'];

fprintf(fid3,'%.15f,%.15f,%.15f\n',dat3);

m1Pts = [];
m2Pts = [];

mChirpPts = [];
mRatPts = [];

j = 1;
k = 1;

while j < 900
   
        m1Pts(k) = m12(j,1);
        mChirpPts(k) = mChirpRat(j,1);
        
        j = j + 30;
        k = k + 1;
        
end

for i = 1:30
   
    m2Pts(i) = m12(i,2);
    mRatPts(i) = mChirpRat(i,2);
    
end

fprintf(fid1,'%.15f,%.15f\n',[m1Pts;m2Pts]);
fprintf(fid2,'%.15f,%.15f\n',[mChirpPts;mRatPts]);
fprintf(fid11,'%.15f\n',P12');
fprintf(fid22,'%.15f\n',PChirpRat');

fclose(fid1);
fclose(fid11);
fclose(fid2);
fclose(fid22);
fclose(fid3);

