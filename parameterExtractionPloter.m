%Loads the data from the text file into the workspace
fid=fopen('RawData500k.txt');

%Extracts the data from the file.
dataFile = fscanf(fid,'%g,%g,%g\n',[3 Inf]);
dataFile = dataFile';
mRatio = dataFile(:,1);
mChirp = dataFile(:,2);
distance = dataFile(:,3);

%calculate credible intervals for variables

%size of credible interval (e.g. 0.90 = 90%)
intervalSize = 0.90;

%assumes same no of datapoints for all parameters
numDatapoints = size(dataFile,1);

%sort file of each parameter and find the values at interval boundaries
sortFile = sort(dataFile);

ratioLower = sortFile(int32(((1-intervalSize)/2)*numDatapoints),1);
ratioUpper = sortFile(int32(((1+intervalSize)/2)*numDatapoints),1);

chirpLower = sortFile(int32(((1-intervalSize)/2)*numDatapoints),2);
chirpUpper = sortFile(int32(((1+intervalSize)/2)*numDatapoints),2);

distLower = sortFile(int32(((1-intervalSize)/2)*numDatapoints),3);
distUpper = sortFile(int32(((1+intervalSize)/2)*numDatapoints),3);

%Plots various pairings of variables
% plot(mRatio)
% figure
% plot(mChirp)
% figure
% plot(distance)
[mR,cM] = plotParams(mRatio, 'mass ratio', 'mass ratio', mChirp, 'chirp mass', 'chirp mass / kg');
% plotParams(mChirp, 'chirp mass', 'chirp mass / kg', distance, 'distance', 'distance / m');
% plotParams( distance, 'distance', 'distance / m', mRatio, 'mass ratio', 'mass ratio');
ma = cM*(1+mR)^0.2*mR^-0.6
mb = cM*(1+mR)^0.2*mR^0.4

for n=1:numDatapoints
    currentChirp = mChirp(n);
    currentRatio = mRatio(n);
    massList(n,1) = currentChirp*(1+currentRatio)^0.2*currentRatio^-0.6;
    massList(n,2) = currentChirp*(1+currentRatio)^0.2*currentRatio^0.4;
end
massListSort = sort(massList);

maLower = massListSort(int32(((1-intervalSize)/2)*numDatapoints),1);
maUpper = massListSort(int32(((1+intervalSize)/2)*numDatapoints),1);

mbLower = massListSort(int32(((1-intervalSize)/2)*numDatapoints),2);
mbUpper = massListSort(int32(((1+intervalSize)/2)*numDatapoints),2);

fclose(fid);