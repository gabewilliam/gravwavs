%Loads the data from the text file into the workspace
fid=fopen('RawData.txt');

%Extracts the data from the file.
dataFile = fscanf(fid,'%g,%g,%g\n',[3 Inf]);
dataFile = dataFile';
mRatio = dataFile(:,1);
mChirp = dataFile(:,2);
distance = dataFile(:,3);

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
fclose(fid);