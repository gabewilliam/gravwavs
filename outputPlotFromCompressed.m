% get length of signal - could just be hardcoded
load('NonAdaptiveSignal.dat')
N=length(NonAdaptiveSignal)

%%

% imports from compresstest.gz to an array X that contains all the data
% I have no idea how this works

% Build a stream chain that reads, decompresses and decodes the file into lines
fileStr = javaObject('java.io.FileInputStream', 'compresstest.gz');
inflatedStr = javaObject('java.util.zip.GZIPInputStream', fileStr);
charStr = javaObject('java.io.InputStreamReader', inflatedStr);
lines = javaObject('java.io.BufferedReader', charStr);

% If you know the size in advance you can preallocate the arrays instead
% of just stating the types to allow vcat to succeed
X = zeros(20,N);
curL = lines.readLine();
for j=1:20
%while ischar(curL) % on EOF, readLine returns null, which becomes [] (type double)
    % Parse a single line from the file
    curX = sscanf(char(curL),'%f,', [1 Inf]);
    % Append new line results 
    for iCol=1:length(X)
        X(j,iCol)= curX(iCol);
    end
    curL = lines.readLine();
    java.lang.Runtime.getRuntime.gc;
    %j=j+1;
end
lines.close(); % Don't forget this or the file will remain open!
