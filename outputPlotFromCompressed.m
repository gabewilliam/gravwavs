load('NonAdaptiveSignal.dat')
N=length(NonAdaptiveSignal)

%%

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

%%
data = testout4;

%%
n=1;
%nmax=size(testout)
figure
m = 1; 
for n = 1:2:20;

    t = X(n,:);
    a = X(n+1,:);
    
    subplot(2,5,m);

    scatter(t,a,'.');
    
    m = m+1;
end
%%
figure
m=1;
for n = 51:2:100;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 101:2:150;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 151:2:200;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 201:2:250;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 251:2:300;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 301:2:350;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 351:2:400;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 401:2:450;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end

figure
m=1;
for n = 451:2:462;

    t = testout(n,:);
    a = testout(n+1,:);
    
    subplot(5,5,m);
  
    scatter(t,a,'.');
    
    m = m+1;
end
%%
for n = 1:2:5;

    t = testout3(n,:);
    a = testout3(n+1,:);
    
    subplot(3,4,m);
  
    plot(a);
    
    m = m+1;
end

for n = 1:2:5;

    t = testout4(n,:);
    a = testout4(n+1,:);
    
    subplot(3,4,m);
  
    plot(a);
    
    m = m+1;
end

t=testout5(1,:);
a=testout5(2,:);

subplot(3,4,m)
plot(a)
m=m+1;