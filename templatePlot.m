figure(1)

fid = fopen('../data/templates.dat');

n=1;

while(fscanf(fid,'%f,%f\n',[1, 2]))

line = fgetl(fid);
t = sscanf(line,'%f,', [1, Inf]);

line = fgetl(fid);
at = sscanf(line,'%f,', [1, Inf]);

subplot(9, 9, n)
plot(t, at)

n = n+1;

end

fclose(fid);

figure(2)

fid = fopen('templatesF.dat');

n=1;

while(fscanf(fid,'%f,%f\n',[1, 2]))

line = fgetl(fid);
t = sscanf(line,'%f,', [1, Inf]);

line = fgetl(fid);
at = sscanf(line,'%f,', [1, Inf]);

subplot(9, 9, n)
plot(t, at)

n = n+1;

end

fclose(fid);