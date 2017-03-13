    
for m = 1:1:3
    
    if(m==1)
        fid = fopen('../data2/temp_30_26.csv');
    end
    
    if(m==2)
        fid = fopen('../data2/conout.dat');  
    end
    
    if(m==3)
        fid = fopen('../data2/output.dat');
    end

    P = [,];
    InDepVar = [,];
    DepVar = [,];

    N = 0;

    line = fgetl(fid);

    dvmax = 0;

    while ischar(line)

      p = sscanf(line,'%f,', [1, 2]);

      %Because kg is clearly the most natural units for black holes
      p = p/1.9891e30;
      p = int64(p);

      P = [P;p];

      line = fgetl(fid);
      InDepVar = [InDepVar; sscanf(line,'%f,',[1,inf])];

      line = fgetl(fid);
      DepVar = [DepVar; sscanf(line,'%f,',[1,inf])];

      if(max(DepVar)>dvmax)
          dvmax = max(DepVar);
      end

      line = fgetl(fid);

      N = N+1;
    end

    fclose(fid);

    figure(m)

    P1range = max(P(:,1)) - min(P(:,1));
    P2range = max(P(:,2)) - min(P(:,2));

    i=0;
    n=1;
    pn=1;

    while n<(N+1)

        if(n>1)
            if(P(n,1)>P(n-1,1))
                i = i+1;
                pn = pn+i; 
            end
        end

        subplot(4, 2, pn)
        plot(InDepVar(n,:), DepVar(n,:))

        if(dvmax ~= 0)
            %ylim([-dvmax, dvmax])
        end

        %ylim([-1e-26, 1e-26])

        n = n+1;
        pn = pn+1;

    end
end