load 'fDomComp.dat'

data = fDomComp;

figure(2)
  
for n = 1:2:161;
    
    t = data(n,:);
    a = data(n+1,:);
    
    plot(t,a);
    hold on
end

figure(3)

m = 1;
  
for n = 1:2:161;

    t = data(n,:);
    a = data(n+1,:);
    
    subplot(9,9,m);
  
    plot(t,a);
    
    m = m+1;
end