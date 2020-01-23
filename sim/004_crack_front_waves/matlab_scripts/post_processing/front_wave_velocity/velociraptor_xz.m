function [velo] = velociraptor(data,rate)

%on donne data qui est time en cs*t et space qui est en x on resort velo
%qui est dx/d(cs*t) ce qui donne v/cs en y et du cs*t en x c est comme
%dx/dt ca donne du x/t en y et du t en x...
%ff: Attention le temps est Cs*t/X (X taille du domaine)
clear velo;

time = data(:,2);
space = data(:,1);

lim = length(space);

range = round(lim*rate/2);

%bool=1;


for i=1:lim
m = 1;

    for j=max(1,i-range):min(lim,i+range)
    
        x(m) = space(j);
        t(m) = time(j);
        m= m+1;
        
    end
    
    p = polyfit(t,x,3);
    
    velo(i,1) = time(i)*time(i)*3*p(1) + time(i)*2*p(2) + p(3);
    %velo(i,2) = time(i);
    velo(i,2) = space(i);
       
    
%     for j=1:m-1
%     
%         xt(j) = t(j)*t(j)*t(j)*p(1) + t(j)*t(j)*p(2) + t(j)*p(3) + p(4);
%                 
%     end

    
    clear x;
    clear t;
    clear xt;
end


end
