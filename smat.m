function [ H ] = smat( h )
%Reverse of operation svec 
    
    %n^2 + n - 2*length(h) = 0

    n = (sqrt(1+8*length(h))-1)/2;
    H = zeros(n,n);
    
    index = 1;s2 = 1/sqrt(2);
    for i = 1:n
        nElements = n-i;
        Scale = s2*eye(nElements+1);Scale(1,1) = 1;
        
        H(i:end,i) = Scale*h(index:index + nElements);

        H(i,i+1:end) = H(i+1:end,i);
        
        index = index + nElements+1;
    end
end

