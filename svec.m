function [ h ] = svec( H )
%Convert symetric matrix nxn into vector of size n(n+1)/2, with isometry preservation
    n = size(H,1);
    h = zeros(n*(n+1)/2,1);
    index = 1;s2 = sqrt(2);
    for i = 1:n
        nElements = n-i;
        Scale = s2*eye(nElements+1);Scale(1,1) = 1;
        
        h(index:index + nElements) = Scale*H(i:end,i);

        index = index + nElements+1;
    end
end

