function [B] = convfilt(A,k)
%Function to compute convolution as a filter

A = double(A);
szA = size(A);
szk = size(k);

pRowCol = (length(k)-1)/2; %number of rows and columns to be padded with
A = padarray(A,[pRowCol pRowCol],'symmetric');
B = zeros(szA);

% for fully overlap only
for x=1:size(B,1)   %rows
    for y=1:size(B,2)  %columns
        xlim = x + szk(1) - 1;
        ylim = y + szk(2) - 1;
        B(x,y) = sum(sum(A(x:xlim,y:ylim).*k));
    end
end

end