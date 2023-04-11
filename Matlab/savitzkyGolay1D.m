function vOUT = savitzkyGolay1D(vIN,order,interval)
%%
% Info:
% This function does 1D Savitzky-Golay smoothing
%
% input:
% vIN      = vector, input vector for smoothing, 
%            if vIN is a matrix, function operates on first column
% order    = int, order of polynomial
% interval = int, size of the window (1.must be a odd number; 2. smaller than size(vIN)) 
% output:
% vOUT     = vector, smoothed data of vIN
% 
% Reference: 
% Abraham Savitzky and Marcel J. E. Golay, 'Smoothing and
% differentiation of data by simplified least sqaure procedures',
% Analytical Chemistry, Vol. 36, No. 8, page 1627-1639, July 1964
%
%%
%vIN = vIN';
vSize = size(vIN);
vSize = vSize(1,1);
vOUT = zeros(vSize,1);

% Compute the Savitzky- Golay coefficients matrix c
% order
n = order;
% moving window
t = -(interval-1)/2:1:(interval-1)/2;

for nn = 1:interval
    for mm = 1:n+1
        A(nn,mm) =t(nn)^(mm-1);
    end
end

AT = A';
AT_A = AT*A;

F = eye(interval);    
c = zeros(interval,interval);
for nn = 1:interval
    CC = AT_A\(AT*F(:,nn));
    for ss = 1:interval
        c(ss,nn) = 0;
        for mm = 1:n+1
            c(ss,nn) = c(ss,nn) + CC(mm)*t(ss)^(mm-1);
        end
    end
end

length = (interval-1)/2;
for p = 1:length
    vOUT(p) = 0;
    for nn = 1:interval
        vOUT(p) = vOUT(p) +c(p,nn)*vIN(nn);
    end
end
for p = length+1:vSize-length
    vOUT(p) = 0;
    for nn = 1:interval
        vOUT(p) = vOUT(p) +c(length+1,nn)*vIN(p+t(nn));
    end
end
row = length+2;
for p = vSize-length+1:vSize
    vOUT(p) =0;
    for nn = 1:interval
        vOUT(p) = vOUT(p) + c(row,nn)*vIN(vSize-(interval-nn));
    end
    row = row+1;
end

vOUT = vOUT';
