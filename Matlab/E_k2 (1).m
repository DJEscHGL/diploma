%%-----------------------------------------------------------------------%%
% This program loads the csv data file(s), and compute the longitudinal and transverse
% structure functions of velocities in all directions, 1-d velocity energy spectrum,
% as well as the longitudinal and transverse correlation functions of velocities.
%
% [E11,E22,E33,R11,R22,R33,DLL,DNN,DMM,k3,rl] = E_k2(fileName,i0,i_last,a,C2,R1,R2)
%
% inputs:
%       fileName: the csv data file name (or set of files) (including the path to folder)
%       i0: the first number of scv file in set (integer)
%       i_last: the last number of scv file in set (integer)
%       a: normalization parameter for the position of the curve with a slope of 5/3 (K41 teory)
%       on the energy spectrum.
%       C2: universal constant (normally C2 =2.12, suggested by Sreenivasan)
%       R1,R2: range bounded by row offsets R1 and R2 in csv data file loading
%
% outputs:
%       E11: 1-d longitudinal velocity energy spectrum of i component
%       E22: 1-d wall-normal velocity energy spectrum of i component
%       E33: 1-d transverse velocity energy spectrum of i component
%       k3:  is wavenumber vector in the i deriction
%       rl:  is spece-separation vector
%
% Last edit: Jul 20, 2020
% Dimitry Ivanov in Heat- and Masstransfer Institute NASB
%%-----------------------------------------------------------------------%%
function [E11,E22,E33,R11,R22,R33,DLL,DNN,DMM,k3,rl] = E_k2(fileName,i0,i_last,a,C2,R1,R2)

firstFile = strcat(fileName, num2str(i0), '.csv');
Arr = csvread(firstFile,R1,13,[R1,13,R2,13]);
LnArr = length(Arr);

% Get array for spline (4x time more)
Splinestep = 0.25;
xx = 1:Splinestep:LnArr;
x = 1:1:LnArr;
LnArr = length(xx);

% Use next highest power of 2 greater than or equal to row (length(velocities)) to calculate FFT.
nfft = 2^nextpow2(LnArr);

%the magnitude of the fft will be symmetric, such that the first (hp/2) points are unique,
%and the rest are symmetrically redundant. Calculate the number of unique points.
NumUniquePts = ceil(nfft/2);

% image count
imageCount = i_last-i0 + 1;

% initialize the arrays of kinetic energy function that will hold data
E11 = zeros(NumUniquePts,1);
E22 = zeros(NumUniquePts,1);
E33 = zeros(NumUniquePts,1);
R11 = zeros(LnArr,1);
R22 = zeros(LnArr,1);
R33 = zeros(LnArr,1);
DLL = zeros(LnArr,1);
DNN = zeros(LnArr,1);
DMM = zeros(LnArr,1);

for i=i0:i_last
    Fname = strcat(fileName, num2str(i), '.csv');
    U_prime = csvread(Fname,R1,9,[R1,9,R2,9]);
    V_prime = csvread(Fname,R1,10,[R1,10,R2,10]);
    W_prime = csvread(Fname,R1,11,[R1,11,R2,11]);
    U_prime = (spline(x,U_prime,xx))';
    V_prime = (spline(x,V_prime,xx))';
    W_prime = (spline(x,W_prime,xx))';

    U = csvread(Fname,R1,3,[R1,3,R2,3]);
    V = csvread(Fname,R1,4,[R1,4,R2,4]);
    W = csvread(Fname,R1,5,[R1,5,R2,5]);
    U = (spline(x,U,xx))';
    V = (spline(x,V,xx))';
    W = (spline(x,W,xx))';

    Fu1 = fft(U_prime,nfft);
    Fhu1 = Fu1(1:end/2);
    Eu1 = (abs(Fhu1).^2)*2;
    Fu2 = fft(V_prime,nfft);
    Fhu2 = Fu2(1:end/2);
    Eu2 = (abs(Fhu2).^2)*2;
    Fu3 = fft(W_prime,nfft);
    Fhu3 = Fu3(1:end/2);
    Eu3 = (abs(Fhu3).^2)*2;

    E11 = E11 + Eu1/imageCount;
    E22 = E22 + Eu2/imageCount;
    E33 = E33 + Eu3/imageCount;

    r11 = zeros(LnArr,1);
    r22 = zeros(LnArr,1);
    r33 = zeros(LnArr,1);
    dLL = zeros(LnArr,1);
    dNN = zeros(LnArr,1);
    dMM = zeros(LnArr,1);
    for j = 1 : LnArr
        for k = 1 : (LnArr-1)
            if (j + k) <= LnArr
                rU = U_prime(j,1) * U_prime(j+k,1);
                r11(k+1,1) = r11(k+1,1) + rU/(LnArr-k);
                rV = V_prime(j,1) * V_prime(j+k,1);
                r22(k+1,1) = r22(k+1,1) + rV/(LnArr-k);
                rW = W_prime(j,1) * W_prime(j+k,1);
                r33(k+1,1) = r33(k+1,1) + rW/(LnArr-k);
                dU = (U(j,1) - U(j+k,1))^2;
                dLL(k,1) = dLL(k,1) + dU/(LnArr-k);
                dV = (V(j,1) - V(j+k,1))^2;
                dNN(k,1) = dNN(k,1) + dV/(LnArr-k);
                dW = (W(j,1) - W(j+k,1))^2;
                dMM(k,1) = dMM(k,1) + dW/(LnArr-k);
            end
%             if (j - k) >= 1
%                 rU = U_prime(j,1) * U_prime(j-k,1);
%                 r11(k+1,1) = r11(k+1,1) + rU/(LnArr-k);
%                 rV = V_prime(j,1) * V_prime(j-k,1);
%                 r22(k+1,1) = r22(k+1,1) + rV/(LnArr-k);
%                 rW = W_prime(j,1) * W_prime(j-k,1);
%                 r33(k+1,1) = r33(k+1,1) + rW/(LnArr-k);
%             end
        end
    end
    r11(1,1) = U_prime(1,1)^2;
    r22(1,1) = V_prime(1,1)^2;
    r33(1,1) = W_prime(1,1)^2;
    R11  = R11 + r11/imageCount;
    R22  = R22 + r22/imageCount;
    R33  = R33 + r33/imageCount;
    DLL  = DLL + dLL/imageCount;
    DNN  = DNN + dNN/imageCount;
    DMM  = DMM + dMM/imageCount;
end
Norm1 = R11(1,1);
R11 = R11/Norm1;

Norm2 = R22(1,1);
R22 = R22/Norm2;

Norm3 = R33(1,1);
R33 = R33/Norm3;

Lz = abs(Arr(end,1)-Arr(1,1));
k3=2*pi*linspace(0,1,NumUniquePts)/Lz;

loglog(k3,E11);
grid on
xlabel('Wave numbers')
ylabel('Power Spectral Density')
title('Kinetic Energy Spectrum')
hold on
loglog(k3,E22);
hold on
loglog(k3,E33);

E_Kol = a * k3.^(-5/3);
loglog(k3, E_Kol)
legend('E11', 'E22', 'E33', '-5/3 slope line')
hold off

figure(2)
rl = (Lz/LnArr):(Lz/LnArr):Lz;
plot(rl,R11);
grid on
xlabel('r_{i}')
ylabel('R_{ij}(r_{3i})')
title('Two-point velocity correlation function')
hold on
plot(rl,R22);
hold on
plot(rl,R33);
legend('R11', 'R22', 'R33')
hold off

DLLpl = (3*DLL/4/C2).^(1.5);
for n = 1:length(DLLpl)
   DLLpl(n) = DLLpl(n)/rl(n);
end

figure(3)
plot(rl,DLLpl);
grid on
xlabel('r_{i}')
ylabel('(3*D_{LL}(r_{i})/4/C2)^{3/2}/r_{i}')
title('The second order transverse velocity structure function')

DNNpl = (3*DNN/4/C2).^(1.5);
for n = 1:length(DNNpl)
   DNNpl(n) = DNNpl(n)/rl(n);
end

figure(4)
plot(rl,DNNpl);
grid on
xlabel('r_{i}')
ylabel('(3*D_{NN}(r_{i})/4/C2)^{3/2}/r_{i}')
title('The second order wall-normal velocity structure function')

DMMpl = (DMM/C2).^(1.5);
for n = 1:length(DMMpl)
   DMMpl(n) = DMMpl(n)/rl(n);
end

figure(5)
plot(rl,DMMpl);
grid on
xlabel('r_{i}')
ylabel('(D_{MM}(r_{i})/C2)^{3/2}/r_{i}')
title('The second order longitudinal velocity structure function')
