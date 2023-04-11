%%-----------------------------------------------------------------------%%
% This program loads the csv data file(s), and compute 1-d velocity energy spectrum
%
% [E11,E22,E33,E,R11,R22,R33,DLL,DNN,DMM,k3,rl] = E_k(fileName,i0,i_last,a,C2,R1,R2)
%
% inputs:
%       fileName: the csv data file name (or set of files) (including the path to folder)
%       Lij: is integral length scale in appropriate direction
%
% outputs:
%       Eii: 1-d velocity energy spectrum of i component
%       Eii: 1-d premultiplied streamwise velosity energy spectrum
%       k3:  is wavenumber vector in the i deriction
%       
%
% Last edit: Oct 14, 2020
% Dimitry Ivanov in Heat- and Masstransfer Institute NASB
%%-----------------------------------------------------------------------%%
function [E11,E22,E33,PreE11,k3,lambda] = PremultSpec(U,V,W,Lz,flag)

LnArr = length(U(1,:));
i_last = length(U(:,1));

Umeanx = mean(U);
Umeany = mean(V);
Umeanz = mean(W);

for j = 1:size(U,1)
    UPrimex(j,:) = U(j,:) - Umeanx(1,:);
    UPrimey(j,:) = V(j,:) - Umeany(1,:);
    UPrimez(j,:) = W(j,:) - Umeanz(1,:);
end
clear j

% Use next highest power of 2 greater than or equal to row (length(velocities)) to calculate FFT.
nfft = 2^nextpow2(LnArr); 
NumUniquePts = ceil(nfft/2);
E11 = zeros(NumUniquePts,1);
E22 = zeros(NumUniquePts,1);
E33 = zeros(NumUniquePts,1);


%h = waitbar(0,'Initializing waitbar...');
for i=1:i_last    
    Fu1 = fft(UPrimex(i,:)',nfft);
    Fhu1 = Fu1(1:end/2);
    Eu1 = (abs(Fhu1).^2)*2;   
    E11 = E11 + Eu1/i_last;
    
    Fu2 = fft(UPrimey(i,:)',nfft);
    Fhu2 = Fu2(1:end/2);
    Eu2 = (abs(Fhu2).^2)*2;   
    E22 = E22 + Eu2/i_last;
    
    Fu3 = fft(UPrimez(i,:)',nfft);
    Fhu3 = Fu3(1:end/2);
    Eu3 = (abs(Fhu3).^2)*2;   
    E33 = E33 + Eu3/i_last;
        
    perc = round((i/i_last)*100);
    %waitbar((i/i_last),h,sprintf('progress in calculating... %d%%',perc))    
end
%close(h)

%Lz = abs(Arr(1,end)-Arr(1,1));
k3=2*pi*linspace(0,1,NumUniquePts)/Lz;
lambda = 1./linspace(0,1,NumUniquePts)/Lz;
PreE11 = times(E11',k3);
PreE22 = times(E22',k3);
PreE33 = times(E33',k3);

if flag == 1
    
    loglog(k3(2:end),E11(2:end));
    grid on
    xlabel('Wavenumber')
    ylabel('Power Spectral Density')
    title('Premultiplied energy spectrum')
    hold on
    loglog(k3(2:end),E22(2:end));
    hold on
    loglog(k3(2:end),E33(2:end));
    legend('E11', 'E22', 'E33')
    hold off
    
    figure(2)
    loglog(lambda(2:end),PreE11(2:end));
    grid on
    xlabel('Wavelength')
    ylabel('Power Spectral Density')
    title('Premultiplied Energy Spectrum')
    hold on
    loglog(lambda(2:end),PreE22(2:end));
    hold on
    loglog(lambda(2:end),PreE33(2:end));
    hold off
    
end