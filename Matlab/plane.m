function [x,e11,e22,e33,preE11,K3,Lambda,vOUT] = plane(dirName,dirNameLast,FileNames)

mcount = 1;
for m = dirName:0.001:dirNameLast
    
    dirNameFull = num2str(m,6)
    
    fullFileNames = strcat('./',dirNameFull,'/',FileNames);
    
    F = dlmread(fullFileNames,' ',2,0);
    count = 1;
    for i=1:length(F)
        if (F(i,1) > 0.865) && (F(i,1) < 0.920) && (F(i,3) > 0.2) && (F(i,3) < 0.255)
            cut(count,:) = F(i,:);
            count = count + 1;
        end
    end
    A{mcount} = sortrows(cut,1);
    
    count = 0;
    l = 1;
    for j=1:length(A{mcount})-1
        if A{mcount}(j,1) < A{mcount}(j+1,1)
            l = l + 1;
            count = 0;
            clear N
            continue
        end
        count = count + 1;
        N(count,:) = A{mcount}(j,:);
        N_cell{l,mcount} = N;
    end
    mcount = mcount + 1;
    clear i j l
end

jcount = 1;
for j = 1:2:size(N_cell,1)
% Here only odd dependences are taken into account, because coincidence in wavenumbers alternates through one.
  for i = 1:size(N_cell,2)
    U(i,:) = N_cell{j,i}(:,4);
    V(i,:) = N_cell{j,i}(:,5);
    W(i,:) = N_cell{j,i}(:,6);
  end

  Lz = abs(N_cell{j,1}(end,3) - N_cell{j,1}(1,3));
  
  [E11,E22,E33,PreE11,k3,lambda] = PremultSpec(U,V,W,Lz,0);
  
  if j==size(N_cell,1)-1
    [E11,E22,E33,PreE11,k3,lambda] = PremultSpec(U,V,W,Lz,1);   
  end
  
  e11(:,jcount) = E11; 
  e22(:,jcount) = E22; 
  e33(:,jcount) = E33; 
  preE11(jcount,:) = PreE11; 
  K3(jcount,:) = k3; 
  Lambda(jcount,:) = lambda;
  x(jcount) = N_cell{j,1}(1,1); 
  jcount = jcount + 1;
  clear E11 E22 E33 PreE11 k3 lambda U V W
end

Lambda = Lambda';
preE11 = preE11';
preE11(1,:) = [];
Lambda(1,:) = [];

for l=1:size(preE11,2)
    %vOUT(l,:) = savitzkyGolay1D(preE11(:,l),5,5);
    vOUT(l,:) = medfiltOne(preE11(:,l),7);
end
vOUT(vOUT<0) = 1e-7;
vOUT = vOUT'; 

figure(3)
surf(x,log(Lambda),log(vOUT))
colorbar
figure(4)
imagesc(x,log(Lambda(:,1)),log(vOUT))
colorbar

