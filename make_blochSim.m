function magnetization = make_blochSim(rf,b1,g,dt,dp,mask)

df = zeros(size(dp,1),1);
mode = 0;
sens = b1(:);
nRF = size(rf,1);
nG = size(g,1);
rfFull = [zeros(nG-nRF,1); rf];
[mxss1_mt,myss1_mt,mzss1_mt] =  blochSim_mex_multiThread(rfFull,sens,g,dt,df,dp,mode);

mx = zeros(size(mask));
my = zeros(size(mask));
mz = zeros(size(mask));
mx(~mask==0) = mxss1_mt;
my(~mask==0) = myss1_mt;
mz(~mask==0) = mzss1_mt;
mxy = abs(mx+1i*my);

magnetization = struct('mx',mx,'my',my,'mz',mz,'mxy',mxy);

figure
slice = round(size(mxy,3)/2);
imagesc(mxy(:,:,slice))
title Excitation Pattern
axis image;colorbar
