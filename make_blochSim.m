function [mxss1_mt,myss1_mt,mzss1_mt] = make_blochSim(rf,b1,g,dt,Positions)
posVecX = squeeze(Positions(1,:,1,1));
posVecY = squeeze(Positions(2,1,:,1));
posVecZ = squeeze(Positions(3,1,1,:));
[X,Y,Z] = meshgrid(posVecX,posVecY,posVecZ);
dp = [X(:) Y(:) Z(:)];
df = zeros(size(dp,1),1);
mode = 0;
sens = b1(:);
nRF = size(rf,1);
nG = size(g,1);
rfFull = [zeros(nG-nRF,1); rf];
[mxss1_mt,myss1_mt,mzss1_mt] =  blochSim_mex_multiThread(rfFull,sens,g,dt,df,dp,mode);

mxss1_mt = reshape(mxss1_mt,size(X));
myss1_mt = reshape(myss1_mt,size(X));
mzss1_mt = reshape(mzss1_mt,size(X));
