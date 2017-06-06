function [mxss1_mt,myss1_mt,mzss1_mt] = make_blochSim(rf,sens,g,tp,df,dp,Positions)
posVecX = Positions(1,:,1,1);
posVecY = Positions(2,1,:,1);
posVecZ = Positions(3,1,1,:);
mode = 0;

[mxss1_mt,myss1_mt,mzss1_mt] =  blochSim_mex_multiThread(rf,sens,g,tp,df,dp,mode);

mxss1_mt = reshape(mxss1_mt,size(X));
myss1_mt = reshape(myss1_mt,size(X));
mzss1_mt = reshape(mzss1_mt,size(X));
