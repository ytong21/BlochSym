%[mx,my,mz] = blochSim_mex_multiThread(complex(real(b1),imag(b1)),sens,gr,tpvec,df,dxyz,0);

%[2:54]  
posVec = -10:0.2:10;
[X,Y,Z] = meshgrid(posVec,posVec,posVec);
dp = [X(:) Y(:) Z(:)];
%%
tp = 0.0001; %100 us
x = -0.002:tp:0.002;
deltaF = 3000;
b1 = sin(deltaF*x*pi)./(deltaF*x*pi);
b1(21) = 1;
sens = complex(400,zeros(size(dp,1),1));
gr = 400*[ones(numel(x),1) zeros(numel(x),2)];
df = zeros(size(dp,1),1);
mode = 0;

tic
[mxss1_mt,myss1_mt,mzss1_mt] =  blochSim_mex_multiThread(complex(b1.',0),sens,gr,tp,df,dp,mode);

multiThread = toc;

mxss1_mt = reshape(mxss1_mt,size(X));
myss1_mt = reshape(myss1_mt,size(X));
mzss1_mt = reshape(mzss1_mt,size(X));
mxyss1_mt = abs(mxss1_mt+1i*myss1_mt);