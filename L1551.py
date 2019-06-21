using HDF5,PyPlot,PyCall,supermodule,FITSIO,Images,GTA,PCA,StatsBase,module_hy
@pyimport numpy as np
@pyimport mpl_toolkits.axes_grid1 as axgrid

f=FITS("13CO_Xs_v0.1_g7.5.fits")
f1=h5open("/L1551_planck2018_10arcmin.h5")
f2=h5open("Planck_offset.h5");

d=read(f[1]);
nx,ny,nv=size(d)

CRPIX3=read_header(f[1])["CRPIX3"];
CRVAL3=read_header(f[1])["CRVAL3"];
CDELT3=read_header(f[1])["CDELT3"];

peak=217;#13co:217,12co:217
I,C,Ch=getICCh(d,CRVAL3,CDELT3,CRPIX3,peak);

#planck
dn=16;
U=read(f1,"U");
Q=read(f1,"Q");
off=read(f2,"L1551");
phi=.5*atan2(-U,Q);
phi1=phi-pi/2;
phi1=phi1-(pi/2-off*pi/180);
phi1[phi1.<-pi/2]+=pi;

Gaussian=[0.32,0.32];
MW1=0.0;
MW2=2.0;
cutoff=-1000000000000;
#function re_rotation(I::Mat,Gaussian,MW1,MW2,cutoff,percentile,b::Mat,dn)
inach,bnach,Ch2=VGT(Ch, Gaussian, MW1,MW2,11, -phi1, dn);
inach=inach+pi/2;
print(string(AM(inach,bnach)))

#########################HRO
theta=abs(inach-bnach-pi/2)*180/pi;
nxx,nyy=size(theta);
n=round((findmax(theta)[1]-180)/180,0)
for k in 1:n
	for i in 1:nxx, j in 1:nyy
		if (theta[i,j]<0)
			theta[i,j]+=180;
		elseif(theta[i,j]>180)
			theta[i,j]-=180;
		end
	end
end

theta[theta.>90]=180-theta[theta.>90];

figure()
PyPlot.plt[:hist](theta[:], bins=100, range=[0,90], normed = true, histtype = "step",label="Yue",color="r");
xlim(0,90)

y=fit(Histogram,abs(theta[!isnan(theta)]),0:0.9:90)
x=[0:0.9:90]
Gaus(x,p)=p[1]*exp(-(x-p[2]).^2/p[3])+p[4];
fit_t=curve_fit(Gaus,x[2:end],y.weights,[maximum(y.weights),88,1.0,1.0]);
z=Gaus(x, fit_t.param)

figure()
bar(x[1:end-1],y.weights*100/sum(y.weights),color="lime",label="Relative angle: L 1551");
plot(x,z*100/sum(z),"--",color="k",label="Gaussian fitting");
ylabel("Percent[%]")
xlabel("Relative angle[degree]")
title("L 1551")
legend()
xlim(0,90)
ylim(0,6)
vlines(fit_t.param[2]-1,0,maximum(z*100/sum(z)),color="k")

nxx,nyy=size(theta);
n=0;
for i in 1:nxx,j in 1:nyy
 if(~isnan(theta[i,j]))
  n+=1
 end
end
error=std(theta[~isnan(theta)])/sqrt(n)

theta=abs(inach-bnach)*180/pi;
nxx,nyy=size(theta);
n=round((findmax(theta)[1]-180)/180,0)
for k in 1:n
	for i in 1:nxx, j in 1:nyy
		if (theta[i,j]<0)
			theta[i,j]+=180;
		elseif(theta[i,j]>180)
			theta[i,j]-=180;
		end
	end
end


theta[theta.>90]=180-theta[theta.>90];

figure()
ax = gca()
im=imshow(Ch3', cmap="Blues",origin="xy")
#colorbar()
contour(Ch3',alpha=0.1);
xticks([0,50,100,150,200,250,300],["68.2","68.1","68.0","67.9","67.8","67.7","67.6"]);
yticks([0,50,100,150,200,250,300,350],["17.8","17.9","18.0","18.1","18.2","18.3","18.4","18.5"]);
title("L1551_Planck")
xlabel("R.A.(J2000)[degree]");
ylabel("DEC.(J2000)[degree]");
title("L1551_Channel_13CO_AM="*string(round(AM(inach,bnach),2)));

nnx,nny=size(inach);
for i in 1:nnx,j in 1:nny
 if(~isnan(inach[i,j]))
  xc=Yd[i,j];
  yc=Xd[i,j];
  if(theta[i,j]<=10)
   plot(xc,yc,color="red","o");
  elseif(theta[i,j]>10 && theta[i,j]<20)
   plot(xc,yc,color="orange","o");
  elseif(theta[i,j]>20 && theta[i,j]<30)
   plot(xc,yc,color="yellow","o");
  elseif(theta[i,j]>30 && theta[i,j]<40)
   plot(xc,yc,color="green","o"); 
  elseif(theta[i,j]>40 && theta[i,j]<50)
   plot(xc,yc,color="lime","o"); 
  elseif(theta[i,j]>50 && theta[i,j]<60)
   plot(xc,yc,color="aqua","o"); 
  elseif(theta[i,j]>60 && theta[i,j]<70)
   plot(xc,yc,color="c","o");
  elseif(theta[i,j]>70 && theta[i,j]<80)
   plot(xc,yc,color="b","o");
  elseif(theta[i,j]>80 && theta[i,j]<90)
   plot(xc,yc,color="purple","o"); 
  end
 end
end

tick_params(direction="in")
divider = axgrid.make_axes_locatable(ax)
cax = divider[:append_axes]("right", size="5%", pad=0.00)
colorbar(im, cax=cax)
##########################error bar##########
am=2*cos(inach-bnach).^2-1;
nxx,nyy=size(am);
n=0;
for i in 1:nxx,j in 1:nyy
 if(~isnan(am[i,j]))
  n+=1
 end
end
error=std(am[~isnan(am)])/sqrt(n)

Ma13[isnan(inach)]=NaN;
Ma13[Ma13.>=10]=NaN;
n=0;
for i in 1:nxx,j in 1:nyy
 if(~isnan(Ma13[i,j]))
  n+=1
 end
end
error=std(Ma13[~isnan(Ma13)])/sqrt(n)

#######################################
Chi=similar(Ch); 
Chi[isnan(Ch)]=1;
Ch[Chi.==1]=minimum(Ch);
Ch2=Images.imfilter_gaussian(Ch,[.32,.32]);
Ch2x=sban2d_intensity(Ch2,dn);
chx,chy=sobel_conv_2d(Ch2);
cha=atan(chy./chx);

dis=dispersion(cha,dn);
nxx,nyy=size(dis);
f2=h5open("/mnt/f/survey2/L1551/L1551_Ma.h5");
Ma13=read(f2,"13CO");

Xd,Yd=np.meshgrid(div(dn,2)+1:dn:div(dn,2)+ny,div(dn,2)+1:dn:div(dn,2)+nx);
Ch3=Images.imfilter_gaussian(Ch,[6,6]);
figure()
ax = gca()
im=imshow(Ch3', cmap="Blues",origin="xy")
#colorbar()
contour(Ch3',alpha=0.1);
xticks([0,50,100,150,200,250,300],["68.2","68.1","68.0","67.9","67.8","67.7","67.6"]);
yticks([0,50,100,150,200,250,300,350],["17.8","17.9","18.0","18.1","18.2","18.3","18.4","18.5"]);
title("L1551_Planck")
xlabel("R.A.(J2000)[degree]");
ylabel("DEC.(J2000)[degree]");
title("L1551_Channel_13CO_AM="*string(AM(inach,bnach)));

nnx,nny=size(inach);
for i in 1:nnx,j in 1:nny
 if(~isnan(inach[i,j]))
  theta=-inach[i,j];
  xc=Yd[i,j];
  yc=Xd[i,j];
  l=6;
  xe=xc+l*cos(theta);
  ye=yc+l*sin(theta);
  xi=xc-l*cos(theta);
  yi=yc-l*sin(theta);

  theta1=theta+4*log(dis[i,j]);
  ld=3;
  xe1=xc+ld*cos(theta1);
  ye1=yc+ld*sin(theta1);
  xi1=xc-ld*cos(theta1);
  yi1=yc-ld*sin(theta1);
  if(xi1<xe1)
   x1=collect(xi1:0.1:xe1);
   y1=x1*tan(theta1)+(yi1*xe1-ye1*xi1)/(xe1-xi1);
  elseif(xi1>xe1)  
   x1=collect(xe1:0.1:xi1);
   y1=x1*tan(theta1)+(yi1*xe1-ye1*xi1)/(xe1-xi1);
  elseif(xi1==xe1) 
   y1=collect(yi1:0.1:ye1);
   x1=xc*zeros(size(y1));
  end 
  plot(x1,y1,color="k");
  
  theta2=theta-4*log(dis[i,j]);
  xe2=xc+ld*cos(theta2);
  ye2=yc+ld*sin(theta2);
  xi2=xc-ld*cos(theta2);
  yi2=yc-ld*sin(theta2);
  if(theta2<pi/2)
   x2=collect(xi2:0.1:xe2);
   y2=x2*tan(theta2)+(yi2*xe2-ye2*xi2)/(xe2-xi2);
  elseif(theta2>pi/2)  
   x2=collect(xe2:0.1:xi2);
   y2=x2*tan(theta2)+(yi2*xe2-ye2*xi2)/(xe2-xi2);
  elseif(theta2==pi/2)
   y2=collect(yi2:0.1:ye2);
   x2=xc*zeros(size(y2));
  end 
  plot(x2,y2,color="k");

  if(xi1<xi2)
   x3=collect(xi1:0.1:xi2);
   y3=x3*(yi2-yi1)/(xi2-xi1)+(yi2*xi1-yi1*xi2)/(xi1-xi2);
  elseif(xi1>xi2)  
   x3=collect(xi2:0.1:xi1);
   y3=x3*(yi2-yi1)/(xi2-xi1)+(yi2*xi1-yi1*xi2)/(xi1-xi2);
  elseif(xi1==xi2)
   if(yi1<yi2)
    y3=collect(yi1:0.1:yi2);
    x3=xc*zeros(size(y3));
   else
    y3=collect(yi2:0.1:yi1);
    x3=xc*zeros(size(y3));
   end 
  end
  plot(x3,y3,color="k");

  if(xe1<xe2)
   x4=collect(xe1:0.1:xe2);
   y4=x4*(ye2-ye1)/(xe2-xe1)+(ye2*xe1-ye1*xe2)/(xe1-xe2);
  elseif(xe1>xe2)  
   x4=collect(xe2:0.1:xe1);
   y4=x4*(ye2-ye1)/(xe2-xe1)+(ye2*xe1-ye1*xe2)/(xe1-xe2);
  elseif(xe1==xe2)
   if(ye1<ye2)
    y4=collect(ye1:0.1:ye2);
    x4=xc*zeros(size(ye));
   else
    y4=collect(ye2:0.1:ye1);
    x4=xc*zeros(size(y4));
   end 
  end
  plot(x4,y4,color="k"); 
   if(theta<pi/2)
   x=collect(xi:0.1:xe);
   y=x*tan(theta)+(yi*xe-ye*xi)/(xe-xi);
  elseif(theta>pi/2)  
   x=collect(xe:0.1:xi);
   y=x*tan(theta)+(yi*xe-ye*xi)/(xe-xi);
  elseif(theta==pi/2)
   y=collect(yi:0.1:ye);
   x=xc*zeros(size(y));
  end 
  if(1/Ma13[i,j]<0.25)
   plot(x,y,color="yellow");
  elseif(1/Ma13[i,j]>0.25 && 1/Ma13[i,j]<0.5)
   plot(x,y,color="orange");
  elseif(1/Ma13[i,j]>0.5 && 1/Ma13[i,j]<0.75)
   plot(x,y,color="hotpink");
  elseif(1/Ma13[i,j]>0.75 && 1/Ma13[i,j]<1.0)
   plot(x,y,color="red"); 
  elseif(1/Ma13[i,j]>1.0)
   plot(x,y,color="maroon"); 
  end
 end
end
tick_params(direction="in")
divider = axgrid.make_axes_locatable(ax)
cax = divider[:append_axes]("right", size="5%", pad=0.00)
colorbar(im, cax=cax)

