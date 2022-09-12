clear
close all

date = 20220908;
file = 1018;

DataFolder = ['C:\Users\co2la\Documents\Data\' num2str(date)];

load(fullfile(DataFolder, [num2str(file) '.mat']));
 
Z = Surface*1000; %hight : mm to um




Xsize = size(Z,1); %size of X axis
Ysize = size(Z,2); %size of Y axis
[X,Y] = meshgrid(-(Xsize-1)/2:(Xsize-1)/2);
xdata = zeros(size(X,1),size(Y,2),2);
pixel=0.329; % Pixel FOV in um
xdata(:,:,1) = X*pixel;%in um
xdata(:,:,2) = Y*pixel;%in um


half_width = 20; %range for fitting (um)
widthpx = round(half_width/pixel); %range for fitting (pixel)
[Xfit,Yfit] = meshgrid(-widthpx:widthpx);


d2 = (Xfit).^2 + (Yfit).^2 <= widthpx.^2;
D2 = double(d2);
D2(D2<1)=nan;


xdata_fit = zeros(widthpx*2+1,widthpx*2+1,2);

xdata_fit(:,:,1) = Xfit*pixel.*D2;%in um exclude outside
xdata_fit(:,:,2) = Yfit*pixel.*D2;%in um
Z_fit = Z((Xsize-1)/2-widthpx:(Xsize-1)/2+widthpx,(Xsize-1)/2-widthpx:(Xsize-1)/2+widthpx).*D2;


% initian fit parameters and limits for Gaussian
 xg0 = [1            ,0   ,500         ,0   ,500         , 0, 10]; %xg0 = [Amp,x0,wx,y0,wy,angle,fi]; Inital guess parameters for gaussian fit
 lbg = [0.01*xg0(1) ,-50 ,0.01*xg0(3) ,-50 ,0.01*xg0(5) , -pi(), 0.01*xg0(7)]; %lower boundry for gaussian fit
 ubg = [100*xg0(1)   ,50  ,5*xg0(3)    ,50  , 5*xg0(5)   ,pi() ,1000*xg0(7)]; %upper boundry for gaussian fit

% initian fit parameters and limits for Hemisphere
 xh0 = [200           ,0   ,0   ,1, 0]; %xh0 = [radius,x0,y0,bg]; Inital guess parameters for Hemisphere fit
 lbh = [0.001*xh0(1)  ,-50 ,-50 ,-0.001*xh0(4) ,-pi()]; %lower boundry for Hemisphere fit
 ubh = [10*xh0(1)     ,50  , 50 , 0.001*xh0(4),pi()]; %upper boundry for Hemisphere fit
 
 
 Zdata = Z_fit(~isnan(Z_fit));

[xg,resnormg,residualg,exitflagg] = lsqcurvefit(@D2GaussFunctionRot,xg0,xdata_fit, Zdata,lbg,ubg);
% 
% [xh,resnormh,residualh,exitflagh] = lsqcurvefit(@D2HemisphereFunctionRot,xh0,xdata_fit,Zdata,lbh,ubh);



RfGx = (xg(3)^2)/xg(1);%radius calculated from gaussian fit
RfGy = (xg(5)^2)/xg(1);%radius calculated from gaussian fit
varg = var(residualg(:),1); %calculate distribution
% varh = var(residualh(:),1); %calculate distribution

difg = Z_fit - D2GaussFunctionRot2(xg,xdata_fit);
% difh = Z_fit - H2(xh,xdata_fit);

ang = xg(6)*180/pi();
ratio = xg(3)/xg(5);
fprintf(1, 'Eccentricity is %.3f   angle(deg): %.3f\n', ratio,  ang)


figure(1)
subplot(1,2,1)
% C = del2(Z);
mesh(xdata(:,:,1),xdata(:,:,2),Z) %plot data
% alpha(0.6) 
% hold on
% surface(xdata_fit(:,:,1),xdata_fit(:,:,2),D2GaussFunctionRot2(xg,xdata_fit),'EdgeColor','none') %plot fit
% axis([-(half_width+50) (half_width+50) -(half_width+50) (half_width+50) -0.5 4 ])
% alpha(0.8)  
% hold off
daspect([1 1 0.05])
title('Measured surface')
%txtg = ['G = -' num2str(xg(1)) 'exp(-(x-' num2str(xg(2)) ')^2/2*' num2str(xg(3)) '+' num2str(xg(4)) ')']

% text(80, 0, -1,txtg)



subplot(1,2,2)

surface(xdata_fit(:,:,1),xdata_fit(:,:,2),D2GaussFunctionRot2(xg,xdata_fit),'EdgeColor','none') %plot fit
axis([-(half_width+50) (half_width+50) -(half_width+50) (half_width+50) -0.5 4 ])
alpha(0.8)  
view(3)

daspect([1 1 0.05])
strg2 = {'Gaussian fit' 'fitting area = ±' num2str(half_width)  'Equivalent RoC for x = ' num2str(RfGx) 'Equivalent RoC for y = ' num2str(RfGy) 'center : x = ' num2str(xg(2)) 'center : y = ' num2str(xg(4)) 'Eccentricity ' num2str(ratio)};
title(strg2)




% figure(2)
% % subplot(1,2,2)
% %C = del2(Z);
% mesh(xdata(:,:,1),xdata(:,:,2),Z,C) %plot data
% alpha(0.6) 
% hold on
% surface(xdata_fit(:,:,1),xdata_fit(:,:,2),D2HemisphereFunctionRot2(xh,xdata_fit),'EdgeColor','none') %plot fit
% axis([-(half_width+50) (half_width+50) -(half_width+50) (half_width+50) -0.5 4 ])
% alpha(0.8)  
% hold off
% daspect([1 1 0.05])
% title('Hemisphere fit')
% %txth = ['H = -1*sqrt((' num2str(xh(1)) '^2-(x-' num2str(xh(2)) ')^2-(x-' num2str(xh(3)) ')).^2) +' num2str(xh(1)) ' + ' num2str(xh(4))]
% % -1*sqrt(x(1).^2-(xdata(:,:,1)-x(2)).^2-(xdata(:,:,2)-x(3)).^2) +x(1) + x(4); %Hemisphere
% txth = {'RoC = ' num2str(xh(1)) 'Distribution = ' num2str(varh) 'center : x = ' num2str(xh(2)) 'center : y = ' num2str(xh(3))};
% text(80, 0, -1,txth)


% %%
 figure(3)
%  (Xsize-1)/2-widthpx:(Xsize-1)/2+widthpx
  figlength = length(difg);
%  left = round(figlength/2);
%  right = figlength-left
 xa = linspace(-half_width,half_width,figlength);
 ya = linspace(-half_width,half_width,figlength);
%  subplot(1, 2, 1)
 % surf(residualg)

 imagesc(xa,ya,difg);
 axis equal tight
 colorbar
 
 strg = ['Residual for Gaussian fit : ' num2str(varg)];
 title(strg)

xg




% %  G2 = @(x,xdata) -1*x(1)*exp(   -((xdata(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdata(:,:,2)-y0rot).^2/(2*x(5)^2) )    )+x(6); %Gaussian
% %  H2 = @(x,xdata) -1*sqrt(x(1).^2-(xdata(:,:,1)-x0rot).^2-(xdata(:,:,2)-y0rot).^2) +x(1) + x(4); %Hemisphere


%  G = @(x,xdata) -1*x(1)*exp(   -((xdata(:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,2)-x(4)).^2/(2*x(5)^2) )    )+x(6); %Gaussian
%  H = @(x,xdata) -1*sqrt(x(1).^2-(xdata(:,1)-x(2)).^2-(xdata(:,2)-x(3)).^2) +x(1) + x(4); %Hemisphere

function Fg = D2GaussFunctionRot(x,xdata)
xaxisrot= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
yaxisrot= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
x0rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(4)*cos(x(6));

xaxisrot = xaxisrot (~isnan (xaxisrot));
yaxisrot = yaxisrot (~isnan (yaxisrot));
L = length(xaxisrot);

xyaxis = zeros(L,2);
xyaxis(:,1) = xaxisrot;
xyaxis(:,2) = yaxisrot;

 Fg = -x(1)*exp(   -((xyaxis(:,1)-x0rot).^2/(2*x(3)^2) + (xyaxis(:,2)-y0rot).^2/(2*x(5)^2) )    )+x(7);
end
 
% function Fh = D2HemisphereFunctionRot(x,xdata)
% xaxisrot= xdata(:,:,1)*cos(x(5)) - xdata(:,:,2)*sin(x(5));
% yaxisrot= xdata(:,:,1)*sin(x(5)) + xdata(:,:,2)*cos(x(5));
% x0rot = x(2)*cos(x(5)) - x(4)*sin(x(5));
% y0rot = x(2)*sin(x(5)) + x(4)*cos(x(5));
% 
% xaxisrot = xaxisrot (~isnan (xaxisrot));
% yaxisrot = yaxisrot (~isnan (yaxisrot));
% L = length(xaxisrot);
% 
% xyaxis = zeros(L,2);
% xyaxis(:,1) = xaxisrot;
% xyaxis(:,2) = yaxisrot;
% Fh = -1*sqrt(x(1).^2-(xyaxis(:,1)-x0rot).^2-(xyaxis(:,1)-y0rot).^2) +x(1) + x(4);
% end

function Fg2 = D2GaussFunctionRot2(x,xdata)
xdatarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
x0rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(4)*cos(x(6));
Fg2 = -x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    )+x(7);
end
%  
% function Fh2 = D2HemisphereFunctionRot2(x,xdata)
% xdatarot(:,:,1)= xdata(:,:,1)*cos(x(5)) - xdata(:,:,2)*sin(x(5));
% xdatarot(:,:,2)= xdata(:,:,1)*sin(x(5)) + xdata(:,:,2)*cos(x(5));
% x0rot = x(2)*cos(x(5)) - x(4)*sin(x(5));
% y0rot = x(2)*sin(x(5)) + x(4)*cos(x(5));
% %  H = @(x,xdata) -1*sqrt(x(1).^2-(xdata(:,1)-x(2)).^2-(xdata(:,2)-x(3)).^2) +x(1) + x(4); %Hemisphere
% Fh2 = -1*sqrt(x(1).^2-(xdatarot(:,:,1)-x0rot).^2-(xdatarot(:,:,1)-y0rot).^2) +x(1) + x(4);
% end