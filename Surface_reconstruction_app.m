function Surface_reconstruction_app

%% Loop for data acquisition
lambd = 532*10^(-6);
start = AtCube.getPosition_z(); %start position in mm
stepsize = lambd/8; %stepsize in mm
n=6; %number of frames
range = (n-1)*stepsize; %range in mm
z_list = linspace(start,start+range,n);
    
for k = 1:6 %loop to move sample in z
    c = z_list(k);
    AtCube.move_z(c);
    z=AtCube.getPosition_z();
    Filename = [num2str(c,'%.8f') 'B.tif'];

for i=1:10 %number of frames to acquire and avearge
[err, M, w, h] = uEye_camera(1); % fetch image
M(M<0)=M(M<0)+256;
M = reshape(M, 3, w*h);
I = reshape(M, 3, w, h);
I = permute(I, [3 2 1])/255;
% IG(:,:,i)=I(:,:,2); % take only the green channel
Im_BG = 0.7*I(:,:,2)+0.3*I(:,:,3);
% ImBG_autoCont = rescale(Im_BG,0,256);
end
    MeanI(:,:,k) = mean(Im_BG,3);   
end
Im_opt_cont = im2double(rescale(MeanI,0,256));
%% Crop image around fibre
sig = 5;
FilterSize = 21;
IR_filt=im2double(imgaussfilt(Im_opt_cont,sig,'FilterSize',FilterSize)); %% Filter
image_size = 500; %Choose cropping area
T = 548-(image_size/2);
B = 548+(image_size/2);
L = 968-(image_size/2);
R = 968+(image_size/2);
IR_crop = IR_filt(T:B,L:R,:); %crop and isolate fibre face by making background 0

%% Create phase map
I1=IR_crop(:,:,1);
I2=IR_crop(:,:,2);
I3=IR_crop(:,:,3);
I4=IR_crop(:,:,4);
I5=IR_crop(:,:,5);
I6=IR_crop(:,:,6);
phi2 = -atan(((3*I2-4*I4+I6)./(I1-4*I3+3*I5))); %calculate phase map

%% Unwrap and create image of surface
Frame = ones(size(phi2,1),size(phi2,2));
[PU,PC,N_unwrap,t_unwrap]=CPULSI(2*phi2,Frame,100,0.0001,300,250,false);
phase_unwrapped=0.5*PU;%Restore correct phase & remove background

%% 3D plot
pixel=0.329; % Pixel FOV in um
[X,Y]=meshgrid(1:image_size+1,1:image_size+1);   % Generate 2D meshgrid full
x=X*0.001*pixel; %in um                   
y=Y*0.001*pixel; %in um
surface=(phase_unwrapped./(2.*pi)).*(1000*lambd); % convert phase to height in mm

s = figure;
surf(x,y,surface) % Show distance calibrated results
daspect([1 1 20])
shading interp
% axis equal
colorbar;
xlabel('mm')
ylabel('mm')
zlabel('um')
