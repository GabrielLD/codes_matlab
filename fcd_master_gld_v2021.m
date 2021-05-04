%% fcd_master _ gld 2021
% basé sur le code de Sanders Wildeman https://github.com/swildeman/fcd


clear all;
close all;
%% LOAD FILES to test 
path = 'F:\Gabriel\Schlieren\Reynolds_ridge\'
date = '20201126'
data_folder = '20201126_wildeman_h10mm_Qps0p3mLmin_c0p2molL_TTAB'
datafiles = strcat(path, date, '\', data_folder, '\', 'image_sequence', '\', '*tif');
files=dir(datafiles)
%% Paramètres de la manip

nair = 1; % air optic index
nliquid = 1.33; % liquid optic index
alpha = 1-nair/nliquid; 
dmm = 6 %mm width of seringe injection measured with a pied a coulisse
dpx = 12 % px width of seringe injection on the reference image
fx = dmm/dpx % ratio mm/px


hp = 130 % mm distance between checkerboard pattern and interface
hp = hp/fx % px;
H = 1000/fx % px ; distance between interface and camera
%% Test du script sur deux images
Iref = imread(files(1).name); % images de reference
% on prend la premiere en evitant le '.' et '..'
figure(1)
imshow(Iref)
Image_number = 800; 
Idef = imread(files(Image_number).name);
figure(2)
imshow(Idef)
% FCD master code de Sanders
% convert images to double to prevent rounding errors
Iref = double(Iref);
Idef = double(Idef);
% get two independent carrier peaks from reference image
[kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf);
% extract carrier signals from reference image and store them for later use
krad = sqrt(sum((kr-ku).^2))/2;
fIref = fft2(Iref);
cr = getcarrier(fIref, kr, krad);
cu = getcarrier(fIref, ku, krad);
% get displacement field and height profile
tic
fIdef = fft2(Idef);
[u,v] = fcd_dispfield(fIdef,cr,cu);
% Here is where the code differs from Wildeman code
hstarinverse = (1/(alpha*hp)+1/H); % 1/hstar = 1/(\alpha*hp)+1/H
u_rescale = u*hstarinverse; % u = u*1/hstar
v_rescale = v*hstarinverse; % v = v*1/hstar
% Then we integrate the gradient field to get h the deformation of the
% interface
h = fftinvgrad(-u_rescale,-v_rescale);
h_mm = h*fx;

toc

% display results
close all;
%set(0,'DefaultFigureWindowStyle','docked')
figure(1)
imshowpair(Iref,Idef,'montage')
%

figure(2)
d = size(h)
X = linspace(0,max(d),max(d))
Y = linspace(0,min(d),min(d))
imagesc(X*fx,Y*fx, h_mm);
xlabel('x (mm)')
ylabel('y (mm)')
cbar =colorbar()
cbar.Label.String = '\zeta (mm)'
caxis([-.1,.1])
set(gca, 'DataAspectRatio', [1 1 1])
%
%%
figure(3)

subplot(2,1,1)
U_norm = sqrt(u_rescale.*u_rescale+v_rescale.*v_rescale); % je calcule la norme du gradiet de h
imagesc(X*fx,Y*fx,U_norm)
cbar1=colorbar()
%caxis([0, 1.1])
xlabel('x (mm)')
ylabel('y (mm)')
cbar1.Label.String = 'delta  / h*'
set(gca, 'DataAspectRatio', [1 1 1])
xlim([100, 300])
ylim([50, 210])
caxis([0, 0.035])
subplot(2,1,2)
imagesc(X*fx,Y*fx, h_mm);
cbar1=colorbar()
%caxis([0, 1.1])
xlabel('x (mm)')
ylabel('y (mm)')
cbar = colorbar()
cbar.Label.String = '\zeta (mm)'
caxis([-.15, .15])
set(gca, 'DataAspectRatio', [1 1 1])
xlim([100, 300])
ylim([50, 210])

%% process full data folder
M = {}
M.name = data_folder;
M.fx = fx;
M.hp = 150;
M.hp_unit = 'mm'
M.H = 1000;
M.Hunit = 'mm'
amount = length(files)-1
for i=750:850
    files(i).name
    Idef = imread(files(i).name);
    % get two independent carrier peaks from reference image
    [kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf);
    % extract carrier signals from reference image and store them for later use
    krad = sqrt(sum((kr-ku).^2))/2;
    fIref = fft2(Iref);
    cr = getcarrier(fIref, kr, krad);
    cu = getcarrier(fIref, ku, krad);
    % get displacement field and height profile
    tic
    fIdef = fft2(Idef);
    [u,v] = fcd_dispfield(fIdef,cr,cu);
    % Here is where the code differs from Wildeman code
    hstarinverse = (1/(alpha*hp)+1/H); % 1/hstar = 1/(\alpha*hp)+1/H
    u_rescale = u*hstarinverse; % u = u*1/hstar
    v_rescale = v*hstarinverse; % v = v*1/hstar
    % Then we integrate the gradient field to get h the deformation of the
    % interface
    h = fftinvgrad(-u_rescale,-v_rescale);
    h_mm = h*fx;
    M(i).u_rescale=u_rescale
    M(i).v_rescale=u_rescale
    M(i).h_mm = h_mm
end

%% Movie frame 
figure(1)

amount = length(files)-1
for i=750:850
    files(i).name
    Idef = imread(files(i).name);
    % get two independent carrier peaks from reference image
    [kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf);
    % extract carrier signals from reference image and store them for later use
    krad = sqrt(sum((kr-ku).^2))/2;
    fIref = fft2(Iref);
    cr = getcarrier(fIref, kr, krad);
    cu = getcarrier(fIref, ku, krad);
    % get displacement field and height profile
    tic
    fIdef = fft2(Idef);
    [u,v] = fcd_dispfield(fIdef,cr,cu);
    % Here is where the code differs from Wildeman code
    hstarinverse = (1/(alpha*hp)+1/H); % 1/hstar = 1/(\alpha*hp)+1/H
    u_rescale = u*hstarinverse; % u = u*1/hstar
    v_rescale = v*hstarinverse; % v = v*1/hstar
    % Then we integrate the gradient field to get h the deformation of the
    % interface
    h = fftinvgrad(-u_rescale,-v_rescale);
    h_mm = h*fx;
    

    
    %subplot(2,1,1)
    U_norm = sqrt(u_rescale.*u_rescale+v_rescale.*v_rescale); % je calcule la norme du gradiet de h
    imagesc(X*fx,Y*fx,U_norm)
    cbar1=colorbar()
    %caxis([0, 1.1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    cbar1.Label.String = '\delta r / h*'
    set(gca, 'DataAspectRatio', [1 1 1])
    xlim([100, 300])
    ylim([50, 210])
    caxis([0, 0.04])
    
%     %subplot(2,1,2)
%     imagesc(X*fx,Y*fx, h_mm);
%     cbar1=colorbar()
%     %caxis([0, 1.1])
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     cbar = colorbar()
%     cbar.Label.String = '\zeta (mm)'
%     caxis([-.15, .15])
%     set(gca, 'DataAspectRatio', [1 1 1])
%     xlim([100, 300])
%     ylim([50, 210])
    export_fig(['solo_gradient_Wildemanmethod_h10mm_Qps0p3mLmin_TTAB_c0p2molL_n',num2str(i),'.jpg'])
end

