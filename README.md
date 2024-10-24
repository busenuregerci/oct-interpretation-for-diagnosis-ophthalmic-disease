%%
% <http://www.mathworks.com MathWorks> 
clear 
close all
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 

cd 'D:\BUSE\Desktop\oct';
matfiles = dir('*.mat'); % Adding files
N = length(matfiles);
N = 1;
for i = 1:N;
    load(matfiles(i).name)
    figure

    h = image(images(:,:,2));
    Y = get(h,'CData');

    net=denoisingNetwork('DnCNN');
    Y=denoiseImage(Y, net);
    figure, montage(Y, 'DisplayRange',[]);

    %First Method
   K = medfilt2(Y,[5 5]);
   [J,noise_out] = wiener2(Y,[5 5]);
   J=imdiffusefilt(J);

    wavelength=4*pi;
    orientation = [-pi/2:pi/180:pi/2];
    g = gabor(wavelength,orientation);
    [outMag,outPhase] = imgaborfilt(J,g);
    outSize = size(outMag);
    outMag = reshape(outMag,[outSize(1:2),1,outSize(3)]);
    figure, montage(outMag,'DisplayRange',[]);
    title('Montage of gabor magnitude output images.');
    
    outSize = size(outPhase);
    outPhase = reshape(outPhase,[outSize(1:2),1,outSize(3)]);
    figure, montage(outPhase,'DisplayRange',[]);
    title('Montage of gabor phase output images.');
    
    figure
    imshow(real(g(1,1).SpatialKernel))
    title('Real part of Spatial Kernel');
    
    
     for x=1:size(outMag,1)
        for y=1:size(outMag,2)
     for k=1:size(outMag,4)       
        R(x,y)=max(outMag(x,y,1,k));
     end    
        end
     end


   

    figure
    Rcap=log(R);
    imagesc(Rcap)
    [counts,x] = imhist(Rcap,256);
    T = otsuthresh(counts);
    figure
    imagesc(Rcap>T);


    %Second Method

     [LoD,HiD] = wfilters('haar','d');
    [cA,cH,cV,cD] = dwt2(Y,LoD,HiD,'mode','symh');
    subplot(2,2,1)
    imagesc(cA)
    colormap gray
    title('Approximation')
    subplot(2,2,2)
    imagesc(cH)
    colormap gray
    title('Horizontal')
    subplot(2,2,3)
    imagesc(cV)
    colormap gray
    title('Vertical')
    subplot(2,2,4)
    imagesc(cD)
    colormap gray
    title('Diagonal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Idiffusion1 = imdiffusefilt(cA);
 figure,imagesc(Idiffusion1 );
 title('diffusion filter of aproximation')

 Idiffusion2 = imdiffusefilt(cH);
 figure,imagesc(Idiffusion2 );
 title('diffusion filter of horizontal')

 Idiffusion3 = imdiffusefilt(cV);
 figure,imagesc(Idiffusion3 );
 title('diffusion filter of vertical')

 Idiffusion4 = imdiffusefilt(cD);
 figure,imagesc(Idiffusion4 );
 title('diffusion filter of diagonal')

   
    
    I1=uint16(Idiffusion1);
    [L,centers]=imsegkmeans(I1,3);
    B1=labeloverlay(I1,L);
    figure,imagesc(B1);
    title('segmentation of diffusion filter of aproximation')

    I2=uint16(Idiffusion2);
    [L,centers]=imsegkmeans(I2,3);
    B2=labeloverlay(I2,L);
    figure,imagesc(B2);
    title('segmentation of diffusion filter of horizontal')

    I3=uint16(Idiffusion3);
    [L,centers]=imsegkmeans(I3,3);
    B3=labeloverlay(I3,L);
    figure,imagesc(B3);
    title('segmentation of diffusion filter of vertical')

    I4=uint16(Idiffusion4);
    [L,centers]=imsegkmeans(I4,3);
    B4=labeloverlay(I4,L);
    figure,imagesc(B4);
    title('segmentation of diffusion filter of diagonal')

end
