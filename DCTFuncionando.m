close all
clear all
x = double(imread('boat.tiff')); %reading image of 512x512 pixels
figure, imshow(x/255); %displaying image (Figura 1)

y = x; %storing image in another variable

%Creating Watermark

a = zeros(512,512); %the size of the matriz depends on the size of the image that im reading
a(100:250,100:350) = 1;

figure, imshow(a); %displaying watermark (Figura 2)

save m.dat a -ascii %saving matriz in a file

%Perform WaterMarking

%RGB content of the image
x1 = x(:,:,1);


%Fazer a dct em cada bloco e dps substituir

%Perform dct on each RGB component and storing in a different variable
dx1 = dct2(x1); dx11 = dx1 ;


%----------------------------------------------- COLOCAR MARCA D'ÁGUA

load m.dat %binary mark for watermarking. Loading the matriz from the file that I created
g = 100; %coeeficient of watermark's strength (Cox call it alpha parameter)
%As the size of g gets bigger, the watermarking gets stronger, but if this value is too big, the image will be deteriorated

[rm, cm] = size(m); %rows and columns of the matriz 300x500

dx1(1:rm,1:cm) = dx1(1:rm,1:cm) + g*m; %adding the watermark to the image by adding the coefficient g to the image

figure, imshow(dx1); %displaying each component of the image after watermarking (Figura 3)

%------------------------------------------------ INVERSA

%Perform idct on each RGB component and storing in a different variable

%Também fazer a inversa bloco à bloco

y1 = idct2(dx1);

y(:,:,1) = y1;

figure, imshow(y1); %Inverse in each component (Figura 6)


figure('Name','Watermarked Image'), imshow(y/255); %displaying the watermarked image (Figura 9)

figure; imshow(abs(y-x)*100); %comparing the original image with the watermarked image( Figura 10)


z = y; %storing the watermarked image in another variable
[r, c, s] = size(z); %rows, columns and channels of the watermarked image

%De-watermarking
% Clean image (known mask)

y = z;
%Perform dct on each RGB component and storing in a different variable
dy1 = dct2(y(:,:,1));

%Subtracting the watermark from the RGB channels
dy1(1:rm,1:cm) = dy1(1:rm,1:cm) - g*m;

%Perform idct on each RGB component to get the original image
y11 = idct2(dy1);


l = double(imread('boat3.tiff'));

%combining the RGB channels to get the original image (yy)
yy(:,:,1) = y11;

%displaying the original image
figure, imshow(yy/255);
%comparison showing all black image for no difference b/w yy and x

%figure, imshow(abs(yy-l)*10000);
figure, imshow(abs(yy-x)*10000);



%figure, imshowpair(y1,x,"diff")



pctImage = imabsdiff(yy, x) ./ double(yy);
meanPct = mean2(pctImage); %difference percentual


%Analysing differences
%Z = imabsdiff(y1,x);

%figure, imshowpair(y1,x,"diff")

%figure, imshowpair(y1, imread("boat2.tiff"),"diff");






