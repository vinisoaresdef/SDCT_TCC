
close all
clear all
Im = imread('boat.tiff'); %reading image of 512x512 pixels
figure('Name','Imagem Original'), imshow(Im); %displaying image (Figura 1)

y = Im; %storing image in another variable

a = zeros(512,512); %the size of the matriz depends on the size of the image that im reading
%a(100:250,100:350) = 1;
a(100:300,100:400) = 1;

figure('Name','MarcaDágua'), imshow(a); %displaying watermark (Figura 2)

save m.dat a -ascii %saving matriz in a file

%----------------------------------------------

%Perform SDCT

%RGB content of the image
%Im = Im(:,:,1);

%Perform dct on each RGB component and storing in a different variable
%sdctx1 = dct2(Im); dx11 = dx1 ;

Im = imread('boat.tiff');  %ler a imagem de entrada
N  = 8;        % tamanho dos blocos
if (size(Im,3)==3)
        %Transformar a imagem RGB em cinza
    Im = rgb2gray(Im);     
end

f1 = figure(1);
f1 = figure('Name', 'Imagem principal');
imagesc(Im);

tam   = size(Im);
i_aux = (tam(1)/N) - 1;
j_aux = (tam(2)/N) - 1;

%Número de ângulos
num_ang  = (N^2 - N)/2;

%A imagem possui valores inteiros, convertemos para double
Im = double(Im);

%Ângulos de rotação (aleatório)
theta    = 2*pi*rand(num_ang,1);

%Matriz de rotação:
M_sdct   = matriz_sdct(theta, N, num_ang);
%Blk --- Bloco
for i = 0 : i_aux
    i_ant = i*N + 1;
    i_dps = i_ant + N - 1;
             
    for j = 0 : j_aux
        j_ant = j*N + 1;
        j_dps = j_ant + N -1;
        
        Blk = Im(i_ant:i_dps, j_ant:j_dps); 
        
        I_dct       = dct2(Blk);
        I_dct_vec   = vetorizar(I_dct);

%Produto da matriz de rotação com a dct vetorizada
        I_sdct_vec  = M_sdct*I_dct_vec;

%resultado é I_sdct
        I_sdct(i_ant:i_dps, j_ant:j_dps)= vetorTOmatriz(I_sdct_vec, N);

    end
end

% imshow(I_sdct);
% imagesc(Im)
% imagesc(I_sdct);

dx1 = I_sdct; dx11 = dx1;

%-----------------------------------------------

load m.dat %binary mark for watermarking. Loading the matriz from the file that I created
g = 10; %coeeficient of watermark's strength (Cox call it alpha parameter)
%As the size of g gets bigger, the watermarking gets stronger, but if this value is too big, the image will be deteriorated

[rm, cm] = size(m); %rows and columns of the matriz 512x512

dx1(1:rm,1:cm) = dx1(1:rm,1:cm) + g*m; %adding the watermark to the image by adding the coefficient g to the image

figure('Name','Componente com Marca dágua'), imshow(dx1); %displaying each component of the image after watermarking (Figura 3)

%Se precisa da imagem original: método cego
%Se não precisa da imagem original: método não cego.


%-----------------------------------------------Início da Inversa
%-------------ESSA PARTE PRECISA RECEBER O VALOR DE DX1 E NÃO DE I_SDCT
%JUSTAMENTE PQ DX1 JÁ POSSUI A MARCA D'ÁGUA INSERIDA. FAZER A INVERSA SÓ
%VAI MOSTRAR O VALOR INICIAL SE N FOR TROCADO PRA DX1.

for i = 0 : i_aux
    i_ant = i*N + 1;
    i_dps = i_ant + N - 1;
             
    for j = 0 : j_aux
        j_ant = j*N + 1;
        j_dps = j_ant + N -1;
        
        Blk   = dx1(i_ant:i_dps, j_ant:j_dps); %troquei I_SDCT por dx1
        
        I_sdct_vec = vetorizar(Blk);
%transposta vezes a I SDCT vetorizada
        I_sdct_vec_inversa = M_sdct'*I_sdct_vec;
        
        I_sdct_blk_inversa = vetorTOmatriz(I_sdct_vec_inversa, N);
        
        I_sdct_inversa(i_ant:i_dps, j_ant:j_dps)= idct2(I_sdct_blk_inversa);

    end
end

%Métricas----------------------------------
%Erro médio quadrático: (uint8 - Transformar em inteiro)
MSE  = immse(uint8(Im), uint8(I_sdct_inversa));
%Índice de similaridade estrutural entre a original e a inversa:
SSIM = ssim(Im, I_sdct_inversa);
%-------------------------------------------

figure(2);
figure('Name','Inversa já com marca dágua'), imagesc(I_sdct_inversa);

figure('Name','Inversa já com marca dágua em preto e branco'), imshow(I_sdct_inversa/255);

figure('Name','Imagem Binária'), imshow(I_sdct_inversa);


%---------------------------------------------Fim da Inversa

pctImage = imabsdiff(Im, I_sdct_inversa) ./ double(Im);
meanPct = mean2(pctImage); %difference percentual

y1 = I_sdct_inversa;

y(:,:,1) = y1; %assigning the inverse on an image component

%%%figure, imshow(y);

figure('Name','Diferenças entre Im e y1'), imshowpair(y1,Im,"diff")




%========================================================================================
%========================================================================================
%========================================================================================
%==============================DEWATERMARKING============================================
%========================================================================================
%========================================================================================

%Imagem já deve estar com a marca d'água a partir daqui. Nos próximos
%passos deve ser feito a retirada da marca d'água.













































%------------------------------Funções

function M = matriz_sdct(theta, N, num_ang)
M = eye(N^2, N^2);  %matriz identidade para os coeff's que n�o s�o rotacionados
c = ordem_varredura_2D(N);

for i = 1 : num_ang
        
    M(c(i,1), c(i,1)) = cos(theta(i));
	M(c(i,1), c(i,2)) = sin(theta(i));
	
		
	M(c(i,2), c(i,1)) = -1*sin(theta(i));
	M(c(i,2), c(i,2)) = cos(theta(i));

		
		
end


end

function [c] = ordem_varredura_2D(blk_size)

a=linspace(1,blk_size^2,blk_size^2);
idx=reshape(a,blk_size,blk_size);
idx2=idx';
c1=[];
c2=[];
for i=2:blk_size
    c1=zigzag(tril(idx,-1));% zigzag ordering
    c1=c1(c1>0);
    c2=zigzag(tril(idx2,-1));
    c2=c2(c2>0);
%     c1=[c1 idx(i,1:i-1)];% row ordering
%     c2=[c2 idx2(i,1:i-1)];
end
c= [c1', c2'];

end

function out=zigzag(in)
% Zig-zag scanning
% This function is used to rearrange a matrix of any size into a 1-D array
% by implementing the ZIG-ZAG SCANNING procedure.
% IN specifies the input matrix of any size
% OUT is the resulting zig-zag scanned (1-D) vector
% having length equal to the total number of elements in the 2-D input matrix
%
% As an example,
% IN = [1	2	6	7
%		3	5	8	11
%		4	9	10	12];
% OUT = ZIGZAG(IN)
% OUT=
%	1     2     3     4     5     6     7     8     9    10    11    12

%
%
% Oluwadamilola (Damie) Martins Ogunbiyi
% University of Maryland, College Park
% Department of Electrical and Computer Engineering
% Communications and Signal Processing
% 22-March-2010
% Copyright 2009-2010 Black Ace of Diamonds.

[num_rows, num_cols]=size(in);

% Initialise the output vector
out=zeros(1,num_rows*num_cols);

cur_row=1;	cur_col=1;	cur_index=1;

% First element
%out(1)=in(1,1);

while cur_row<=num_rows & cur_col<=num_cols
	if cur_row==1 & mod(cur_row+cur_col,2)==0 & cur_col~=num_cols
		out(cur_index)=in(cur_row,cur_col);
		cur_col=cur_col+1;							%move right at the top
		cur_index=cur_index+1;
		
	elseif cur_row==num_rows & mod(cur_row+cur_col,2)~=0 & cur_col~=num_cols
		out(cur_index)=in(cur_row,cur_col);
		cur_col=cur_col+1;							%move right at the bottom
		cur_index=cur_index+1;
		
	elseif cur_col==1 & mod(cur_row+cur_col,2)~=0 & cur_row~=num_rows
		out(cur_index)=in(cur_row,cur_col);
		cur_row=cur_row+1;							%move down at the left
		cur_index=cur_index+1;
		
	elseif cur_col==num_cols & mod(cur_row+cur_col,2)==0 & cur_row~=num_rows
		out(cur_index)=in(cur_row,cur_col);
		cur_row=cur_row+1;							%move down at the right
		cur_index=cur_index+1;
		
	elseif cur_col~=1 & cur_row~=num_rows & mod(cur_row+cur_col,2)~=0
		out(cur_index)=in(cur_row,cur_col);
		cur_row=cur_row+1;		cur_col=cur_col-1;	%move diagonally left down
		cur_index=cur_index+1;
		
	elseif cur_row~=1 & cur_col~=num_cols & mod(cur_row+cur_col,2)==0
		out(cur_index)=in(cur_row,cur_col);
		cur_row=cur_row-1;		cur_col=cur_col+1;	%move diagonally right up
		cur_index=cur_index+1;
		
	elseif cur_row==num_rows & cur_col==num_cols	%obtain the bottom right element
        out(end)=in(end);							%end of the operation
		break										%terminate the operation
    end
end
end

function C = vetorizar( coeff )
%--------------------------------------------------------------
%  Funcao que vetoriza os coeficientes 
%  >> converte a matriz de coeficientes 
%     em um vetor de coeficentes    
%--------------------------------------------------------------


[row,column] = size(coeff);
aux = row*column;
C = zeros(aux, 1);

cont = 0;

for k = 1 : row
    for l = 1 : column
            cont = cont + 1;            
            C(cont) = coeff(k,l);
    end
end

end

function [ C ] = vetorTOmatriz( coeff, tam )

C  = zeros(tam,tam); 
cont = 0;

for k = 1 : tam
    for l = 1 : tam   
            cont = cont + 1; 
          
            %matriz
            C(k,l) = coeff(cont);
    end
end

end
