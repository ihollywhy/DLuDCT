function Detect_caption(img_path)
if (nargin < 1)
    img_path = 'caption.jpg';
end   
im_read = imread(img_path);
I = rgb2gray(im_read);
[m, n] = size(I);
m = floor((m-0.1)/8); n = floor((n-0.1)/8);
I = imresize(I, [8*m+8, 8*n+8]);
I = im2double(I);

%------------------------------------------
%                                   8x8 DCT
%------------------------------------------
T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(I,[8 8],dct);
mask = [0 0 0 1 1 1 1 1; 
        0 0 0 1 1 1 1 1; 
        0 0 0 2 1 1 1 1;
        3 3 2 2 2 1 1 1;
        3 3 3 2 2 2 2 1;
        3 3 3 3 2 2 2 2;
        3 3 3 3 2 2 2 2;
        3 3 3 3 3 2 2 2];
 
%------------------------------------------
%               8x8 calc the energy by mask
%------------------------------------------
fun_energy = @(block_struct) calc_energy(block_struct.data, mask);
B_enegy= blockproc(B,[8 8],fun_energy);
energy = zeros(m+1, n+1, 3);
for i = 1 : m+1
    for j = 1 : n+1
        energy(i, j, :) = B_enegy((i-1)*3+1, (j-1)*3+1 : j*3);
    end
end

e_para = 0.95;                  %one of parameter control the threshold
energy1 = energy(:,:,1);
sort_energy = reshape(energy1, [], 1);
sort_energy = sort(sort_energy);
length = size(sort_energy,1)
mean_energy1 = sort_energy(floor(length*e_para));     %core parameter
energy1(energy1>mean_energy1) = 1;
energy1(energy1<mean_energy1) = 0;

energy2 = energy(:,:,2);
sort_energy = reshape(energy2, [], 1);
sort_energy = sort(sort_energy);
length = size(sort_energy,1)
mean_energy3 = sort_energy(floor(length*e_para));     %core parameter
energy2(energy2>mean_energy3) = 1;
energy2(energy2<mean_energy3) = 0;

energy3 = energy(:,:,3);
sort_energy = reshape(energy3, [], 1);
sort_energy = sort(sort_energy);
length = size(sort_energy,1)
mean_energy3 = sort_energy(floor(length*e_para));     %core parameter
energy3(energy3>mean_energy3) = 1;
energy3(energy3<mean_energy3) = 0;

caption = bitand(int8(energy1), int8(energy3));
caption = bitand(int8(caption), int8(energy2));
imshow(caption,[]);
% figure();
% imshow(energy2,[]);
se1=strel('square',2);   %create a square area to erode and dilate 
A2 = imerode(caption,se1);
% figure();
% imshow(A2,[]);
se1=strel('square',3);
A2 = imdilate(A2,se1);
% figure();
% imshow(A2,[]);
[ccl, n] = bwlabel(A2,4);
dets = [];
for i = 1:n
    tmp = ccl ==i;
    [y, x] = find(tmp > 0);
    size_x = max(x) - min(x);
    size_y = max(y) - min(y);
    if(size_x * size_y > 2)
        dets = [dets ; [min(x),max(x),min(y),max(y)]];
    end
end
dets = dets*8;
%------------------------------------------
%                                image show
%------------------------------------------
figure();
imshow(ccl,[]);
figure();
imshow(im_read, []);
hold on;
for i = 1:size(dets,1)
    rectangle('Position',[dets(i,1), dets(i,3), dets(i,2) - dets(i,1), dets(i,4) - dets(i,3)],...
        'Curvature',[0,0],'EdgeColor','r');
end
end

function energy = calc_energy(data, mask)
    energy = zeros(3,3);
    energy(1,1) = sum(abs(data(mask == 1)));
    energy(1,2) = sum(abs(data(mask == 2)));
    energy(1,3) = sum(abs(data(mask == 3)));
end