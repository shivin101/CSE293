% addpath(genpath('/data/gal/ImageGraphs'));
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\';
dir = '\M1-2';
nam = join([root,dir,'\KA0013_180115_001_001_summed_50.tif']);
title_str = 'KA0013_180115_001_001';

info = imfinfo(nam);
T = length(info);   % number of frames
d1 = info.Height;   % height of the image
d2 = info.Width;    % width of the image
Ysiz = [d1, d2, T]';
A_full_quant = bigread2(nam, 1, 8000);
A_full_quant= double(A_full_quant);
%
region = preprocess_data(A_full_quant);
%
figure;
subplot(121);imagesc(max(A_full_quant,[],3));axis image
subplot(122);imagesc(max(region,[],3)); axis image;
figure
subplot(121);imagesc(mean(A_full_quant,3)); axis image;
subplot(122);imagesc(mean(region,3)); axis image;
%%
% if 1
%     for i =1:1000
%         subplot(121);imagesc(A_full_quant(:,:,i),[0 200]);axis image
%         subplot(122);imagesc(region(:,:,i),[0 200]);axis image;title(num2str(i));pause(0.1)
%     end
%     
%     figure(2);
%     inds = [75 150 259 365];
%     for i =1:4
%         subplot(2,4,i);
%         imagesc(A_full_quant(:,:,inds(i)),[0 150]);axis image;colormap gray
%         subplot(2,4,i+4);
%         imagesc(region(:,:,inds(i)),[0 150]);axis image
%     end
% end