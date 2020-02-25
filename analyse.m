% File to analyse tiff stacks 
% This is the begining of the code analysis
clc;
clear all;
close all;
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\';
dir = '\M1-2';
img_name = join([root,dir,'\KA0013_180106_001_001_summed_50.tif']);
tiff_info = imfinfo(img_name); % return tiff structure, one element per image
tiff_stack = imread(img_name, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(img_name, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

size_vol = 50; % The size of the volume that is averaged to downsize the 
                % image for computation
avg_stack = [];
for k=1 : size_vol : size(tiff_info,1)-size_vol
    temp_avg = mean(tiff_stack(:,:,k:k+size_vol),3);
    avg_stack = cat(3,avg_stack,temp_avg);
end

size_vol = 50; % The size of the volume that is averaged to downsize the 
                % image for computation
max_stack = [];
for k=1 : size_vol : size(tiff_info,1)-size_vol
    temp_max = max(tiff_stack(:,:,k:k+size_vol),[],3);
    max_stack = cat(3,max_stack,temp_max);
end

save(join([root,dir,'\modified_stack1.mat']),'avg_stack','max_stack');
% for ii = 1 : size(avg_stack, 3)
%     imwrite(avg_stack(:,:,ii) , join([root,dir,'\avg_stack.tif']) , 'WriteMode' , 'append') ;
% end