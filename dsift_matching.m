
% File to analyse tiff stacks 
% This is the begining of the code analysis
clc;
clear all;
ts = sym('ts');
dir_root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\M1-2\';
dir1 = '';
dir2 = '';
dir3 = '181216';
modified_stack1 = join([dir_root,dir1,'\modified_stack1.mat']);
load(modified_stack1);
avg_stackl= avg_stack;
max_stack1= max_stack;

modified_stack2 = join([dir_root,dir2,'\modified_stack2.mat']);
load(modified_stack2);
avg_stack2= avg_stack;
max_stack2= max_stack;
% tiff_info = imfinfo(img_name); % return tiff structure, one element per image
% tiff_stack2 = imread(img_name, 1) ; % read in first image
% %concatenate each successive tiff to tiff_stack
% for ii = 2 : size(tiff_info, 1)
%     temp_tiff = imread(img_name, ii);
%     tiff_stack2 = cat(3 , tiff_stack2, temp_tiff);
% end

% % % SIFT PART
binSize = 16;
magnif = 8 ;
I1 = normalize(single(max_stack1(:,:,3)),2);
I1 = linearize_val(I1);
I2 = normalize(single(max_stack2(:,:,12)),2);
I2 = linearize_val(I2);
% I2 = imresize(I2,size(I1));
% I1 = vl_imsmooth(single(I1), sqrt((binSize/magnif)^2 - .25)) ;
% I2 = vl_imsmooth(single(I2), sqrt((binSize/magnif)^2 - .25)) ;
sampling_freq = 300;
[f1, d1] = vl_dsift(I1, 'size', binSize) ;
[f2, d2] = vl_dsift(I2, 'size', binSize) ;
f1 = f1(:,sampling_freq:sampling_freq:end);
d1 = d1(:,sampling_freq:sampling_freq:end);
f2 = f2(:,sampling_freq:sampling_freq:end);
d2 = d2(:,sampling_freq:sampling_freq:end);
[matches, scores] = vl_ubcmatch(d1, d2,1.3) ;
feature_loc_1 = f1(:,matches(1,:));
feature_loc_2 = f2(:,matches(2,:));
figure; ax = axes;
showMatchedFeatures(I1,I2,feature_loc_1',feature_loc_2','montage','Parent',ax);