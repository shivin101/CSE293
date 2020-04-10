% File to analyse tiff stacks 
% This is the begining of the code analysis
clc;
clear all;
close all;

%Root directory where the directories to be processed are stored
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\Process_data\';
save_root  = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\';

%First two files are .. and .
dir_list = dir(root);
dir_list = dir_list(3:end);


for j = 1:size(dir_list)
    val = dir_list(j);
    image_files = dir([val.folder,'\',val.name,'\*.tif']);
    roi_files = dir([val.folder,'\',val.name,'\*.roi']);
    
    save_dir = [save_root,'\',val.name];
    if isdir(save_dir)
        continue;
    else
        mkdir(save_dir);
        for i = 1:size(image_files)
            img_file = image_files(i);
            roi_file = roi_files(i);
            img_name = [img_file.folder,'\',img_file.name];
            roi_name = [roi_file.folder,'\',roi_file.name];
            ROI = load(roi_name,'-mat');
            disp(img_name);
            ROI = ROI.polygon.ROI;
            tiff_info = imfinfo(img_name); % return tiff structure, one element per image
    %         tiff_stack = imread(img_name, 1) ; % read in first image
            tiff_stack = bigread2(img_name, 1);
    %         concatenate each successive tiff to tiff_stack
    %         for ii = 2 : size(tiff_info, 1)
    %             temp_tiff = bigread2(img_name, ii);
    %             tiff_stack = cat(3 , tiff_stack, temp_tiff);
    %         end

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

            save(join([save_dir,'\modified_stack',num2str(i),'.mat']),'avg_stack','max_stack','ROI');
        end
    end
end

clear all;
close all;
% for ii = 1 : size(avg_stack, 3)
%     imwrite(avg_stack(:,:,ii) , join([root,dir,'\avg_stack.tif']) , 'WriteMode' , 'append') ;
% end