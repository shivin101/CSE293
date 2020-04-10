clc;
clear all;
close all;

tic

ids=[1,4];
denoise = 0;
root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\';
dir_list = dir(root);
dir_list = dir_list(3:end);

dir_idx =10;
val = dir_list(dir_idx);
mod_dir = [val.folder,'\',val.name];
%ROI ={1,size(ids,2)};

num_sessions = dir([mod_dir,'\','*.mat']);

file_list = {};
%Hard coding some naming conventions in this stage of the process
lenId = size(ids,2);
for i = 1:size(ids,2)
    file_list(i) = {[mod_dir,'\','modified_stack',num2str(ids(i)),'.mat']};
    file(i) = file_list(i);
    load(file_list{i});
    
    %Create an array of ROI from preprocessed arrays
%     filename = root+file(i);
%     ret = importdata(filename);
    Roi(i) = {ROI};
%     bg_roi = ret.bgROI;
    
    %Create an array of Image array
    
    % imagesc(max(max_stack,[],3));
    img(i) = {max_stack};
    if denoise
        for j =1:size(img{i},3)
            curr_img = max_stack;
            den_img = preprocess_data(double(curr_img(:,:,j)));
            curr_img(:,:,j)=den_img;
        end
        img(i)={curr_img};
    end
    % hold on;
end

%Compare 2 at a time for now
for i =1:lenId
    for j=i+1:lenId
        img1 = img{i};
        img2 = img{j};
        ROI1 = Roi{i};
        ROI2 = Roi{j};
        [matches1,matches2]=get_corresp(img1,ROI1,img2,ROI2);
        I1 = max(img1,[],3);
        I2 = max(img2,[],3);
        
        idx = (lenId-1)*(i-1)+j-1;
%         P(idx) = subtightplot(lenId-1,lenId-1,idx);
%         showMatchedFeatures(I1,I2,matches1,matches2);
        
        tform  = estimateGeometricTransform(matches1,matches2,'affine','MaxDistance',1.5,'Confidence',99,'MaxNumTrials',2000);
        outputView = imref2d(size(I1));
        figure;
        ax=axes;
        showMatchedFeatures(I1*20,I2*20,matches1,matches2,'montage','Parent',ax);
        Ir = imwarp(I1,tform,'OutputView',outputView);
        figure;
        ax=axes;
        showMatchedFeatures(Ir*20,I2*20,matches1(1,:),matches2(1,:));
        fprintf('%d-%d\n',i,j);
        pause;
        close all;
    end
end

time_el = toc



