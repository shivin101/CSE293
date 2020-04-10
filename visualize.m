clc;
clear all;
close all;

tic

ids=[1,2,3,4,5,6,7];

root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\New_data\';
dir_list = dir(root);
dir_list = dir_list(3:end);

dir_idx = 3;
val = dir_list(dir_idx);
mod_dir = [val.folder,'\',val.name];
ROI ={1,size(ids,2)};

num_sessions = dir([mod_dir,'\','*.mat']);

file_list = {};
%Hard coding some naming conventions in this stage of the process
for i = 3:3
    file_list(i) = {[mod_dir,'\','modified_stack',num2str(ids(i)),'.mat']};
    file(i) = file_list(i)
    load(file_list{i});
    
    %Create an array of ROI from preprocessed arrays
%     filename = root+file(i);
%     ret = importdata(filename);
    Roi(i) = {ROI};
%     bg_roi = ret.bgROI;
    
    %Create an array of Image array
    
    % imagesc(max(max_stack,[],3));
%     img(i) = {uint8(max(max_stack,[],3)/3)};
    img(i) = {uint8(max_stack/3)};
    % hold on;
    f = figure;
    roi = Roi{i}
%     implay(img{i});
    imshow(max(img{i},[],3));
    hold on;
    for j = 1:size(Roi{i},2)
       x= roi{j}(:,1);
       y= roi{j}(:,2);
       k  = convhull(x,y);
       plot(x(k),y(k ),'g','LineWidth',1.5);
       img = bwconvhull()
    end
    hold off 
    saveas(f,['../img',num2str(i),'.jpg']);
     
  


end
