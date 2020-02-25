clc;
clear all;
close all;

tic

root = 'C:\Users\shivi\OneDrive\Documents\Courses\CSE293\M1-2\';
addpath('C:\Users\shivi\OneDrive\Documents\Courses\CSE293\RSC\code\subtightplot');

% IDS to be compared are placed here[replace by some kind of file fxn]
ids = [2,6,9];
ids = sort(ids);
lenId = size(ids,2);
figure1=figure('Position', [0, 0, 1200, 1200]);
P = [];
%Compare 2 at a time for now
for i =1:lenId
    for j=i+1:lenId
        id1 = ids(i);
        id2 = ids(j);
        
        file = 'KA0013_1801'+string(num2str(id1+5,'%02d'))+'_001_001_summed_50.roi';
        filename = root+file;
        ret = importdata(filename);
        ROI1 = ret.ROI;
        bg_roi = ret.bgROI;

        filename = root+string('modified_stack'+string(id1)+'.mat');
        load(filename);
        % imagesc(max(max_stack,[],3));
        img1 = max_stack;
        % hold on;

        file = 'KA0013_1801'+string(num2str(id2+5,'%02d'))+'_001_001_summed_50.roi';

        filename = root+file;
        ret = importdata(filename);
        ROI2 = ret.ROI;
        bg_roi = ret.bgROI;

        filename = root+string('modified_stack'+string(id2)+'.mat');
        load(filename);
        % imagesc(max(max_stack,[],3));
        img2 = max_stack;
        [matches1,matches2]=get_corresp(img1,ROI1,img2,ROI2);
        I1 = max(img1,[],3);
        I2 = max(img2,[],3);
        
        idx = (lenId-1)*(i-1)+j-1;
        P(idx) = subtightplot(lenId-1,lenId-1,idx);
%         showMatchedFeatures(I1,I2,matches1,matches2);
        
        tform  = estimateGeometricTransform(matches1,matches2,'affine','MaxDistance',0.5,'Confidence',99.5);
        outputView = imref2d(size(I1));
        Ir = imwarp(I1,tform,'OutputView',outputView);
        showMatchedFeatures(Ir,I2,matches1(1,:),matches2(1,:));
    end
end

time_el = toc



