function [tform]=getTform(pMatch,views)
    
    view1 = views(1);
    view2 = views(2);
    idx1 = pMatch.matchInfo.match(1,pMatch.X);
    idx2 = pMatch.matchInfo.match(2,pMatch.X);
    feat1 = double(view1.frame(1:2,idx1))';
    feat2 = double(view2.frame(1:2,idx2))';
    img1 = view1.img;
    img2 = view2.img;
%     % imgInput = mat2gray(imgInput);
%    
%     % figure; ax = axes;
%     % showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'Parent',ax);
%     % 
    tform  = estimateGeometricTransform(feat1,feat2,'affine','MaxDistance',0.2,'Confidence',99);
    showMatch=true;
    if showMatch
        
        outputView = imref2d(size(img1));
        Ir = imwarp(img1,tform,'OutputView',outputView);
        % [H corrPtIdx] = findHomography(matchedPoints1',matchedPoints2');
        % tform = affine2d(round(H'));
        % [I21,I2_ref] = imwarp(I2,tform); % reproject img2
        % I21 = imtranslate(I2,-H(1:2,3)');
        figure;
        showMatchedFeatures(img1,img2,feat1,feat2);
        figure;
        showMatchedFeatures(Ir,img2,feat1(1,:),feat2(1,:));
        pause;
        close all;
    end
end
