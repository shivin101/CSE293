function tiff_writer(outputFileName,mov)

% function Y = tiff_writer(name,mov)
%
% Code to write tiff stack.
%
% 2015 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs and set variables

if nargin < 2
    error('Need to give a name and a file to write!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write images

for kk=1:length(mov(1, 1, :))
   imwrite(mov(:, :, kk), outputFileName, 'WriteMode', 'append','compression','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%