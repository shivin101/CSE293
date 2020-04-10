function write_TPM_movie(mov, save_name)
    
% function write_TPM_movie(name,mov)
%
% Write the 3-D array in mov to a file. Can writ to either a tiff stack, a
% fits file, or a matlab '.mat' file. The mode is automatically chosen
% based on the extension of the file-name passed to the function
% (save_name). The inputs to this function are
%     - mov       - 3D array of data to save (mov(:,:,kk) is the kk^th
%                   frame of the movie)
%     - save_name - String containing the file-name to save the data as
%
% 2017 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

if ~ischar(save_name)
    error('Need to provide a string as the file-name (second input)!')     % Check that the file-name is a string
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Data

switch save_name(end-3:end)                                                % Check the file type that the save-name string implies
    case 'fits'
        fitswrite(C_poi, [data_path,save_name(1:end-4),'_C.fits']);        % Write as a fits file OR
    case '.tif'
        tiff_writer(save_name,mov);                                        % Write as a tif file
    case '.mat'
        save(save_name,'mov','-v7.3')                                      % Write as a mat file
    otherwise
        error('Unknown file type. Can only save to .fits, .mat, or .tif file-types!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%