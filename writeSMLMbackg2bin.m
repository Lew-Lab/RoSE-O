function writeSMLMbackg2bin(source_dir, image_stack_name, backg_stack_name)
%writeSMLMbackg2bin writes image_stack_name and backg_stack_name data into
%binary format with the same name at source_dir
%INPUT
%source_dir: directory containing required data
%image_stack_name: name of the image stack
%backg_stack_name: name of the background stack or image (in case you have only one background image)

%% load the data

S_img = load(fullfile(source_dir, image_stack_name));
S_backg = load(fullfile(source_dir, backg_stack_name));

%% save data in binary format

% image stack
fieldName = fieldnames(S_img);
fileID = fopen(fullfile(source_dir, 'SMLM_img.bin'), 'w');
fwrite(fileID, getfield(S_img, fieldName{1}), 'single');
fclose(fileID);

% background
fieldName = fieldnames(S_backg);
fileID = fopen(fullfile(source_dir, 'backg.bin'), 'w');
fwrite(fileID, getfield(S_backg, fieldName{1}), 'single');
fclose(fileID);
