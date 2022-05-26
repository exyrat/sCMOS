%%%% sCMOS Calibration following Huang et al 2013, Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms (Bewersdorf)
%%%% Code written 04/21/2022, A. Raterink
%
%
%
%

%% User Inputs

cam_size = 256;         % camera size (NxN pixels)
dark_frames = 250;      % # of dark images for offset
var_seqs = 3;          % # of dark sequences for variance (need M^2 frames total, but this can only handle ~20,000 at a time so you need to make L substacks of <20,000 frames)
var_frames = 1617;      % # of frames per dark sequence for variance
illum_levels = 4;       % # of illumination intensity levels for gain calibration
illum_frames = 2000;   % # of frames per intensity ( not actually necessary, just used for simulation. But should be <20,000 )

offset_darks_fn = "vars1.tif"; % filename of dark sequence for offset calibration
var_darks_fn = "vars";          % filename prefix for the L dark sequences for variance calibration. e.g. "vars1.tiff", "vars2.tiff", ..., "varsL.tiff". Must be TIFF for now
illum_fn = "illum";             % filename prefix for the K illumination sequences for gain estimation. e.g. "ilum1.tiff", "illum2.tiff", ..., "illumK.tiff". Must be TIFF for now

%
%
%
%
%

%% Offset Calibration

%dark_ims = round(20*rand(cam_size,cam_size,dark_frames));   % sim dark images
%t = Tiff(offset_darks_fn,'r');
dark_ims = FastTiff('darks.tif'); % load dark sequence for offset calc
offset = mean(dark_ims,3); % calculate offsets

dark_var = var(dark_ims,[],3);
dark_std = std(dark_ims,[],3);

%% Variance Calibration

varsum = 0;
for i = 1:var_seqs
    %var_ims = round(20*rand(cam_size,cam_size,var_frames));     % sim ith dark sequence
    
    var_ims = FastTiff(var_darks_fn + num2str(i) + ".tif");  % load ith dark sequence which has <20,000 frames
    offset_stack = repmat(offset, [1 1 length(var_ims(1,1,:))]);
    varsum = varsum+mean(var_ims.^2 - offset_stack.^2,3);   % takes temporal mean and adds to sum of means
end

var_map = varsum/var_seqs;                                         % final mean over L to get variance map

%% Gain Estimation

A = [];
B = [];
for k = 1:illum_levels
    illum_im = round(20*rand(cam_size,cam_size,illum_frames))+20*10^k;          % sim kth illumination sequence

    %illum_im = imread(illum_fn + num2str(k) + ".tiff"); % load kth illumination sequence
    Dk = mean(illum_im,3);                              % mean ADU per pixel
    vark = var(illum_im,[],3);                          % mean variance per pixel

    % subtract dark variance and offset and concatenate to matrix 
    A = cat(3,A,vark-var_map);                          % A = {(var1-var_map), ... , (varK-var_map)}                             
    B = cat(3,B,Dk-offset);                             % B = {(D1-offset), ... , (Dk-offset)} 

end

g = zeros([cam_size cam_size]);
for i=1:cam_size           % doing this iteratively because can't take transpose of N-D matrix (FIX ?)
    for j = 1:cam_size
        Aij = A(i,j,:); Aij = Aij(:)';                  % Aij (and Bij) is column vector of pixel i,j across K illumination sequences
        Bij = B(i,j,:); Bij = Bij(:)';                  
        g(i,j) = inv(Bij*Bij')*Bij*Aij';                % gain estimate at pixel i,j
    end
end

%% Save the Maps

save('offset_map.mat','offset');
save('variance_map','var_map');
save('gain_map.mat','g');

%% Plot the Maps and Distributions

subplot(3,1,1); imshow(offset,[]); title('Offset Map')
subplot(3,1,2); imshow(var_map,[]); title('Variance Map')
subplot(3,1,3); imshow(g,[]); title('Gain Map')

figure;
subplot(3,1,1); histogram(offset); xlabel('ADU Offset'); ylabel('Counts');
subplot(3,1,2); histogram(var_map); xlabel('ADU Variance'); ylabel('Counts');
subplot(3,1,3); histogram(g); xlabel('Gain'); ylabel('Counts');


%% Tif reader

function data = FastTiff(filename)
warning('off','all') % Suppress all the tiff warnings
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(I,J,K);
data(:,:,1)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(:,:,n) = tstack.read();
end
warning('on','all')
end