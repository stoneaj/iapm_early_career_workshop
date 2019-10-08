% demonstration of diffusion weighted imaging analysis for IAPM Early Careers Workshop 2019
% example data is publically available from https://www.fmrib.ox.ac.uk/primers/intro_primer/ExBox2/IntroBox2.html
% Alan Stone, TCD, 19/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop'; % NOTE: Change this to directory on your machine

% add nifti toolbox to path
addpath(genpath([tutorial_directory filesep 'downloads' filesep 'NIfTI_20140122']))

% read diffusion imaging data into matlab
S_dwi = load_nii([tutorial_directory filesep 'data' filesep 'dwi' filesep 'ExBox2' filesep 'data.nii.gz']);

% b=0 is first volume (4th dimension)
S_b0 = S_dwi.img(:,:,:,1);

% average b-weighted images to get diffusion trace
S_b1500 = mean(S_dwi.img(:,:,:,2:end),4);

% diffusion weighting used
b = 1500; % [s/mm2]

% ADC calculation
ADC = -1./b .* log(S_b1500./S_b0);

% set min and max viewing window
min = 0;
max = 700;

% view b=0 ... first volume (4th dimension) is diffusion weighting
figure('name','DWI','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
subplot(3,3,1), imshow(S_b0(:,:,end/2), 'displayrange', [min max])
hold on, title('S_{b0} - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,2), imshow(imrotate(squeeze(S_b0(:,end/2,:)),90), 'displayrange', [min max])
hold on, title('S_{b0} - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,3), imshow(imrotate(squeeze(S_b0(end/2,:,:)),90), 'displayrange', [min max])
hold on, title('S_{b0} - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

% view b=1500 ... first volume (4th dimension) is diffusion weighting
subplot(3,3,4), imshow(S_b1500(:,:,end/2), 'displayrange', [min max])
hold on, title('S_{b1500} - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,5), imshow(imrotate(squeeze(S_b1500(:,end/2,:)),90), 'displayrange', [min max])
hold on, title('S_{b1500} - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,6), imshow(imrotate(squeeze(S_b1500(end/2,:,:)),90), 'displayrange', [min max])
hold on, title('S_{b1500} - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

% view b=1500 ... first volume (4th dimension) is diffusion weighting
subplot(3,3,7), imshow(ADC(:,:,end/2), 'displayrange', [0 0.002])
hold on, title('ADC - Axial'), c = colorbar; c.Label.String = '[mm2/s]';
subplot(3,3,8), imshow(imrotate(squeeze(ADC(:,end/2,:)),90), 'displayrange', [0 0.002])
hold on, title('ADC - Coronal'), c = colorbar; c.Label.String = '[mm2/s]';
subplot(3,3,9), imshow(imrotate(squeeze(ADC(end/2,:,:)),90), 'displayrange', [0 0.002])
hold on, title('ADC - Sagittal'), c = colorbar; c.Label.String = '[mm2/s]';
