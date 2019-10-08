% demonstration of Arterial Spin Labelling analysis for IAPM Early Careers Workshop 2019
% Alan Stone, TCD, 19/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop';

% add nifti toolbox to path
addpath(genpath([tutorial_directory filesep 'downloads' filesep 'NIfTI_20140122']))

% read asl data into matlab
asl_dataset = load_nii([tutorial_directory filesep 'data' filesep 'asl' filesep 'singlePLDpcASL' filesep 'asltc.nii.gz']);
S_asl = asl_dataset.img;

% mean brain signal time-course
S_asl_global_timecourse = squeeze(mean(mean(mean(S_asl,1),2),3));
figure, plot(S_asl_global_timecourse)
title('ASL tag-control series (mean brain signal)')
xlabel('Brain Volume')
ylabel('[A.U.]')

% tag-control subtraction and average t-c pairs
dS_asl = mean(S_asl(:,:,:,2:2:end) - S_asl(:,:,:,1:2:end),4);
mask = dS_asl > 0.3; % crude brain mask

% Calculate calibration image
pd_dataset = load_nii([tutorial_directory filesep 'data' filesep 'asl' filesep 'singlePLDpcASL' filesep 'aslcalib.nii.gz']); % read proton density weighted image
S_pd = double(pd_dataset.img);
t1gm = 1.3; % [s] t1 grey matter for proton density correction
tr = 4.8; % [s] sequence repition time
S_pdc = S_pd./(1-exp(-tr./t1gm)); % signal intensity of proton density weighted image corrected for T1

% CBF map [ml/100g/min]
% constants taken from Alsop et al. (2015). Magnetic Resonance in Medicine, 73(1), 102–116. https://doi.org/10.1002/mrm.25197
lambda = 0.9; % blood brain partition coefficent
t1b = 1.650; % [s] @ 3T
alpha = 0.85; % [%] labeling efficiency for PCASL
tau = 1.8; % [s] labelling duration
convfact = 6000; % factor converts units from mL/g/s to mL/100 g/min

% slice-timing correction
slice_timing_gap = 0.0452; % [s] timing between consecutive slice acquisitions
PLD = 1.8; % [s] Post-labelling delay
nslices = size(dS_asl,3); % number of slices
% slice-timing correction for PLD
PLDcorr = PLD .* ones(size(dS_asl));
for id = 1:nslices
    PLDcorr(:,:,id) = PLDcorr(:,:,id) + (slice_timing_gap * id);
end

% CBF quantification
CBF = (convfact .* lambda .* dS_asl .* exp(PLDcorr/t1b)) ./ (2 .* alpha .* t1b .* S_pdc .* (1-exp(-tau./t1b))); % [ml/100g/min]
% crude brain mask
CBF = CBF .* mask;

% set min and max viewing window
min = 0;
max = 600;

% view ASL volume
figure('name','ASL','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
subplot(3,3,1), imshow(S_asl(:,:,end/2,1), 'displayrange', [min max])
hold on, title('S_{b0} - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,2), imshow(imrotate(squeeze(S_asl(:,end/2,:,1)),90), 'displayrange', [min max])
hold on, title('S_{b0} - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,3), imshow(imrotate(squeeze(S_asl(end/2,:,:,1)),90), 'displayrange', [min max])
hold on, title('S_{b0} - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

% averaged tag-control subtraction
subplot(3,3,4), imshow(dS_asl(:,:,end/2), 'displayrange', [min/100 max/100])
hold on, title('\Delta ASL - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,5), imshow(imrotate(squeeze(dS_asl(:,end/2,:)),90), 'displayrange', [min/100 max/100])
hold on, title('\Delta ASL - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
subplot(3,3,6), imshow(imrotate(squeeze(dS_asl(end/2,:,:)),90), 'displayrange', [min/100 max/100])
hold on, title('\Delta ASL - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

% view CBF map
c1 = subplot(3,3,7); imshow(CBF(:,:,end/2), 'displayrange', [0 100]), colormap(c1, 'Jet')
hold on, title('CBF - Axial'), c = colorbar; c.Label.String = '[ml/100g/min]';
c2 = subplot(3,3,8); imshow(imrotate(squeeze(CBF(:,end/2,:)),90), 'displayrange', [0 100]), colormap(c2, 'Jet')
hold on, title('CBF - Coronal'), c = colorbar; c.Label.String = '[ml/100g/min]';
c3 = subplot(3,3,9); imshow(imrotate(squeeze(CBF(end/2,:,:)),90), 'displayrange', [0 100]), colormap(c3, 'Jet')
hold on, title('CBF - Sagittal'), c = colorbar; c.Label.String = '[ml/100g/min]';
