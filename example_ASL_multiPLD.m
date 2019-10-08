% demonstration of diffusion weighted imaging analysis for IAPM Early Careers Workshop 2019
% example data is publically available from https://www.fmrib.ox.ac.uk/primers/intro_primer/ExBox2/IntroBox2.html
% Alan Stone, TCD, 19/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop';

% add nifti toolbox to path
addpath(genpath([tutorial_directory filesep 'downloads' filesep 'NIfTI_20140122']))

% read asl data into matlab
asl_dataset = load_nii([tutorial_directory filesep 'data' filesep 'asl' filesep 'multiPLDpcASL' filesep 'asltc.nii.gz']);
S_asl = double(asl_dataset.img);

% mean brain signal time-course
S_asl_global_timecourse = squeeze(mean(mean(mean(S_asl,1),2),3));
figure, plot(S_asl_global_timecourse)
title('ASL tag-control series (mean brain signal)')
xlabel('Brain Volume')
ylabel('[A.U.]')

% tag-control subtraction
dS_asl_tmp = S_asl(:,:,:,2:2:end) - S_asl(:,:,:,1:2:end);

% block average multiple TIs
dS_asl = cat(4, mean(dS_asl_tmp(:,:,:,1:2:8),4), ...
                mean(dS_asl_tmp(:,:,:,9:2:16),4), ...
                mean(dS_asl_tmp(:,:,:,17:2:24),4), ...
                mean(dS_asl_tmp(:,:,:,25:2:32),4), ...
                mean(dS_asl_tmp(:,:,:,33:2:40),4), ...
                mean(dS_asl_tmp(:,:,:,41:2:48),4));

%% Fit Model
% lsqcurvefit params
x0 = [60 1]; %initial guess for cbf & aat
lb = [0, 0]; % vector of lower bounds
ub = [300,10]; % vector of upper bounds
options = optimset('display','off');

% constants taken from Alsop et al. (2015). Magnetic Resonance in Medicine, 73(1), 102–116. https://doi.org/10.1002/mrm.25197
input.lambda = 0.9; % blood brain partition coefficent
input.t1b = 1.650; % [s] @ 3T
input.alpha = 0.85; % [%] labeling efficiency for PCASL
input.tau = 1.4; % [s] labelling duration
input.convfact = 6000; % factor converts units from mL/g/s to mL/100 g/min
input.m0b = 800; % longitudinal magnetisation of blood
input.t1d = 1.3;
input.alpha = 0.85; % [%] labeling efficiency for PCASL


% slice-timing correction parameters
slice_timing_gap = 0.0452; % [s] timing between consecutive slice acquisitions
PLD = [0.25 0.5 0.75 1.0 1.25 1.5]; % [s] Post-labelling delay

% Fit full Kinetic Model to measure CBF
% dti_z = 0.047; % [s] incremental increase in ti from slice to slice
CBF = zeros(size(dS_asl(:,:,:,1)));
AAT = zeros(size(dS_asl(:,:,:,1)));
resnorm = zeros(size(dS_asl(:,:,:,1)));
%residual = zeros(x,y,z);

[x,y,z] = size(dS_asl(:,:,:,1));
fprintf('GO MAKE A CUP OF TEA THIS WILL TAKE A FEW MINS !!! \n')
tic
for xID = 1:x
    for yID = 1:y
        for zID = z/2%1:z

            fprintf(1,'%.0f,%.0f,%.0f \n' ,xID,yID,zID);
            input.PLD = PLD + (slice_timing_gap * (zID-1)) + input.tau; %inversion time
            [x, resnorm, residual] = lsqcurvefit('ASL_buxton_model',x0,input,...
                                        squeeze(dS_asl(xID,yID,zID,:))',lb,ub, options);
            CBF(xID,yID,zID) = x(1);
            AAT(xID,yID,zID) = x(2);
            resnorm(xID,yID,zID) = resnorm;
            %residual(xID,yID,zID) = residual;

        end
    end
end
toc

% % Calculate calibration image
% pd_dataset = load_nii([tutorial_directory filesep 'data' filesep 'asl' filesep 'singlePLDpcASL' filesep 'aslcalib.nii.gz']); % read proton density weighted image
% S_pd = double(pd_dataset.img);
% t1gm = 1.3; % [s] t1 grey matter for proton density correction
% tr = 4.8; % [s] sequence repition time
% S_pdc = S_pd./(1-exp(-tr./t1gm)); % signal intensity of proton density weighted image corrected for T1
%
% % CBF map [ml/100g/min]
% % constants taken from Alsop et al. (2015). Magnetic Resonance in Medicine, 73(1), 102–116. https://doi.org/10.1002/mrm.25197
% lambda = 0.9; % blood brain partition coefficent
% t1b = 1.650; % [s] @ 3T
% alpha = 0.85; % [%] labeling efficiency for PCASL
% tau = 1.4; % [s] labelling duration
% convfact = 6000; % factor converts units from mL/g/s to mL/100 g/min
%
% % slice-timing correction
% slice_timing_gap = 0.0452; % [s] timing between consecutive slice acquisitions
% PLD = PLDs (s) - 0.25, 0.5, 0.75, 1.0, 1.254, 1.5.; % [s] Post-labelling delay
% nslices = size(dS_asl,3); % number of slices
% % slice-timing correction for PLD
% PLDcorr = PLD .* ones(size(dS_asl));
% for id = 1:nslices
%     PLDcorr(:,:,id) = PLDcorr(:,:,id) + (slice_timing_gap * id);
% end
%
% % CBF quantification
% CBF = (convfact .* lambda .* dS_asl .* exp(PLDcorr/t1b)) ./ (2 .* alpha .* t1b .* S_pdc .* (1-exp(-tau./t1b))); % [ml/100g/min]
% % crude brain mask
% CBF = CBF .* mask;

% set min and max viewing window
min = 0;
max = 500;

% view ASL volume
figure('name','ASL','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
countA = 1;
countB = 1;
for id = 1:3:18
    subplot(6,3,id), imshow(S_asl(:,:,end/2,countA), 'displayrange', [min max])
    hold on, title(sprintf('ASL PLD_%s - Axial',num2str(countB))), c = colorbar; c.Label.String = '[A.U.]';
    subplot(6,3,id+1), imshow(imrotate(squeeze(S_asl(:,end/2,:,countA)),90), 'displayrange', [min max])
    hold on, title(sprintf('ASL PLD_%s - Coronal',num2str(countB))), c = colorbar; c.Label.String = '[A.U.]';
    subplot(6,3,id+2), imshow(imrotate(squeeze(S_asl(end/2,:,:,countA)),90), 'displayrange', [min max])
    hold on, title(sprintf('ASL PLD_%s - Sagittal',num2str(countB))), c = colorbar; c.Label.String = '[A.U.]';
    countA = countA + 16;
    countB = countB + 1;
end

% view averaged tag-control subtraction
figure('name','\Delta ASL','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
countA = 1;
countB = 1;
for id = 1:3:18
    subplot(6,3,id), imshow(dS_asl(:,:,end/2,countA), 'displayrange', [0 10])
    hold on, title(sprintf('dASL PLD_%s - Axial',num2str(countA))), c = colorbar; c.Label.String = '[A.U.]';
    subplot(6,3,id+1), imshow(imrotate(squeeze(dS_asl(:,end/2,:,countA)),90), 'displayrange', [0 10])
    hold on, title(sprintf('dASL PLD_%s - Coronal',num2str(countA))), c = colorbar; c.Label.String = '[A.U.]';
    subplot(6,3,id+2), imshow(imrotate(squeeze(dS_asl(end/2,:,:,countA)),90), 'displayrange', [0 10])
    hold on, title(sprintf('dASL PLD_%s - dSagittal',num2str(countA))), c = colorbar; c.Label.String = '[A.U.]';
    countA = countA + 1;
end

% view CBF & AAT map
figure('name','CBF & AAT','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
c1 = subplot(2,3,1); imshow(CBF(:,:,end/2), 'displayrange', [0 100]), colormap(c1, 'Jet')
hold on, title('CBF - Axial'), c = colorbar; c.Label.String = '[ml/100g/min]';
c2 = subplot(2,3,2); imshow(imrotate(squeeze(CBF(:,end/2,:)),90), 'displayrange', [0 100]), colormap(c2, 'Jet')
hold on, title('CBF - Coronal'), c = colorbar; c.Label.String = '[ml/100g/min]';
c3 = subplot(2,3,3); imshow(imrotate(squeeze(CBF(end/2,:,:)),90), 'displayrange', [0 100]), colormap(c3, 'Jet')
hold on, title('CBF - Sagittal'), c = colorbar; c.Label.String = '[ml/100g/min]';

c1 = subplot(2,3,1); imshow(AAT(:,:,end/2), 'displayrange', [0 3]), colormap(c1, 'Jet')
hold on, title('AAT - Axial'), c = colorbar; c.Label.String = '[ml/100g/min]';
c2 = subplot(2,3,2); imshow(imrotate(squeeze(AAT(:,end/2,:)),90), 'displayrange', [0 3]), colormap(c2, 'Jet')
hold on, title('AAT - Coronal'), c = colorbar; c.Label.String = '[ml/100g/min]';
c3 = subplot(2,3,3); imshow(imrotate(squeeze(AAT(end/2,:,:)),90), 'displayrange', [0 3]), colormap(c3, 'Jet')
hold on, title('AAT - Sagittal'), c = colorbar; c.Label.String = '[ml/100g/min]';
