% demonstration of diffusion weighted imaging analysis for IAPM Early Careers Workshop 2019
% Data is from:
% Blockley, N. P., Harkin, J. W., Stone, A. J., & Bulte, D. P. (2016).
% Data acquired to investigate new approaches to cerebrovascular reactivity mapping using MRI.
% University of Oxford.
% Alan Stone, TCD, 19/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop';

% read asl data into matlab
load([tutorial_directory filesep 'data' filesep 'bold_cvr' filesep 'bold_cvr.mat']);
bold = bold_cvr;

% read end tidal breath-by-breath respiratory data
end_tidal = load([tutorial_directory filesep 'data' filesep 'bold_cvr' filesep 'sub-01_end_tidal.txt']);
scan_start_id = find(abs(end_tidal(:,1) - 76.877) < eps); % find index of start time = 76.877 mins
scan_end_id = find(abs(end_tidal(:,1) - 85.847) < eps); % find index of start time = 85.875 mins

% extract temporary ET-CO2
time_tmp =  (end_tidal(scan_start_id:scan_end_id,1) - end_tidal(scan_start_id,1)) .* 60; % [secs] normalise time
et_co2_tmp = end_tidal(scan_start_id:scan_end_id,2); % end tidal co2

% mean brain signal time-course
bold_global = squeeze(mean(mean(mean(bold,1),2),3));
tr = 6.1; % [sec] TR
timei = (1:size(bold,4)) .* tr; % repetition time

% co2 interpolated to TR
et_co2i_tmp = interp1(time_tmp,et_co2_tmp,timei,'linear','extrap');

% normalise ET-CO2 & BOLD signals to overlay
et_co2i_tmp_norm = (et_co2i_tmp - min(et_co2i_tmp))./(max(et_co2i_tmp) - min(et_co2i_tmp));
bold_global_norm = (bold_global - min(bold_global))./(max(bold_global) - min(bold_global));

% calculate lung to brain delay
[acor lag] = xcorr(bold_global_norm,et_co2i_tmp_norm);
[~,I] = max(acor);
timelag = lag(I);
fprintf('Lag = %f\n',timelag)

% extract lagged ET-CO2
scan_start_id_lag = scan_start_id - timelag;
scan_end_id_lag = scan_end_id - timelag;
time =  (end_tidal(scan_start_id_lag:scan_end_id_lag,1) - end_tidal(scan_start_id_lag,1)) .* 60; % [secs] normalise time
et_co2 = end_tidal(scan_start_id_lag:scan_end_id_lag,2); % end tidal co2

% lagged co2 interpolated to TR
et_co2i = interp1(time,et_co2,timei,'linear','extrap');

% normalise ET-CO2 & BOLD signals to overlay
et_co2i_norm = (et_co2i - min(et_co2i))./(max(et_co2i) - min(et_co2i));
bold_global_norm = (bold_global - min(bold_global))./(max(bold_global) - min(bold_global));

% General Linear Model - Global
% MRI signal = β1⋅EtCO2 + β2⋅t + β0 + ε
nvols = length(bold_global); % number of MRI images
x0 = ones(nvols,1); % mean
x1 = et_co2i'; % delta with CO2
x2 = (timei(1):timei(end)/(nvols):timei(end))'; % drift
X = [x0 x1 x2];
b = regress(bold_global,X);
mri_signal_model = (b(2).*et_co2i) + (b(3).*timei) + b(1);

% Calculate Global CVR ... this is meaningless global signal includes voxels outside brain
%                      ... more meaningufl if caclulated over a grey matter mask
% CVR = β1 / (β0 + (min(EtCO2)*β1))
CVR_global = b(2) / (b(1) + (min(et_co2i)*b(2)));
% CVR_global = b(2) / (b(1) + (prctile(et_co2i,10)*b(2)));
fprintf('Global CVR: %f %%dBOLD/mmHg-CO2 \n', CVR_global.*100)

% Voxel-wise fitting
[x,y,z] = size(bold(:,:,:,1));
CVR = zeros(x,y,z);

for xID = 1:x
    for yID = 1:y
        for zID = 1:z

            % fprintf(1,'%.0f,%.0f,%.0f \n' ,xID,yID,zID);
            b = regress(squeeze(bold(xID,yID,zID,:)),X);
            CVR(xID,yID,zID) = b(2) / (b(1) + (min(et_co2i)*b(2))) .* 100;
            % CVR(xID,yID,zID) = b(2) / (b(1) + (prctile(et_co2i,10)*b(2)));

        end
    end
end

% make crude mask
mask = bold(:,:,:,1) > 100;
CVR_brain = mask .* CVR;

% view end-tidal CO2
figure, hold on
subplot(3,2,1), plot(time_tmp,et_co2_tmp)
title('End-tidal CO_2')
xlabel('Time [s]')
ylabel('P_{ET|CO2} [mmHg]'), grid on

% view global brain signal
subplot(3,2,2), plot(timei, bold_global)
title('BOLD-CVR (mean brain signal)')
xlabel('Time [s]')
ylabel('[A.U.]'), grid on

% view normalised signals overlayed
subplot(3,2,3),
h1 = plot(timei, et_co2i_tmp_norm); hold on
h2 = plot(timei, bold_global_norm);
title('Normalised signals')
xlabel('Time [s]'), legend([h1 h2], {'ET-CO2_{norm}','BOLD_{norm}'}), grid on

% view normalised signals overlayed
subplot(3,2,4),
h1 = plot(timei, et_co2i_norm); hold on
h2 = plot(timei, bold_global_norm);
title('Normalised signals (Lagged)')
xlabel('Time [s]'), legend([h1 h2], {'ET-CO2_{norm|lag}','BOLD_{norm}'}), grid on

% view normalised signals overlayed
subplot(3,2,5),
h1 = plot(timei, mri_signal_model); hold on
h2 = plot(timei, bold_global);
title('General Linear Model')
xlabel('Time [s]'), legend([h1 h2], {'MRI Signal (GLM)','BOLD_{norm}'}), grid on


figure('name','CVR','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
subplot(2,3,1); imshow(bold(:,:,16,1), 'displayrange', [0 1000])
hold on, title('BOLD - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(2,3,2); imshow(imrotate(squeeze(bold(:,end/2,:,1)),90), 'displayrange', [0 1000])
hold on, title('BOLD - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
c3 = subplot(2,3,3); imshow(imrotate(squeeze(bold(end/2,:,:,1)),90), 'displayrange', [0 1000])
hold on, title('BOLD - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

c1 = subplot(2,3,4); imshow(CVR_brain(:,:,16), 'displayrange', [0 0.6]), colormap(c1, 'Parula')
hold on, title('CVR - Axial'), c = colorbar; c.Label.String = '[%\DeltaBOLD/mmHg_{CO2}]';
c2 = subplot(2,3,5); imshow(imrotate(squeeze(CVR_brain(:,end/2,:)),90), 'displayrange', [0 0.6]), colormap(c2, 'Parula')
hold on, title('CVR - Coronal'), c = colorbar; c.Label.String = '[%\DeltaBOLD/mmHg_{CO2}]';
c3 = subplot(2,3,6); imshow(imrotate(squeeze(CVR_brain(end/2,:,:)),90), 'displayrange', [0 0.6]), colormap(c3, 'Parula')
hold on, title('CVR - Sagittal'), c = colorbar; c.Label.String = '[%\DeltaBOLD/mmHg_{CO2}]';
