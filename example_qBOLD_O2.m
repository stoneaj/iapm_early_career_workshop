% demonstration of qBOLD imaging for mapping R2' & cerebral O2 for IAPM Early Careers Workshop 2019
% Data is from:
% Stone, A. J., & Blockley, N. P. (2016).
% Data acquired to demonstrate a streamlined approach to mapping and quantifying brain oxygenation using quantitative BOLD.
% University of Oxford.
% Alan Stone, TCD, 23/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop';

% read asl data into matlab
load([tutorial_directory filesep 'data' filesep 'qbold' filesep 'qbold.mat']);

% mean brain signal time-course
qbold_global = squeeze(mean(mean(mean(qbold,1),2),3));
tau = -28:4:64; % [ms] spin-echo offset
tau = tau.*10^-3; % Convert tau to seconds

% view global signal
figure, hold on
plot(tau, qbold_global)
title('Global qBOLD')
xlabel('\tau [ms]')
ylabel('[a.u.]'), grid on

% Y
ln_Sase = log(qbold);
ln_Sase(isnan(ln_Sase)) = 0;
ln_Sase(isinf(ln_Sase)) = 0;

Tc = 0.015; % cutoff time for monoexponential regime [s]
tau_lineID = find(tau > Tc); % tau's to be used for R2' fitting
s0_id = find(tau == 0);
tau_lineID = [s0_id tau_lineID];
w = 1./tau(1,tau_lineID)'; % weightings for lscov
w(1) = w(2);
[x,y,z,t] = size(qbold);
p = zeros(x,y,z,3);
stdp = zeros(x,y,z,3);

% Fit qBOLD model to voxel-wise data
for xID = 1:x
    for yID = 1:y
        for zID = 1:z
            %fprintf('xID%d ; yID%d ; zID%d \n',xID,yID,zID)
            X = [ones(length(tau(1,tau_lineID)'),1) -tau(1,tau_lineID)' ones(length(tau(1,tau_lineID)'),1)];
            X(1,3) = 0;
            Y = squeeze(ln_Sase(xID,yID,zID,tau_lineID));
            % least squares solution to the linear system of equations
            [p(xID,yID,zID,:),stdp(xID,yID,zID,:),~,S(xID,yID,zID,:,:)]= lscov(X,Y,w);
        end
    end
end

% Constants
Hct = 0.34; % hct ratio in small vessels
dChi0 = 0.264*10^-6; % ppm, sus difference between fully oxy & deoxy rbc's
gamma = 2.675*10^4; % rads/(secs.Gauss) gyromagnetic ratio
B0 = 3*10^4; %Gauss, Field strength
phi = 1.34; % mlO2/gHb
k = 0.03; % conversion factor Hct (% rbc's in blood) to Hb (iron-containing molecules in the rbc's used to transport o2)

% Calculate Physiological Parameters
r2p = p(:,:,:,2); % reversible transverse relaxation rate
c = p(:,:,:,1);
dbv = p(:,:,:,3); % deoxygenated blood volume
oef = r2p./(dbv.*gamma.*(4./3).*pi.*dChi0.*Hct.*B0); % Oxygen Extraction Fraction
dhb = r2p./(dbv.*gamma.*(4./3).*pi.*dChi0.*B0.*k); % DeOxyheamoglbin concentration

oef(isinf(oef)) = 0;
dhb(isinf(dhb)) = 0;

% Error calculation dHb
cov_r2p_dbv = S(:,:,:,3,2);
std_dhb = abs(r2p./dbv).*sqrt( (stdp(:,:,:,2)./r2p).^2 + (stdp(:,:,:,3)./dbv).^2 - 2.*(cov_r2p_dbv./(r2p.*dbv)) );
std_dhb = std_dhb .* (1 ./ (gamma.*(4./3).*pi.*dChi0.*B0.*k));

% make crude mask
mask = qbold(:,:,:,s0_id) > 100;
figure, imshow(mask(:,:,end/2),'displayrange',[])
r2p = r2p .* mask;
dbv = dbv .* mask;
oef = oef .* mask;

% view qBOLD maps
figure('name','qBOLD','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings

subplot(4,3,1), imshow(qbold(:,:,end/2,s0_id), 'displayrange', [0 200])
hold on, title('Spin Echo - Axial'), c = colorbar; c.Label.String = '[A.U.]';
subplot(4,3,2), imshow(imrotate(squeeze(qbold(:,end/2,:,s0_id)),90), 'displayrange', [0 200])
hold on, title('Spin Echo - Coronal'), c = colorbar; c.Label.String = '[A.U.]';
subplot(4,3,3), imshow(imrotate(squeeze(qbold(end/2,:,:,s0_id)),90), 'displayrange', [0 200])
hold on, title('Spin Echo - Sagittal'), c = colorbar; c.Label.String = '[A.U.]';

c1 = subplot(4,3,4); imshow(r2p(:,:,end/2), 'displayrange', [0 10]), colormap(c1, 'Parula')
hold on, title('R_2'' - Axial'), c = colorbar; c.Label.String = '[s-1]';
c2 = subplot(4,3,5); imshow(imrotate(squeeze(r2p(:,end/2,:)),90), 'displayrange', [0 10]), colormap(c2, 'Parula')
hold on, title('R_2'' - Coronal'), c = colorbar; c.Label.String = '[s-1]';
c3 = subplot(4,3,6); imshow(imrotate(squeeze(r2p(end/2,:,:)),90), 'displayrange', [0 10]), colormap(c3, 'Parula')
hold on, title('R_2'' - Sagittal'), c = colorbar; c.Label.String = '[s-1]';

c1 = subplot(4,3,7); imshow(dbv(:,:,end/2), 'displayrange', [0 0.12]), colormap(c1, 'Parula')
hold on, title('DBV - Axial'), c = colorbar; c.Label.String = '[%]';
c2 = subplot(4,3,8); imshow(imrotate(squeeze(dbv(:,end/2,:)),90), 'displayrange', [0 0.12]), colormap(c2, 'Parula')
hold on, title('DBV - Coronal'), c = colorbar; c.Label.String = '[%]';
c3 = subplot(4,3,9); imshow(imrotate(squeeze(dbv(end/2,:,:)),90), 'displayrange', [0 0.12]), colormap(c3, 'Parula')
hold on, title('DBV - Sagittal'), c = colorbar; c.Label.String = '[%]';

c1 = subplot(4,3,10); imshow(oef(:,:,end/2), 'displayrange', [0 1]), colormap(c1, 'Parula')
hold on, title('OEF - Axial'), c = colorbar; c.Label.String = '[%]';
c2 = subplot(4,3,11); imshow(imrotate(squeeze(oef(:,end/2,:)),90), 'displayrange', [0 1]), colormap(c2, 'Parula')
hold on, title('OEF - Coronal'), c = colorbar; c.Label.String = '[%]';
c3 = subplot(4,3,12); imshow(imrotate(squeeze(oef(end/2,:,:)),90), 'displayrange', [0 1]), colormap(c3, 'Parula')
hold on, title('OEF - Sagittal'), c = colorbar; c.Label.String = '[%]';

%
% c1 = subplot(3,3,7); imshow(CBF(:,:,end/2), 'displayrange', [0 100]), colormap(c1, 'Jet')
% hold on, title('CBF - Axial'), c = colorbar; c.Label.String = '[ml/100g/min]';
% c2 = subplot(3,3,8); imshow(imrotate(squeeze(CBF(:,end/2,:)),90), 'displayrange', [0 100]), colormap(c2, 'Jet')
% hold on, title('CBF - Coronal'), c = colorbar; c.Label.String = '[ml/100g/min]';
% c3 = subplot(3,3,9); imshow(imrotate(squeeze(CBF(end/2,:,:)),90), 'displayrange', [0 100]), colormap(c3, 'Jet')
% hold on, title('CBF - Sagittal'), c = colorbar; c.Label.String = '[ml/100g/min]';
