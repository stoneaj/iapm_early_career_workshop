% demonstration of TRUST for measuring cerebral oxygenation for IAPM Early Careers Workshop 2019
% Data is from:
% My own brain!
% Alan Stone, TCD, 23/09/2019

% directory of tutorial toolbox
tutorial_directory = '/Users/astone/albin/iapm_early_career_workshop';

% read asl data into matlab
load([tutorial_directory filesep 'data' filesep 'trust' filesep 'trust.mat']);
trust = squeeze(trust);

% subtract tag-controls
dS_trust = trust(:,:,1:2:end) - trust(:,:,2:2:end);

% average eTE's
dS_trust_mean = cat(3, mean(dS_trust(:,:,1:5:end),3), mean(dS_trust(:,:,2:5:end),3),...
                       mean(dS_trust(:,:,3:5:end),3), mean(dS_trust(:,:,4:5:end),3),...
                       mean(dS_trust(:,:,5:5:end),3));

% create superior saggital sinus mask
dS_trust_eTE1 = dS_trust_mean(:,:,1);
[val ind] = sort(dS_trust_eTE1(:),'descend');
mask_SSS = dS_trust_eTE1 >= val(3);

% View
figure('name','TRUST - Tag-Control difference','NumberTitle','off'), set(gcf,'color','w'), hold on % set figure name & other settings
for id = 1:5
    subplot(1,6,id), imshow(dS_trust_mean(:,:,id), 'displayrange', [0 300])
    hold on, title(sprintf('eTE_{%1.0f}',id))
end

subplot(1,6,6), imshow(mask_SSS(:,:,1), 'displayrange', []); title('SSS Mask')

% extract SSS signal for each eTE
eTE = [0,40,80,120,160]; % [ms] effective echo time

for id = 1:length(eTE)
    dS(id) = mean(nonzeros(mask_SSS .* dS_trust_mean(:,:,id)));
end

% fit exponential to SSS signal
f = fit(eTE',dS','exp1');

figure('name','TRUST O2','NumberTitle','off'), set(gcf,'color','w'), hold on
plot(f,eTE,dS,'o')
xlabel('eTE [ms]')
ylabel('MR signal in Saggital Sinus [a.u.]')

% calculate T2 of blood using fitted exponential
t1b = 1624; % [ms] @ 3T ... t1 of blood assumed to be 1624ms
t2b = 1/((1/t1b)-(f.b));

% Convert to Venous Oxygenation
% constants taken from "Lu et al. MRM 67:42 (2012)"
%  these are values specifically for tau_cpmg = 10 ms
%  Hct between
a1 = -13.5; a2 = 80.2; a3 = -75.9; b1 = -0.5; b2 = 3.4; c1 = 247.4; % [s-1]
hct = 0.42; % the ratio of the volume of red blood cells to the total volume of blood

A = a1 + a2*hct + a3*hct^2;
B = b1*hct + b2*hct^2;
C = c1*hct*(1 - hct);

r = roots([C B A-(1/(t2b/1000))]);
x = r(r>=0);

if length(x) == 1
    Y = 1-x;
else
    fprintf('Root of quadratic equation problematic')
end

% print results
fprintf('TE used [ms] ... \n')
disp(eTE)
fprintf('T2b = %0.1f [ms] \n',t2b)
fprintf('Assuming Hct of 0.42 \n')
fprintf('Y = %0.3f [%%] \n',Y)
