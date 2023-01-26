% Generate Fig. 4

% requires BrainSpace
% (https://brainspace.readthedocs.io/en/latest/index.html)

% requires function smooth2a
% https://uk.mathworks.com/matlabcentral/fileexchange/23287-smooth2a

clear;clc;close all;

% load HCP data

% ddir   = '~/Dropbox/HCP/rest/';
% subs   = {'rest_hcp100307_ts';...
%     'rest_hcp100408_ts';...
%     'rest_hcp101107_ts';...
%     'rest_hcp101309_ts';...
%     'rest_hcp101915_ts';...
%     'rest_hcp103111_ts';...
%     'rest_hcp103414_ts';...
%     'rest_hcp103818_ts';...
%     'rest_hcp105014_ts';...
%     'rest_hcp105115_ts';...
%     'rest_hcp106016_ts';...
%     'rest_hcp108828_ts';...
%     'rest_hcp110411_ts';...
%     'rest_hcp111312_ts';...
%     'rest_hcp111716_ts';...
%     'rest_hcp113922_ts';...
%     'rest_hcp114419_ts';...
%     'rest_hcp115320_ts';...
%     'rest_hcp116524_ts';...
%     'rest_hcp117122_ts';...
%     'rest_hcp118528_ts';...
%     'rest_hcp118730_ts';...
%     'rest_hcp118932_ts';...
%     'rest_hcp120111_ts';...
%     'rest_hcp122317_ts';...
%     'rest_hcp122620_ts';...
%     'rest_hcp123117_ts';...
%     'rest_hcp124422_ts';...
%     'rest_hcp125525_ts'};
% z      = [];
% for ii = 1:numel(subs)
%     disp(['loading data, ' num2str(round(ii*100/numel(subs))) '% complete'])
%     zt = load([ddir subs{ii}]);
%     z  = cat(3,z,zt);
% end
% save('HCP_data.mat','z')

load('HCP_data.mat','z')
states = 10;
grpin = 10;

for ii = 1:size(z,3)

    disp(['generating Fig. 4, ' num2str(round(ii*100/size(z,3))) '% complete'])

    z1 = z(:,:,ii);

    for jj = 1:size(z1,2)-grpin+1

        % first probability distribution
        z2 = z1(:,jj:jj+grpin-1);
        n  = histcounts(z2(:),states);
        N  = sum(n);
        pn = n/N;

        for kk = 1:size(z1,2)-grpin+1

            % second probability distribution
            z3 = z1(:,kk:kk+grpin-1);
            m  = histcounts(z3(:),states);
            M  = sum(m);
            pm = m/M;

            % KL divergence
            KLt = pn.*log2(pn./pm);
            KLt(find(isinf(KLt))) = nan;
            KL(jj,kk,ii)          = nansum(KLt);

            x1      = find(pn>pm);
            x2      = find(pn<pm);

            nums(jj,kk,ii) = numel(x1)/numel(x2);

            % Selection entropy
            SE1(jj,kk,ii) = nansum(pm(x1).*log2(pn(x1)./pm(x1)-1)-pn(x1).*log2(1-pm(x1)./pn(x1)));
            SE2(jj,kk,ii) = nansum(pn(x2).*log2(pm(x2)./pn(x2)-1)-pm(x2).*log2(1-pn(x2)./pm(x2)));

        end
    end
end

KLm  = nanmean(KL,3);
SE1m = nanmean(SE1,3);
SE2m = nanmean(SE2,3);

KLm(logical(eye(size(KLm))))   = nan;
SE1m(logical(eye(size(SE1m)))) = nan;
SE2m(logical(eye(size(SE2m)))) = nan;

SEm = SE1m + SE2m;

smthwin = 10;

figure
subplot(2,2,1)
KLm2 = smooth2a(KLm,smthwin);
imagesc2(KLm2)
colormap gray
colorbar
title('KL')

subplot(2,2,2)
SEm2 = smooth2a(SEm,smthwin);
imagesc2(SEm2)
colormap gray
colorbar
title('SE')

N = 2^4;
t = 2^4;

subplot(2,2,3)
datkl = nanmean(KLm2);
datkl = smoothdat(datkl,N,t);
plot(datkl);
title('KL')

subplot(2,2,4)
datse = nanmean(SEm2);
datse = smoothdat(datse,N,t);
plot(datse);
title('SE')
axis tight

nums = 10;
c = 1;
for ii = 1:nums

    disp([num2str(ii) ', ' num2str(c:c+100/nums-1)])
    
    dat = ones(1,100);
    dat(c:c+100/nums-1) = 0;
    addpath(genpath('~/Dropbox/BrainSpace/matlab/'))
    [surf_lh, surf_rh] = load_conte69();
    labeling           = load_parcellation('schaefer',100);
    conn_matices       = load_group_fc('schaefer',100);
    plot_hemispheres(dat', {surf_lh,surf_rh},'parcellation',labeling.schaefer_100,'views','ms');
    c = c + 100/nums;
end

function B = smoothdat(B,N,t)
for ii = 1:N
    B = movmean(B,t);
end
end
