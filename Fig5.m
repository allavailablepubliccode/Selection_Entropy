% Generate Fig. 5

% requires function smooth2a
% https://uk.mathworks.com/matlabcentral/fileexchange/23287-smooth2a

clear;clc;close all;

% load HCP data

% gr = importdata('/Users/erik/Library/CloudStorage/Dropbox/Margulies/1');
% [a,b] = sort(gr);
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
%     zt = zt(:,b);
%     zt = detrend(zt,1,50:50:1200);
%     z  = cat(3,z,zt);
% end
% z_1 = z(:,:,1:15);
% z_2 = z(:,:,16:end);
% save('HCP_data_1.mat','z_1')
% save('HCP_data_2.mat','z_2')

load('HCP_data_1.mat')
load('HCP_data_2.mat')
z = cat(3,z_1,z_2);
clear z_1 z_2
states = 100;

for ii = 1:size(z,3)

    disp(['generating Fig. 5, ' num2str(round(ii*100/size(z,3))) '% complete'])

    z1 = z(:,:,ii);

    for jj = 1:size(z1,2)

        % first probability distribution
        z2 = z1(:,jj);
        n  = histcounts(z2(:),states);
        N  = sum(n);
        pn = n/N;

        for kk = 1:size(z1,2)

            % second probability distribution
            z3 = z1(:,kk);
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

smthwin = 2^2;

figure
subplot(3,2,3)
KLm2 = smooth2a(KLm,smthwin);
imagesc2(KLm2)
colormap gray
colorbar
title('KL')

subplot(3,2,4)
SEm2 = smooth2a(SEm,smthwin);
imagesc2(SEm2)
colormap gray
colorbar
title('SE')

for ii = 1:99
    a_se_m(ii) = nanmean(diag(SEm,ii));
    a_kl_m(ii) = nanmean(diag(KLm,ii));
    b_se_m(ii) = nanmean(diag(SEm,-ii));
    b_kl_m(ii) = nanmean(diag(KLm,-ii));

    a_se_s(ii) = nanstd(diag(SEm,ii))/sqrt(99-ii+1);
    a_kl_s(ii) = nanstd(diag(KLm,ii))/sqrt(99-ii+1);
    b_se_s(ii) = nanstd(diag(SEm,-ii))/sqrt(99-ii+1);
    b_kl_s(ii) = nanstd(diag(KLm,-ii))/sqrt(99-ii+1);
end
c_se_m = mean([a_se_m;b_se_m]);
c_kl_m = mean([a_kl_m;b_kl_m]);

c_se_s = mean([a_se_s;b_se_s]);
c_kl_s = mean([a_kl_s;b_kl_s]);

subplot(3,2,6)
errorbar(c_se_m,c_se_s)
title('se')
subplot(3,2,5)
errorbar(c_kl_m,c_kl_s)
title('kl')

subplot(3,2,1)
dat = squeeze(z(:,1,1));
dat = normalize(dat,'range');
plot(dat)
title('raw','k')

subplot(3,2,2)
[bincounts,edg] = histcounts(dat,states);
bincounts = bincounts/sum(bincounts);
histogram('BinCounts', bincounts, 'BinEdges', edg);

% requires BrainSpace
% (https://brainspace.readthedocs.io/en/latest/index.html)

% nums = 4;
% c = 1;
% for ii = 1:nums
% 
%     disp([num2str(ii) ', ' num2str(c:c+100/nums-1)])
%     
%     dat = ones(1,100);
%     dat(b(c:c+100/nums-1)) = 0;
%     addpath(genpath('~/Dropbox/BrainSpace/matlab/'))
%     [surf_lh, surf_rh] = load_conte69();
%     labeling           = load_parcellation('schaefer',100);
%     conn_matices       = load_group_fc('schaefer',100);
%     plot_hemispheres(dat', {surf_lh,surf_rh},'parcellation',labeling.schaefer_100,'views','ms');
%     c = c + 100/nums;
% 
% end
