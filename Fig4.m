% Generate Fig. 4

% requires BrainSpace
% (https://brainspace.readthedocs.io/en/latest/index.html)

clear;clc;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ddir   = '~/Dropbox/HCP/rest/';
% gdir   = '~/Dropbox/Schaeffer/';
% gr     = 3;
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
% g      = load([gdir '/' num2str(gr)]);
% z      = [];
% for ii = 1:numel(subs)
%     ii
%     zt = load([ddir subs{ii}]);
%     z  = cat(3,z,zt);
% end
% save('zdat_gr3.mat','z')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('zdat_gr3.mat','z')

states = 100;

for ii = 1:size(z,3)
    ii

    z1 = z(:,:,ii);

    for jj = 1:size(z1,2)

        z2 = z1(:,jj);
        n  = histcounts(z2(:),states);
        N  = sum(n);
        pn = n/N;

        for kk = 1:size(z1,2)

            z3 = z1(:,kk);
            m  = histcounts(z3(:),states);
            M  = sum(m);
            pm = m/M;

            KLt = pn.*log2(pn./pm);
            KLt(find(isinf(KLt))) = nan;
            KL(jj,kk,ii)          = nansum(KLt);

            x1      = find(pn>pm);
            x2      = find(pn<pm);

            nums(jj,kk,ii) = numel(x1)/numel(x2);

            SE1(jj,kk,ii) = nansum(pm(x1).*log2(pn(x1)./pm(x1)-1)-pn(x1).*log2(1-pm(x1)./pn(x1)));
            SE2(jj,kk,ii) = nansum(pn(x2).*log2(pm(x2)./pn(x2)-1)-pm(x2).*log2(1-pn(x2)./pm(x2)));

        end
    end
end
% save('SEdat_gr3.mat','KL','SE1','SE2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('SEdat_gr3.mat')

KLm  = nanmean(KL,3);
SE1m = nanmean(SE1,3);
SE2m = nanmean(SE2,3);

KLm(logical(eye(size(KLm))))   = nan;
SE1m(logical(eye(size(SE1m)))) = nan;
SE2m(logical(eye(size(SE2m)))) = nan;

SEm = SE1m;

figure
imagesc2(KLm)
colormap gray
colorbar
title('KL')

figure
imagesc2(SEm)
colormap gray
colorbar
title('SE')

N = 5;

figure;
dat = nanmean(KLm);
dat = smoothlub(dat,N);
[r_kl,p_kl] = corr((1:numel(dat))',dat');
plot(dat);
title('KL')

addpath(genpath('~/Dropbox/BrainSpace/matlab/'))
[surf_lh, surf_rh] = load_conte69();
labeling           = load_parcellation('schaefer',100);
conn_matices       = load_group_fc('schaefer',100);
plot_hemispheres(dat', {surf_lh,surf_rh},'parcellation',labeling.schaefer_100,'views','lmisap');
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'~/Desktop/kl_0_90.eps');
% close all

figure;
dat = nanmean(SEm);
dat = smoothlub(dat,N);
[r_se,p_se] = corr((1:numel(dat))',dat');
plot(dat);
title('se')
axis tight

addpath(genpath('~/Dropbox/BrainSpace/matlab/'))
[surf_lh, surf_rh] = load_conte69();
labeling           = load_parcellation('schaefer',100);
conn_matices       = load_group_fc('schaefer',100);
plot_hemispheres(dat', {surf_lh,surf_rh},'parcellation',labeling.schaefer_100,'views','lmisap');
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'~/Desktop/90_0.eps');
% close all
% [0 90], [-90 0], [90, 0]

figure
dat = squeeze(z(:,1,2));
dat = normalize(dat,'range');
plot(dat)
title('raw','k')

figure
[bincounts,edg] = histcounts(dat,100);
bincounts = bincounts/sum(bincounts);
histogram('BinCounts', bincounts, 'BinEdges', edg);

% c = 1;
% for ii = 1:4
%     ii
%     dat = ones(1,100);
%     dat(c:c+24) = 0;
%     addpath(genpath('~/Dropbox/BrainSpace/matlab/'))
%     [surf_lh, surf_rh] = load_conte69();
%     labeling           = load_parcellation('schaefer',100);
%     conn_matices       = load_group_fc('schaefer',100);
%     plot_hemispheres(dat', {surf_lh,surf_rh},'parcellation',labeling.schaefer_100,'views','ms');
% %     set(gcf,'Renderer','Painter')
% %     hgexport(gcf,['~/Desktop/lub' num2str(ii) '.eps']);
% %     close all
%     c = c + 25;
% end
function B = smoothlub(B,N)
for ii = 1:N
    B = movmean(B,4);
end
end
