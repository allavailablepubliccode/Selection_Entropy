% Generate Fig. 3

clear;clc;close all;

sigr        = 1:0.0001:2;   % range of standard deviations
states      = 100;          % number of bins
x           = -10:0.001:10; % range for gaussians 

for ii      = 1:numel(sigr)

    if mod(ii,100)==0
        disp([num2str(ii) '/' num2str(numel(sigr))])
    end

    % standard deviations for the two gaussians
    sig1    = sigr(ii);
    sig2    = sigr(end-ii+1);

    % generate probability distributions
    pn      = normpdf(x,0,sig1);
    pn      = histcounts(pn,states);
    pn      = pn/sum(pn);

    pm      = normpdf(x,0,sig2);
    pm      = histcounts(pm,states);
    pm      = pm/sum(pm);

    KL(ii)  = sum(pn.*log2(pn./pm));    % KL divergence

    % find values that give real values for selection entropy
    x1      = find(pn>pm);
    x2      = find(pn<pm);
    nums(ii) = numel(x1)/numel(x2);

    % calculate the two flavours of selection entropy
    % i.e. select N from M, and M from N
    SE1(ii) = sum(pm(x1).*log2(pn(x1)./pm(x1)-1)-pn(x1).*log2(1-pm(x1)./pn(x1)));
    SE2(ii) = sum(pn(x2).*log2(pm(x2)./pn(x2)-1)-pm(x2).*log2(1-pn(x2)./pm(x2)));

end
% save('lub.mat','KL','SE1','SE2')
% 
% load('lub.mat')

% plot results
figure
plot3(sigr-1.5, zeros(size(KL)), normalize(KL,'range'), 'k','Linewidth',1.5)
hold on
plot3(sigr-1.5, zeros(size(KL)), normalize(SE1,'range'), 'k','Linewidth',5)
hold off
axis tight
grid on
grid minor
view(43,56)
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'~/Desktop/lub1.eps');
% close all

x(x<-6)=[];
x(x>6)=[];

bb = 200;
bb = [ones(bb,1) ; zeros(bb,1)];
bb = repmat(bb,round(size(x,2)/numel(bb)),1);
bb(numel(bb)+1) = 0;

figure
hold on
c = -1;
for ii      = [1, floor(numel(sigr)/2), numel(sigr)]
c = c + 0.5;

    sig1    = sigr(ii);
    sig2    = sigr(end-ii+1);

    pn      = normpdf(x,0,sig1);
    pm      = normpdf(x,0,sig2);

    pn      = pn.*bb';

    pn(pn==0) = nan;
    pm(pm==0) = nan;

    plot3(c*ones(size(pn)), x, pn, 'k.')
    plot3(c*ones(size(pn)), x, pm, 'k.')

end
hold off
view(43,56)
axis tight
grid on 
grid minor
% set(gcf,'Renderer','Painter')
% hgexport(gcf,'~/Desktop/lub2.eps');
% close all
