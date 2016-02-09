% Decoding the position and computing the level of precision with which 
% expression levels encode it.
% PNAS Positional information, in bits. Eq 11
% Get the value for the relative position
clear all
close all
load('wt_130104_Kni_Kr_Gt_Hb_AP.mat')

FigFolder='figs';

% data = 
% 
% 1x243 struct array with fields:
% 
%     index
%     orient
%     dist
%     age
%     genotype
%     Kni
%     Kr
%     Gt
%     Hb
%     AP
%
%   I think all these have genotype = 1, which is the wild type
%   (this is part of a series of experiments with mutants, so 
%   they need to keep track ...)
%
ages = [data.age];
side = [data.orient];
rawg1 = vertcat(data.Hb);
rawg2 = vertcat(data.Kr);
rawg3 = vertcat(data.Gt);
rawg4 = vertcat(data.Kni);
lengths = [data.AP];

%
%   all of these embryos are in nuclear cycle 14, but things evolve with
%   time; we need to pick out those embryos in a relatively small window.
%   the clock is provided by the progress of the cellularization membrane
%   ...
%   also, you want to look only along one side of the embryo
idx = (ages>38 & ages<48) &  side==1;
Nem = sum(idx);
%   impressively, we seem to have 100 embryos here (should check with
%   Mariela that this makes sense)
%
%   we should normalize the expression levels of the genes
%   I'll look just at one, and you can play with the rest

g=struct('Hb',[],'Kr',[],'Gt',[],'Kni',[]);

g.Hb = rawg1(idx,:);
g.Kr = rawg2(idx,:);
g.Gt = rawg3(idx,:);
g.Kni= rawg4(idx,:);


gNames = fieldnames(g);
for loopIndex = 1:numel(gNames) 
    tmp=g.(gNames{loopIndex});
    offset1 = min(nanmean(tmp));
    range1 = max(nanmean(tmp))-min(nanmean(tmp));
    g.(gNames{loopIndex}) = (tmp-offset1)/range1;
end


LL = lengths(idx);
% %   nanstd(LL)/nanmean(LL) = 0.0293
% %   so we have +/- 3% fluctuations in length
% 
% %   the index for the array g1 is length in FRACTIONAL units of the egg
% %   length, so from [1:1000]/1000.  In generl we tend to look only in the
% %   middle 80% (101:900) just because special things happen at the ends,
% %   and because imaging at the ends can be problematic
% 
% %   But we could also make an x-axis in pixels (shoudl check units!)
% %   for each embryo
% 

XX = LL'*[1:1000]/1000;
% %   we could put these absolute positions into bins
% %   which are a bit bigger than the pixels
nbin=5;
xx = ceil(XX/nbin);

% 
% %   along this absolute axis, we could compute a mean and covariance of the
% %   genes expression level given the space
% 

temporary=struct('Hb',[],'Kr',[],'Gt',[],'Kni',[]);
absolute=repmat(struct('Cov',[],'mu',[]),1,1000/nbin);
 for n=min(min(xx)):max(max(xx));
    [ii,jj] = find(xx==n);  
    for loopIndex = 1:numel(gNames) 
        tmp=g.(gNames{loopIndex});  
         samples=[];
        for k=1:length(ii);
            samples = [samples tmp(ii(k),jj(k))];
        end
        temporary.(gNames{loopIndex})=samples;
    end
    C=[ temporary.(gNames{1});temporary.(gNames{2});temporary.(gNames{3});temporary.(gNames{4})];
    absolute(n).Cov=nancov(C');
    absolute(n).mu=nanmean(C'); 
 end

 
% %Check how many data for each n
% 
% for n=1:max(max(xx));
%     valx(n)=length(find(xx==n));
% end
% 
% %   can also do this along the scaled axes

yy= ones(Nem,1)*ceil([1:1000]/nbin);

temporary=struct('Hb',[],'Kr',[],'Gt',[],'Kni',[]);
relative=repmat(struct('Cov',[],'mu',[]),1,1000/nbin);
 for n=min(min(yy)):max(max(yy));
    [ii,jj] = find(yy==n);  
    for loopIndex = 1:numel(gNames) 
        tmp=g.(gNames{loopIndex});  
         samples=[];
        for k=1:length(ii);
            samples = [samples tmp(ii(k),jj(k))];
        end
        temporary.(gNames{loopIndex})=samples;
    end
    C=[temporary.(gNames{1});temporary.(gNames{2});temporary.(gNames{3});temporary.(gNames{4})];
    relative(n).Cov=nancov(C');
    relative(n).mu=nanmean(C'); 
 end
 
 
  
% Compute the precision of the absolute coding position
 media=[];
 for i=1:numel(absolute)
    media(i,:)=absolute(i).mu(:);
 end
 media(i+1,:)=media(i,:);
 gradmu=[];
 for i=1:numel(gNames)
     gradmu(:,i)=diff(media(:,i));
 end
 
 for i=1:numel(absolute)
    invsigma_abs(i)=gradmu(i,:)*inv(absolute(i).Cov)*gradmu(i,:)'; 
 end

 
 
% Compute the precision of the relative coding position
 media=[];
 for i=1:numel(relative)
    media(i,:)=relative(i).mu(:);
 end
 media(i+1,:)=media(i,:);
 gradmu=[];
 for i=1:numel(gNames)
     gradmu(:,i)=diff(media(:,i));
 end
 invsigma_rel=[];
 for i=1:numel(relative)
    invsigma_rel(i)=gradmu(i,:)*inv(relative(i).Cov)*gradmu(i,:)'; 
 end

 
 %Compute positional error for each gene
for n=1:4
    for i=1:numel(relative)
        invSinglesigma_rel(i,n)=gradmu(i,n)/relative(i).Cov(n,n)*gradmu(i,n);
    end
    
    for i=1:numel(relative)
        invSinglesigma_abs(i,n)=gradmu(i,n)/absolute(i).Cov(n,n)*gradmu(i,n);
    end
end
 colors=colormap(summer(4));
 len=1000/nbin;
 figure(1)
 plot([1:len]/len,1./sqrt(invsigma_rel)/1000,'r')
 hold on
 plot([1:numel(absolute)]/len,1./sqrt(invsigma_abs)/1000,'b')
 plot([1:numel(absolute)]/len,0.01*ones(numel(absolute)),'k')
 plot([1:len]/len,1./sqrt(invSinglesigma_rel(:,1))/1000,'color',colors(1,:),'linestyle','--')
 plot([1:len]/len,1./sqrt(invSinglesigma_abs(:,1))/1000,'color',colors(3,:),'linestyle','--')
 %plot([1:len-1]/len,1./sqrt(invSinglesigma_rel(:,2))/1000,'color',colors(2,:),'linestyle','--')
 % plot([1:len-1]/len,1./sqrt(invSinglesigma_rel(:,3))/1000,'color',colors(3,:),'linestyle','--')
 %  plot([1:len-1]/len,1./sqrt(invSinglesigma_rel(:,4))/1000,'color',colors(4,:),'linestyle','--')
 xlabel('x/L')
 ylabel('sigma/L')
 set(gca,'FontSize',16,'Box','Off','TickDir','Out');
 axis([0. 1.2 0 .1])

 legend('Relative','Absolute','Hb rel','Hb abs')



