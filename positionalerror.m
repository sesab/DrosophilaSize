% Decoding the position and computing the level of precision with which 
% expression levels encode it.
% PNAS Positional information, in bits. Eq 11
% Get the value for the relative position

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
 
 gradmu=[];
 for i=1:numel(gNames)
     gradmu(:,i)=diff(media(:,i));
 end
 
 for i=2:numel(absolute)
    invsigma_abs(i)=gradmu(i-1,:)*inv(absolute(i).Cov)*gradmu(i-1,:)'; 
 end

 
 
% Compute the precision of the relative coding position
 media=[];
 for i=1:numel(relative)
    media(i,:)=relative(i).mu(:);
 end
 
 gradmu=[];
 for i=1:numel(gNames)
     gradmu(:,i)=diff(media(:,i));
 end
 
 for i=2:numel(relative)
    invsigma_rel(i)=gradmu(i-1,:)*inv(relative(i).Cov)*gradmu(i-1,:)'; 
 end

 figure(1)
 plot([1:200]/200,1./sqrt(invsigma_rel)/1000,'r')
 hold on
 plot([1:numel(absolute)]/200,1./sqrt(invsigma_abs)/1000,'k')
 xlabel('x/L')
 ylabel('sigma/L')
 set(gca,'FontSize',16,'Box','Off','TickDir','Out');
 axis([0. 1.2 0 .011])
 axis square
 legend('Relative','Absolute')



% for n=1:max(max(yy));
%     valy(n)=length(find(yy==n));
% end
% 
% 
% figure(1)
% 
% plot(nanmeang1_abs,nanvarg1_abs,'-',nanmeang1_rel,nanvarg1_rel,'-')
% 
% %figure(3)
% %v=nbin*[1:n];
% %plot(v,nanvarg1_rel./nanmeang1_rel(1:n))
% %hold on
% %plot(v,nanvarg1_abs(1:n)./nanmeang1_abs(1:n),'r')
% %hold off
% %figure(4)
% %v=mean(XX(:,[1:5:1000]));
% %plot(v,nanmeang1_rel,'b-',[v v],[nanmeang1_rel-nanstdg1_rel nanmeang1_rel+nanstdg1_rel],'b-')
% 
% %   this is very encouraging, let's look with error bars
% 
% for kk=1:20;
%     list = randperm(Nem);
%     list = list(1:round(Nem/2));
%     test = g1(list,:);
%     xtest = xx(list,:);
%     ytest = yy(list,:);
%     for n=1:max(max(xx));
%         [ii,jj] = find(xtest==n);
%         samples = [];
%         for k=1:length(ii);
%             samples = [samples test(ii(k),jj(k))];
%         end
%         Mg1_abs(kk,n) = nanmean(samples);
%         Vg1_abs(kk,n) = nanvar(samples);
%     end
%     for n=1:max(max(yy));
%         [ii,jj] = find(ytest==n);
%         samples = [];
%         for k=1:length(ii);
%             samples = [samples test(ii(k),jj(k))];
%         end
%         Mg1_rel(kk,n) = nanmean(samples);
%         Vg1_rel(kk,n) = nanvar(samples);
%     end
% end
% 
% figure(2)
% for n=1:length(Mg1_abs);
%     am = nanmean(Mg1_abs(:,n));
%     bm = nanstd(Mg1_abs(:,n));
%     av = nanmean(Vg1_abs(:,n));
%     bv = nanstd(Vg1_abs(:,n));
%     plot([am-bm am+bm],[av av],'b-',[am am],[av-bv av+bv],'b-')
%     hold on
% end
% for n=1:length(Mg1_rel);
%     cm = nanmean(Mg1_rel(:,n));
%     dm = nanstd(Mg1_rel(:,n));
%     cv = nanmean(Vg1_rel(:,n));
%     dv = nanstd(Vg1_rel(:,n));
%     plot([cm-dm cm+dm],[cv cv],'r-',[cm cm],[cv-dv cv+dv],'r-')
%     hold on
% end
% hold off
% label=sprintf('mean of %s expression level',genename);
% xlabel(label)
% label=sprintf('variance of %s expression level',genename);
% ylabel(label)
% set(gca,'FontSize',16,'Box','Off','TickDir','Out');
% limit=nanmax(cd+dv);
% axis([-0.1 1.1 0 .1])
% axis square
% 
% filename=fullfile('figs',sprintf('Fig%s.pdf',genename));
% print('-dpdf','-r200',filename); 
% 
% 



