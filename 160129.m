load('wt_130104_Kni_Kr_Gt_Hb_AP.mat')


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

g1 = rawg1(idx,:);
offset1 = min(mean(g1));
range1 = max(mean(g1))-min(mean(g1));
g1 = (g1-offset1)/range1;

%   we also need the lengths of the embryos
LL = lengths(idx);
%   std(LL)/mean(LL) = 0.0293
%   so we have +/- 3% fluctuations in length

%   the index for the array g1 is length in FRACTIONAL units of the egg
%   length, so from [1:1000]/1000.  In general we tend to look only in the
%   middle 80% (101:900) just because special things happen at the ends,
%   and because imaging at the ends can be problematic

%   But we could also make an x-axis in pixels (shoudl check units!)
%   for each embryo

XX = LL'*[1:1000]/1000;
%   we could put these absolute positions into bins
%   which are a bit bigger than the pixels
xx = round(XX/5);

%   along this absolute axis, we could compute a mean and variance of the
%   Hb expression level, at least as a start ...

for n=1:max(max(xx));
    [ii,jj] = find(xx==n);
    samples = [];
    for k=1:length(ii);
        samples = [samples g1(ii(k),jj(k))];
    end
    meang1_abs(n) = mean(samples);
    varg1_abs(n) = var(samples);
end



%   can also do this along the scaled axes
yy = ones(Nem,1)*round([1:1000]/5);
for n=1:max(max(yy));
    [ii,jj] = find(yy==n);
    samples = [];
    for k=1:length(ii);
        samples = [samples g1(ii(k),jj(k))];
    end
    meang1_rel(n) = mean(samples);
    varg1_rel(n) = var(samples);
end

figure(1)
plot(meang1_abs,varg1_abs,'-',meang1_rel,varg1_rel,'-')

%   this is very encouraging, let's look with error bars

for kk=1:20;
    list = randperm(Nem);
    list = list(1:round(Nem/2));
    test = g1(list,:);
    xtest = xx(list,:);
    ytest = yy(list,:);
    for n=1:max(max(xx));
        [ii,jj] = find(xtest==n);
        samples = [];
        for k=1:length(ii);
            samples = [samples test(ii(k),jj(k))];
        end
        Mg1_abs(kk,n) = mean(samples);
        Vg1_abs(kk,n) = var(samples);
    end
    for n=1:max(max(yy));
        [ii,jj] = find(ytest==n);
        samples = [];
        for k=1:length(ii);
            samples = [samples test(ii(k),jj(k))];
        end
        Mg1_rel(kk,n) = mean(samples);
        Vg1_rel(kk,n) = var(samples);
    end
end

figure(2)
for n=1:length(Mg1_abs);
    am = mean(Mg1_abs(:,n));
    bm = std(Mg1_abs(:,n));
    av = mean(Vg1_abs(:,n));
    bv = std(Vg1_abs(:,n));
    plot([am-bm am+bm],[av av],'b-',[am am],[av-bv av+bv],'b-')
    hold on
end
for n=1:length(Mg1_rel);
    cm = mean(Mg1_rel(:,n));
    dm = std(Mg1_rel(:,n));
    cv = mean(Vg1_rel(:,n));
    dv = std(Vg1_rel(:,n));
    plot([cm-dm cm+dm],[cv cv],'r-',[cm cm],[cv-dv cv+dv],'r-')
    hold on
end
hold off
xlabel('mean Hb expression level')
ylabel('variance of Hb expression level')
set(gca,'FontSize',16,'Box','Off','TickDir','Out');
axis([-0.1 1.1 0 0.06])
axis square

