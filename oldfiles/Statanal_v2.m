%% final stat analysis
% writeen by PA, Jan. 2019

%% Load the data
load('WholeComp3.mat')

% import the valid directions
ecc = EC;
ecc(ecc == 0) = nan;
%% Collision rate
% Plot averaged + individuals for final resutls
figure(100); subplot(1,2,1);hold on
for i = 1:18
    plot(1:3,[NCc(2,1,i) NCc(1,1,i) NCc(3,1,i)],':k','LineWidth',0.5)
end
tt = zeros(6,2);
for i = 1:6
    tt(i,:) = mean(NCc(i,1,:),3);
    tt(i,2) = std(NCc(i,1,:));
end
errorbar(1:3,[ tt(2,1) tt(1,1) tt(3,1)],[ tt(2,2) tt(1,2) tt(3,2)],':k','LineWidth',2.5)
xlabel('Head Rotations (deg)'); ylabel('Number of collisions')
set(gca,'XTick',1:3)
set(gca, 'XTickLabels',{'HR = CCW' , 'HR = 0' , 'HR = CW'})
xlim([0.5,3.5])

figure(100); subplot(1,2,2);hold on
for i = 1:18
    plot(1:3,[NCc(2+3,1,i) NCc(1+3,1,i) NCc(3+3,1,i)],'-k','LineWidth',0.5)
end
tt = zeros(6,2);
for i = 1:6
    tt(i,:) = mean(NCc(i,1,:),3);
    tt(i,2) = std(NCc(i,1,:));
end
errorbar(1:3,[ tt(2+3,1) tt(1+3,1) tt(3+3,1)],[ tt(2+3,2) tt(1+3,2) tt(3+3,2)],':k','LineWidth',2.5)
xlabel('Head Rotations (deg)'); ylabel('Number of collisions')
set(gca,'XTick',1:3)
set(gca, 'XTickLabels',{'HR = CCW' , 'HR = 0' , 'HR = CW'})
xlim([0.5,3.5])

% stat analysis
[h1,p1,~,stats1] = ttest(NCc(4,1,:),NCc(5,1,:));
[h2,p2,~,stats2] = ttest(NCc(4,1,:),NCc(6,1,:));
[h3,p3,~,stats3] = ttest(NCc(4,1,:),mean(NCc(5:6,1,:),1));

% repeated measure anova
varNames = cell(6,1);
for i = 1 : 6
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
HRs = {'HR0';'HRL';'HRR';'HR0';'HRL';'HRR'}; % head roll conditions
VFs = {'Y';'Y';'Y';'N';'N';'N'}; % Visual feedback conditions

factorNames = {'HRs','VisualFeedback'};
within = table(HRs, VFs, 'VariableNames', factorNames);
cols = reshape(NCc(:,1,:),6,18);
t = array2table(cols', 'VariableNames',varNames);
rmcol = fitrm(t,'V1-V6~1','WithinDesign',within); % overal fit
[ranovatbcol] = ranova(rmcol,  'WithinModel','HRs*VisualFeedback');
mrm1 = multcompare(rmcol,'HRs','By','VisualFeedback','ComparisonType','bonferroni');

figname = 'Collisions2-averaged and individuals'; 
savefig(gcf,figname)
saveas(gcf,figname,'tiff')
close all

%% Variability of the movement
% Plot averages and individuals (Final results)    
% averaged for HRs
mvar = zeros(2,6);
for i = 1 : 6
    mvar(1,i) = nanmean(reshape(ecc.*variability(i,:,:,:),1,[]));
    mvar(2,i) = nanstd(reshape(ecc.*variability(i,:,:,:),1,[]));
end
mvars = zeros(18,6);    
for i = 1 : 6
    mvars(:,i) = nanmean(reshape(ecc.*variability(i,:,:,:),[],18));        
end

for c = 1 : 2
    figure(100); subplot(1,2,c);hold on
    for i= 1 : 18
        plot(1:3,[mvars(i,2+3*(c-1)) mvars(i,1+3*(c-1)) mvars(i,3+3*(c-1))],':k','LineWidth',0.5)
    end
    errorbar(1:3,[mvar(1,2+3*(c-1)) mvar(1,1+3*(c-1)) mvar(1,3+3*(c-1))],...
        [ mvar(2,2+3*(c-1)) mvar(2,1+3*(c-1)) mvar(2,3+3*(c-1))],'-k','LineWidth',2.5)
    xlabel('Head Rotations (deg)'); ylabel('Movement variability')
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'HR = CCW' , 'HR = 0' , 'HR = CW'})
    xlim([0.5 3.5])
end 

% stat analysis
% simple ttest, just to check
[h1,p1,~,stats1] = ttest(mvars(:,4),mvars(:,5));
[h2,p2,~,stats2] = ttest(mvars(:,4),mvars(:,6));
[h3,p3,~,stats3] = ttest(mvars(:,4),mean(mvars(:,5:6),2));

% repeated measure anova
varNames = cell(6,1);
for i = 1 : 6
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
HRs = {'HR0';'HRL';'HRR';'HR0';'HRL';'HRR'}; % head roll conditions
VFs = {'Y';'Y';'Y';'N';'N';'N'}; % Visual feedback conditions

factorNames = {'HRs','VisualFeedback'};
within = table(HRs, VFs, 'VariableNames', factorNames);
t = array2table(mvars, 'VariableNames',varNames);
rmvar = fitrm(t,'V1-V6~1','WithinDesign',within); % overal fit
[ranovatbvar] = ranova(rmvar,  'WithinModel','HRs*VisualFeedback');
mrm2 = multcompare(rmvar,'HRs','By','VisualFeedback','ComparisonType','bonferroni');

figname = 'Variability2-averaged and individuals'; 
savefig(gcf,figname) 
saveas(gcf,figname,'tiff')
close all 

%% Percentage of rightward MVs
pr = zeros(6,5,18);
for i = 1:6
    for j = 1 : 5            
        r = reshape(RD(i+0*3,j,1,:),1,18);r = flipud(r);
        l = reshape(RD(i+0*3,j,2,:),1,18);l = flipud(l);
        pr(i,j,:) = r.*100./(r+l);         
    end           
end 

figure(100); subplot(2,1,1);hold on   
xa = 5:-1:1;
for i = [2 1 3]
    r = reshape(RD(i,:,1,:),5,18);r = flipud(r);
    l = reshape(RD(i,:,2,:),5,18);l = flipud(l);
    pr = r.*100./(r+l); 
    for j = 1 : 5
        scatter(repmat(j,18,1),pr(j,:),'MarkerEdgeColor',colors{i})
    end
    plot(1:5,nanmean(pr,2),'Color',colors{i},'LineStyle',':','LineWidth',2.5); 
end  

ylabel('Percentage of Rightward MVs') 
set(gca,'XTick',fliplr(xa))
xlim([0,6])

figure(100); subplot(2,1,2);hold on   
for i = [2 1 3]
    r = reshape(RD(i+3,:,1,:),5,18);r = flipud(r);
    l = reshape(RD(i+3,:,2,:),5,18);l = flipud(l);
    pr = r.*100./(r+l); 
    for j = 1 : 5
        scatter(repmat(j,18,1),pr(j,:),'MarkerEdgeColor',colors{i})
    end
    plot(1:5,nanmean(pr,2),'Color',colors{i},'LineStyle','-','LineWidth',2.5); 
end  

ylabel('Percentage of Rightward MVs') 
set(gca,'XTick',fliplr(xa))
set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
xlim([0,6])

% repeated measure anova
varNames = cell(6*5,1);
for i = 1 : 6*5
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end

 % Create a table reflecting the within subject factors
HRs = cell(6*5,1); % head roll conditions
VFs = cell(6*5,1); % Visual feedback conditions
OPs = cell(6*5,1); % Obstacle Positions

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'HR0'; c1 = repmat(c1,10,1); HRs(1:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRL'; c1 = repmat(c1,10,1); HRs(2:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRR'; c1 = repmat(c1,10,1); HRs(3:3:end,1) = c1; 
ind = [ (1:6:6*5) (2:6:6*5) (3:6:6*5)];
c1 = cell(1,1); c1{1} = 'Y'; c1 = repmat(c1,15,1); VFs(ind,1) = c1;
c1 = cell(1,1); c1{1} = 'N'; c1 = repmat(c1,15,1); VFs(ind+3,1) = c1;
for i = 1 : 5
    o = strcat('O',num2str(i));
    c1 = cell(1,1); c1{1} = o; c1 = repmat(c1,6,1); OPs((i-1)*6+1:i*6,1) = c1;
end
pr2 = reshape(pr, 6*5,18);

factorNames = {'HRs','VisualFeedback','ObstaclePos'};
within = table(HRs, VFs,OPs, 'VariableNames', factorNames);
t = array2table(pr2', 'VariableNames',varNames);
rmpr = fitrm(t,'V1-V30~1','WithinDesign',within); % overal fit
[ranovatbpr] = ranova(rmpr,  'WithinModel','HRs*VisualFeedback*ObstaclePos');
mrm3 = multcompare(rmpr,'HRs','By','VisualFeedback','ComparisonType','bonferroni');

%% save results

save('stats_finale.mat','mrm1','ranovatbcol','rmcol','mrm2','ranovatbvar','rmvar','mrm3','ranovatbpr','rmpr')