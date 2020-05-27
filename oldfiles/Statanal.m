%% Obstacle Avoidacne Data analysis
% Ploting setup 3, version2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'};
obstaclepos = [22.4750   11.0000
               20.6750   11.0000
               18.8750   11.0000
               17.0750   11.0000
               15.2750   11.0000
               ];
colors = { 'k' ; 'g' ; 'r'};
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
init = 0;
% Choose which section to run!
createdata = 0;
statanal = 1;


%% Pooling
if createdata == 1
    load('WholeComp.mat','trajc')
    % 1080 = 6*5*2*18 %conditions(3HR*2VF)*obstaclepos*right/left*subjects
    trajxstats = reshape(trajc(:,1,:,1,:,:,:),200,1080); 
    trajxstats = trajxstats';
    trajystats = reshape(trajc(:,2,:,1,:,:,:),200,1080); 
    trajystats = trajystats';
    group = cell(1,4); % condition , obstacle#, right/left, subjectID
    group{1,1} = repmat((1:6)',1080/6,1);
    group{1,2} = repmat(repmat(reshape((repmat((1:5),6,1)),[],1),2,1),18,1);
    group{1,3} = repmat(reshape((repmat((0:1),6*5,1)),[],1),18,1);
    group{1,4} = reshape((repmat((1:18),6*5*2,1)),[],1);
    
    %% test the data arrangement
    for i = [1 2 3]
        for j = 5
            figure(j); hold on
            for d = 1 : 2 
                indx = find(group{1}==i & group{2}==j & group{3}==d-1 &group{4}==14);
                tempx = trajxstats(indx,:) - 1*18.875;
                ys = trajystats(indx,:); 
                if ~isnan(ys)
                    plot(tempx,ys,'Color',colors{i-0*3},'Linewidth',2)
                    scatter([obstaclepos(j,1)-18.8750 0 0],[30-obstaclepos(j,2) 9 29],'ro')  
                end
            end
        end       
    end
    %% Save the data
    txt = ' group = condition, obstacles, direction, subjects';
    save('StatAnal.mat','trajxstats','trajystats','group','txt');    

end

%% Functional statistical analysis


if statanal == 1
    load('StatAnal.mat')
    % fanovan to detect the general head roll effect!
    % indicate the conditions you want to compare
    indx1 = find(group{3} == 1 &...
                group{2} == 3 & ...
                group{1} == 4& group{4} ~= 14);
    indx2 = find(group{3} == 1 &...
                group{2} == 3 & ...
                group{1} == 6 & group{4} ~= 14);
    indx3 = find(group{3} == 1 &...
                group{2} == 3 & ...
                group{1} == 6 & group{4} ~= 14);
    indx = [indx1; indx2];
    
    newGroupx = {group{1}(indx) group{4}(indx)};
    xRs = trajxstats(indx,:);
    %continuous = [0 0];
    %[y,allgrps,nanrow,varinfo] = removenans(xRs,newGroupx,continuous);
    [px,corrPx,tx,statsx] = fanovan(xRs, newGroupx,'model','full',...
                                    'random' , length(newGroupx),...
                                    'varnames',{'conditions' 'subjects'});
    figure; hold on
    plot((1:200)/2,corrPx(1,:),'linewidth',2);
    line([0 100],[0.05 0.05],'color','r','linewidth',2);
    ylim([0 1]);
    xlabel('Percent Y Distance');
    ylabel('p-value of X deviation');
                                
    % evaluate the interaction effect
    % basically we are subtracting the baseline (HR=0)from each HR in both
    % side and compare them together!
    indx1 = find(group{3} == 0 &...
                group{2} == 5 & ...
                group{1} == 4 & group{4}~=14);
    indx2 = find(group{3} == 0 &...
                group{2} == 5 & ...
                group{1} == 6 & group{4}~=14 );
    indx11 = find(group{3} == 1 &...
                group{2} == 5 & ...
                group{1} == 4 & group{4}~=14);
    indx22 = find(group{3} == 1 &...
                  group{2} == 3 & ...
                  group{1} == 6 & group{4}~=14);
    indx = [indx2; indx22 ];

    newGroupx = {group{3}(indx) group{4}(indx)};
%     xRs = trajxstats(indx,:);
    xRs = [trajxstats(indx2,:)  - 1*trajxstats(indx1,:);...
           trajxstats(indx22,:) - 1*trajxstats(indx11,:)];
    
    [px,corrPx,tx,statsx] = fanovan(xRs, newGroupx,'model','full',...
                                    'random' , length(newGroupx),...
                                    'varnames',{'conditions' 'subjects'});
    
    figure; hold on
    plot((1:200)/2,corrPx(1,:),'linewidth',2);
    line([0 100],[0.05 0.05],'color','r','linewidth',2);
    ylim([0 1]);
    xlabel('Percent Y Distance');
    ylabel('p-value of X deviation');
    close all
    
    % check the factorial anova idea! Not working really!
%     indx1 = find(group{3} == 0 &...
%                 group{2} == 3 & ...
%                 group{1} == 4 & group{4}~=14);
%     indx2 = find(group{3} == 0 &...
%                 group{2} == 3 & ...
%                 group{1} == 5 & group{4}~=14 );
%     indx3 = find(group{3} == 0 &...
%                 group{2} == 3 & ...
%                 group{1} == 6 & group{4}~=14);
%     indx = [indx1; indx2 ;indx3];
%     yRs = nanmean(trajystats(indx,:));
%     yRss = repmat(yRs,[length(indx),1]);yRss = yRss'; yRss = yRss(:);
%     cond = group{1}(indx); conds = repmat(cond,[200,1]); conds = conds';conds = conds(:);
%     subs = group{4}(indx); subss = repmat(subs,[200,1]); subss = subss';subss = subss(:);
%     newGroupx = {yRss, conds, subss};
%     xRs = trajxstats(indx,:);
%     xRss = xRs'; xRss = xRss(:);
%     [px,tablex,statsx] = anovan(xRss, newGroupx,'model','interaction',...
%                                     'random' , length(newGroupx),...
%                                     'varnames',{'Ys' 'conditions' 'subjects'});
end
%%
close all

% -----------------------------------
function [y,allgrps,nanrow,varinfo] = removenans(y,group,continuous)

% Find NaNs among response and group arrays

%FANOVA edit
%change so n and nanrow reflects the number of rows (use size)
n = size(y,1);
nanrow = isnan(mean(y,2));
ng = length(group);
for j=1:ng
   gj = group{j};
   if (size(gj,1) == 1), gj = gj(:); end
   if (size(gj,1) ~= n)
      error('stats:anovan:InputSizeMismatch',...
            'Group variable %d must have %d elements.',j,n);
   end
   if (ischar(gj)), gj = cellstr(gj); end
   if ~isvector(gj)
       error('stats:anovan:BadGroup',...
             'Grouping variable must be a vector or a character array.');
   end

   group{j} = gj;
   if (isnumeric(gj))
      nanrow = (nanrow | isnan(gj));
   elseif isa(gj,'categorical')
      nanrow = (nanrow | isundefined(gj));
   else
      nanrow = (nanrow | strcmp(gj,''));
   end
end

% Remove rows with NaN anywhere
%FANOVA edit
%change to remove the entire function (row) for trials with nans, also have
%n reflect functional dims
y(nanrow,:) = [];
n = size(y,1);
ng = length(group);
dfvar = zeros(ng,1);
allgrps = zeros(n, ng);
allnames = cell(ng,1);

% Get arrays describing the groups
for j=1:ng
   gj = group{j};
   gj(nanrow,:) = [];
   group{j} = gj;
   if continuous(j)
      dfvar(j) = 1;
      allgrps(:,j) = gj;
   else
      [gij,gnj] = grp2idx(gj);
      nlevels = size(gnj,1);
      dfvar(j) = nlevels - 1;
      allgrps(:,j) = gij;
      allnames{j} = gnj;
   end
end

% The df and allnames information does not yet reflect nesting.  We will
% fix up allnames for nested variables later, and not use df for them.

varinfo.df = dfvar;
varinfo.allnames = allnames;
end
