%% Analysis of rotation and expansion biases
% Written by PA, Nov. 2018
% Updated by PA, Jan. 2019

%% Match the relevant obstacles together
% mirrored obstacles are matched to have both rightward and leftward MVs (3L&3R, 1&5, 2&4)

% Initialization
alfabeta = zeros(3,2,2,18); %  3 pairs of obstacles (left and right), 2 parameters, #subjects
y = 19; % point of interest along y-axis
c = 1; % with/without visual feedback (0=> with, 1=> without)
valid = 1; % check if there is enough data on both right and left direction
obs = [3 3 ; 5 1 ; 4 2];
% Load the data
load('WholeComp3.mat')
j = 1;
for i = 1 : 18
    for j = 1:3 
        if j == 1
            ec = zeros(3,2);
            for t = 1 : 3
                rrd = 100*(reshape(RD(t+c*3,3,:,i),[],1)./sum(reshape(RD(t+c*3,3,:,i),[],1)));
                ec(t,:) = rrd > 10;            
            end
            if sum(ec(:)==0) > 0
                valid = 0; 
            else
                valid = 1;
            end
        else
            valid = 1;
        end
        if valid == 1
            % find pointers 
            p1 = find(trajc(:,2,1+c*3,1,obs(j,1),1,i) == y | trajc(:,2,1+c*3,1,obs(j,1),1,i) > y,1); 
            p1p = find(trajc(:,2,1+c*3,1,obs(j,2),2,i) == y | trajc(:,2,1+c*3,1,obs(j,2),2,i) > y,1); 
            p3 = find(trajc(:,2,3+c*3,1,obs(j,1),1,i) == y | trajc(:,2,3+c*3,1,obs(j,1),1,i) > y,1);
            p3p = find(trajc(:,2,3+c*3,1,obs(j,2),2,i) == y | trajc(:,2,3+c*3,1,obs(j,2),2,i) > y,1);
            p33 = find(trajc(:,2,2+c*3,1,obs(j,1),1,i) == y | trajc(:,2,2+c*3,1,obs(j,1),1,i) > y,1);
            p33p = find(trajc(:,2,2+c*3,1,obs(j,2),2,i) == y | trajc(:,2,2+c*3,1,obs(j,2),2,i) > y,1);

            % calculate deltaX for different directions and conditions
            dxrhr = trajc(p3,1,3+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
            dxlhr = trajc(p3p,1,3+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);
            dxrhl = trajc(p3,1,2+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
            dxlhl = trajc(p3p,1,2+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);

            % calculate alfa and beta
            alfabeta(j,1,1,i) = (+(dxlhr-dxlhl))/2; % alfa (Left)
            alfabeta(j,1,2,i) = (-(dxlhr+dxlhl))/2; % beta (Left)
            alfabeta(j,2,1,i) = (+(dxrhr-dxrhl))/2; % alfa (Right)
            alfabeta(j,2,2,i) = (+(dxrhr+dxrhl))/2; % beta (Right)
        end
    end    
end


%% Plot the result
% sor the results first
alfabetap = zeros(6,2,18);
s = [ 2 3 1 1 3 2;1 1 1 2 2 2]';
for j = 1 : 6
    alfabetap(j,1,:) = alfabeta(s(j,1),s(j,2),1,:);
    alfabetap(j,2,:) = alfabeta(s(j,1),s(j,2),2,:);
end
ms = zeros(6,2);
xs = 0:3:17;
for j = 1 : 6
    ddata = reshape(alfabetap(j,2,:),1,18);
    figure(100); hold on; scatter(repmat(xs(j),1,18),ddata)
    ms(j,1) = nanmean(ddata); ms(j,2) = nanstd(ddata);
end
figure(100); hold on; plot(xs,ms(:,1))

ylabel('Expansion biases') 
ylim([-2.5 2.5])
set(gca,'XTick',xs)
set(gca, 'XTickLabels',{'Most Rightward','Right','Center-L','Center-R','Left','Most Leftward'})
xlim([0-3,18])

%% save the results
save('RotExp.mat','alfabeta','alfabetap')

%% stat analysis
% ttest
[h1,p1,~,stats1] = ttest(reshape(alfabetap(1,2,:),1,18));




% Anova (it doesn't have any meaning though!!!!)
% check if the expansion is consistently greater than zero for all
% participants and osbtacle positions
% repeated measure anova

varNames = cell(6*1,1);
for i = 1 : 6*1
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end

 % Create a table reflecting the within subject factors
% HRs = cell(6*5,1); % head roll conditions
% VFs = cell(6*5,1); % Visual feedback conditions
%OPs = cell(6*1,1); % Obstacle Positions
OPs={'ML'; 'L';'C-L'; 'C-R'; 'R'; 'MR'};
% Assiging the values to the parameters based on the data sorting
% c1 = cell(1,1); c1{1} = 'HR0'; c1 = repmat(c1,10,1); HRs(1:3:end,1) = c1;
% c1 = cell(1,1); c1{1} = 'HRL'; c1 = repmat(c1,10,1); HRs(2:3:end,1) = c1;
% c1 = cell(1,1); c1{1} = 'HRR'; c1 = repmat(c1,10,1); HRs(3:3:end,1) = c1; 
% ind = [ (1:6:6*5) (2:6:6*5) (3:6:6*5)];
% c1 = cell(1,1); c1{1} = 'Y'; c1 = repmat(c1,15,1); VFs(ind,1) = c1;
% c1 = cell(1,1); c1{1} = 'N'; c1 = repmat(c1,15,1); VFs(ind+3,1) = c1;
% for i = 1 : 5
%     o = strcat('O',num2str(i));
%     c1 = cell(1,1); c1{1} = o; c1 = repmat(c1,6,1); OPs((i-1)*6+1:i*6,1) = c1;
% end
betas = reshape(alfabetap(:,2,:), 6 , 18);

factorNames = {'ObstaclePos'};
within = table(OPs, 'VariableNames', factorNames);
t = array2table(betas', 'VariableNames',varNames);
rmpr = fitrm(t,'V1-V6~1','WithinDesign',within); % overal fit
[ranovatbpr] = ranova(rmpr,  'WithinModel','ObstaclePos');
mrm3 = multcompare(rmpr,'ObstaclePos','ComparisonType','bonferroni');


