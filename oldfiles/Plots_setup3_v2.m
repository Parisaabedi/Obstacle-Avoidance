%% Obstacle Avoidacne Data analysis
% Ploting setup 3, version2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'};
colors = { 'k' ; 'g' ; 'r'};
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
obstaclepos = [22.4750   11.0000
               20.6750   11.0000
               18.8750   11.0000
               17.0750   11.0000
               15.2750   11.0000
               ];
%obstaclepos = flipud(obstaclepos);
init = 0;
pooled = 0;
poolploting = 1;
plotting = 0;
trajectories = 0;

%% Pooling
if pooled == 1
    % individual participants
    % Initialization (condition, mean/std, target, right(0)/left(1),subjects)
    ima = zeros(6,2,5,2,18); % initilal movement angle    
    mcurv = zeros(6,2,5,2,18); % movement curvature
    DpassO = zeros(6,2,5,2,18); % Distance passing the Obstacle
    DfromO = zeros(6,2,5,2,18); % Distance from Obstacle
    MaxVel = zeros(6,2,5,2,18); % Max Velocity
    TMaxVel = zeros(6,2,5,2,18); % time of max verlocity wrt mvonset
    MaxAcc = zeros(6,2,5,2,18); % Max Acceleration
    TMaxAcc = zeros(6,2,5,2,18); %time of max acceleration wrt mvonset
    trajc = zeros(200,2,6,2,5,2,18); % Normalized trajectories
    RT = zeros(6,2,5,2,18); % Reaction time
    MD = zeros(6,2,5,2,18); % Movement Duration
    RD = zeros(6,2,5,18); % Right MV. Direction (Cond, right(0)/left(1), targets, subs)
    NCd = zeros(6,5,2,18);
    NC = zeros(6,18); % number of collisions
    NR = zeros(6,18); % number of repeats
    NMT = zeros(6,18); % number of missed target
    
    subs = 1:18;     
    for s = 1:18        
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp','_',sID,'.mat');
        load(filename)       
        for i = 1 : 6
            temp = params(1:Ns(i),:,i);
            ttraj = ntraj(:,1:Ns(i),:,:);
            rep = (temp(:,13)-temp(:,12)>1);
            for j = 1 : 2
                for t = 1 : 5
                    % no collision, no repeat, and either right to left -ward
                    % directions, good trials
                    qnc = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t; 
                    
                    ima(i,1,t,j,s) = mean(temp(qnc,18));
                    ima(i,2,t,j,s) = std(temp(qnc,18));
                    mcurv(i,1,t,j,s) = mean(temp(qnc,19));
                    mcurv(i,2,t,j,s) = std(temp(qnc,19));
                    DpassO(i,1,t,j,s) = mean(temp(qnc,21));
                    DpassO(i,2,t,j,s) = std(temp(qnc,21));
                    DfromO(i,1,t,j,s) = mean(temp(qnc,21));
                    DfromO(i,2,t,j,s) = std(temp(qnc,21));
                    MaxVel(i,1,t,j,s) = mean(temp(qnc,14));
                    MaxVel(i,2,t,j,s) = std(temp(qnc,14));
                    TMaxVel(i,1,t,j,s) = mean(temp(qnc,15)-temp(qnc,12));
                    TMaxVel(i,2,t,j,s) = std(temp(qnc,15)-temp(qnc,12));
                    MaxAcc(i,1,t,j,s) = mean(temp(qnc,16));
                    MaxAcc(i,2,t,j,s) = std(temp(qnc,16));
                    TMaxAcc(i,1,t,j,s) = mean(temp(qnc,17)-temp(qnc,12));
                    TMaxAcc(i,2,t,j,s) = std(temp(qnc,17)-temp(qnc,12));
                    if sum(qnc) ~= 0
                        trajc(:,1,i,1,t,j,s) = nanmean(ttraj(:,qnc,1,i),2); % mean xtraj
                        trajc(:,1,i,2,t,j,s) = nanstd(ttraj(:,qnc,1,i),0,2); % std xtraj
                        trajc(:,2,i,1,t,j,s) = nanmean(ttraj(:,qnc,2,i),2); % mean ytraj
                        trajc(:,2,i,2,t,j,s) = nanstd(ttraj(:,qnc,2,i),0,2); % std ytraj                        
                    else
                        trajc(:,1,i,1,t,j,s) = zeros(200,1)*nan; % mean xtraj
                        trajc(:,1,i,2,t,j,s) = zeros(200,1)*nan; % std xtraj
                        trajc(:,2,i,1,t,j,s) = zeros(200,1)*nan; % mean ytraj
                        trajc(:,2,i,2,t,j,s) = zeros(200,1)*nan; % std ytraj                        
                    end
                    RT(i,1,t,j,s) = mean(temp(qnc,12)-temp(qnc,11));
                    RT(i,2,t,j,s) = std(temp(qnc,12)-temp(qnc,11));
                    MD(i,1,t,j,s) = mean(temp(qnc,13)-temp(qnc,12));
                    MD(i,2,t,j,s) = std(temp(qnc,13)-temp(qnc,12));                      
                    qrc = temp(:,6) == 0 &temp(:,7) == (j-1) & temp(:,8) == t;
                    RD(i,j,t,s) = sum(qrc);
                    qncc = temp(:,5) == 1 & temp(:,7) == (j-1) & temp(:,8) == t;
                    NCd(i,t,j,s) = sum(qncc);
                end               
            end             
            NC(i,s) = sum(temp(:,5));
            NR(i,s) = sum(temp(:,4));             
            NMT(i,s) = sum(temp(:,10));
        end
    end
    save('WholeComp.mat','ima','mcurv','DpassO','DfromO','MaxVel','TMaxVel',...
         'MaxAcc','TMaxAcc','trajc','RT','MD','NC','NR','NMT','RD', 'NCd')
end
%%
if poolploting == 1
    %%
    load('WholeComp.mat')
    %%
    imac = zeros(6,4,5,2); % initilal movement angle    
    mcurvc = zeros(6,4,5,2); % movement curvature
    DpassOc = zeros(6,4,5,2); % Distance passing the Obstacle
    DfromOc = zeros(6,4,5,2); % Distance from Obstacle
    MaxVelc = zeros(6,4,5,2); % Max Velocity
    TMaxVelc = zeros(6,4,5,2); % time of max verlocity wrt mvonset
    MaxAccc = zeros(6,4,5,2); % Max Acceleration
    TMaxAccc = zeros(6,4,5,2); %time of max acceleration wrt mvonset
    trajcc = zeros(200,2,6,4,5,2); % Normalized trajectories(samples,x/y, cond,mean/std,obst, right/left)
    RTc = zeros(6,4,5,2); % Reaction time
    MDc = zeros(6,4,5,2); % Movement Duration
    NCdc = zeros(6,2,5,2); % cond, mean/std, obst, right/left
    %NCc = zeros(6,1,18); % number of collisions(condition,mean/std, subject)
    %NRc = zeros(6,1,18); % number of repeats(condition,mean/std, subject)
    %NMTc = zeros(6,1,18); % number of missed target(condition,mean/std, subject)
    
    % Standardize trajectories
    variability = zeros(6,5,2,18);% cond,tar,dir,sub
    biases = zeros(6,5,2,18); % cond, tar, dir, sub
    signs = zeros(6,5,2,18); % cond, tar, dir, sub
    % exclusion criteria (based on the number of MVs for each direction <10%)
    EC = zeros(6,5,2,18);
    PRD = EC;  % percentage of rightward MVs
    for s = 1:18
        for i = 1 : 6
            for j = 1:5
                for d = 1 : 2
                    tempxb = trajc(:,1,i,1,j,d,s);
                    if i <5
                        uEb=trajc(:,1,1,1,j,d,s);                        
                    else
                        uEb= trajc(:,1,4,1,j,d,s);                        
                    end
                    tempx = trajc(:,1,i,1,j,d,s);
                    ys = trajc(:,2,i,1,j,d,s);
                    uE=tempx + trajc(:,1,i,2,j,d,s);
                    lE=tempx - trajc(:,1,i,2,j,d,s);
                    xP=[lE;flipud(uE)];
                    yP=[ys;flipud(ys)];
                    variability(i,j,d,s) = polyarea(xP,yP);
                    xPb=[tempxb;flipud(uEb)];
                    yPb=[ys;flipud(ys)];
                    si = tempxb - uEb; % detect the sign
                    biases(i,j,d,s) = polyarea(xPb,yPb);
                    if sum(si) ~= 0
                        signs(i,j,d,s)  = (sum(si)/abs(sum(si)));
                    elseif isnan(sum(si))
                        signs(i,j,d,s) = nan;
                        
                    end
                end
                rrd = 100*reshape(RD(i,:,j,s),[],1)./sum(reshape(RD(i,:,j,s),[],1));
                t = rrd > 10;
                EC(i,j,:,s) = t;
                PRD(i,j,:,s) = rrd;
            end
        end
    end
    variabilityc = zeros(6,2,5,2); % cond, mean/std, obst, right/left
    for i = 1 : 6
        for d = 1:2
            for t = 1 : 5
                imac(i,1,t,d) = nanmean(reshape(ima(i,1,t,d,:),18,1));
                imac(i,2,t,d) = nanstd(reshape(ima(i,1,t,d,:),18,1));
                imac(i,3,t,d) = nanmean(reshape(ima(i,2,t,d,:),18,1));
                imac(i,4,t,d) = nanstd(reshape(ima(i,2,t,d,:),18,1));
                
                mcurvc(i,1,t,d) = nanmean(reshape(mcurv(i,1,t,d,:),18,1));
                mcurvc(i,2,t,d) = nanstd(reshape(mcurv(i,1,t,d,:),18,1));
                mcurvc(i,3,t,d) = nanmean(reshape(mcurv(i,2,t,d,:),18,1));
                mcurvc(i,4,t,d) = nanstd(reshape(mcurv(i,2,t,d,:),18,1));
                
                DpassOc(i,1,t,d) = nanmean(reshape(DpassO(i,1,t,d,:),18,1));
                DpassOc(i,2,t,d) = nanstd(reshape(DpassO(i,1,t,d,:),18,1));
                DpassOc(i,3,t,d) = nanmean(reshape(DpassO(i,2,t,d,:),18,1));
                DpassOc(i,4,t,d) = nanstd(reshape(DpassO(i,2,t,d,:),18,1));
                
                DfromOc(i,1,t,d) = nanmean(reshape(DfromO(i,1,t,d,:),18,1));
                DfromOc(i,2,t,d) = nanstd(reshape(DfromO(i,1,t,d,:),18,1));
                DfromOc(i,3,t,d) = nanmean(reshape(DfromO(i,2,t,d,:),18,1));
                DfromOc(i,4,t,d) = nanstd(reshape(DfromO(i,2,t,d,:),18,1));
                
                MaxVelc(i,1,t,d) = nanmean(reshape(MaxVel(i,1,t,d,:),18,1));
                MaxVelc(i,2,t,d) = nanstd(reshape(MaxVel(i,1,t,d,:),18,1));
                MaxVelc(i,3,t,d) = nanmean(reshape(MaxVel(i,2,t,d,:),18,1));
                MaxVelc(i,4,t,d) = nanstd(reshape(MaxVel(i,2,t,d,:),18,1));
                
                TMaxVelc(i,1,t,d) = nanmean(reshape(TMaxVel(i,1,t,d,:),18,1));
                TMaxVelc(i,2,t,d) = nanstd(reshape(TMaxVel(i,1,t,d,:),18,1));
                TMaxVelc(i,3,t,d) = nanmean(reshape(TMaxVel(i,2,t,d,:),18,1));
                TMaxVelc(i,4,t,d) = nanstd(reshape(TMaxVel(i,2,t,d,:),18,1));
                
                MaxAccc(i,1,t,d) = nanmean(reshape(MaxAcc(i,1,t,d,:),18,1));
                MaxAccc(i,2,t,d) = nanstd(reshape(MaxAcc(i,1,t,d,:),18,1));
                MaxAccc(i,3,t,d) = nanmean(reshape(MaxAcc(i,2,t,d,:),18,1));
                MaxAccc(i,4,t,d) = nanstd(reshape(MaxAcc(i,2,t,d,:),18,1));
                
                TMaxAccc(i,1,t,d) = nanmean(reshape(TMaxAcc(i,1,t,d,:),18,1));
                TMaxAccc(i,2,t,d) = nanstd(reshape(TMaxAcc(i,1,t,d,:),18,1));
                TMaxAccc(i,3,t,d) = nanmean(reshape(TMaxAcc(i,2,t,d,:),18,1));
                TMaxAccc(i,4,t,d) = nanstd(reshape(TMaxAcc(i,2,t,d,:),18,1));
                
                trajcc(:,1,i,1,t,d) = nanmean(reshape(trajc(:,1,i,1,t,d,:),[200,18]),2);
                trajcc(:,1,i,2,t,d) = nanstd(reshape(trajc(:,1,i,1,t,d,:),[200,18]),0,2);
                trajcc(:,1,i,3,t,d) = nanmean(reshape(trajc(:,1,i,2,t,d,:),[200,18]),2);
                trajcc(:,1,i,4,t,d)= nanstd(reshape(trajc(:,1,i,2,t,d,:),[200,18]),0,2);
                
                trajcc(:,2,i,1,t,d) = nanmean(reshape(trajc(:,2,i,1,t,d,:),[200,18]),2);
                trajcc(:,2,i,2,t,d) = nanstd(reshape(trajc(:,2,i,1,t,d,:),[200,18]),0,2);
                trajcc(:,2,i,3,t,d) = nanmean(reshape(trajc(:,2,i,2,t,d,:),[200,18]),2);
                trajcc(:,2,i,4,t,d)= nanstd(reshape(trajc(:,2,i,2,t,d,:),[200,18]),0,2);
                
                RTc(i,1,t,d) = nanmean(reshape(RT(i,1,t,d,:),18,1));
                RTc(i,2,t,d) = nanstd(reshape(RT(i,1,t,d,:),18,1));
                RTc(i,3,t,d) = nanmean(reshape(RT(i,2,t,d,:),18,1));
                RTc(i,4,t,d) = nanstd(reshape(RT(i,2,t,d,:),18,1));
                
                MDc(i,1,t,d) = nanmean(reshape(MD(i,1,t,d,:),18,1));
                MDc(i,2,t,d) = nanstd(reshape(MD(i,1,t,d,:),18,1));
                MDc(i,3,t,d) = nanmean(reshape(MD(i,2,t,d,:),18,1));
                MDc(i,4,t,d) = nanstd(reshape(MD(i,2,t,d,:),18,1));  
                
                NCdc(i,1,t,d) = nanmean(reshape(NCd(i,t,d,:),18,1));
                NCdc(i,2,t,d) = nanstd(reshape(NCd(i,t,d,:),18,1));
                
                variabilityc(i,1,t,d) = nanmean(reshape(variability(i,t,d,:),18,1));
                variabilityc(i,2,t,d) = nanstd(reshape(variability(i,t,d,:),18,1));
            end
        end
        %NCc(i,1,:) = sum(reshape(NC(i,:,:,:),[10,18]));
        %NCc(i,2,:) = std(reshape(NC(i,:,:,:),[10,18]));
        
        %NRc(i,1,:) = sum(reshape(NR(i,:,:,:),[10,18]));
        %NRc(i,2,:) = std(reshape(NR(i,:,:,:),[10,18]));
        
        %NMTc(i,1,:) = sum(reshape(NMT(i,:,:,:),[10,18]));
        %NMTc(i,2,:) = std(reshape(NMT(i,:,:,:),[10,18]));
    end
    
    %% Plot Initial movement angles
    % for each obstacle
    figure(12); subplot(2,1,1);hold on  
    %xa = 5:-1:1;            
    h1 = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(imac(i,1,:,1),5,1));
        t2 = flipud(reshape(imac(i,2,:,1),5,1));
        h1(1,i)=errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i});
    end  

    figure(12); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(imac(i,1,:,2),5,1));
        t2 = flipud(reshape(imac(i,2,:,2),5,1));
        errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i});
    end
    
    figure(12); subplot(2,1,1);hold on
    xa = 5:-1:1;            
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(imac(i+3,1,:,1),5,1));
        t2 = flipud(reshape(imac(i+3,2,:,1),5,1));
        h(1,i)=errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end 
    
    
    legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
        {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
         'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
    legend('Location','southeast')    
    ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    title('Initial movement angle modulation by varying head roll - Rightward MVs')
    
    figure(12); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(imac(i+3,1,:,2),5,1));
        t2 = flipud(reshape(imac(i+3,2,:,2),5,1));
        errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end
    xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Initial movement angle modulation by varying head roll - Leftward MVs')

    figname = strcat('IMA-OBs-Whole'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

    %% Plot curvature
    % for each obstacle
    figure(22); subplot(2,1,1);hold on  
    %xa = 5:-1:1;            
    h1 = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(mcurvc(i,1,:,1),5,1));
        t2 = flipud(reshape(mcurvc(i,2,:,1),5,1));
        h1(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end  

    figure(22); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(mcurvc(i,1,:,2),5,1));
        t2 = flipud(reshape(mcurvc(i,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end
    
    figure(22); subplot(2,1,1);hold on
    xa = 5:-1:1;            
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(mcurvc(i+3,1,:,1),5,1));
        t2 = flipud(reshape(mcurvc(i+3,2,:,1),5,1));
        h(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end     
    legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
        {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
         'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
    legend('Location','southeast') ; 
    ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    title('Maximum curvature modulation by varying head roll - Rightward MVs')
    
    figure(22); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(mcurvc(i+3,1,:,2),5,1));
        t2 = flipud(reshape(mcurvc(i+3,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end
    xlabel('Obstacle position'); ylabel('Maximum cruvature (cm)'); 
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Maximum curvature modulation by varying head roll - Leftward MVs')

    figname = strcat('MC-OBs-Whole'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

    %% Plot distance passing the obstacle
    % for each obstacle
    figure(32); subplot(2,1,1);hold on  
    %xa = 5:-1:1;            
    h1 = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(DpassOc(i,1,:,1),5,1));
        t2 = flipud(reshape(DpassOc(i,2,:,1),5,1));
        h1(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end  

    figure(32); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(DpassOc(i,1,:,2),5,1));
        t2 = flipud(reshape(DpassOc(i,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end
    
    figure(32); subplot(2,1,1);hold on
    xa = 5:-1:1;            
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(DpassOc(i+3,1,:,1),5,1));
        t2 = flipud(reshape(DpassOc(i+3,2,:,1),5,1));
        h(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end 
    
    
    legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
        {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
         'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
    legend('Location','southeast') ; 
    ylabel('Distance passing Obstacle (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    title('Distance passing Obstacle modulation by varying head roll - Rightward MVs')
    
    figure(32); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(DpassOc(i+3,1,:,2),5,1));
        t2 = flipud(reshape(DpassOc(i+3,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end
    xlabel('Obstacle position'); ylabel('Distance passing Obstacle (cm)'); 
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Distance passing Obstacle modulation by varying head roll - Leftward MVs')

    figname = 'DPO-OBs-Whole'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

%% Plot distance from the obstacle
    % for each obstacle
    figure(42); subplot(2,1,1);hold on  
    %xa = 5:-1:1;            
    h1 = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(DfromOc(i,1,:,1),5,1));
        t2 = flipud(reshape(DfromOc(i,2,:,1),5,1));
        h1(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end  

    figure(42); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(DfromOc(i,1,:,2),5,1));
        t2 = flipud(reshape(DfromOc(i,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i});
    end
    
    figure(42); subplot(2,1,1);hold on
    xa = 5:-1:1;            
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = flipud(reshape(DfromOc(i+3,1,:,1),5,1));
        t2 = flipud(reshape(DfromOc(i+3,2,:,1),5,1));
        h(1,i)=errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end 
        
    legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
        {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
         'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
    legend('Location','southeast') ; 
    ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    title('Distance from Obstacle modulation by varying head roll - Rightward MVs')
    
    figure(42); subplot(2,1,2);hold on
    for i = [2 1 3]
        t1 = flipud(reshape(DfromOc(i+3,1,:,2),5,1));
        t2 = flipud(reshape(DfromOc(i+3,2,:,2),5,1));
        errorbar(1:5,abs(t1).*1,t2./sqrt(18),'Color',colors{i},'LineStyle',':');
    end
    xlabel('Obstacle position'); ylabel('Distance from Obstacle (cm)'); 
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Distance from Obstacle modulation by varying head roll - Leftward MVs')

    figname = 'DfO-OBs-Whole'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

    %% Plot Number of collisition
    figure(41); hold on%subplot(1,2,1); hold on
    HR = [0 -30 30];
    for i = [2 1 3]
        errorbar(HR(i), mean(NC(i,:)),std(NC(i,:))./sqrt(18),'Color',colors{i},'Marker','x')        
    end
    tt = mean(NC,2);
    plot(HR([2 1 3]),tt([2 1 3]),'-b')
    %plot(repmat([-30;0;30],1,18),NC([2 1 3],:),'--')
    xlabel('Head Rotations (deg)'); ylabel('Obstacle collision rate')
    title('With visual feedback')
    figure(42);hold on %subplot(1,2,2); hold on
    for i = [2 1 3]
        errorbar(HR(i), mean(NC(i+3,:)),std(NC(i+3,:))./sqrt(18),'Color',colors{i},'Marker','x')
    end
    tt = mean(NC,2);
    plot(HR([2 1 3]),tt([5 4 6]),':b')
    %plot(repmat([-30;0;30],1,18),NC([5 4 6],:),'--')
    xlabel('Head Rotations (deg)'); ylabel('Obstacle collision rate')
    title('No visual feedback')
    %supertitle('Obstacle collision modulation by varying head roll')

    figname = 'NCs-OBs-Whole'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %% Plot Reaction times
    figure(61); hold on%subplot(1,2,1); hold on
    HR = [0 -30 30];
    mrt = zeros(3,1); srt = zeros(3,1);
    for i = [2 1 3]
        mrt(i) = nanmean(1000*reshape([RTc(i+3*1,1,5:-1:3,1) RTc(i+3*1,1,3:-1:1,2)],[],1));
        srt(i) = nanmean(1000*reshape([RTc(i+3*1,2,5:-1:3,1) RTc(i+3*1,2,3:-1:1,2)],[],1));
        errorbar(HR(i), mrt(i),srt(i)/sqrt(18),'Color',colors{i},'Marker','x')        
    end
    
    plot(HR([2 1 3]),mrt([2 1 3]),':b')
    %plot(repmat([-30;0;30],1,18),NC([2 1 3],:),'--')
    xlabel('Head Rotations (deg)'); ylabel('Reaction time (ms)')
    title('Modulation of Reaction time with varying head roll')

    figname = 'RTs-OBs-Whole'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %% Plot Movement durations
    figure(71); hold on%subplot(1,2,1); hold on
    HR = [0 -30 30];
    mmd = zeros(3,1); smd = zeros(3,1);
    for i = [2 1 3]
        mmd(i) = nanmean(1000*reshape([MDc(i+3*0,1,5:-1:3,1) MDc(i+3*0,1,3:-1:1,2)],[],1));
        smd(i) = nanmean(1000*reshape([MDc(i+3*0,2,5:-1:3,1) MDc(i+3*0,2,3:-1:1,2)],[],1));
        errorbar(HR(i), mmd(i),smd(i),'Color',colors{i},'Marker','x')        
    end
    
    plot(HR([2 1 3]),mmd([2 1 3]),'-b')
    %plot(repmat([-30;0;30],1,18),NC([2 1 3],:),'--')
    xlabel('Head Rotations (deg)'); ylabel('Reaction time (ms)')
    title('Modulation of Movement Durations with varying head roll')

    figname = 'MDs-OBs-Whole'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %% Plot Maximum velocity
    figure(112); hold on%subplot(1,2,1); hold on
    HR = [0 -30 30];
    mmd = zeros(3,1); smd = zeros(3,1);
    for i = [2 1 3]
        mmd(i) = nanmean(1000*reshape([MaxVelc(i+3*1,1,5:-1:3,1) MaxVelc(i+3*1,1,3:-1:1,2)],[],1));
        smd(i) = nanmean(1000*reshape([MaxVelc(i+3*1,2,5:-1:3,1) MaxVelc(i+3*1,2,3:-1:1,2)],[],1));
        errorbar(HR(i), mmd(i),smd(i)/sqrt(18),'Color',colors{i},'Marker','x')        
    end
    
    plot(HR([2 1 3]),mmd([2 1 3]),':b')
    %plot(repmat([-30;0;30],1,18),NC([2 1 3],:),'--')
    xlabel('Head Rotations (deg)'); ylabel('Reaction time (ms)')
    title('Modulation of Maximum Velocity with varying head roll')

    figname = 'MDs-OBs-Whole'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    %% Plot Percentage of rightward movements
    figure(7); hold on               
    rl = zeros(3,5,2);
    h = zeros(1,3);
    pr = zeros(3,5);
    for i = 3
        r = reshape(RD(i,1,:,:),5,18);r = flipud(r);
        l = reshape(RD(i,2,:,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        h(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i},'LineStyle','-'); 
        for ss = 1 : 18
            %scatter(1:5, pr(:,ss),'MarkerEdgeColor',colors{i})
            plot(1:5, pr(:,ss),'Color',colors{i},'LineStyle',':');
        end
    end        
    %legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'},'Location','best')
    %legend('Location','southeast') 
    %figure(6); hold on
    rl = zeros(3,5,2);
    h1 = zeros(1,3);
    pr = zeros(3,5);
    figure(6); hold on
    for i = [2 1 3]
        r = reshape(RD(i+3,1,:,:),5,18);r = flipud(r);
        l = reshape(RD(i+3,2,:,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        h1(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i},'LIneStyle',':'); 
%         for ss = 1 : 18
%             scatter(1:5, pr(:,ss),'MarkerEdgeColor',colors{i})
%             %plot(1:5, pr(:,ss),'MarkerEdgeColor',colors{i},'Marker',markers{ss},'LineStyle',':');
%         end
    end  
%     legend([h1(1,1) h1(1,2) h1(1,3)],{'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF'},'Location','best')
%     legend([h1(1,1) h1(1,2) h1(1,3) h(1,1) h(1,2) h(1,3)],...
%         {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
%         'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'},'Location','best')
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs')    

    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    %title('Movement Directions')
    xlim([0,6])
    
    figname = strcat('LR-MDir-Whole');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    %% Plot the averaged trajectories
    figure(51); hold on
    color = ['-K';  '-g'; '-r'];
    %obstaclepos = flipud(D.obstacle_posi);
    for i = [1 2 3]+0*3
        for j = 1:5
            figure(j); hold on
            for d = 1 : 2                
                tempx = trajcc(:,1,i,1,j,d) - 1*18.875;
                ys = trajcc(:,2,i,1,j,d);
                uE=tempx + trajcc(:,1,i,2,j,d);
                lE=tempx - trajcc(:,1,i,2,j,d);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(i-0*3,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tempx,ys,color(i-0*3,:),'Linewidth',2)
                scatter([obstaclepos(j,1)-18.8750 0 0],[30-obstaclepos(j,2) 9 29],'ro')  
            end
        end       
    end
    
    
end

%%
if plotting == 1
    subs = 1:18;     
    for s = 1        
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp','_',sID,'.mat');
        load(filename)
        % individual participants
        % Initialization (condition, mean/std, target, left/right)
        ima = zeros(6,2,5,2); % initilal movement angle    
        mcurv = zeros(6,2,5,2); % movement curvature
        DpassO = zeros(6,2,5,2); % Distance passing the Obstacle
        DfromO = zeros(6,2,5,2); % Distance from Obstacle
        MaxVel = zeros(6,2,5,2); % Max Velocity
        TMaxVel = zeros(6,2,5,2); % time of max verlocity wrt mvonset
        MaxAcc = zeros(6,2,5,2); % Max Acceleration
        TMaxAcc = zeros(6,2,5,2); %time of max acceleration wrt mvonset
        RT = zeros(6,2,5,2); % Reaction time
        MD = zeros(6,2,5,2); % Movement Duration
        NC = zeros(6,5,2); % number of collisions
        NR = zeros(6,5,2); % number of repeats
        NMT = zeros(6,5,2); % number of missed target

        for i = 1 : 6
            temp = params(1:Ns(i),:,i);
            rep = (temp(:,13)-temp(:,12)>1);
            for j = 1 : 2
                for t = 1 : 5
                    % no collision, no repeat, and either right to left -ward
                    % directions, good trials
                    qnc = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t; 
                    
                    ima(i,1,t,j) = mean(temp(qnc,18));
                    ima(i,2,t,j) = std(temp(qnc,18));
                    mcurv(i,1,t,j) = mean(temp(qnc,19));
                    mcurv(i,2,t,j) = std(temp(qnc,19));
                    DpassO(i,1,t,j) = mean(temp(qnc,21));
                    DpassO(i,2,t,j) = std(temp(qnc,21));
                    DfromO(i,1,t,j) = mean(temp(qnc,21));
                    DfromO(i,2,t,j) = std(temp(qnc,21));
                    MaxVel(i,1,t,j) = mean(temp(qnc,14));
                    MaxVel(i,2,t,j) = std(temp(qnc,14));
                    TMaxVel(i,1,t,j) = mean(temp(qnc,15)-temp(qnc,12));
                    TMaxVel(i,2,t,j) = std(temp(qnc,15)-temp(qnc,12));
                    MaxAcc(i,1,t,j) = mean(temp(qnc,16));
                    MaxAcc(i,2,t,j) = std(temp(qnc,16));
                    TMaxAcc(i,1,t,j) = mean(temp(qnc,17)-temp(qnc,12));
                    TMaxAcc(i,2,t,j) = std(temp(qnc,17)-temp(qnc,12));
                    RT(i,1,t,j) = mean(temp(qnc,12)-temp(qnc,11));
                    RT(i,2,t,j) = std(temp(qnc,12)-temp(qnc,11));
                    MD(i,1,t,j) = mean(temp(qnc,13)-temp(qnc,12));
                    MD(i,2,t,j) = std(temp(qnc,13)-temp(qnc,12));
                    qc = rep == 0 & temp(:,5) == 1 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t; 
                    NC(i,t,j) = sum(qc);
                    qr =  temp(:,5) == 0  &...
                          temp(:,7) == (j-1) & temp(:,8) == t; 
                    NR(i,t,j) = sum(qr);
                    qmt =  temp(:,10) == 1 &...
                          temp(:,7) == (j-1) & temp(:,8) == t; 
                    NR(i,t,j) = sum(qmt);                    
                end
            end
        end
        %% Plot Initial movement angles
        % for each obstacle
        figure(12); subplot(2,1,1);hold on  
        xa = 5:-1:1;            
        h1 = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(ima(i,1,:,1),5,1));
            t2 = flipud(reshape(ima(i,2,:,1),5,1));
            h1(1,i)=errorbar(1:5,abs(t1-90),t2,'Color',colors{i});

        end  
        %legend([h1(1,1) h1(1,2) h1(1,3)],{'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
        %legend('Location','southeast')    
%         ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
%         plot(xa,[nan nan nan nan nan])
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'','','','',''}))
%         xlim([0,6])
%         title('Initial movement angle modulation by varying head roll - Rightward MVs')

        figure(12); subplot(2,1,2);hold on
        for i = [2 1 3]
            t1 = flipud(reshape(ima(i,1,:,2),5,1));
            t2 = flipud(reshape(ima(i,2,:,2),5,1));
            errorbar(1:5,abs(t1-90),t2,'Color',colors{i});

        end              
%         xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
%         plot(xa,[nan nan nan nan nan])
%         xlim([0,6])
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
%         title('Initial movement angle modulation by varying head roll - Leftward MVs')

%         figname = strcat(sID,'-IMA-OBs'); 
%         savefig(gcf,figname) 
%         saveas(gcf,figname,'tiff')
                
        figure(12); subplot(2,1,1);hold on  
        xa = 5:-1:1;            
        h = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(ima(i+3,1,:,1),5,1));
            t2 = flipud(reshape(ima(i+3,2,:,1),5,1));
            h(1,i)=errorbar(1:5,abs(t1-90),t2,'Color',colors{i},'LineStyle',':');

        end  
        legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
            {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
             'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
        legend('Location','southeast')    
        ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
        plot(xa,[nan nan nan nan nan])
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'','','','',''}))
        xlim([0,6])
        title('Initial movement angle modulation by varying head roll - Rightward MVs')

        figure(12); subplot(2,1,2);hold on
        for i = [2 1 3]
            t1 = flipud(reshape(ima(i+3,1,:,2),5,1));
            t2 = flipud(reshape(ima(i+3,2,:,2),5,1));
            errorbar(1:5,abs(t1-90),t2,'Color',colors{i},'LineStyle',':');

        end              
        xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
        plot(xa,[nan nan nan nan nan])
        xlim([0,6])
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
        title('Initial movement angle modulation by varying head roll - Leftward MVs')

        figname = strcat(sID,'-IMA-OBs'); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')

        %%
        % Plot curvature
        % for each obstacle
        figure(22); subplot(2,1,1);hold on 
        xa = 5:-1:1;
        h1 = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(mcurv(i,1,:,1,c),5,1));
            t2 = flipud(reshape(mcurv(i,2,:,1,c),5,1));
            h1(1,i)=errorbar(1:5,t1,t2,'Color',colors{i});                
        end  
%         legend([h1(1,1) h1(1,2) h1(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
%         legend('Location','southeast')           
%         ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
%         plot(xa,[nan nan nan nan nan])
%         xlim([0,6]);
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'','','','',''}))
%         title('Maximum movement curvature modulation by varying head roll - Rightward MVs')

        figure(22); subplot(2,1,2);hold on
        for i = 1:3
            t1 = flipud(reshape(mcurv(i,1,:,2,c),5,1));
            t2 = flipud(reshape(mcurv(i,2,:,2,c),5,1));
            errorbar(1:5,t1,t2,':','Color',colors{i});                
        end                
%         xlabel('Obstacle position'); ylabel('Maximum curvature (cm)')
%         plot(xa,[nan nan nan nan nan])
%         plot(xa,[nan nan nan nan nan])
%         xlim([0,6])
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
%         title('Maximum movement curvature modulation by varying head roll - Leftward MVs')
            
        figure(22); subplot(2,1,1);hold on 
        xa = 5:-1:1;
        h1 = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(mcurv(i,1,:,1,c),5,1));
            t2 = flipud(reshape(mcurv(i,2,:,1,c),5,1));
            h1(1,i)=errorbar(1:5,t1,t2,':','Color',colors{i});                
        end  
        legend([h1(1,1) h1(1,2) h1(1,3) h(1,1) h(1,2) h(1,3)],...
            {'HR = 30CCW - F','HR = 0 - F','HR = 30CW - F',...
            'HR = 30CCW - NF','HR = 0 - NF','HR = 30CW - NF'})
        legend('Location','southeast')           
        ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
        plot(xa,[nan nan nan nan nan])
        xlim([0,6]);
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'','','','',''}))
        title('Maximum movement curvature modulation by varying head roll - Rightward MVs')

        figure(22); subplot(2,1,2);hold on
        for i = [2 1 3]
            t1 = flipud(reshape(mcurv(i+3,1,:,2),5,1));
            t2 = flipud(reshape(mcurv(i+3,2,:,2),5,1));
            errorbar(1:5,t1,t2,'Color',colors{i});                
        end                
        xlabel('Obstacle position'); ylabel('Maximum curvature (cm)')
        plot(xa,[nan nan nan nan nan])
        plot(xa,[nan nan nan nan nan])
        xlim([0,6])
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
        title('Maximum movement curvature modulation by varying head roll - Leftward MVs')
        
        figname = strcat(sID,cond{c,1},'-MC-OBS');
        savefig(gcf,figname)
        saveas(gcf,figname,'tiff')

        %%
        % Plot path distance from obstacle            
        % for each obstacle
        figure(32); subplot(2,1,1);hold on       
        xa = 5:-1:1;
        h = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(DfromO(i,1,:,1,c),5,1));
            t2 = flipud(reshape(DfromO(i,2,:,1,c),5,1));
            h(1,i)=errorbar(1:5,abs(t1),t2,'Color',colors{i});                
        end
%         legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
%         legend('Location','southeast')           
%         ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
%         plot(xa,[nan nan nan nan nan])
%         xlim([0,6]);
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'','','','',''})) 
%         title('Distance form obstacle modulation by varying head roll - Rightawrd MVs')

        figure(32); subplot(2,1,2);hold on
        for i = [2 1 3]
            t1 = flipud(reshape(DfromO(i,1,:,2,c),5,1));
            t2 = flipud(reshape(DfromO(i,2,:,2,c),5,1));
            errorbar(1:5,abs(t1),t2,'Color',colors{i});                
        end
%         xlabel('Obstacle position');ylabel('Distance from Obstacle'); % xlabel('Obstacle position');
%         plot(xa,[nan nan nan nan nan])
%         xlim([1,6]);
%         set(gca,'XTick',fliplr(xa))
%         set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))            
%         title('Distance form obstacle modulation by varying head roll - Leftward MVs')            
%         
        figure(32); subplot(2,1,1);hold on       
        xa = 5:-1:1;
        h1 = zeros(1,3);
        for i = [2 1 3]
            t1 = flipud(reshape(DfromO(i+3,1,:,1),5,1));
            t2 = flipud(reshape(DfromO(i+3,2,:,1),5,1));
            h1(1,i)=errorbar(1:5,abs(t1),t2,'Color',colors{i});                
        end
        legend([h1(1,1) h1(1,2) h1(1,3) h(1,1) h(1,2) h(1,3)],...
            {'HR = 30CCW - NF','HR = 0 - NF','HR = 30CW - NF',...
            'HR = 30CCW - F','HR = 0 - F','HR = 30CW - F'})
        legend('Location','southeast')           
        ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
        plot(xa,[nan nan nan nan nan])
        xlim([0,6]);
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'','','','',''})) 
        title('Distance form obstacle modulation by varying head roll - Rightawrd MVs')

        figure(32); subplot(2,1,2);hold on
        for i = [2 1 3]
            t1 = flipud(reshape(DfromO(i+3,1,:,2),5,1));
            t2 = flipud(reshape(DfromO(i+3,2,:,2),5,1));
            errorbar(1:5,abs(t1),t2,'Color',colors{i});                
        end
        xlabel('Obstacle position');ylabel('Distance from Obstacle'); % xlabel('Obstacle position');
        plot(xa,[nan nan nan nan nan])
        xlim([1,6]);
        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))            
        title('Distance form obstacle modulation by varying head roll - Leftward MVs')            

        figname = strcat(sID,cond{c,1},'-pdfo-OBS');
        savefig(gcf,figname)
        saveas(gcf,figname,'tiff')

        %%
        % Plot Number of collisition
        figure(4); subplot(1,2,1); hold on
        HR = [-30 0 30];
        for i = 1:3
            bar(HR(i), sum(NC(i,:,1,c)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
        end
        xlabel('Head Rotations (deg)'); ylabel('Obstacle collision - Rightward MVs')
        figure(4); subplot(1,2,2); hold on
        for i = 1:3
            bar(HR(i), sum(NC(i,:,2,c)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
        end
        xlabel('Head Rotations (deg)'); ylabel('Obstacle collision - Leftward MVs')
        supertitle('Obstacle collision modulation by varying head roll')

%             figname = strcat(sID,cond{c,1},'-NC');
%             savefig(gcf,figname)
%             saveas(gcf,figname,'tiff')
% 
%             % for each obstacle
%             figure(42); hold on
%             xa = [250 150 50 -50 -150];
%             xad = [5.5 15.5 25.5];
%             colorss = [204 255 204; 204 204 204; 255 102 102]./255;
%             h = zeros(2,3);
%             for i = 1:3
%                 t1 = reshape(NC(i,:,1,c),5,1);
%                 t2 = reshape(NC(i,:,2,c),5,1);
%                 h(1,i)=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
%                 h(2,i)=bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
%             end  
%             legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
%                          'HR = 0-R','HR = 0-L',...
%                          'HR = 30CW-R','HR = 30CW-L',...
%                     })
%             legend('Location','northeastoutside')    
%             ylabel('Nomber of collisions');  xlabel('Obstacle position');
%             plot(xa,[nan nan nan nan nan])
%             set(gca,'XTick',fliplr(xa))
%             set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
%             title('Variation of collision number by varying head roll')    
% 
% 
%             figname = strcat(sID,cond{c,1},'-NC-OBS');
%             savefig(gcf,figname)
%             saveas(gcf,figname,'tiff')
% 
        %%   
        % Plot Number of right vs left -ward movements 
        figure(5); hold on               
        rl = zeros(3,5,2);
        h = zeros(1,3);
        pr = zeros(3,5);
        for i = 1:3
            for j = 1 : 5
                qnc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 0 & Comp{i}(:,8) == j ; % no collision, no repeat
                rl(i,j,1) = sum(Comp{i}(qnc,7)==0);            
                rl(i,j,2) = sum(Comp{i}(qnc,7)==1); 
                pr(i,j) = rl(i,j,1)*100/(rl(i,j,1)+rl(i,j,2));
            end 

            h(1,i)=plot(1:5,fliplr(pr(i,:)),'Color',colors{i});                
        end        
        legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
        legend('Location','southeast') 
        xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs')    

        set(gca,'XTick',fliplr(xa))
        set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
        title('Modulation of movement direction by rolling the head')
        xlim([0,6])
        figname = strcat(sID,cond{c,1},'LR-MDir');
        savefig(gcf,figname)
        saveas(gcf,figname,'tiff')
        close all

    end
end
%%
close all
%%
if trajectories == 1
    colorss = [204 255 204; 204 204 204; 255 102 102]./255;
    files = dir;
    trajx = 1000*ones(1000,300,3);
    trajy = 1000*ones(1000,300,3);
    time = 1000*ones(1000,300,3);
    l = zeros(300,3);
    o = l;
    td = o;
    tn = td;
    for j = 1 : length(files)
        k = strfind(files(j).name,'_marked');
        color = [ '-g'; '-k'; '-r'];
%         figure(10); hold on; % row trajectories        
        for i = 1 : 3
            figure(50); hold on
            % load the data        
            if ~isempty(k)
                k1 = strfind(files(j).name,['9' cond1{i}]);
                if ~isempty(k1)   
                    load(files(j).name)
                    load('CompF_s9.mat')
                    tti = 1;
                    for tt = 1 : N
                        tt
%                         ontime = D.handON{1,tt};
                        if param(tt,24) == 0 && param(tt,25) == 0 && param(tt,27) == 0 && Comp{i,1}(tt,9)==0 % select good trials                            
                            on = find(D.t{tt,1}(:) > D.handON{1,tt},1);
                            off = size(D.handlepos{tt,1},1);
                            trajx(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:end,1);
                            trajy(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:end,2);
%                             time(1:off-on+1,tti,i) = D.t{tt,1}(on:end,2);
%                             figure(10); 
%                             H=plot(trajx(1:off-on+1,tti,i),trajy(1:off-on+1,tti,i),color(i,:));
%                             H.Color(4) = 0.2;
                            l(tti,i) = off-on+1;                            
                            tti = tti + 1;
                            figure(5*i);hold on
                            plot(D.handscreen{1,tt}(on:end,1),D.handscreen{1,tt}(on:end,2),'Color',colorss(i,:))
                        end
                    end 
                    trialN = 1:N;
                    q = param(:,24) == 0 & param(:,25) == 0 & param(:,27) == 0 & Comp{i,1}(:,9)==0; % select good trials
                    o(1:tti-1,i) = D.obstacleorder(q);
                    td(1:tti-1,i) = D.trajdir(q);
                    tn(1:tti-1,i) = trialN(q);
                end
            end
        end        
    end
    q = l ~= 0;
%     s = min(l(q));
%     s = floor(s); % normalization factor
    s = max(l(q)) + 10; % fixed normalization for the pooled data
    mtrajx = zeros(s,5,2,3); % mean trajectories along x-axis
    strajx = zeros(s,5,2,3); % strandard deviation of trajectories along x-axis
    figure(5); hold on
    ys = 9.5:(29-9.5)/(s-1):29; 
    ys = ys';
%     color = [ '-g'; '-k'; '-r'];
    for i = 1:3
        [mtrajx(1:s,:,:,i),strajx(1:s,:,:,i),scheck] = trajnormalization2_2(trajx(:,:,i),trajy(:,:,i),l(:,i),o(:,i),td(:,i),tn(:,i),s);         
        for j = 1:5
            %figure(5);subplot(1,5,j); hold on
            figure(j);hold on
            if sum(strajx(:,j,1,i))~=0
                tttempt = mtrajx(:,j,1,i) - 18.875;
                uE=tttempt+strajx(:,j,1,i);
                lE=tttempt-strajx(:,j,1,i);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(i,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tttempt,ys,color(i,:),'Linewidth',2)
%                 shadedErrorBar(mtrajx(:,j,1,i),ys,strajx(:,j,1,i),'lineprops',color(i,:),'patchSaturation',0.33)            
            else
                %plot(mtrajx(:,j,1,i),ys,color(i,:))
            end
            if sum(strajx(:,j,2,i))~=0
                tttempt = mtrajx(:,j,2,i) - 18.875;
                uE=tttempt+strajx(:,j,2,i);
                lE=tttempt-strajx(:,j,2,i);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(i,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tttempt,ys,color(i,:),'Linewidth',2)
%                 shadedErrorBar(mtrajx(:,j,2,i),ys,strajx(:,j,2,i),'lineprops',color(i,:),'patchSaturation',0.33)            
            else
                %plot(mtrajx(:,j,2,i),ys,color(i,:))
            end            
        end        
    end  
    obstaclepos = flipud(D.obstacle_posi);    
    for i = 1 : 5
        figname = 's9F-trajnormal';
        %figure(5);subplot(1,5,i);set(gca,'xlim',[-20 20])
        figure(i);set(gca,'xlim',[-15 15])
        scatter([obstaclepos(i,1)-18.8750 0 0],[30-obstaclepos(i,2) 9 29],'ro')    
        fign = num2str(i);
        figname = [figname fign];
        savefig(gcf,figname)
        saveas(gcf,figname,'tiff')
    end      
    save('trajFnormal-s9F.mat','mtrajx','strajx')

end
close all