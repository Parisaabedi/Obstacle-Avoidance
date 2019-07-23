%% Obstacle Avoidacne Data analysis
% Ploting setup 
% Written/editted by PA (July 2019)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'};
color = [ '-k'; '-g'; '-b'];
colors = { 'k' ; 'g' ; 'b'};
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
obstaclepos = [22.4750   11.0000
               20.6750   11.0000
               18.8750   11.0000
               17.0750   11.0000
               15.2750   11.0000
               ];
%obstaclepos = flipud(obstaclepos);
sn = 2000;
init = 0;

%% Different analysis
pooled = 0; %% Pooling the data
poolploting = 1; %% Ploting the pooled data
plotting = 0; %% Ploting different conditions
trajectories = 0; %% Plotting trajectories

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
    trajc = zeros(sn+1,2,6,2,5,2,18); % Normalized trajectories
    velc = zeros(sn+1,2,6,2,5,2,18); % Normalized velocity/acc profiles
    RTs = zeros(6,9,5,2,18); % Reaction time (mean&median%std,extrap, vel thre, acc thre) 
    RTs2 = zeros(6,6,5,2,18); % Reaction time (mean&median over averaged traj,extrap, vel thre, acc thre) 
    MDs = zeros(6,9,5,2,18); % Movement Duration (mean&median%std,extrap, vel thre, acc thre)     
    MDs2 = zeros(6,6,5,2,18); % Movement Duration (mean&median over averaged traj,extrap, vel thre, acc thre)
    RD = zeros(6,5,2,18); % Right MV. Direction (Cond, targets, right(0)/left(1), subs)
    CDir = zeros(6,5,2,18); % Changed Directions (Cond, targets, right(0)/left(1), subs)
    NCd = zeros(6,5,2,18); % number of collisions (for each conditoins, targets,....)
    URT = zeros(6,5,6,2,18); % Unacceptable RTs (for each conditoins, targets,different measurement methods,....)    
    NR = zeros(6,5,2,18); % number of repeats
    NMT = zeros(6,5,2,18); % number of missed target
    
    subs = 1:18;     
    for s = 1:18        
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp','_',sID,'.mat');
        load(filename)       
        for i = 1 : 6
            temp = params(1:Ns(i),:,i);
            ttraj = ntraj(:,1:Ns(i),:,:);
            rep = (temp(:,13)-0*temp(:,12)>1);
            for j = 1 : 2
                for t = 1 : 5
                    % no collision, no repeat, and either right to left -ward
                    % directions, good trials, and no changed direction (measured)
                    qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0; 
                    
                    ima(i,1,t,j,s) = mean(temp(qnc,18));
                    ima(i,2,t,j,s) = std(temp(qnc,18));
                    mcurv(i,1,t,j,s) = mean(temp(qnc,19));
                    mcurv(i,2,t,j,s) = std(temp(qnc,19));
                    DpassO(i,1,t,j,s) = mean(temp(qnc,21));
                    DpassO(i,2,t,j,s) = std(temp(qnc,21));
                    DfromO(i,1,t,j,s) = mean(temp(qnc,21));
                    DfromO(i,2,t,j,s) = std(temp(qnc,21));
                    MaxVel(i,1,t,j,s) = nanmedian(temp(qnc,14));
                    MaxVel(i,2,t,j,s) = nanstd(temp(qnc,14));
                    TMaxVel(i,1,t,j,s) = mean(temp(qnc,15));% wrt target ONset
                    TMaxVel(i,2,t,j,s) = std(temp(qnc,15));% wrt target ONset
                    MaxAcc(i,1,t,j,s) = mean(temp(qnc,16));
                    MaxAcc(i,2,t,j,s) = std(temp(qnc,16));
                    TMaxAcc(i,1,t,j,s) = mean(temp(qnc,17));% wrt target ONset
                    TMaxAcc(i,2,t,j,s) = std(temp(qnc,17));% wrt target ONset
                    if sum(qnc) ~= 0
                        trajc(:,1,i,1,t,j,s) = nanmedian(ttraj(:,qnc,1,i),2); % mean xtraj
                        trajc(:,1,i,2,t,j,s) = nanstd(ttraj(:,qnc,1,i),0,2); % std xtraj
                        trajc(:,2,i,1,t,j,s) = nanmedian(ttraj(:,qnc,2,i),2); % mean ytraj
                        trajc(:,2,i,2,t,j,s) = nanstd(ttraj(:,qnc,2,i),0,2); % std ytraj
                        
                        % velocity
                        velc(:,1,i,1,t,j,s) = nanmean(nvelacc(:,qnc,1,i),2); % mean velocity
                        velc(:,1,i,2,t,j,s) = nanstd(nvelacc(:,qnc,1,i),0,2); % std velocity
                        velc(:,2,i,1,t,j,s) = nanmean(nvelacc(:,qnc,2,i),2); % mean acceleration
                        velc(:,2,i,2,t,j,s) = nanstd(nvelacc(:,qnc,2,i),0,2); % std acceleration
                        
                        % RT calculations
                        % Velocity Calculation
                        win = 0.005;
                        sfr = sn/nanmean(temp(qnc,13));
                        handXv = GRAIdiff(1/sfr,win,trajc(:,1,i,1,t,j,s));
                        handYv = GRAIdiff(1/sfr,win,trajc(:,2,i,1,t,j,s));
                        vel = sqrt(handXv.^2 + handYv .^ 2);
                        
                        % Acceleration Calculation
                        handXa = GRAIdiff(1/sfr,win,handXv);
                        handYa = GRAIdiff(1/sfr,win,handYv);
                        acc = sqrt(handXa.^2 + handYa .^ 2);
                        
                        % Threshold method
                        indv = find(vel > 0.1 * max(vel));
                        posv = GRAIgroup(indv,2,4);
                        posonv = indv(posv(1));
                        RTs2(i,2,t,j,s) = mean(temp(qnc,11)+posonv*(1/sfr));
                        if RTs2(i,2,t,j,s) < 0.1
                            URT(i,t,5,j,s) = URT(i,t,5,j,s) + 1;
                        end
                        inda = find(acc > 0.1 * max(acc));
                        posa = GRAIgroup(inda,2,5);
                        posona = inda(posa(1));
                        RTs2(i,3,t,j,s) = mean(temp(qnc,11)+posona*(1/sfr));
                        if RTs2(i,3,t,j,s) < 0.1
                            URT(i,t,6,j,s) = URT(i,t,6,j,s) + 1;
                        end
                        
                        % Extrapolation method (Based Smeets and Brenner, 2018) 
                        ts = 0:1/sfr:mean(temp(qnc,13));
                        Rarrv = find(ts >= 0, 1 ):find(ts <= 0+0.15, 1, 'last' );
                        x1 = [0 1 5]; y1 = repmat(nanmean(vel(Rarrv)),3,1); % the baseline calculation 
                        y2 = [0.25*max(vel) 0.75* max(vel)];
                        point1 = find(vel >= 0.25*max(vel));point2 = find(vel >= 0.75*max(vel));
                        x2 = [ts(point1(1)) ts(point2(1))];
                        % draw the line between two points (y = m*x+n)
                        m = (y2(2)-y2(1))/(x2(2)-x2(1)); n = y2(1)-m*x2(1);
                        x2p = [0 1 5];y2p=m.*x2p + n;
                        [tx,ty] = polyxpoly(x1,y1,x2p,y2p); % The intersection of the two lines
                        
                        if isempty(tx)
                            RTs2(i,1,t,j,s) = -1;
                            URT(i,t,4,j,s) = URT(i,t,4,j,s) + 1;
                        else
                            RTs2(i,1,t,j,s) = tx;
                        end
                        MDs2(i,1,t,j,s) = nanmean(temp(qnc,13));%
                        MDs2(i,2,t,j,s) = nanmean(temp(qnc,13));%
                        MDs2(i,3,t,j,s) = nanmean(temp(qnc,13));%
                    else
                        trajc(:,1,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean xtraj
                        trajc(:,1,i,2,t,j,s) = zeros(sn+1,1)*nan; % std xtraj
                        trajc(:,2,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean ytraj
                        trajc(:,2,i,2,t,j,s) = zeros(sn+1,1)*nan; % std ytraj 
                        
                        velc(:,1,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean vel
                        velc(:,1,i,2,t,j,s) = zeros(sn+1,1)*nan; % std vel
                        velc(:,2,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean acc
                        velc(:,2,i,2,t,j,s) = zeros(sn+1,1)*nan; % std acc  
                    end
                    
                    % Extrapolation method
                    qnc1 = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t &...
                          temp(:,24) == 0 & temp(:,12)> 0.1;
                    RTs(i,1,t,j,s) = mean(temp(qnc1,12));
                    RTs(i,2,t,j,s) = std(temp(qnc1,12));
                    RTs(i,3,t,j,s) = median(temp(qnc1,12));
                    URT(i,t,1,j,s) = sum(qnc1);
                    
                    % Vel Threshold
                    qnc2 = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t &...
                          temp(:,24) == 0 & temp(:,22)> 0.1;
                    RTs(i,4,t,j,s) = mean(temp(qnc,22));
                    RTs(i,5,t,j,s) = std(temp(qnc,22));
                    RTs(i,6,t,j,s) = median(temp(qnc,22)); 
                    URT(i,t,2,j,s) = sum(qnc2);
                    
                    % Acc Threshold
                    qnc3 = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                          temp(:,7) == (j-1) & temp(:,8) == t &...
                          temp(:,24) == 0 & temp(:,23)> 0.1;
                    RTs(i,7,t,j,s) = mean(temp(qnc,23));
                    RTs(i,8,t,j,s) = std(temp(qnc,23));
                    RTs(i,9,t,j,s) = median(temp(qnc,23));
                    URT(i,t,3,j,s) = sum(qnc3);
                    
                    % Extrapolation method
                    MDs(i,1,t,j,s) = nanmean(temp(qnc1,13));
                    MDs(i,2,t,j,s) = std(temp(qnc1,13));
                    MDs(i,3,t,j,s) = nanmedian(temp(qnc1,13));
                    
                    % Vel Threshold
                    MDs(i,4,t,j,s) = nanmean(temp(qnc2,13));
                    MDs(i,5,t,j,s) = std(temp(qnc2,13));
                    MDs(i,6,t,j,s) = nanmedian(temp(qnc2,13));
                    
                    % Acc Threshold
                    MDs(i,7,t,j,s) = nanmean(temp(qnc3,13));
                    MDs(i,8,t,j,s) = std(temp(qnc3,13));
                    MDs(i,9,t,j,s) = nanmedian(temp(qnc3,13));
                                      
                    qrc = temp(:,6) == 0 & temp(:,7) == (j-1) & temp(:,8) == t; % Rightward direction
                    RD(i,t,j,s) = sum(qrc);
                    
                    qncc = temp(:,5) == 1 & temp(:,7) == (j-1) & temp(:,8) == t;% Collision rate
                    NCd(i,t,j,s) = sum(qncc);
                    
                    qcdir = temp(:,5) == 0 & temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 1;% Changed Direction
                    CDir(i,t,j,s) = sum(qcdir);
                    
                    qcdl = temp(:,5) == 0 & temp(:,7) == (j-1) & temp(:,8) == t;% Changed Direction
                    
                    NR(i,t,j,s) = sum(temp(qcdl,4));             
                    NMT(i,t,j,s) = sum(temp(qcdl,10));
                end 
                
            end             
            
        end
    end
    save('WholeComp.mat','ima','mcurv','DpassO','DfromO','MaxVel','TMaxVel',...
         'MaxAcc','TMaxAcc','trajc','RTs','MDs','NCd','NR','NMT','RD',...,
         'CDir','RTs2','MDs2','velc')
    
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
    TMaxVelc = zeros(6,4,5,2); % time of max verlocity wrt target onset
    MaxAccc = zeros(6,4,5,2); % Max Acceleration
    TMaxAccc = zeros(6,4,5,2); %time of max acceleration wrt target onset
    trajcc = zeros(sn+1,2,6,4,5,2); % Normalized trajectories(samples,x/y, cond,mean/std,obst, right/left)
    velcc = zeros(sn+1,2,6,4,5,2); % Normalized vel/acc profiles(samples,vel/acc, cond,mean/std,obst, right/left)
    RTc = zeros(6,2,9,5,2); % Reaction time
    MDc = zeros(6,2,9,5,2); % Movement Duration
    RTc2 = zeros(6,2,6,5,2); % Reaction time
    MDc2 = zeros(6,2,6,5,2); % Movement Duration
    NCdc = zeros(6,2,5,2); % Numver os collisions (cond, mean/std, obst, right/left)
    NCc = zeros(6,2,18); % number of collisions(condition,mean/std, subject)
    NRc = zeros(6,2,18); % number of repeats(condition,mean/std, subject)
    NMTc = zeros(6,2,18); % number of missed target(condition,mean/std, subject)
    
    % Standardize trajectories
    variability = zeros(6,5,2,18);% cond,tar,dir,sub
    
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
                end
                rrd = 100*reshape(RD(i,j,:,s),[],1)./sum(reshape(RD(i,j,:,s),[],1));
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
                ecc = reshape(EC(i,t,d,:),18,1);
                ecc(ecc == 0) = nan;
                imac(i,1,t,d) = nanmean(ecc.*reshape(ima(i,1,t,d,:),18,1));
                imac(i,2,t,d) = nanstd(ecc.*reshape(ima(i,1,t,d,:),18,1));
                imac(i,3,t,d) = nanmean(ecc.*reshape(ima(i,2,t,d,:),18,1));
                imac(i,4,t,d) = nanstd(ecc.*reshape(ima(i,2,t,d,:),18,1));
                
                mcurvc(i,1,t,d) = nanmean(ecc.*reshape(mcurv(i,1,t,d,:),18,1));
                mcurvc(i,2,t,d) = nanstd(ecc.*reshape(mcurv(i,1,t,d,:),18,1));
                mcurvc(i,3,t,d) = nanmean(ecc.*reshape(mcurv(i,2,t,d,:),18,1));
                mcurvc(i,4,t,d) = nanstd(ecc.*reshape(mcurv(i,2,t,d,:),18,1));
                
                DpassOc(i,1,t,d) = nanmean(ecc.*reshape(DpassO(i,1,t,d,:),18,1));
                DpassOc(i,2,t,d) = nanstd(ecc.*reshape(DpassO(i,1,t,d,:),18,1));
                DpassOc(i,3,t,d) = nanmean(ecc.*reshape(DpassO(i,2,t,d,:),18,1));
                DpassOc(i,4,t,d) = nanstd(ecc.*reshape(DpassO(i,2,t,d,:),18,1));
                
                DfromOc(i,1,t,d) = nanmean(ecc.*reshape(DfromO(i,1,t,d,:),18,1));
                DfromOc(i,2,t,d) = nanstd(ecc.*reshape(DfromO(i,1,t,d,:),18,1));
                DfromOc(i,3,t,d) = nanmean(ecc.*reshape(DfromO(i,2,t,d,:),18,1));
                DfromOc(i,4,t,d) = nanstd(ecc.*reshape(DfromO(i,2,t,d,:),18,1));
                
                MaxVelc(i,1,t,d) = nanmean(ecc.*reshape(MaxVel(i,1,t,d,:),18,1));
                MaxVelc(i,2,t,d) = nanstd(ecc.*reshape(MaxVel(i,1,t,d,:),18,1));
                MaxVelc(i,3,t,d) = nanmean(ecc.*reshape(MaxVel(i,2,t,d,:),18,1));
                MaxVelc(i,4,t,d) = nanstd(ecc.*reshape(MaxVel(i,2,t,d,:),18,1));
                
                TMaxVelc(i,1,t,d) = nanmean(ecc.*reshape(TMaxVel(i,1,t,d,:),18,1));
                TMaxVelc(i,2,t,d) = nanstd(ecc.*reshape(TMaxVel(i,1,t,d,:),18,1));
                TMaxVelc(i,3,t,d) = nanmean(ecc.*reshape(TMaxVel(i,2,t,d,:),18,1));
                TMaxVelc(i,4,t,d) = nanstd(ecc.*reshape(TMaxVel(i,2,t,d,:),18,1));
                
                MaxAccc(i,1,t,d) = nanmean(ecc.*reshape(MaxAcc(i,1,t,d,:),18,1));
                MaxAccc(i,2,t,d) = nanstd(ecc.*reshape(MaxAcc(i,1,t,d,:),18,1));
                MaxAccc(i,3,t,d) = nanmean(ecc.*reshape(MaxAcc(i,2,t,d,:),18,1));
                MaxAccc(i,4,t,d) = nanstd(ecc.*reshape(MaxAcc(i,2,t,d,:),18,1));
                
                TMaxAccc(i,1,t,d) = nanmean(ecc.*reshape(TMaxAcc(i,1,t,d,:),18,1));
                TMaxAccc(i,2,t,d) = nanstd(ecc.*reshape(TMaxAcc(i,1,t,d,:),18,1));
                TMaxAccc(i,3,t,d) = nanmean(ecc.*reshape(TMaxAcc(i,2,t,d,:),18,1));
                TMaxAccc(i,4,t,d) = nanstd(ecc.*reshape(TMaxAcc(i,2,t,d,:),18,1));
                
                ecctraj = repmat(ecc',sn+1,1);
                trajcc(:,1,i,1,t,d) = nanmean(ecctraj.*reshape(trajc(:,1,i,1,t,d,:),[sn+1,18]),2);
                trajcc(:,1,i,2,t,d) = nanstd(ecctraj.*reshape(trajc(:,1,i,1,t,d,:),[sn+1,18]),0,2);
                trajcc(:,1,i,3,t,d) = nanmean(ecctraj.*reshape(trajc(:,1,i,2,t,d,:),[sn+1,18]),2);
                trajcc(:,1,i,4,t,d)= nanstd(ecctraj.*reshape(trajc(:,1,i,2,t,d,:),[sn+1,18]),0,2);
                
                trajcc(:,2,i,1,t,d) = nanmean(ecctraj.*reshape(trajc(:,2,i,1,t,d,:),[sn+1,18]),2);
                trajcc(:,2,i,2,t,d) = nanstd(ecctraj.*reshape(trajc(:,2,i,1,t,d,:),[sn+1,18]),0,2);
                trajcc(:,2,i,3,t,d) = nanmean(ecctraj.*reshape(trajc(:,2,i,2,t,d,:),[sn+1,18]),2);
                trajcc(:,2,i,4,t,d)= nanstd(ecctraj.*reshape(trajc(:,2,i,2,t,d,:),[sn+1,18]),0,2);
                
                ecctraj = repmat(ecc',sn+1,1);
                velcc(:,1,i,1,t,d) = nanmean(ecctraj.*reshape(velc(:,1,i,1,t,d,:),[sn+1,18]),2);
                velcc(:,1,i,2,t,d) = nanstd(ecctraj.*reshape(velc(:,1,i,1,t,d,:),[sn+1,18]),0,2);
                velcc(:,1,i,3,t,d) = nanmean(ecctraj.*reshape(velc(:,1,i,2,t,d,:),[sn+1,18]),2);
                velcc(:,1,i,4,t,d)= nanstd(ecctraj.*reshape(velc(:,1,i,2,t,d,:),[sn+1,18]),0,2);
                
                velcc(:,2,i,1,t,d) = nanmean(ecctraj.*reshape(velc(:,2,i,1,t,d,:),[sn+1,18]),2);
                velcc(:,2,i,2,t,d) = nanstd(ecctraj.*reshape(velc(:,2,i,1,t,d,:),[sn+1,18]),0,2);
                velcc(:,2,i,3,t,d) = nanmean(ecctraj.*reshape(velc(:,2,i,2,t,d,:),[sn+1,18]),2);
                velcc(:,2,i,4,t,d)= nanstd(ecctraj.*reshape(velc(:,2,i,2,t,d,:),[sn+1,18]),0,2);
                
                for rtc = 1 : 6
                    RTc(i,1,rtc,t,d) = nanmean(ecc.*reshape(RTs(i,rtc,t,d,:),18,1));
                    RTc(i,2,rtc,t,d) = nanstd(ecc.*reshape(RTs(i,rtc,t,d,:),18,1));

                    MDc(i,1,rtc,t,d) = nanmean(ecc.*reshape(MDs(i,rtc,t,d,:),18,1));
                    MDc(i,2,rtc,t,d) = nanstd(ecc.*reshape(MDs(i,rtc,t,d,:),18,1));
                    
                    RTc2(i,1,rtc,t,d) = nanmean(ecc.*reshape(RTs2(i,rtc,t,d,:),18,1));
                    RTc2(i,2,rtc,t,d) = nanstd(ecc.*reshape(RTs2(i,rtc,t,d,:),18,1));

                    MDc2(i,1,rtc,t,d) = nanmean(ecc.*reshape(MDs2(i,rtc,t,d,:),18,1));
                    MDc2(i,2,rtc,t,d) = nanstd(ecc.*reshape(MDs2(i,rtc,t,d,:),18,1));
                end
                for rtc = 7 : 9
                    RTc(i,1,rtc,t,d) = nanmean(ecc.*reshape(RTs(i,rtc,t,d,:),18,1));
                    RTc(i,2,rtc,t,d) = nanstd(ecc.*reshape(RTs(i,rtc,t,d,:),18,1));

                    MDc(i,1,rtc,t,d) = nanmean(ecc.*reshape(MDs(i,rtc,t,d,:),18,1));
                    MDc(i,2,rtc,t,d) = nanstd(ecc.*reshape(MDs(i,rtc,t,d,:),18,1));
                end
                
                NCdc(i,1,t,d) = nanmean(reshape(NCd(i,t,d,:),18,1));
                NCdc(i,2,t,d) = nanstd(reshape(NCd(i,t,d,:),18,1));
                
                variabilityc(i,1,t,d) = nanmean(ecc.*reshape(variability(i,t,d,:),18,1));
                variabilityc(i,2,t,d) = nanstd(ecc.*reshape(variability(i,t,d,:),18,1));
            end
        end
        NCc(i,1,:) = sum(reshape(NCd(i,:,:,:),[10,18]));
        NCc(i,2,:) = std(reshape(NCd(i,:,:,:),[10,18]));
        
        NRc(i,1,:) = sum(reshape(NR(i,:,:,:),[10,18]));
        NRc(i,2,:) = std(reshape(NR(i,:,:,:),[10,18]));
        
        NMTc(i,1,:) = sum(reshape(NMT(i,:,:,:),[10,18]));
        NMTc(i,2,:) = std(reshape(NMT(i,:,:,:),[10,18]));
    end
    
    % Expansion Biases
    % Initialization
    % 1 - > Most leftward obstacle
    % mirrored obstacles are matched to have both rightward and leftward MVs (3L&3R, 1&5, 2&4)
    alfabeta = zeros(3,2,2,18)*nan; %  6 obstacles (the central can have both directions),left / right, 2 parameters, #subjects
    y = 19; % point of interest along y-axis
    c = 0; % with/without visual feedback (0=> with, 1=> without)
    %valid = 1; % check if there is enough data on both right and left direction
    obs = [3 3 ; 5 1 ; 4 2];
    % Load the data  
    load('WholeComp.mat')    
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
    %             dxrhr = trajc(p3,1,3+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
    %             dxlhr = trajc(p3p,1,3+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);
    %             dxrhl = trajc(p3,1,2+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
    %             dxlhl = trajc(p3p,1,2+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);

                dxrhr = mcurv(3+c*3,1,obs(j,1),1,i) - mcurv(1+c*3,1,obs(j,1),1,i);
                dxlhr = mcurv(3+c*3,1,obs(j,2),2,i) - mcurv(1+c*3,1,obs(j,2),2,i);
                dxrhl = mcurv(2+c*3,1,obs(j,1),1,i) - mcurv(1+c*3,1,obs(j,1),1,i);
                dxlhl = mcurv(2+c*3,1,obs(j,2),2,i) - mcurv(1+c*3,1,obs(j,2),2,i);

                % calculate deltax
                dxaL = (-(dxlhr-dxlhl))/2; % dx alfa (Left)
                dxaR = (+(dxrhr-dxrhl))/2; % dx alfa (Right)
                dxbL = (+(dxlhr+dxlhl))/2; % dx beta (Left)
                dxbR = (+(dxrhr+dxrhl))/2; % dx beta (Right)

                % calculate alfa and beta
                x0r = abs(mcurv(1+c*3,1,obs(j,1),1,i));
                x0l = abs(mcurv(1+c*3,1,obs(j,2),2,i));
                tmpt = abs(trajc(:,1,:,1,:,:,:)-18.75); 
                py0r = find(tmpt(:,1,1+c*3,1,obs(j,1),1,i)==max(tmpt(:,1,1+c*3,1,obs(j,1),1,i)),1);
                py0l = find(tmpt(:,1,1+c*3,1,obs(j,2),2,i)==max(tmpt(:,1,1+c*3,1,obs(j,2),2,i)),1);
                y0r = trajc(py0r,1,1+c*3,2,obs(j,1),1,i);
                y0l = trajc(py0l,1,1+c*3,2,obs(j,2),2,i);
                alfabeta(j,1,1,i) = (dxaL/abs(dxaL))*(atand(y0l/x0l) - atand(y0l/(x0l+abs(dxaL)))); % alfa (Left)
                alfabeta(j,2,1,i) = (dxaR/abs(dxaR))*(atand(y0r/x0r) - atand(y0r/(x0r+abs(dxaR)))); % alfa (Right)
                alfabeta(j,1,2,i) =  100*dxbL/x0l; % beta (Left)            
                alfabeta(j,2,2,i) = 100*dxbR/x0r; % beta (Right)
            end
        end    
        
    end
    
    alfabetapy = zeros(6,2,18); %6 obstacles, 2 direction, 18 subjects
    s = [ 2 3 1 1 3 2;1 1 1 2 2 2]'; % to select the proper value
    for j = 1 : 6
        alfabetapy(j,1,:) = alfabeta(s(j,1),s(j,2),1,:);
        alfabetapy(j,2,:) = alfabeta(s(j,1),s(j,2),2,:);
    end
    save('RotExp_YV_final.mat','alfabeta','alfabetapy')    
    % without visual feedback
    alfabeta = zeros(3,2,2,18)*nan; %  6 obstacles (the central can have both directions),left / right, 2 parameters, #subjects   
    c = 1; % with/without visual feedback (0=> with, 1=> without)    
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
    %             dxrhr = trajc(p3,1,3+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
    %             dxlhr = trajc(p3p,1,3+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);
    %             dxrhl = trajc(p3,1,2+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
    %             dxlhl = trajc(p3p,1,2+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);

                dxrhr = mcurv(3+c*3,1,obs(j,1),1,i) - mcurv(1+c*3,1,obs(j,1),1,i);
                dxlhr = mcurv(3+c*3,1,obs(j,2),2,i) - mcurv(1+c*3,1,obs(j,2),2,i);
                dxrhl = mcurv(2+c*3,1,obs(j,1),1,i) - mcurv(1+c*3,1,obs(j,1),1,i);
                dxlhl = mcurv(2+c*3,1,obs(j,2),2,i) - mcurv(1+c*3,1,obs(j,2),2,i);

                % calculate deltax
                dxaL = (-(dxlhr-dxlhl))/2; % dx alfa (Left)
                dxaR = (+(dxrhr-dxrhl))/2; % dx alfa (Right)
                dxbL = (+(dxlhr+dxlhl))/2; % dx beta (Left)
                dxbR = (+(dxrhr+dxrhl))/2; % dx beta (Right)

                % calculate alfa and beta
                x0r = abs(mcurv(1+c*3,1,obs(j,1),1,i));
                x0l = abs(mcurv(1+c*3,1,obs(j,2),2,i));
                tmpt = abs(trajc(:,1,:,1,:,:,:)-18.75); 
                py0r = find(tmpt(:,1,1+c*3,1,obs(j,1),1,i)==max(tmpt(:,1,1+c*3,1,obs(j,1),1,i)),1);
                py0l = find(tmpt(:,1,1+c*3,1,obs(j,2),2,i)==max(tmpt(:,1,1+c*3,1,obs(j,2),2,i)),1);
                y0r = trajc(py0r,1,1+c*3,2,obs(j,1),1,i);
                y0l = trajc(py0l,1,1+c*3,2,obs(j,2),2,i);
                alfabeta(j,1,1,i) = (dxaL/abs(dxaL))*(atand(y0l/x0l) - atand(y0l/(x0l+abs(dxaL)))); % alfa (Left)
                alfabeta(j,2,1,i) = (dxaR/abs(dxaR))*(atand(y0r/x0r) - atand(y0r/(x0r+abs(dxaR)))); % alfa (Right)
                alfabeta(j,1,2,i) =  100*dxbL/x0l; % beta (Left)            
                alfabeta(j,2,2,i) = 100*dxbR/x0r; % beta (Right)
            end
        end    
        
    end
    
    alfabetap = zeros(6,2,18); %6 obstacles, 2 direction, 18 subjects
    s = [ 2 3 1 1 3 2;1 1 1 2 2 2]'; % to select the proper value
    for j = 1 : 6
        alfabetap(j,1,:) = alfabeta(s(j,1),s(j,2),1,:);
        alfabetap(j,2,:) = alfabeta(s(j,1),s(j,2),2,:);
    end
    save('RotExp_NV_final.mat','alfabeta','alfabetap') 
    
    
    %% Plot Variabilities
    % Plot averages and individuals (Final results)    
    % averaged for HRs
    mvars = zeros(18,6);    
    for i = 1 : 6
        mvars(:,i) = nanmean(reshape(EC.*variability(i,:,:,:),[],18));        
    end
    mvar = nanmean(mvars);
    mvar = [mvar; nanstd(mvars)];
    for c = 1 : 2
        figure(1); subplot(1,2,c);hold on
        for i= 1 : 18
            plot(1:3,mvars(i,[2 1 3]+(c-1)*3,1),':k','LineWidth',0.5)
        end
        errorbar(1:3,mvar(1,[2 1 3]+(c-1)*3),...
            mvar(2,[2 1 3]+(c-1)*3) ,'-k','LineWidth',2.5)
        xlabel('Head Rotations (deg)'); ylabel('Movement variability')
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'30 CCW' , '0' , '30 CW'})
        xlim([0.5 3.5])
        ylim([14 60])
    end     
    
    figname = 'Variability-averaged and individuals'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    
    
    %% Plot percentage of rightward movements
     % with visual feedback     
    figure(2); subplot(2,1,1);hold on  
    xa = 5:-1:1;
    for i = [2 1 3]
        r = reshape(RD(i,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        for ss = 1 : 18
            scatter(1:5, pr(:,ss),'MarkerEdgeColor',colors{i})
        end
        plot(1:5,nanmean(pr'),'Color',colors{i},'LineStyle',':','LineWidth',2.5);        
        
    end   
    ylabel('Percentage of Rightward MVs') 
    set(gca,'XTick',fliplr(xa))
    xlim([0,6])
    % without visual feedback
    figure(2); subplot(2,1,2);hold on   
    for i = [2 1 3]
        r = reshape(RD(i+3,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i+3,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        for ss = 1 : 18
            scatter(1:5, pr(:,ss),'MarkerEdgeColor',colors{i})
        end
        plot(1:5,nanmean(pr'),'Color',colors{i},'LineStyle',':','LineWidth',2.5);
        
    end   
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs') 
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    xlim([0,6])
    
    figname = 'RightwardMVsPer-averaged and individuals'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    
    %% Plot Collision Rate
    %Collision-Rate (Number of collision/Number of MVs per each direction)    
    tt1 = zeros(6,2,18);  
    temp = NCd./RD;
    for j = 1 :6
        for i = 1:18
            tt1(j,1,i) = nanmean(reshape(temp(j,:,:,i),[],1));
            tt1(j,2,i) = nanstd(reshape(temp(j,:,:,i),[],1));
        end   
    end 

    figure(3); subplot(1,2,1);hold on
    for i = 1:18
        plot(1:3,tt1([2 1 3],1,i),':k','LineWidth',0.5)
    end
    tt = zeros(6,2);
    for i = 1:6
        tt(i,1) = nanmean(tt1(i,1,:),3);
        tt(i,2) = nanmean(tt1(i,1,:),3);
    end
    errorbar(1:3,tt([2 1 3],1),tt([2 1 3],2),':k','LineWidth',2.5)
    xlabel('Head Rotations (deg)'); ylabel('Collision rate')
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'30CCW',' 0' , '30CW'})
    xlim([0.5,3.5])
    
    figure(3); subplot(1,2,2);hold on
    for i = 1:18
        plot(1:3,tt1([2 1 3]+3,1,i),':k','LineWidth',0.5)
    end
    
    errorbar(1:3,tt([2 1 3]+3,1),tt([2 1 3]+3,2)./sqrt(18),'-k','LineWidth',2.5)
    xlabel('Head Rotations (deg)'); ylabel('Collision rate')
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'Hr = 30CCW','HR = 0' , 'HR = 30CW'})
    xlim([0.5,3.5])

    figname = 'Collisions-averaged and individuals3'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %% Plot expansion Biases
    % mirrored obstacles are matched to have both rightward and leftward MVs (3L&3R, 1&5, 2&4)
    load('RotExp_NV_final.mat') % Without Visual feedback
    load('RotExp_YV_final.mat') % With visual feedback
    
    % sort the results 
    expansiony = reshape(alfabetapy(:,2,:),6,18); % With VF
    expansionn = reshape(alfabetap(:,2,:),6,18); % Without VF
    expansion = [expansiony;expansionn];expansion = expansion';
    expansionc = [ nanmean(expansion(:,1:3),2) nanmean(expansion(:,4:6),2)...
                   nanmean(expansion(:,7:9),2) nanmean(expansion(:,10:12),2)];
               
    rotationy = reshape(alfabetapy(:,1,:),6,18);
    rotationn = reshape(alfabetap(:,1,:),6,18);
    rotation = [rotationy;rotationn];rotation = rotation';
    rotationc = [ nanmean(rotation(:,1:3),2) nanmean(rotation(:,4:6),2)...
                  nanmean(rotation(:,7:9),2) nanmean(rotation(:,10:12),2)];
    % Plotting
    figure(4); subplot(1,2,1); hold on 
    bar([0.5 2.5 1.5 3.5],nanmean(rotationc))
    errorbar([0.5 2.5 1.5 3.5], nanmean(rotationc),nanstd(rotationc)./sqrt(18),'.')
    xlabel('Movement Direction'); ylabel('Rotational biases (deg)')
    set(gca,'XTick',[0.5 1.5 2.5 3.5])
    set(gca,'XTickLabels',{'L-VF' , 'L-NVF' , 'R-VF' , 'R-NVF'})
    
    figure(4); subplot(1,2,2); hold on 
    bar([0.5 2.5 1.5 3.5],nanmean(expansionc))
    errorbar([0.5 2.5 1.5 3.5], nanmean(expansionc),nanstd(expansionc)./sqrt(18),'.')
    xlabel('Movement Direction'); ylabel('Expansion biases (%)')
    set(gca,'XTick',[0.5 1.5 2.5 3.5])
    set(gca,'XTickLabels',{'L-VF' , 'L-NVF' , 'R-VF' , 'R-NVF'})
    
    figname = 'RotationExpansion-averaged and individuals3'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %% Plot Trajectories
    color = [ '-k'; '-g'; '-b'];
    valid = [ nan 1 ; nan 1 ; 1 1; 1 nan; 1 nan];
    obstaclepos = 18.875-[22.4750 20.6750 18.8750 17.0750 15.2750 ];
    c = 0; % with visual feedback
    for jj = 1:5 % obstacles
        figure(5); subplot(1,5,jj); hold on
        for ii = 1:3 % HR conditions
            for d = 1:2
                tempx = valid(jj,d)*trajcc(:,1,ii+c*3,1,jj,d) - 1*18.875;
                ys = valid(jj,d)*trajcc(:,2,ii+c*3,1,jj,d);
                uE=tempx + valid(jj,d)*trajcc(:,1,ii+c*3,2,jj,d)./sqrt(18);
                lE=tempx - valid(jj,d)*trajcc(:,1,ii+c*3,2,jj,d)./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tempx,ys,color(ii,:),'Linewidth',2, 'LineStyle',':')
                scatter(obstaclepos(6-jj),19,'ro')
                xlim([-5 5])
            end
        end
    end
    figname = 'MTraj2-VF-selectedobstacles'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    
    c = 1; % without visual feedback
    for jj = 1:5 % obstacles
        figure(6); subplot(1,5,jj); hold on
        for ii = 1:3 % HR conditions
            for d = 1:2
                tempx = valid(jj,d)*trajcc(:,1,ii+c*3,1,jj,d) - 1*18.875;
                ys = valid(jj,d)*trajcc(:,2,ii+c*3,1,jj,d);
                uE=tempx + valid(jj,d)*trajcc(:,1,ii+c*3,2,jj,d)./sqrt(18);
                lE=tempx - valid(jj,d)*trajcc(:,1,ii+c*3,2,jj,d)./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tempx,ys,color(ii,:),'Linewidth',2, 'LineStyle',':')
                scatter(obstaclepos(6-jj),19,'ro')
                xlim([-5 5])
            end
        end
    end
    figname = 'MTraj2-NVF-selectedobstacles'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    
end
   