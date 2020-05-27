%% Obstacle Avoidacne Data analysis
% Ploting setup 3, version3
% updated by PA (Nov. 2018)
% Some new plotting is adapted based on new analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'};
color = [ '-k'; '-g'; '-r'];
colors = { 'k' ; 'g' ; 'r'};
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
pooled = 1;
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
    NC = zeros(6,5,2,18); % number of collisions
    NR = zeros(6,5,2,18); % number of repeats
    NMT = zeros(6,5,2,18); % number of missed target
    
    subs = 1:18;     
    for s = 1:18        
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp2','_',sID,'.mat');
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
                    MaxVel(i,1,t,j,s) = mean(temp(qnc,14));
                    MaxVel(i,2,t,j,s) = std(temp(qnc,14));
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
                        MDs2(i,1,t,j,s) = nanmean(temp(qnc,13));%- RTs2(i,1,t,j,s);
                        MDs2(i,2,t,j,s) = nanmean(temp(qnc,13));%- RTs2(i,2,t,j,s);
                        MDs2(i,3,t,j,s) = nanmean(temp(qnc,13));%- RTs2(i,3,t,j,s);
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
                    MDs(i,1,t,j,s) = nanmean(temp(qnc1,13));%-temp(qnc1,12));
                    MDs(i,2,t,j,s) = std(temp(qnc1,13));%-temp(qnc1,12));
                    MDs(i,3,t,j,s) = nanmedian(temp(qnc1,13));%-temp(qnc1,12));
                    % Vel Threshold
                    MDs(i,4,t,j,s) = nanmean(temp(qnc2,13));%-temp(qnc2,22));
                    MDs(i,5,t,j,s) = std(temp(qnc2,13));%-temp(qnc2,22));
                    MDs(i,6,t,j,s) = nanmedian(temp(qnc2,13));%-temp(qnc2,22));                    
                    % Acc Threshold
                    MDs(i,7,t,j,s) = nanmean(temp(qnc3,13));%-temp(qnc3,23));
                    MDs(i,8,t,j,s) = std(temp(qnc3,13));%-temp(qnc3,23));
                    MDs(i,9,t,j,s) = nanmedian(temp(qnc3,13));%-temp(qnc3,23));
                                      
                    qrc = temp(:,6) == 0 & temp(:,7) == (j-1) & temp(:,8) == t; % Rightward direction
                    RD(i,t,j,s) = sum(qrc);
                    
                    qncc = temp(:,5) == 1 & temp(:,7) == (j-1) & temp(:,8) == t;% Collision rate
                    NCd(i,t,j,s) = sum(qncc);
                    
                    qcdir = temp(:,5) == 0 & temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 1;% Changed Direction
                    CDir(i,t,j,s) = sum(qcdir);
                    
                    qcdl = temp(:,5) == 0 & temp(:,7) == (j-1) & temp(:,8) == t;% Changed Direction
                    NC(i,t,j,s) = sum(temp(qcdl,5));% not true!!!
                    NR(i,t,j,s) = sum(temp(qcdl,4));             
                    NMT(i,t,j,s) = sum(temp(qcdl,10));
                end 
                
            end             
            %NC(i,s) = sum(temp(:,5));
            %NR(i,s) = sum(temp(:,4));             
            %NMT(i,s) = sum(temp(:,10));
        end
    end
    save('WholeComp3.mat','ima','mcurv','DpassO','DfromO','MaxVel','TMaxVel',...
         'MaxAcc','TMaxAcc','trajc','RTs','MDs','NC','NR','NMT','RD', 'NCd',...,
         'CDir','RTs2','MDs2','velc')
     % wholecomp3 contains the median values!
end
%%
if poolploting == 1
    %%
    load('WholeComp3.mat')
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
                    si = tempxb - uEb; % detect the sign(rightside of baseline = positive)
                    biases(i,j,d,s) = polyarea(xPb,yPb);
                    if sum(si) ~= 0
                        signs(i,j,d,s)  = (sum(si)/abs(sum(si)));
                    elseif isnan(sum(si))
                        signs(i,j,d,s) = nan;                        
                    end
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
    
    %% Plot Trajectories
    % for the averaged ones
    for jj = 1:5
        for ii = 1:3
            figure(jj); hold on
            for d = 1 : 2                
                tempx = trajcc(:,1,ii+0*3,1,jj,d) - 1*18.875;
                ys = trajcc(:,2,ii+0*3,1,jj,d);
                uE=tempx + trajcc(:,1,ii+0*3,2,jj,d);%./sqrt(18);
                lE=tempx - trajcc(:,1,ii+0*3,2,jj,d);%./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tempx,ys,color(ii,:),'Linewidth',2, 'LineStyle','-')
                scatter([obstaclepos(jj,1)-18.8750 0 0],[30-obstaclepos(jj,2) 9 29],'ro') 
            end
        end 
        figname = strcat(['MTraj-VF-obstacle# ',num2str(jj)]); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')        
    end
    close all
    for jj = 1:5
        for ii = 1:3
            figure(jj); hold on
            for d = 1 : 2                
                tempx = trajcc(:,1,ii+1*3,1,jj,d) - 1*18.875;
                ys = trajcc(:,2,ii+1*3,1,jj,d);
                uE=tempx + trajcc(:,1,ii+1*3,2,jj,d);%./sqrt(18);
                lE=tempx - trajcc(:,1,ii+1*3,2,jj,d);%./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                patch(xP,yP,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(tempx,ys,color(ii,:),'Linewidth',2, 'LineStyle',':')
                scatter([obstaclepos(jj,1)-18.8750 0 0],[30-obstaclepos(jj,2) 9 29],'ro') 
            end
        end 
        figname = strcat(['MTraj-NVF-obstacle# ',num2str(jj)]); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')        
    end
    close all
    % for individual subjects
    for jj = 3
        for ii = 1:3
            for ss = 1 : 18
                figure(jj);subplot(6,3,ss); hold on
                for d = 1 : 2 
                    ecc = 1;
                    if EC(ii+1*3,jj,d,ss)== 0
                        ecc = nan;
                    end
                    ecc = repmat(ecc,sn+1,1);
                    tempx = ecc.*trajc(:,1,ii+1*3,1,jj,d,ss) - 1*18.875;
                    ys = ecc.*trajc(:,2,ii+1*3,1,jj,d,ss);
                    uE=tempx + ecc.*trajc(:,1,ii+1*3,2,jj,d,ss);%./sqrt(18);
                    lE=tempx - ecc.*trajc(:,1,ii+1*3,2,jj,d,ss);%./sqrt(18);
                    xP=[lE;flipud(uE)];
                    yP=[ys;flipud(ys)];
                    patch(xP,yP,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                    plot(tempx,ys,color(ii,:),'Linewidth',2, 'LineStyle',':')
                    scatter([obstaclepos(jj,1)-18.8750 0 0],[30-obstaclepos(jj,2) 9 29],'ro') 
                end
            end
        end 
%         figname = strcat(['MTraj-indviduals-NVF-obstacle# ',num2str(jj)]); 
%         savefig(gcf,figname) 
%         saveas(gcf,figname,'tiff')        
    end
    close all
    
    %% Plot variabilities
    % averaged for HRs
    mvar = zeros(2,6);
    for i = 1 : 6
        mvar(1,i) = nanmean(reshape(variability(i,:,:,:),1,[]));
        mvar(2,i) = nanstd(reshape(variability(i,:,:,:),1,[]));
    end
    HR = [ 0 -30 30];
    for i = [2 1 3]
        figure(1);hold on
        errorbar(HR(i),mvar(1,i+0*3),mvar(2,i+0*3),'Color',colors{i})
        errorbar(HR(i),mvar(1,i+1*3),mvar(2,i+1*3),'Color',colors{i})
    end
    plot([-30 0 30], mvar(1,[5-1*3 4-1*3 6-1*3]),'-')
    plot([-30 0 30], mvar(1,[5-0*3 4-0*3 6-0*3]),':')
    
    ylabel('Variability (cm^2)');  xlabel('Head Rotation');    
    set(gca,'XTick',-30:30:30)
    set(gca, 'XTickLabels',{'30 CCW','0','30 CW',})
    xlim([-35,35])
    
    figname = strcat('Variability-averaged'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all    
    
    % for individuals for HRs
    mvars = zeros(18,6);    
    for i = 1 : 6
        mvars(:,i) = nanmean(reshape(variability(i,:,:,:),[],18));        
    end
    HR = [ 0 -30 30];
    
    % with visual feedback
    for i = [2 1 3]
        figure(1);hold on
        scatter(repmat(HR(i),18,1),mvars(:,i+0*3))
    end
    for i = 1 : 18
        figure(1); hold on
        plot(HR([2, 1, 3]),mvars(i,[5-3,4-3,6-3]),'-')
    end
    ylabel('Variability (cm^2)');  xlabel('Head Rotation');    
    set(gca,'XTick',-30:30:30)
    set(gca, 'XTickLabels',{'30 CCW','0','30 CW',})
    xlim([-35,35])
    
    figname = strcat('Variability-ind-VF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all   
    
    % without visual feedback
    for i = [2 1 3]
        figure(1);hold on
        scatter(repmat(HR(i),18,1),mvars(:,i+1*3))
    end
    for i = 1 : 18
        figure(1); hold on
        plot(HR([2, 1, 3]),mvars(i,[5,4,6]),':')
    end
    ylabel('Variability (cm^2)');  xlabel('Head Rotation');    
    set(gca,'XTick',-30:30:30)
    set(gca, 'XTickLabels',{'30 CCW','0','30 CW',})
    xlim([-35,35])
    
    figname = strcat('Variability-ind-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all 
    
    %% Plot Biases
    msbiases = zeros(6,5,2,2); % cond,target,right/left, mean/std
    mbiases = zeros(6,5,2,2); % cond,target,right/left, mean/std
    for i = 1 : 6
        for t = 1 : 5
            for d = 1:2
                msbiases(i,t,d,1) = nanmean(signs(i,t,d,:).*biases(i,t,d,:));
                msbiases(i,t,d,2) = nanstd(signs(i,t,d,:).*biases(i,t,d,:));
                mbiases(i,t,d,1) = nanmean(abs(biases(i,t,d,:)));
                mbiases(i,t,d,2) = nanstd(abs(biases(i,t,d,:)));
            end
        end
    end
    
    % plot the biases for the majority directions - averaged    
    figure(1); hold on
    xa = 6:-1:1;            
    h = zeros(1,3);
    for i = [2 3]
        t1 = fliplr([msbiases(i+1*3,5:-1:3,1,1) msbiases(i+1*3,3:-1:1,2,1)]);% 1 => most rightward
        t2 = fliplr([msbiases(i+1*3,5:-1:3,1,1) msbiases(i+1*3,3:-1:1,2,1)]);
        h(1,i)=errorbar(1:6,t1,t2,'Color',colors{i},'LineStyle',':');
    end 
    legend([h(1,2) h(1,3) ],...
        {'HR = 30CCW','HR = 30CW'})
    legend('Location','best') ; 
    ylabel('Biases'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan nan])
    set(gca,'XTick',1:6)
    set(gca, 'XTickLabels',{'Most leftward','Leftward','Center-rightward',...
                            'Center-leftward','Rightward','Most rightward'})
    xlim([0,7])
    
    figname = strcat('Biases-averaged-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all    
    
    % plot the biases for the majority directions - individuals  
    for ss = 1 : 18
        figure(1); subplot(6,3,ss); hold on
        xa = 6:-1:1; 
        for i = [2 3]
            t1 = fliplr([signs(i+1*3,5:-1:3,1,ss).*biases(i+1*3,5:-1:3,1,ss)...
                signs(i+1*3,3:-1:1,2,ss).*biases(i+1*3,3:-1:1,2,ss)]);% 1 => most rightward
            plot(1:6,t1,'Color',colors{i},'LineStyle',':')
        end         
        set(gca,'XTick',1:6)        
        xlim([0,7])
    end
    
    figname = strcat('Biases-ind-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all 
    
    %% Plot percentage of rightward MVs
    % averaged
    figure(1); hold on    
    % with visual feedback
    h = zeros(1,3);    
    for i = [2 1 3]
        r = reshape(RD(i,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        h(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i},'LineStyle','-');
    end   
    % without visual feedback
    h1 = zeros(1,3);   
    figure(1); hold on
    for i = [2 1 3]
        r = reshape(RD(i+3,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i+3,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        h1(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i},'LIneStyle',':'); 

    end  
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs') 
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    xlim([0,6])
    
    figname = strcat('LR-MDir-averaged');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % individuals
    for i = [2 1 3]
        r = reshape(RD(i+0*3,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i+0*3,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l); 
        %h(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i},'LineStyle','-'); 
        for ss = 1 : 18
            figure(1); subplot(6,3,ss);hold on  
            %scatter(1:5, pr(:,ss),'MarkerEdgeColor',colors{i})
            plot(1:5, pr(:,ss),'Color',colors{i},'LineStyle',':');
        end
    end  
    
    figname = strcat('LR-MDir-indv');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    %% Plot Initial movement angles
%     % for each obstacle (averaged)
%     figure(1); subplot(2,1,1);hold on  
%     %xa = 5:-1:1;            
%     h1 = zeros(1,3);
%     for i = [2 1 3]
%         t1 = flipud(reshape(imac(i,1,:,1),5,1));
%         t2 = flipud(reshape(imac(i,2,:,1),5,1));
%         h1(1,i)=errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i});
%     end  
% 
%     figure(1); subplot(2,1,2);hold on
%     for i = [2 1 3]
%         t1 = flipud(reshape(imac(i,1,:,2),5,1));
%         t2 = flipud(reshape(imac(i,2,:,2),5,1));
%         errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i});
%     end
%     
%     figure(1); subplot(2,1,1);hold on
%     xa = 5:-1:1;            
%     h = zeros(1,3);
%     for i = [2 1 3]
%         t1 = flipud(reshape(imac(i+3,1,:,1),5,1));
%         t2 = flipud(reshape(imac(i+3,2,:,1),5,1));
%         h(1,i)=errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i},'LineStyle',':');
%     end 
%     
%     
%     legend([h(1,1) h(1,2) h(1,3) h1(1,1) h1(1,2) h1(1,3)],...
%         {'HR = 30CCW-NF','HR = 0-NF','HR = 30CW-NF',...
%          'HR = 30CCW-F','HR = 0-F','HR = 30CW-F'})
%     legend('Location','southeast')    
%     ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
%     plot(xa,[nan nan nan nan nan])
%     set(gca,'XTick',fliplr(xa))
%     set(gca, 'XTickLabels',fliplr({'','','','',''}))
%     xlim([0,6])
%     title('Initial movement angle modulation by varying head roll - Rightward MVs')
%     
%     figure(1); subplot(2,1,2);hold on
%     for i = [2 1 3]
%         t1 = flipud(reshape(imac(i+3,1,:,2),5,1));
%         t2 = flipud(reshape(imac(i+3,2,:,2),5,1));
%         errorbar(1:5,abs(t1-90),t2./sqrt(18),'Color',colors{i},'LineStyle',':');
%     end
%     xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
%     plot(xa,[nan nan nan nan nan])
%     xlim([0,6])
%     set(gca,'XTick',fliplr(xa))
%     set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
%     title('Initial movement angle modulation by varying head roll - Leftward MVs')
% 
%     figname = strcat('IMA-OBs-Whole'); 
%     savefig(gcf,figname) 
%     saveas(gcf,figname,'tiff')
%     close all

    %% Plot curvature
    % Based on the majority (averaged)
    figure(22); hold on  
    %xa = 5:-1:1;            
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = [flipud(reshape(mcurvc(i+1*3,1,3:5,1),3,1)); flipud(reshape(mcurvc(i+1*3,1,1:3,2),3,1))];
        t2 = [flipud(reshape(mcurvc(i+1*3,2,3:5,1),3,1)); flipud(reshape(mcurvc(i+1*3,2,1:3,2),3,1))];
        h(1,i)=errorbar(1:6,abs(t1),t2,'Color',colors{i});
    end     
    legend([h(1,1) h(1,2) h(1,3) ],...
        {'HR = 30CCW','HR = 0','HR = 30CW'})
    legend('Location','best') ; 
    ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',1:6)
    set(gca, 'XTickLabels',fliplr({'Most leftward','Leftward','Center-rightward',...
                                   'Center-leftward','Rightward','Most rightward'}))
    xlim([0,7])
    title('Maximum curvature modulation by varying head roll')
    
    figname = strcat('MC-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all
    
    % Based on the majority (individuals)    
    for ss = 1 : 18
        figure(22);subplot(6,3,ss); hold on 
        for i = [2 1 3]
            t1 = [flipud(reshape(mcurv(i+1*3,1,3:5,1,ss),3,1)); flipud(reshape(mcurv(i+1*3,1,1:3,2,ss),3,1))];
            t2 = [flipud(reshape(mcurv(i+1*3,2,3:5,1,ss),3,1)); flipud(reshape(mcurv(i+1*3,2,1:3,2,ss),3,1))];
            errorbar(1:6,abs(t1),t2,'Color',colors{i});
        end   
        plot(xa,[nan nan nan nan nan nan])
        set(gca,'XTick',1:6)
        set(gca, 'XTickLabels',fliplr({'ML','L','CR','CL','R','MR'}))
        xlim([0,7])
    end
    
    figname = strcat('MC-ind-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

    %% Plot distance passing the obstacle
    % Based on majority (averaged)
    figure(32); hold on           
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = [flipud(reshape(DpassOc(i+1*3,1,3:5,1),3,1)); flipud(reshape(DpassOc(i+1*3,1,1:3,2),3,1))];
        t2 = [flipud(reshape(DpassOc(i+1*3,2,3:5,1),3,1)); flipud(reshape(DpassOc(i+1*3,2,1:3,2),3,1))];
        h(1,i)=errorbar(1:6,abs(t1),t2,'Color',colors{i});
    end  
    legend([h(1,1) h(1,2) h(1,3) ],...
        {'HR = 30CCW','HR = 0','HR = 30CW'})
    legend('Location','best') ; 
    ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',1:6)
    set(gca, 'XTickLabels',fliplr({'Most leftward','Leftward','Center-rightward',...
                                   'Center-leftward','Rightward','Most rightward'}))
    xlim([0,7])
    title('Distance passing Obstacle modulation by varying head roll')
    
    figname = 'DPO-NVF'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all
    
    % Based on the majority (individuals)    
    for ss = 1 : 18
        figure(22);subplot(6,3,ss); hold on 
        for i = [2 1 3]
            t1 = [flipud(reshape(DpassO(i+1*3,1,3:5,1,ss),3,1)); flipud(reshape(DpassO(i+1*3,1,1:3,2,ss),3,1))];
            t2 = [flipud(reshape(DpassO(i+1*3,2,3:5,1,ss),3,1)); flipud(reshape(DpassO(i+1*3,2,1:3,2,ss),3,1))];
            errorbar(1:6,abs(t1),t2,'Color',colors{i});
        end   
        plot(xa,[nan nan nan nan nan nan])
        set(gca,'XTick',1:6)
        set(gca, 'XTickLabels',fliplr({'ML','L','CR','CL','R','MR'}))
        xlim([0,7])
    end
    
    figname = strcat('DPO-ind-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all    

%% Plot distance from the obstacle
    figure(1); hold on           
    h = zeros(1,3);
    for i = [2 1 3]
        t1 = [flipud(reshape(DfromOc(i+1*3,1,3:5,1),3,1)); flipud(reshape(DfromOc(i+1*3,1,1:3,2),3,1))];
        t2 = [flipud(reshape(DfromOc(i+1*3,2,3:5,1),3,1)); flipud(reshape(DfromOc(i+1*3,2,1:3,2),3,1))];
        h(1,i)=errorbar(1:6,abs(t1),t2,'Color',colors{i});
    end  
    legend([h(1,1) h(1,2) h(1,3) ],...
        {'HR = 30CCW','HR = 0','HR = 30CW'})
    legend('Location','best') ; 
    ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',1:6)
    set(gca, 'XTickLabels',fliplr({'Most leftward','Leftward','Center-rightward',...
                                   'Center-leftward','Rightward','Most rightward'}))
    xlim([0,7])
    title('Distance from Obstacle modulation by varying head roll')
    
    figname = 'DfO-NVF'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all
    
    % Based on the majority (individuals)    
    for ss = 1 : 18
        figure(1);subplot(6,3,ss); hold on 
        for i = [2 1 3]
            t1 = [flipud(reshape(DfromO(i+1*3,1,3:5,1,ss),3,1)); flipud(reshape(DfromO(i+1*3,1,1:3,2,ss),3,1))];
            t2 = [flipud(reshape(DfromO(i+1*3,2,3:5,1,ss),3,1)); flipud(reshape(DfromO(i+1*3,2,1:3,2,ss),3,1))];
            errorbar(1:6,abs(t1),t2,'Color',colors{i});
        end   
        plot(xa,[nan nan nan nan nan nan])
        set(gca,'XTick',1:6)
        set(gca, 'XTickLabels',fliplr({'ML','L','CR','CL','R','MR'}))
        xlim([0,7])
    end
    
    figname = strcat('DfO-ind-NVF'); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all

    %% Plot Number of collisition
    % the average for different HRs
    figure(1); hold on
    HR = [0 -30 30];
    tt = zeros(3,1);
    for i = [2 1 3]
        errorbar(HR(i), mean(reshape(NCc(i,1,:),[],1)),mean(reshape(NCc(i,2,:),[],1)),'Color',colors{i},'Marker','x')        
        tt(i) = mean(reshape(NCc(i,1,:),[],1));
    end
    plot(HR([2 1 3]),tt([2 1 3]),'-b')
    for i = [2 1 3]
        errorbar(HR(i), mean(reshape(NCc(i+3,1,:),[],1)),mean(reshape(NCc(i+3,2,:),[],1)),'Color',colors{i},'Marker','x')        
        tt(i) = mean(reshape(NCc(i+3,1,:),[],1));
    end
    plot(HR([2 1 3]),tt([2 1 3]),':b')
    xlabel('Head Rotations (deg)'); ylabel('Number of collisions')    

    figname = 'NCs-averaged-HRs'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % the average for different HRs    
    HR = [0 -30 30];
    for ss = 1 : 18
        figure(1); subplot(6,3,ss); hold on
        tt = zeros(3,1);
        for i = [2 1 3]
            errorbar(HR(i), NCc(i,1,ss),NCc(i,2,ss),'Color',colors{i},'Marker','x')        
            tt(i) = NCc(i,1,ss);
        end
        plot(HR([2 1 3]),tt([2 1 3]),'-b')
    end
   
    for ss = 1 : 18
        figure(1); subplot(6,3,ss); hold on
        tt = zeros(3,1);
        for i = [2 1 3]
            errorbar(HR(i), NCc(i+3,1,ss),NCc(i+3,2,ss),'Color',colors{i},'Marker','x')        
            tt(i) = NCc(i+3,1,ss);
        end
        plot(HR([2 1 3]),tt([2 1 3]),':b')
    end
    
    figname = 'NCs-ind-HRs'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % Collision-Rate (Number of collision/Number of MVs per each direction)
    % the average for different obstacles (Rightward vs. Leftward)
    
    % with visual feedback
    figure(1);subplot(2,1,1); hold on
    tt = zeros(5,2);
    t = 5:-1:1;    
    for j = 1 :3
        for i = 1:5
            m = nanmean(reshape(NCd(j,i,1,:)./RD(j,i,1,:),1,18));
            s = nanstd(reshape(NCd(j,i,1,:)./RD(j,i,1,:),1,18));
            errorbar(t(i), m,s,'-','Color',colors{j},'Marker','x')%rightward        
            tt(i,1) = m;
            m = nanmean(reshape(NCd(j,i,2,:)./RD(j,i,2,:),1,18));
            s = nanstd(reshape(NCd(j,i,2,:)./RD(j,i,2,:),1,18));
            errorbar(t(i), m,s,':','Color',colors{j},'Marker','x')% leftward        
            tt(i,2) = m;
        end
        plot(1:5,tt(5:-1:1,1),'-','Color',colors{j})
        plot(1:5,tt(5:-1:1,2),':','Color',colors{j})
    end
    set(gca,'XTick',1:5)
    xlim([0,6])
    set(gca, 'XTickLabels',fliplr({' ',' ',' ',' ',' '}))
    title('With visual feedback')
    % no visual feedback
    figure(1);subplot(2,1,2); hold on
    tt = zeros(5,2);
    t = 5:-1:1;    
    for j = 1 :3
        for i = 1:5
            m = nanmean(reshape(NCd(j+3,i,1,:)./RD(j+3,i,1,:),1,18));
            s = nanstd(reshape(NCd(j+3,i,1,:)./RD(j+3,i,1,:),1,18));
            errorbar(t(i), m,s,'-','Color',colors{j},'Marker','x')        
            tt(i,1) = m;
            m = nanmean(reshape(NCd(j+3,i,2,:)./RD(j+3,i,2,:),1,18));
            s = nanstd(reshape(NCd(j+3,i,2,:)./RD(j+3,i,2,:),1,18));
            errorbar(t(i), m,s,':','Color',colors{j},'Marker','x')        
            tt(i,2) = m;
        end
        plot(1:5,tt(5:-1:1,1),'-','Color',colors{j})
        plot(1:5,tt(5:-1:1,2),':','Color',colors{j})
    end   
    xlabel('Obstacle position'); ylabel('Number of collisions')  
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',1:5)
    set(gca, 'XTickLabels',fliplr({'MR','R','C','L','ML'}))
    xlim([0,6])
    title('Without visual feedback')

    figname = 'CR-averaged-Obs'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % Individual participants for different obstacles (Rightward vs. Leftward)    
    % with visual feedback
    for ss = 1 : 18
        figure(1);subplot(6,3,ss); hold on
        tt = zeros(5,2);          
        for j = 1 :3
            for i = 1:5
                m = NCd(j,i,1,ss)./RD(j,i,1,ss);                       
                tt(i,1) = m;
                m = NCd(j,i,2,ss)./RD(j,i,2,ss);                        
                tt(i,2) = m;
            end
            plot(1:5,tt(5:-1:1,1),'-','Color',colors{j})
            plot(1:5,tt(5:-1:1,2),':','Color',colors{j})
        end
    end
    
    figname = 'CR-ind-VF-Obs'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % no visual feedback
    for ss = 1 : 18
        figure(1);subplot(6,3,ss); hold on
        tt = zeros(5,2);          
        for j = 1 :3
            for i = 1:5
                m = NCd(j+3,i,1,ss)./RD(j+3,i,1,ss);                       
                tt(i,1) = m;
                m = NCd(j+3,i,2,ss)./RD(j+3,i,2,ss);                        
                tt(i,2) = m;
            end
            plot(1:5,tt(5:-1:1,1),'-','Color',colors{j})
            plot(1:5,tt(5:-1:1,2),':','Color',colors{j})
        end
    end
    
    figname = 'CR-ind-NVF-Obs';  
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    %% Plot the rigthward movement percentage 
    % averaged
    % Plot Number of right vs left -ward movements 
    figure(1); hold on  
    for i = 1:3
        for j = 1 : 5            
            r = reshape(RD(i+0*3,:,1,:),5,18);r = flipud(r);
            l = reshape(RD(i+0*3,:,2,:),5,18);l = flipud(l);
            pr = r.*100./(r+l); 
            h1(1,i)=plot(1:5,nanmean(pr,2),'Color',colors{i},'LIneStyle','-');
        end           
    end 
    for i = 1:3
        for j = 1 : 5            
            r = reshape(RD(i+3,:,1,:),5,18);r = flipud(r);
            l = reshape(RD(i+3,:,2,:),5,18);l = flipud(l);
            pr = r.*100./(r+l); 
            h1(1,i)=plot(1:5,nanmean(pr,2),'Color',colors{i},'LIneStyle',':');
        end           
    end        
%     legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
%     legend('Location','southeast') 
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs')    

    set(gca,'XTick',1:5)
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Modulation of movement direction by rolling the head')
    xlim([0,6])
    
    figname = 'Rper-averaged-mean';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    figure(1); hold on  
    for i = 1:3
        for j = 1 : 5            
            r = reshape(RD(i+0*3,:,1,:),5,18);r = flipud(r);
            l = reshape(RD(i+0*3,:,2,:),5,18);l = flipud(l);
            pr = r.*100./(r+l); 
            h1(1,i)=plot(1:5,nanmedian(pr,2),'Color',colors{i},'LIneStyle','-');
        end           
    end 
    for i = 1:3
        for j = 1 : 5            
            r = reshape(RD(i+3,:,1,:),5,18);r = flipud(r);
            l = reshape(RD(i+3,:,2,:),5,18);l = flipud(l);
            pr = r.*100./(r+l); 
            h1(1,i)=plot(1:5,nanmedian(pr,2),'Color',colors{i},'LIneStyle',':');
        end           
    end   
    
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs') 
    set(gca,'XTick',1:5)
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    title('Modulation of movement direction by rolling the head')
    xlim([0,6])
    
    figname = 'Rper-averaged-meadian';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
    
    % for individuals    
    for i = [2 1 3]
        r = reshape(RD(i+0*3,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i+0*3,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l);          
        for ss = 1 : 18
            figure(1); subplot(6,3,ss);hold on             
            plot(1:5, pr(:,ss),'Color',colors{i},'LineStyle','-');
        end
    end 
    for i = [2 1 3]
        r = reshape(RD(i+1*3,:,1,:),5,18);r = flipud(r);
        l = reshape(RD(i+1*3,:,2,:),5,18);l = flipud(l);
        pr = r.*100./(r+l);          
        for ss = 1 : 18
            figure(1); subplot(6,3,ss);hold on             
            plot(1:5, pr(:,ss),'Color',colors{i},'LineStyle',':');
        end
    end   
    figname = 'Rper-ind';  
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
        
    %% Plot Velocity profiles
    % for the averaged ones
    for jj = 1:5
        for ii = 1:3
            figure(jj); hold on
            for d = 1 : 2                
                tempx = velcc(:,1,ii+1*3,1,jj,d) ;
                ys = 0:2000;ys = ys';
                if d == 2
                    ys = -ys;
                end
                uE=tempx + velcc(:,1,ii+0*3,2,jj,d);%./sqrt(18);
                lE=tempx - velcc(:,1,ii+0*3,2,jj,d);%./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                xpp = xP(~isnan(xP));
                ypp = yP(~isnan(xP));
                patch(ypp,xpp,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(ys,tempx,color(ii,:),'Linewidth',2, 'LineStyle','-')                
            end
        end 
        figname = strcat(['velprof-VF-obstacle# ',num2str(jj)]); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')        
    end
    close all
    for jj = 1:5
        for ii = 1:3
            figure(jj); hold on
            for d = 1 : 2                
                tempx = velcc(:,1,ii+1*3,1,jj,d) ;
                ys = 0:2000;ys = ys';
                if d == 2
                    ys = -ys;
                end
                uE=tempx + velcc(:,1,ii+1*3,2,jj,d);%./sqrt(18);
                lE=tempx - velcc(:,1,ii+1*3,2,jj,d);%./sqrt(18);
                xP=[lE;flipud(uE)];
                yP=[ys;flipud(ys)];
                xpp = xP(~isnan(xP));
                ypp = yP(~isnan(xP));
                patch(ypp,xpp,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
                plot(ys,tempx,color(ii,:),'Linewidth',2, 'LineStyle','-')
            end
        end 
        figname = strcat(['velprof-NVF-obstacle# ',num2str(jj)]); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')        
    end
    close all
    % for individual subjects
%     for jj = 1:5
%         for ii = 1:3
%             for ss = 1 : 18
%                 figure(jj);subplot(6,3,ss); hold on
%                 for d = 1 : 2 
%                     ecc = 1;
%                     if EC(ii+1*3,jj,d,ss)== 0
%                         ecc = nan;
%                     end
%                     ecc = repmat(ecc,sn+1,1);
%                     tempx = ecc.*velc(:,1,ii+1*3,1,jj,d,ss);
%                     ys = ecc.*(0:2000);ys = ys';
%                     if d == 2
%                         ys = -ys;
%                     end
%                     uE=tempx + ecc.*velc(:,1,ii+1*3,2,jj,d,ss);%./sqrt(18);
%                     lE=tempx - ecc.*velc(:,1,ii+1*3,2,jj,d,ss);%./sqrt(18);
%                     xP=[lE;flipud(uE)];
%                     yP=[ys;flipud(ys)];
%                     xpp = xP(~isnan(xP));
%                     ypp = yP(~isnan(xP));
%                     patch(ypp,xpp,color(ii,:),'FaceAlpha',0.3,'EdgeColor','none');
%                     plot(ys,tempx,color(ii,:),'Linewidth',2, 'LineStyle','-')
%                 end
%             end
%         end 
%         figname = strcat(['velprofile-indviduals-NVF-obstacle# ',num2str(jj)]); 
%         savefig(gcf,figname) 
%         saveas(gcf,figname,'tiff')        
%     end
    close all
    
     %% Plot Maximum velocity
     % for individual subjects
    for ss = 1:18
        names = strcat('Comp2_s',num2str(ss),'.mat');        
        rep = (temp(:,13)-temp(:,12)>1);
        for jj = 1:5
            for ii = 1 : 3
                temp = params(:,:,ii);
                figure(jj);subplot(6,3,ss); hold on
                for d = 1  
                    ecc = 1;
                    if EC(ii+1*3,jj,d,ss)== 0
                        ecc = nan;
                    end                    
                    q = rep == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                        temp(:,7) == d & temp(:,8) == jj & temp(:,24) == 0;
                    ecc = repmat(ecc,sum(q),1);
                    tempx = ecc.*temp(q,12);
                    hist(tempx)
                end
            end
        end 
        figname = strcat(['velprofile-indviduals-NVF-obstacle# ',num2str(jj)]); 
        savefig(gcf,figname) 
        saveas(gcf,figname,'tiff')        
    end
    close all
    
%     figure(1); hold on%subplot(1,2,1); hold on
%     HR = [0 -30 30];
%     mmd = zeros(3,1); smd = zeros(3,1);
%     for i = [2 1 3]
%         mmd(i) = nanmean(1000*reshape([MaxVelc(i+3*1,1,5:-1:3,1) MaxVelc(i+3*1,1,3:-1:1,2)],[],1));
%         smd(i) = nanmean(1000*reshape([MaxVelc(i+3*1,2,5:-1:3,1) MaxVelc(i+3*1,2,3:-1:1,2)],[],1));
%         errorbar(HR(i), mmd(i),smd(i),'Color',colors{i},'Marker','x')        
%     end
%     
%     plot(HR([2 1 3]),mmd([2 1 3]),':b')
%     %plot(repmat([-30;0;30],1,18),NC([2 1 3],:),'--')
%     xlabel('Head Rotations (deg)'); ylabel('Peak velocity')
%     title('Modulation of Maximum Velocity with varying head roll')
% 
%     figname = 'Peakvel-averaged'; 
%     savefig(gcf,figname)
%     saveas(gcf,figname,'tiff')
%     close all
    
    %% Plot Reaction times
    
    figure(1); hold on%subplot(1,2,1); hold on
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
    
    % for individual subjects and for selected targets (#3, #1&5, #2&4)
    obs = [3 3 ; 5 1 ; 4 2];
    colnew = [ 'k' ; 'c' ; 'm'];
    for ss = 1:18
        figure(33); subplot(6,3,ss); hold on
        for jj = 1:3
            ys = zeros(2,3);
            for ii = 1 : 3
                ys(1,ii) = mean([RTs(ii+1*3,1+3*0,obs(jj,1),1,ss),RTs(ii+1*3,1+3*0,obs(jj,2),2,ss)]);
                ys(2,ii) = mean([RTs(ii+1*3,3,obs(jj,1),1,ss),RTs(ii+1*3,3,obs(jj,2),2,ss)]);
            end
            figure(33); subplot(6,3,ss);%errorbar(1:3,ys(1,[2,1,3]),ys(2,[2,1,3]),'color',colnew(jj))
            plot(1:3,ys(1,[2,1,3])*1000,'color',colnew(jj))
        end               
    end
    figname = strcat(['RT-indviduals-NVF-obstacle# ',num2str(jj)]); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    close all
    
    
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
    
    % for individual subjects and for selected targets (#3, #1&5, #2&4)
    obs = [3 3 ; 5 1 ; 4 2];
    colnew = [ 'k' ; 'c' ; 'm'];
    for ss = 1:18
        figure(44); subplot(6,3,ss); hold on
        for jj = 1:3
            ys = zeros(2,3);
            for ii = 1 : 3
                ys(1,ii) = mean([MDs(ii+1*3,1,obs(jj,1),1,ss),MDs(ii+1*3,1,obs(jj,2),2,ss)]);
                ys(2,ii) = mean([MDs(ii+1*3,3,obs(jj,1),1,ss),MDs(ii+1*3,3,obs(jj,2),2,ss)]);
            end
            figure(44); subplot(6,3,ss);%errorbar(1:3,ys(1,[2,1,3]),ys(2,[2,1,3]),'color',colnew(jj))
            plot(1:3,ys(1,[2,1,3])*1,'color',colnew(jj))
        end               
    end
    
    figname = 'MDs-OBs-Whole'; 
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
   
    
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
        RTs = zeros(6,2,5,2); % Reaction time
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
                    RTs(i,1,t,j) = mean(temp(qnc,12)-temp(qnc,11));
                    RTs(i,2,t,j) = std(temp(qnc,12)-temp(qnc,11));
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