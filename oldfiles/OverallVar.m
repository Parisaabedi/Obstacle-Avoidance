%% General Variability Analysis
% Written by PA, Jan. 2019
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

%% Pooling
if pooled == 1
    % individual participants
    % Initialization (condition, mean/std, target, right(0)/left(1),subjects)   
    trajc = zeros(sn+1,2,6,2,5,2,18); % Normalized trajectories 
    variabilityoverall = zeros(6,5,2,18);% cond,tar,dir,sub
    imao = zeros(6,2,5,2,18); % initilal movement angle
    subs = 1:18;     
    for s = 1:18        
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp3','_',sID,'.mat');
        load(filename)       
        for i = 1 : 6
            temp = params(1:Ns(i),:,i);
            ttraj = ntraj(:,1:Ns(i),:,:);
            rep = (temp(:,13)-0*temp(:,12)>1);
            for j = 1 : 2
                for t = 1 : 5
                    % Remove bad trials
                    qnc = temp(:,6) == 0 & temp(:,7) == (j-1) & temp(:,8) == t;        
                    
                    if sum(qnc) ~= 0
                        imao(i,1,t,j,s) = nanmean(temp(qnc,18));
                        imao(i,2,t,j,s) = nanstd(temp(qnc,18));
                        pointers = find(qnc ==1);
                        temptraj=zeros(sn+1,2,sum(qnc));
                        for tt = 1 : sum(qnc)
                            tstart = 0;
                            tend = temp(pointers(tt),13);
                            tr = (tend-tstart)/2000; % time points
                            ts = tstart:tr:tend;
                            k = find(ttraj(:,pointers(tt),2,i)>12,1);
                            tend = ts(k);
                            tr = (tend-tstart)/2000; % time points
                            ts = tstart:tr:tend;
                            bpx = straj{pointers(tt),1,i};
                            bpy = straj{pointers(tt),2,i};
                            temptraj(:,1,tt)=ppval(bpx,ts); % xtraj
                            temptraj(:,2,tt)=ppval(bpy,ts); % ytraj  
                        end
                        trajc(:,1,i,1,t,j,s) = nanmean(temptraj(:,1,:),3); % mean xtraj
                        trajc(:,1,i,2,t,j,s) = nanstd(temptraj(:,1,:),0,3); % std xtraj
                        trajc(:,2,i,1,t,j,s) = nanmean(temptraj(:,2,:),3); % mean ytraj
                        trajc(:,2,i,2,t,j,s) = nanstd(temptraj(:,2,:),0,3); % std ytraj
                    else
                        trajc(:,1,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean xtraj
                        trajc(:,1,i,2,t,j,s) = zeros(sn+1,1)*nan; % std xtraj
                        trajc(:,2,i,1,t,j,s) = zeros(sn+1,1)*nan; % mean ytraj
                        trajc(:,2,i,2,t,j,s) = zeros(sn+1,1)*nan; % std ytraj  
                    end
                end 
                
            end
        end
    end
    EC = zeros(6,5,2,18);
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
                    variabilityoverall(i,j,d,s) = polyarea(xP,yP);
                end
                rrd = 100*reshape(RD(i,j,:,s),[],1)./sum(reshape(RD(i,j,:,s),[],1));
                t = rrd > 10;
                EC(i,j,:,s) = t;
            end
        end
    end
    variabilityoc = zeros(6,5,2,18); % cond, mean/std, obst, right/left
    imaoc = zeros(6,5,2,18); % initilal movement angle    
    for i = 1 : 6
        for d = 1:2
            for t = 1 : 5
                ecc = reshape(EC(i,t,d,:),18,1);
                ecc(ecc == 0) = nan;                
                variabilityoc(i,t,d,:) = (ecc.*reshape(variabilityoverall(i,t,d,:),18,1));                
                imaoc(i,t,d,:) = (ecc.*reshape(imao(i,1,t,d,:),18,1));                
            end
        end
       
    end
    save('varoverall.mat','variabilityoverall','imao','variabilityoc','imaoc')
end