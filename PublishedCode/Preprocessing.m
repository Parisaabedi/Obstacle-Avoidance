%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by PAK, October 2018
% This codes go over all the subjects data to calculate the interested
% parameteres for analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for subnum = 1 : 18
    dir1 = 'C:\Users\Paris\Desktop\PhD\Projects\Project #4\Mostupdated\ObstacleAvoidance\WholeData';
    dir2 = strcat(dir1,'\',num2str(subnum));
    cd(dir2)
    %%
    sn = 2000;
    straj = cell(300,2,6); % spline trajectories=trials, x and y, conditions
    ntraj = zeros(sn+1,300,2,6); % normalized trajectories = samples,trialnums,xandy,conditions
    nvelacc = zeros(sn+1,300,2,6); % normalized vel&acc = samples,trialnums,vel/acc,conditions
    params = zeros(300,24,6); % other parameters = trialnum, parameters,conditions
    Ns = zeros(6,1);
    
    cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'}; % different experimental conditions
    % Read changed direction data
    sheet = strcat('s',num2str(subnum));
    num = readtable('C:\Users\Paris\Desktop\PhD\Projects\Project #4\Mostupdated\ObstacleAvoidance\ChangedDecisions.xlsx',...
                    'Sheet',sheet);
    cddata = table2array(num);
    num = readtable('C:\Users\Paris\Desktop\PhD\Projects\Project #4\Mostupdated\ObstacleAvoidance\Cutoffend.xlsx',...
                    'Sheet',sheet);
    coffdata = table2array(num);
    for cnd = 1 : 6    
        sID = num2str(subnum);       
        filename = strcat(sID,cond{cnd,1},'_marked','.mat');
        load(filename)
        N = length(D.trajdir);
        temp = zeros(N,24); % 
        
        %[1->sid,2->trialnum, 3->condition(Fl,F0,FR,NFL,NF0,NFR),4-> repeat,
        % 5-> collision,6-> good/bad; 7->direction; 8->obstaclen; 9->changedDirection, 
        % 10-> missed target, 11->target onset, 12-> movementONset(extrapolation),
        % 13-> movementOFfset, 14 -> maxVel, 15-> timeofMaxVel, 16->Maxacc,
        % 17->TimeofMax, 18->Init.Mv.Angle, 19-> Maxcurvature, 
        % 20->dispssinobs, 21->disfromobtraj, 
        % 22 -> mvONset(velthreshold), 23 -> mvONset(accthreshold), 
        % 24 ->measured changed direction and fast movements

        % Initial values
        temp(:,1) = SID; % subject ID
        temp(:,2) = 1:N; % trial num
        temp(:,3) = cnd; % condition
        temp(:,4) = D.repeat; % trial repeat
        temp(:,5) = D.obstaclehit; % collision => 1
        temp(:,6) = cell2mat(D.good); % good / bad
        temp(:,7) = D.trajdir; % trajectory direction 
        temp(:,8) = D.obstacleorder; % obstacle order
        % add changed direction and missed target trials 
        cds = cddata(:,cnd);
        cds(isnan(cds)) = [];
        if ~isempty(cds)
            temp(cds,9) = 1; % changed direction
        end
        cof = coffdata(:,cnd);
        cof(isnan(cof)) = [];
        if ~isempty(cof)
            temp(cof,10) = 1; % missed target
        end

        obpos = flipud(D.obstacle_posi); % to counter-blance screen mirroring

        for i = 1 : N  
            if temp(i,6) == 0
                tON = D.eventtime{i}(2); % time of target presentation
                temp(i,11) = tON;            
                tEM = D.handOFF{i}; % End of the trial
                Rarr = find(D.t{i} >= tON, 1 ):find(D.t{i} <= tEM, 1, 'last' );
                sr = length(Rarr) / (tEM - tON);
                x = nan(length(Rarr)*2-1,2);
                x(1:2:end,:)= D.handscreen{1,i}(Rarr(1):Rarr(end),:);

                % intropolate to find the missing data points
                xint = inpaint_nans(x,2);
                %figure(1); hold on; plot(xint(:,1),xint(:,2))

                handx = xint(:,1);
                handy = xint(:,2);
                % adapt the time points based on the increased sampling frequency
                clear ti;
                ti(1:2:length(handx)) = D.t{i}(Rarr);
                ti(2:2:end) = ti(1:2:length(handx)-1) + diff(ti(1:2:length(handx)))./2 ;

                %% Functional analysis of the Trajectory
                % fit S-pline on the trajectories
                Rarr2 = find(ti >= tON, 1 ):find(ti <= tEM, 1, 'last' );
                t = ti(Rarr2)-tON;
                Hxy = [handx(Rarr2) handy(Rarr2)];
                q = isnan(Hxy(:,1));
                Hxy(q,:) = [];
                t(q) = [];
                bpx = splinefit(t,Hxy(:,1),10,6);
                bpy = splinefit(t,Hxy(:,2),10,6);
                straj{i,1,cnd} = bpx;
                straj{i,2,cnd} = bpy;

                % sample from the trajectory (Normalization along y-axis)            
                temp(i,13) = t(end);

                tr = (t(end)-t(1))/2000; % time points
                ts = t(1):tr:t(end);
                ys = ppval(bpy,ts); % oversample in time (y-axis)
                xs = ppval(bpx,ts); % oversample in time (x-axis)

                %% Calculate the Velocity and Acceleration profiles
                % low pass filtering
                sfr = length(ys)/(ts(end)-ts(1));
                cutoff = 100;            
                tempx = GRAIautoregfilt(sfr,cutoff,xs);
                tempy = GRAIautoregfilt(sfr,cutoff,ys);
                figure(7);hold on; plot(tempx,tempy)
                ntraj(:,i,1,cnd) = xs; ntraj(:,i,2,cnd) = ys;

                % Velocity Calculation
                win = 0.005;
                handXv = GRAIdiff(1/sfr,win,tempx);
                handYv = GRAIdiff(1/sfr,win,tempy);
                %figure(3); hold on; plot(ts,handXv,'b',ts,handYv,'r')

                % Acceleration Calculation
                handXa = GRAIdiff(1/sfr,win,handXv);
                handYa = GRAIdiff(1/sfr,win,handYv);

                % check for the changed directions
                Rarrcd = find(ts>ts(1)+0.2,1 ):find(ts <= ts(end)-0.2, 1, 'last' );
                if min(handYa(Rarrcd)) < -70
                    temp(i,24) = 1;
                end
                %figure(4); hold on; plot(ts,handXa,'b',ts,handYa,'r')

                %% ONset and OFFset calculations
                vel = sqrt(handXv.^2 + handYv .^ 2);
                nvelacc(:,i,1,cnd) = vel;
                acc = sqrt(handXa.^2 + handYa .^ 2);
                nvelacc(:,i,2,cnd) = acc;

                % Threshold method
                indv = find(vel > 0.1 * max(vel));
                inda = find(acc > 0.1 * max(acc));

                posv = GRAIgroup(indv,2,5);
                posonv = indv(posv(1));
                tonv = ts(posonv); % time of onset based on velocity threshold
                temp(i,22) = tonv;

                posa = GRAIgroup(inda,2,5);
                posona = inda(posa(1));
                tona = ts(posona); % time of onset based on acceleration threshold
                temp(i,23) = tona;

                % Extrapolation method (Based Smeets and Brenner, 2018)            
                Rarrv = find(ts >= 0, 1 ):find(ts <= 0+0.15, 1, 'last' );
                x1 = [0 1 5]; y1 = repmat(nanmean(vel(Rarrv)),3,1); % the baseline calculation 
                y2 = [0.25*max(vel) 0.75* max(vel)];
                point1 = find(vel >= 0.25*max(vel));point2 = find(vel >= 0.75*max(vel));
                x2 = [ts(point1(1)) ts(point2(1))];
                % draw the line between two points (y = m*x+n)
                m = (y2(2)-y2(1))/(x2(2)-x2(1)); n = y2(1)-m*x2(1);
                x2p = [0 1 5];y2p=m.*x2p + n;
                [tx,ty] = polyxpoly(x1,y1,x2p,y2p); % The intersection of the two lines
                tonex = tx;
                if isempty(tx)
                    tonex = -1;
                    %pause
                end
                temp(i,12) = tonex;
                % check visually if the calculations are right
                %figure(10); hold on; plot(ts, vel); plot(x1,y1,x2p,y2p) 
                %figure(11); hold on; plot(ts, acc);
                
                % check again for the fast responses
                if temp(i,12) < 0.1
                     temp(i,24) = 1;
                end

                %% calculate the rest of the parameters
                if isnan(max(vel))
                else
                    temp(i,14)= max(vel); % max velocity
                    pos3 = find(vel == temp(i,14));
                    temp(i,15) = ts(pos3); % time of peak velocity
                    temp(i,16) = max(acc); % max acceleration
                    pos4 = find(acc == temp(i,16));
                    temp(i,17) = ts(pos4); % time of peak acc.

                    % Initial MV angle
                    Rarr3 = find(ts>=0, 1 ):find(ts<=0.1, 1, 'last' );
                    b = robustfit(tempx(Rarr3),tempy(Rarr3)); % fit on trajectory
                    temp(i,18) = CalAngle(90,b(2));% initial movement direction (deg)

                    % Maximum curvature
                    temp(i,19) = max(abs(handx(Rarr2,1) - D.xOrig));

                    % Distance passing the obstacle 
                    if temp(i,5) == 0
                        posdof = find(ys >= 30-obpos(D.obstacleorder(i),2));
                        if isempty(posdof)
                            pause
                        end
                        temp(i,20) = abs(xs(posdof(1)) - obpos(D.obstacleorder(i),1));
                    end

                    % Min distance from the obstacle along the trajectory
                    disfob = sqrt((xs-obpos(D.obstacleorder(i),1)).^2 +...
                                  (ys-(30-obpos(D.obstacleorder(i),2))).^2);
                    temp(i,21) = min(disfob);                
                end
                close all
            end
        end
        params(1:N,:,cnd) = temp;
        Ns(cnd) = N;
        %close all
    end
    
    finalname = strcat('Comp_s', num2str(subnum),'.mat');
    save(finalname,'straj','ntraj','nvelacc','params','Ns')
end