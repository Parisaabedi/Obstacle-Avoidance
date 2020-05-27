% test my new analysis
straj = cell(300,2,6); % spline trajectories=trials, x and y, conditions
ntraj = zeros(200,300,2,6); % normalized trajectories = samples,trialnums,xandy,conditions
params = zeros(300,21,6); % other parameters = trialnum, parameters,conditions
Ns = zeros(6,1);
% [filename,pathname] = uigetfile('*.mat','Please choose the DataSet you want to load...');
% cd(pathname);
% % trn = 1;
% FILES = filename;
% leng = length(filename);
% 
% if exist([filename(1:leng-4) '.mat'],'file') == 2         
%     if contains(filename,'marked')
%         load([filename(1:leng-4) '.mat']);
%     end
% end
cond = {'F0';'FL';'FR';'NF0';'NFL';'NFR'};
% Read changed direction data
sheet = strcat('s',num2str(9));
num = readtable('C:\Users\pabed\OneDrive\Desktop\PHD\ObstacleAvoidance\ChangedDecisions.xlsx',...
                'Sheet',sheet);
cddata = table2array(num);
num = readtable('C:\Users\pabed\OneDrive\Desktop\PHD\ObstacleAvoidance\Cutoffend.xlsx',...
                'Sheet',sheet);
coffdata = table2array(num);
for cnd = 1 : 6  
    sID = '5';
    filename = strcat(sID,cond{cnd,1},'_marked','.mat');
    load(filename)
    N = length(D.trajdir);
    temp = zeros(N,21); % 
    %[1->sid,2->trialnum, 3->condition(Fl,F0,FR,NFL,NF0,NFR),4-> repeat,
    % 5-> collision,6-> good/bad; 7->direction; 8->obstaclen; 9->changedDirection, 
    % 10-> missed target, 11->target onset, 12-> movementONset, 13-> movementOFfset
    % 14 -> maxVel, 15->timeofMaxVel, 16->Maxacc, 17->TimeofMax, 18->Init.Mv.Angle
    % 19-> Maxcurvature, 20->dispssinobs, 21->disfromobtraj

    % Initial values
    temp(:,1) = param(:,1); % subject ID
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
        tON = D.eventtime{i}(2); % time of target presentation
        temp(i,11) = tON;
        %ton = find(D.t{i}(:) >= tON, 1 ); 
        %tEM = D.handOFF{i}; % End of movement (targethit/collision)
        tEM = D.t{i}(end); % End of the trial
        %tem = find(D.t{i}(:) >= tEM, 1 ); 
        Rarr = find(D.t{i} >= tON, 1 ):find(D.t{i} <= tEM, 1, 'last' );
        sr = length(Rarr) / (tEM - tON);
        x = nan(length(Rarr)*2-1,2);
        x(1:2:end,:)= D.handscreen{1,i}(Rarr(1):Rarr(end),:);
        % intropolate to find the missing data points
        xint = inpaint_nans(x,2);
        %figure(1); hold on; plot(xint(:,1),xint(:,2))

        % Low-pass filtering (to remove noise)
        sfr = sr * 2;
        cutoff = 50;
        %xint = [repmat(xint(1,:),3,1) ; xint; repmat(xint(end,:),3,1)];
        tempx = GRAIautoregfilt(sfr,cutoff,xint(:,1));
        tempy = GRAIautoregfilt(sfr,cutoff,xint(:,2));
        handx = tempx;
        handy = tempy;
        % if you are interested in the end-point errors, you should use the
        % original hand points
        %figure(2); hold on; plot(handx,handy,'b',xint(:,1),xint(:,2),'r')

        % Velocity Calculation
        win = 0.005;
        handXv = GRAIdiff(1/sfr,win,handx);
        handYv = GRAIdiff(1/sfr,win,handy);
        clear ti;
        ti(1:2:length(handx)) = D.t{i}(Rarr);
        ti(2:2:end) = ti(1:2:length(handx)-1) + diff(ti(1:2:length(handx)))./2 ;
        %figure(3); hold on; plot(ti,handXv,'b',ti,handYv,'r')

        % Acceleration Calculation
        handXa = GRAIdiff(1/sfr,win,handXv);
        handYa = GRAIdiff(1/sfr,win,handYv);
        %figure(4); hold on; plot(ti,handXa,'b',ti,handYa,'r')

        % ON-OFF calculations
        vel = sqrt(handXv.^2 + handYv .^ 2);
        acc = sqrt(handXa.^2 + handYa .^ 2);
        %figure(5);plot(ti,vel,'b',ti,acc,'r')
        if isnan(max(vel))
            poson = 0;
        else
            indv = find(vel > 0.1 * max(vel));
            inda = find(acc > 0.1 * max(acc));
            posv = GRAIgroup(indv,2,5);
            %posa = GRAIgroup(inda,2,8);
            poson = 0;
            for j = 1 : length(posv)/2
                if sum(acc(indv(posv(j)):indv(posv(j))+4)) > 0.1 * max(acc)
                    poson = indv(posv(j));
                    break
                end
            end
        end
        if poson == 0 % for delayed movements
            tOn = D.handON{i};
        else
            tOn = D.t{i}(Rarr(1)+round(poson/2));
            poson = poson + 2*Rarr(1);
        end          
        
        if abs(tOn - D.handON{i}) > 0.25
            tOn = D.handON{i};
            poson = 2*find(D.t{i}(:)>tOn,1);
            %pause % check if the detected movement onset is too early/late
        end
        if abs(D.handOFF{i}- D.t{i}(end))<0.100
            tOff = D.t{i}(end);
            temp(i,10) = 0;
        else
            temp(i,10) = 1;
            indv = find(vel(poson+50:end) < 0.01 * max(vel));
            %posv = GRAIgroup(indv,2,5);
            %pos2 = indv(posv(3));
            if isempty(indv) || length(indv)==1
                tOff = D.handOFF{i};
            else
                post = GRAIgroup(indv,2,5);                
                if isempty(post)
                    tOff = D.handOFF{i};
                else
                    pos2 = indv(post(1));
                    tOff = D.t{i}(round((50+poson+pos2(1))/2));
                end
            end
        end
        if tOff - tOn < .25
            tOff = D.handOFF{i};
        end
        %figure(5);hold on; plot([tOn tOn],[-1,400],[tOff tOff],[-1, 400])

        temp(i,12) = tOn; % hand Onset
        temp(i,13) = tOff; % hand offset
        Rarr2 = find(ti>=tOn, 1 ):find(ti<=tOff, 1, 'last' );
        if isnan(max(vel))
        else
            temp(i,14)= max(vel(Rarr2)); % max velocity
            pos3 = find(vel(Rarr2) == max(vel(Rarr2)));
            temp(i,15) = D.t{i}(Rarr(1)*0+round((poson+pos3)/2)); % time of peak velocity
            temp(i,16) = max(acc(Rarr2)); % max acceleration
            pos4 = find(acc(Rarr2) == max(acc(Rarr2)));
            temp(i,17) = D.t{i}(Rarr(1)*0+round((poson+pos4)/2)); % time of peak acc.

            % Initial MV angle
            Rarr3 = find(ti>=tOn, 1 ):find(ti<=tOn+0.1, 1, 'last' );
            b = robustfit(xint(Rarr3,1),xint(Rarr3,2)); % fit on trajectory
            temp(i,18) = CalAngle(90,b(2));% initial movement direction (deg)

            % Maximum curvature
            temp(i,19) = max(abs(handx(Rarr2,1) - D.xOrig));

            % Distance passing the obstacle 
            if temp(i,5) == 0
                posdof = find(xint(Rarr2,2) >= 30-obpos(D.obstacleorder(i),2));
                if isempty(posdof)
                    pause
                end
                temp(i,20) = xint(Rarr2(1)+posdof(1),1) - obpos(D.obstacleorder(i),1);
            end

            % Min distance from the obstacle along the trajectory
            disfob = sqrt((xint(Rarr2,1)-obpos(D.obstacleorder(i),1)).^2 +...
                          (xint(Rarr2,2)-(30-obpos(D.obstacleorder(i),2))).^2);
            temp(i,21) = min(disfob);

            % fit S-pline on the trajectories    
            t = ti(Rarr2)-tOn;
            Hxy = [handx(Rarr2) handy(Rarr2)];
            q = isnan(Hxy(:,1));
            Hxy(q,:) = [];
            t(q) = [];
            bpx = splinefit(t,Hxy(:,1),10,6);
            bpy = splinefit(t,Hxy(:,2),10,6);
            straj{i,1,cnd} = bpx;
            straj{i,2,cnd} = bpy;
            %tt = t(1):0.0001:t(end);
            %bpvalues = ppval(bpx,tt);
            %figure(6); hold on; plot(t,Hxy(:,1),'b',tt,bpvalues,':r')

            % sample from the trajectory
            tr = (t(end)-t(1))/2000; % time points
            ts = t(1):tr:t(end);
            ys = ppval(bpy,ts); % oversample in time (y-axis)
            xs = ppval(bpx,ts); % oversample in time (x-axis)

            % Normalize along the path distance
            sn = 200;    
            % find equally distanced points along the path
            path = sqrt((diff(ys)).^2+(diff(xs)).^2);
            dp = sum(path)/sn;
            trajn = zeros(sn, 2); % normalized trajectory
            trajn(1,2) = ys(1); trajn(end,2) = ys(end);
            trajn(1,1) = xs(1); trajn(end,1) = xs(end);
            tn = zeros(sn,1); tn(1) = ts(1);tn(sn) = tn(end); 
            pos = 2;
            po1 = 1;
            while pos < 200
                pt = 0;
                while pt < dp && po1 < length(ts)-1
                    pt = pt + path(po1);
                    po1 = po1 + 1;
                end
                trajn(pos,1) = xs(po1+1);
                trajn(pos,2) = ys(po1+1);
                tn(pos) = ts(po1+1);
                pos = pos + 1;        
            end
            figure(7);hold on; plot(trajn(:,1),trajn(:,2))
            ntraj(:,i,1,cnd) = trajn(:,1); ntraj(:,i,2,cnd) = trajn(:,2); 
        end
    end
    params(1:N,:,cnd) = temp;
    Ns(cnd) = N;
    close all
end

save('Comp_s9.mat','straj','ntraj','params','Ns')