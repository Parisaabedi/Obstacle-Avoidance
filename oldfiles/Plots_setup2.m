%% Obstacle Avoidacne Data analysis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond = { 'Left' ; '0' ; 'Righ'};
colors = { 'g' ; 'k' ; 'r'};
init = 0;
pooled = 0;
poolploting = 0;
plotting = 0;
trajectories = 1;
if init == 1
    % individual participants
    % Initialization (col1 => mean, col2 => standard deviation)
    ima = zeros(3,2); % initilal movement angle    
    mcurv = zeros(3,2); % movement curvature
    DfromO = zeros(3,2); % Distance from Obstacle
    NC = zeros(3,1); % number of collisions 
    % columns : 1-collision, 2-repeat, 3-initial movment angle, 4-movement curvature, 
    % 5-distance from obstacle, 6-movement duration
    Comp = cell(3,1); 
    
    %
    files = dir;

    for j = 1 : length(files)
        k = strfind(files(j).name,'_marked');
        for i = 1 : 3
            % load the data        
            if ~isempty(k)
                k1 = strfind(files(j).name,cond{i});
                if ~isempty(k1)   
                    load(files(j).name)
                    tk = size(param,1) - size(D.repeat,1);
                    if tk > 0
                        D.repeat = [D.repeat ; zeros(tk,1)];     
                    end
                    if strfind(filename,'_marked_marked')
                        m = matfile(filename(1:end-7), 'Writable', true);

                    else
                        m = matfile(filename, 'Writable', true);
                    end   
                    if size(param,2) < 26
                        param = [param D.repeat];
                        m.param = param;                                
                    end
                    q =  param(:,24) == 0; % Select valid trials   
                    Comp{i,1} = zeros(sum(q),6);
                    Comp{i,1}(:,1) = param(q,25); % collisions
                    Comp{i,1}(:,2) = param(q,26); % repeated trail (too long)
                    Comp{i,1}(:,3) = param(q,8);  % Initial movement angle
                    Comp{i,1}(:,4) = param(q,22); % Maximum curvature
                    Comp{i,1}(:,5) = param(q,23); % Distance from obstacle
                    Comp{i,1}(:,6) = param(q,17); % Movement duration

                end
            end

        end
    end
    save('Comp.mat','Comp')
elseif init == 2
    load('Comp.mat')
end

if pooled == 1
    files = dir;
    % columns : 1-collision, 2-repeat, 3-initial movment angle, 4-movement curvature, 
    % 5-distance from obstacle, 6-movement duration, 7-subject ID
    Comps = cell(3,1); 
    mtrajs = zeros(80,5*2,3);
    strajs = zeros(80,5*2,3);
    for j = 1 : length(files)
        k = strfind(files(j).name,'Comp');
        if ~isempty(k)
            load(files(j).name)
            sid = str2double(files(j).name(end-4));
            for i = 1 : 3
                Comp{i,1} = [Comp{i,1} sid*ones(length(Comp{i,1}),1)];
                Comps{i,1} = [Comps{i,1} ;Comp{i,1}];
            end
            
        end
        l = strfind(files(j).name,'traj');
        if ~isempty(l)
            load(files(j).name)
            sid = str2double(files(j).name(end-4));
            for i = 1 : 3                
                mtrajs(:,1+(sid-1)*2:2+(sid-1)*2,i) = mtraj(:,:,i);
                strajs(:,1+(sid-1)*2:2+(sid-1)*2,i) = straj(:,:,i);
            end
            
        end
    end
    save('Comps.mat','Comps')
    save('trajs.mat','mtrajs','strajs')
elseif pooled == 2
    load('Comps.mat')
    load('trajs.mat')
end

if poolploting == 1
        % Plot Initial movement angles
    for i = 1 : 3
        qnc = Comps{i}(:,1) == 0 & Comps{i}(:,2) == 0; % no collision, no repeat
        qc = Comps{i}(:,1) == 1 ;  % collision, no repeat
        ima(i,1) = mean(Comps{i}(:,3));
        ima(i,2) = std(Comps{i}(:,3));
        mcurv(i,1) = mean(Comps{i}(:,4));
        mcurv(i,2) = std(Comps{i}(:,4));
        DfromO(i,1) = mean(Comps{i}(:,5));
        DfromO(i,2) = std(Comps{i}(:,5));
        NC(i) = sum(qc);
    end
    
    figure(1); subplot(1,2,1); hold on
    for i = 1 : 3
        bar(HR(i), abs(ima(i,1)-90),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('initial movement angle')    
    figure(1); subplot(1,2,2);  hold on
    for i = 1: 3
        bar(HR(i), ima(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (deg)')
    supertitle('Initial movement angle modulation by varying head roll')
    figname = 'pooled-IMA';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot curvature
    figure(2); subplot(1,2,1); hold on
    for i = 1:3
        bar(HR(i), mcurv(i,1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('movement curvature')
    supertitle('movement curvature modulation by varying head roll')
    figure(2); subplot(1,2,2); hold on
    for i = 1: 3
        bar(HR(i), mcurv(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (cm)')
    supertitle('movement curvature modulation by varying head roll')
    figname = 'pooled-MC';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot distance from obstacle
    figure(3); subplot(1,2,1); hold on
    for i = 1: 3
        bar(HR(i), abs(DfromO(i,1)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('distance from obstacle')    
    figure(3); subplot(1,2,2); hold on
    for i = 1:3
        bar(HR(i), DfromO(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (cm)')
    supertitle('distance from obstacle modulation by varying head roll')
    figname = 'pooled-dfo';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot Number of collisition
    figure(4); hold on
    for i = 1:3
        bar(HR(i), NC(i),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('obstacle collision')
    title('obstacle collision modulation by varying head roll')
    figname = 'pooled-NC';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    % Plot the averaged trajectories
    figure(5); hold on
    color = [ '-g'; '-k'; '-r'];
    mtraj = zeros(80,2,3);
    straj = zeros(80,2,3);
    for i = 1 : 3
        mtraj(:,1,i) = mean(abs(18.875-mtrajs(:,1:2:10,i)),2);
        mtraj(:,2,i) = mean(mtrajs(:,2:2:10,i),2);   
        straj(:,1,i) = mean(strajs(:,1:2:10,i),2);
        straj(:,2,i) = mean(strajs(:,2:2:10,i),2); 
        figure(5); hold on
        shadedErrorBar(mtraj(:,2,i),mtraj(:,1,i),straj(:,1,i),'lineprops',color(i,:),'patchSaturation',0.33)        
    end 
    set(gca,'View',[90 -90])
    scatter([9 29 19 19 19],[18.875-18.875 18.875-18.875 18.875-18.875 25.975-18.875 11.775-18.875],'ro')
    figname = 'pooled-traj';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
end


if plotting == 1
    % Plot Initial movement angles
    for i = 1 : 3
        qnc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 0; % no collision, no repeat
        qc = Comp{i}(:,1) == 1 ;  % collision, no repeat
        ima(i,1) = mean(Comp{i}(:,3));
        ima(i,2) = std(Comp{i}(:,3))/sqrt(sum(qnc));
        mcurv(i,1) = mean(Comp{i}(:,4));
        mcurv(i,2) = std(Comp{i}(:,4))/sqrt(sum(qnc));
        DfromO(i,1) = mean(Comp{i}(:,5));
        DfromO(i,2) = std(Comp{i}(:,5))/sqrt(sum(qnc));
        NC(i) = sum(qc);
    end
    
    figure(1); subplot(1,2,1); hold on
    for i = 1 : 3
        bar(HR(i), abs(ima(i,1)-90),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('initial movement angle')    
    figure(1); subplot(1,2,2);  hold on
    for i = 1: 3
        bar(HR(i), ima(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (deg)')
    supertitle('Initial movement angle modulation by varying head roll')
    figname = 's2-IMA';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot curvature
    figure(2); subplot(1,2,1); hold on
    for i = 1:3
        bar(HR(i), mcurv(i,1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('movement curvature')
    figure(2); subplot(1,2,2); hold on
    for i = 1: 3
        bar(HR(i), mcurv(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (cm)')
    supertitle('movement curvature modulation by varying head roll')
    figname = 's2-MC';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot distance from obstacle
    figure(3); subplot(1,2,1); hold on
    for i = 1: 3
        bar(HR(i), abs(DfromO(i,1)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('distance from obstacle')
    figure(3); subplot(1,2,2); hold on
    for i = 1:3
        bar(HR(i), DfromO(i,2),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('Standard deviation (cm)')
    supertitle('distance from obstacle modulation by varying head roll')
    figname = 's2-dfo';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    % Plot Number of collisition
    figure(4); hold on
    for i = 1:3
        bar(HR(i), NC(i,1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
    end
    xlabel('Head Rotations'); ylabel('obstacle collision')
    title('obstacle collision modulation by varying head roll')
    figname = 's2-NC';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
end

if trajectories == 1
    files = dir;
    trajx = 1000*ones(1000,300,3);
    trajy = 1000*ones(1000,300,3);
    l = zeros(300,3);
    for j = 1 : length(files)
        k = strfind(files(j).name,'_marked');
        color = [ '-g'; '-k'; '-r'];
%         figure(10); hold on; % row trajectories        
        for i = 1 : 3
            % load the data        
            if ~isempty(k)
                k1 = strfind(files(j).name,cond{i});
                if ~isempty(k1)   
                    load(files(j).name)
                    tti = 1;
                    for tt = 1 : N
%                         ontime = D.handON{1,tt};
                        if param(tt,24) == 0 && param(tt,25) == 0 && param(tt,26) == 0 % select good trials                            
                            on = find(D.t{tt,1}(:) > D.handON{1,tt},1);
                            off = size(D.handlepos{tt,1},1);
                            trajx(1:off-on+1,tti,i) = 18.875-D.handscreen{1,tt}(on:end,1);
                            trajy(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:end,2);
%                             figure(10); 
%                             H=plot(trajx(1:off-on+1,tti,i),trajy(1:off-on+1,tti,i),color(i,:));
%                             H.Color(4) = 0.2;
                            l(tti,i) = off-on+1;
                            tti = tti + 1;
                        end
                    end                   
                end
            end
        end        
    end
    q = l ~= 0;
%     s = min(l(q));
%     s = floor(s); % normalization factor
    s = max(l(q)) + 10; % fixed normalization for the pooled data
    mtraj = zeros(s,2,3); % mean trajectories
    straj = zeros(s,2,3); % strandard deviation of trajectories
    figure(5); hold on
%     color = [ '-g'; '-k'; '-r'];
    for i = 1 : 3
        [mtraj(1:s,:,i),straj(1:s,:,i)] = trajnormalization2(trajx(:,:,i),trajy(:,:,i),l(:,i),s);
        figure(5); hold on
        shadedErrorBar(mtraj(:,2,i),mtraj(:,1,i),straj(:,1,i),'lineprops',color(i,:),'patchSaturation',0.33)
        
    end  
    set(gca,'View',[90 -90])
    scatter([9 29 19 19 19],[18.875-18.875 18.875-18.875 18.875-18.875 25.975-18.875 11.775-18.875],'ro')
    figname = 's4-traj';
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    save('traj-s5.mat','mtraj','straj')

end
close all