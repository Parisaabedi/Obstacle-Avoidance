%% This function is for trajectory analysis
% normalize trajectories
% find the ones with the changed dirction even a small one! 

% Written by PA, June 2018

%% 
% Initialization
HR = [-30 0 30];
cond1 = { 'FL'; 'F0'; 'FR'};
cond2 = {'NFL' ; 'NF0' ; 'NFR'};
color = [ '-g'; '-k'; '-r'];
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
subs = 1:7; subs = [subs 9];
cond = {'F';'NF'};
markers = {'o';'+';'*';'>';'s';'d';'x';'^'};
obstaclepos = flipud([15.2750   11.0000
                      17.0750   11.0000
                      18.8750   11.0000
                      20.6750   11.0000
                      22.4750   11.0000]);

trajectories = 0;
overal = 0; 
oplot = 1;

if trajectories == 1    
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
            %figure(50); hold on
            % load the data        
            if ~isempty(k)
                k1 = strfind(files(j).name,['Marie' cond1{i}]);
                if ~isempty(k1)   
                    load(files(j).name)
                    load('CompF_s1.mat')
                    tti = 1;
                    for tt = 1 : N
                        %tt
%                         ontime = D.handON{1,tt};
                        if param(tt,24) == 0 && param(tt,25) == 0 && param(tt,27) == 0 && Comp{i,1}(tt,9)== 0 % select good trials                            
                            on = find(D.t{tt,1}(:) > D.handON{1,tt},1);
                            off1 = find(D.t{tt,1}(:) < D.handOFF{1,tt});
                            off = off1(end);
                            trajx(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:off,1);
                            trajy(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:off,2);
%                             time(1:off-on+1,tti,i) = D.t{tt,1}(on:end,2);
%                             figure(10); 
%                             H=plot(trajx(1:off-on+1,tti,i),trajy(1:off-on+1,tti,i),color(i,:));
%                             H.Color(4) = 0.2;
                            l(tti,i) = off-on+1;                            
                            tti = tti + 1;                            
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
%     s = max(l(q)) + 10; % fixed normalization for the pooled data
    s = 50;
    mtrajx = zeros(s,5,2,3); % mean trajectories along x-axis
    strajx = zeros(s,5,2,3); % strandard deviation of trajectories along x-axis
    %figure(5); hold on
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
        figname = 's1F-trajnormal';
        %figure(5);subplot(1,5,i);set(gca,'xlim',[-20 20])
        figure(i);%set(gca,'xlim',[-15 15])
        scatter([obstaclepos(i,1)-18.8750 0 0],[30-obstaclepos(i,2) 9 29],'ro')    
        fign = num2str(i);
        figname = [figname fign];
        savefig(gcf,figname)
        saveas(gcf,figname,'tiff')
    end      
    save('trajnormal-s1F.mat','mtrajx','strajx')

end
close all

%% Put the whole data together
if overal == 1
    mtrajxs = zeros(50,5,2,3,2,8);% ypositions,obstacles,R/L,HRs,F/NF,Subjects
    strajxs = zeros(50,5,2,3,2,8);
    for s = 1 : 8 % for all subjects    
        for c = 1 : 2 % for both with and without feedback conditions
            sID = strcat('s',num2str(subs(s)));
            filename = strcat('trajnormal-',sID,cond{c,1},'.mat');
            load(filename)
            mtrajxs(:,:,:,:,c,s) = mtrajx;
            strajxs(:,:,:,:,c,s) = strajx;
        end
    end
    save('trajnormal-overal.mat','mtrajxs','strajxs')
end

%% Plot the pooled data
if oplot == 1
    load('trajnormal-overal.mat')
    s = 50;    
    ys = 9.5:(29-9.5)/(s-1):29; ys = ys';
    for c = 1 : 2
        for i = 1 : 3
            for o = 1 : 5
                figure(o); hold on
                for d = 1:2
                    temp1 = mtrajxs(:,o,d,i,c,:);
                    temp2 = strajxs(:,o,d,i,c,:);
                    temp1 = reshape(temp1,50,8);
                    ms = nanmean(temp1,2);
                    temp2 = reshape(temp2,50,8);
                    ss = nanmean(temp2,2);
                    if sum(ss)~=0
                        tttempt = ms - 18.875;
                        uE=tttempt+ss;
                        lE=tttempt-ss;
                        xP=[lE;flipud(uE)];
                        yP=[ys;flipud(ys)];
                        patch(xP,yP,color(i,:),'FaceAlpha',0.3,'EdgeColor','none');
                        plot(tttempt,ys,color(i,:),'Linewidth',2)       
                    end
%                     for s = 1 : 8
%                         if ~isnan(temp1(:,s))
%                             scatter(temp1(:,s)-18.875,ys,'MarkerEdgeColor',colorss(i,:),'Marker',markers{s});
%                         end
%                     end

                end
            end
        end
        for i = 1 : 5
            figname = 'overal-trajnormal';           
            figure(i);%set(gca,'xlim',[-15 15])
            scatter([obstaclepos(i,1)-18.8750 0 0],[30-obstaclepos(i,2) 9 29],'ro')   
            fign = num2str(i);
            figname = strcat(figname,cond{c,1},fign);            
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')
        end 
        close all
    end
end