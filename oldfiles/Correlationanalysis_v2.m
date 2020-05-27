%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by PAK (Dec. 2018) -updated Jan. 2019
% plotting the correlation between several parameters:
% Max curvature and Reaction time
% Movement time and peak velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcrtcor = zeros(6,5,2,18); % correlation (cond, obstacle, direction,subs)
% for s = 1:18
%     filename = strcat('Comp3_s',num2str(s),'.mat');
%     load(filename)    
%     for cond = 1:6
%         temp = params(1:Ns(cond),:,cond);       
%         for j = 1 : 2
%             for t = 1 : 5
%                 figure(s);subplot(6,5,t+(cond-1)*5); hold on
%                 qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
%                       temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0; % repeats are not allowed
%                 rts = temp(qnc,12);
%                 mcs = temp(qnc,19);
%                 scatter(rts,(-2*(j-1)+1)*mcs)
%                 covtemp = cov(rts,mcs);
%                 mcrtcor(cond,t,j,s) = covtemp(1,2);
%             end
%         end
%     end
%     figname = strcat(['RTMCcorr-s# ',num2str(s)]); 
%     savefig(gcf,figname) 
%     saveas(gcf,figname,'tiff')  
% end
% close all

% mcmtcor = zeros(6,5,2,18); % correlation (cond, obstacle, direction,subs)
% for s = 1:18
%     filename = strcat('Comp3_s',num2str(s),'.mat');
%     load(filename)    
%     for cond = 1:6
%         temp = params(1:Ns(cond),:,cond);       
%         for j = 1 : 2
%             for t = 1 : 5
%                 figure(s);subplot(6,5,t+(cond-1)*5); hold on
%                 qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
%                       temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0; % repeats are not allowed
%                 mts = temp(qnc,13)-temp(qnc,12);
%                 mcs = temp(qnc,19);
%                 scatter(mts,(-2*(j-1)+1)*mcs)
%                 covtemp = cov(mts,mcs);
%                 mcmtcor(cond,t,j,s) = covtemp(1,2);
%             end
%         end
%     end
%     figname = strcat(['MTMCcorr-s# ',num2str(s)]); 
%     savefig(gcf,figname) 
%     saveas(gcf,figname,'tiff')  
% end
% close all
%% Initialization

plotting = 0;
plotting2 = 0;
poolplotting = 1;
poolplotting2 = 0;
vf = 1;
color = {'r', 'b', 'g', 'k', 'm','y'};
if plotting == 1
    %  3 pairs of obstacles, HRs, MT, velocity and MaxCurvature, Mean&STD, #subjects
    mtmc = zeros(3,3,3,2,18); 
    %  6 pairs of obstacles, HRs, MT, velocity and MaxCurvature, Mean&STD, #subjects
    mtmcc = zeros(6,3,3,2,18);%6 pairs of obstacles(1L,2L,3L,3R,4R,5R)
    CR = zeros(3,3,2,18); % 3 pairs of obstacles, HRs, Coll. and Rep., #subjects
    c = 1; % with/without visual feedback (0=> with, 1=> without)
    obs = [ 5 1 ; 4 2; 3 3 ];
    dir = [0 1 ; 0 1 ; 0 1];
    for s = 1 : 18
        filename = strcat('Comp3_s',num2str(s),'.mat');
        load(filename)
        for cond = 1:3 
            temp = params(1:Ns(cond+c*3),:,cond+c*3);
            for j = 1 : 3
                % no collision, no repeat, and either right to left -ward
                % directions, good trials, and no changed direction or predicted movement
                qnc1 = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == dir(j,1) & temp(:,8) == obs(j,1) & temp(:,24) == 0; 
                qnc2 = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == dir(j,2) & temp(:,8) == obs(j,2) & temp(:,24) == 0; 
                mtmc(j,cond,1,1,s) = mean([temp(qnc1,13);temp(qnc2,13)]-[temp(qnc1,12);temp(qnc2,12)]);
                mtmc(j,cond,1,2,s) = std([temp(qnc1,13);temp(qnc2,13)]-[temp(qnc1,12);temp(qnc2,12)]);
                mtmc(j,cond,2,1,s) = mean([temp(qnc1,14);temp(qnc2,14)]);
                mtmc(j,cond,2,2,s) = std([temp(qnc1,14);temp(qnc2,14)]);
                mtmc(j,cond,3,1,s) = mean([temp(qnc1,19);temp(qnc2,19)]);
                mtmc(j,cond,3,2,s) = std([temp(qnc1,19);temp(qnc2,19)]);
                qnc3 = temp(:,8) ==  obs(j,1) | temp(:,8)== obs(j,2);
                CR(j,cond,1,s) = sum(temp(qnc3,5));
                CR(j,cond,2,s) = sum(temp(qnc3,4));
            end
        end    
    end

    for s = 1:18
        figure(1);subplot(6,3,s); hold on % plot movement time
        errorbar(1:3,mtmc(1:3,1,1,1,s),mtmc(1:3,1,1,2,s),'-k') % HR = 0
        k1 = mean(reshape(mtmc(1:3,2:3,1,1,s), 3 ,2),2); %HR !=0 , mean
        k2 = mean(reshape(mtmc(1:3,2:3,1,2,s), 3 ,2),2); %HR !=0 , std
        errorbar(1:3,k1,k2,':b')
        ylim([0,0.8])
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'1&5','2&4','3'})

        figure(2);subplot(6,3,s); hold on % plot maximum velocity
        errorbar(1:3,mtmc(1:3,1,2,1,s),mtmc(1:3,1,2,2,s),'-k') % HR = 0
        k1 = mean(reshape(mtmc(1:3,2:3,2,1,s), 3 ,2),2); %HR !=0 , mean
        k2 = mean(reshape(mtmc(1:3,2:3,2,2,s), 3 ,2),2); %HR !=0 , std
        errorbar(1:3,k1,k2,':b')
        %ylim([0,10])
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'1&5','2&4','3'})        
        
        figure(3);subplot(6,3,s); hold on % plot maximum curvature
        errorbar(1:3,mtmc(1:3,1,3,1,s),mtmc(1:3,1,3,2,s),'-k') % HR = 0
        k1 = mean(reshape(mtmc(1:3,2:3,3,1,s), 3 ,2),2); %HR !=0 , mean
        k2 = mean(reshape(mtmc(1:3,2:3,3,2,s), 3 ,2),2); %HR !=0 , std
        errorbar(1:3,k1,k2,':b')
        ylim([0,10])
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'1&5','2&4','3'})
        
        figure(4);subplot(6,3,s); hold on % plot collisions
        plot(1:3,CR(1:3,1,1,s),'-k') % HR = 0        
        plot(1:3,CR(1:3,2,1,s)+CR(1:3,3,1,s),':b') % HR != 0
        ylim([0,30])
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'1&5','2&4','3'})
        
        figure(5);subplot(6,3,s); hold on % plot repetitions
        plot(1:3,CR(1:3,1,2,s),'-k') % HR = 0        
        plot(1:3,CR(1:3,2,2,s)+CR(1:3,3,2,s),':b') % HR != 0
        ylim([0,100])
        set(gca,'XTick',1:3)
        set(gca, 'XTickLabels',{'1&5','2&4','3'})
    end
    figure(1)
    figname = 'Movement Time - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 

    figure(2)
    figname = 'Maximum Velocity - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    

    figure(3)
    figname = 'Maximum Curvature - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    close all
    
    figure(4)
    figname = 'Collisions - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')    
    
    figure(5)
    figname = 'Repetitions - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff') 
    close all
    
    expansions = [-0.0448 0.5199 0.7062 -0.0317 0.1265 0.1 1.0347 0.4292 0.4228...
              0.1210 0.0114 0.1141 0.7644 0.0325 0.0955 -.2259 -0.3771 -0.0260; ...
              -0.0285 0.268 0.4449 -0.3926 0.1414 0.0525 0.8495 0.3593 0.2670 ...
              -0.0021 0.0878 -0.2776 0.6556 0.0554 -0.1537 -0.0835 -0.3959 0.5169;...
              0 0.4561 0.2741 0 -0.261 0 0 0.0525 0.314 0.4104 0.0780 -0.4163 ...
              0.4165 0 -0.1844 -0.2774 -0.3191 0.0267];
    save('corr.mat','mtmc','expansions','CR')
end

% plot the pooled data (once for average and once for median)
if poolplotting == 1
    %color = {'r', 'b', 'g'};
    load('corr.mat')
    deltas = zeros(3,3,18); % obstacles, params(MT,Vel, MC), subjects
    m = zeros(3,1); % median of time for each obstacle pair
    CRmedianstd = zeros(3,2); % the median and IQR
    % the compensation factor: It is a signed parameters evaluating the
    % amount of compensation for the added noise (Positive -> good compensation)
    compensation = zeros(3,18);
    compensation2 = zeros(3,18);
    failure = sum(CR,3);
%     nbins = 25;
%     for j = 1:3
%         figure(1);subplot(3,3,j+3*0); hold on
%         histogram(failure(j,1,1,:),nbins,'FaceColor','k')
%         histogram(sum(failure(j,2:3,1,:),2),nbins,'FaceColor','r')
%         xlim([0,60])
%         figure(1);subplot(3,3,j+3*1); hold on
%         histogram(CR(j,1,1,:),nbins,'FaceColor','k')
%         histogram(sum(CR(j,2:3,1,:),2),nbins,'FaceColor','r')
%         xlim([0,60])
%         figure(1);subplot(3,3,j+3*2); hold on
%         histogram(CR(j,1,2,:),nbins,'FaceColor','k')
%         histogram(sum(CR(j,2:3,2,:),2),nbins,'FaceColor','r')
%         xlim([0,60])
%     end
    for j = 1:3        
        for p = 1:3
            %figure(p); hold on
            k1 = reshape(mtmc(j,1,p,1,:),1,18);
            k2 = mean([reshape(mtmc(j,2,p,1,:),1,18);reshape(mtmc(j,3,p,1,:),1,18)]);
            if p == 1
                %deltas(j,p,:) = k1;
                m(j) = median(k1);
%             else
%                 deltas(j,p,:) = (k2 - k1);
            end
            deltas(j,p,:) = (k2 - k1)./k1;
            %scatter(deltas(j,p,:),expansions(j,:),'MarkerEdgeColor',color{j})               
        end
        figure(4); hold on
        %scatter(deltas(j,2,:),deltas(j,3,:),'MarkerEdgeColor',color{j});% curveture against velocity
        scatter(deltas(j,2,:),expansions(j,:),'MarkerEdgeColor',color{j});% expansion against velocity
%         for k = 1:length(deltas(j,2,:))
%             text(deltas(j,2,k),deltas(j,3,k),num2str(k))
%         end
        compensation(j,:) = sqrt(deltas(j,2,:).^2 + deltas(j,3,:).^2);
        angle = atan2d(deltas(j,3,:),deltas(j,2,:)); angle = angle(:)';
        compensation(j,(angle < 0 & angle > -90)) = ...
            -1 *compensation(j,(angle < 0 & angle > -90));
        compensation2(j,:) = sqrt(reshape(deltas(j,2,:),1,18).^2 + expansions(j,:).^2);
        angle = atan2d(expansions(j,:),reshape(deltas(j,2,:),1,18)); angle = angle(:)';
        compensation2(j,(angle < 0 & angle > -90)) = ...
            -1 *compensation2(j,(angle < 0 & angle > -90)); 
        
            %         
%         figure(5); hold on
%         errorbar(1:3,mean(deltas(:,1,:),3),std(deltas(:,1,:),3),'-k')
%         errorbar(1:3,mean(deltas(:,1,:),3),std(deltas(:,1,:),3),'-k')
%         figure(p+5);  hold on
%         pointer = deltas(j,1,:) > m(j);
%         k1 = median([reshape(deltas(j,2,pointer),[],1) reshape(deltas(j,2,~pointer),[],1)]);
%         k2 = std([reshape(deltas(j,2,pointer),[],1) reshape(deltas(j,2,~pointer),[],1)]);
%         plot(1:2,k1,'Color',color{j})  
%         figure(p+6);  hold on
%         pointer = deltas(j,1,:) > m(j);
%         k1 = mean([reshape(deltas(j,3,pointer),[],1) reshape(deltas(j,3,~pointer),[],1)]);
%         k2 = std([reshape(deltas(j,3,pointer),[],1) reshape(deltas(j,3,~pointer),[],1)]);
%         plot(1:2,k1,'Color',color{j})  
    end  
    figure(5); hold on
    for ss = 1 : 18
        plot(1:3, compensation2(:,ss),':k','LineWidth',0.5)
    end
    errorbar(1:3, mean(compensation2,2),std(compensation2,0,2),'-k','LineWidth',2.5)
    xlabel('Task difficualty'); ylabel('Compensation')
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'easy' , 'Medium' , 'Hard'})
    xlim([0.5 3.5])    

    figname = 'Compensation - No visual feedback'; 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')
    close all 
end



% Run for individual obstacles (1L,2L,3L,3R,4R,5R)
if plotting2 == 1
    %  6 pairs of obstacles, HRs, MT, velocity and MaxCurvature, Mean&STD, #subjects
    mtmcc = zeros(6,3,3,2,18);%6 pairs of obstacles(1L,2L,3L,3R,4R,5R)
    
    c = 1; % with/without visual feedback (0=> with, 1=> without)
    obs = [ 1 2 3 3 4 5];
    dir = [ 1 1 1 0 0 0];
    for s = 1 : 18
        filename = strcat('Comp3_s',num2str(s),'.mat');
        load(filename)
        for cond = 1:3 
            temp = params(1:Ns(cond+c*3),:,cond+c*3);
            for j = 1 : 6
                % no collision, no repeat, and either right to left -ward
                % directions, good trials, and no changed direction or predicted movement
                qnc1 = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == dir(j) & temp(:,8) == obs(j) & temp(:,24) == 0; 
                
                mtmcc(j,cond,1,1,s) = mean(temp(qnc1,13)-temp(qnc1,12));
                mtmcc(j,cond,1,2,s) = std(temp(qnc1,13)-temp(qnc1,12));
                mtmcc(j,cond,2,1,s) = mean(temp(qnc1,14));
                mtmcc(j,cond,2,2,s) = std(temp(qnc1,14));
                mtmcc(j,cond,3,1,s) = mean(temp(qnc1,19));
                mtmcc(j,cond,3,2,s) = std(temp(qnc1,19));                
            end
        end    
    end

    
    
    expansionsobs = alfabetap;
    save('corrobs.mat','mtmcc','expansionsobs')
    
    color = {'r', 'b', 'g','k','m','y'};
    load('corrobs.mat')
    deltas = zeros(6,3,18); % obstacles, params(MT,Vel, MC), subjects   
    compensation = zeros(6,18);
    compensation2 = zeros(6,18);
    failure = sum(CR,3);
    
     for j = 1:6        
        for p = 1:3
            k1 = reshape(mtmcc(j,1,p,1,:),1,18);
            k2 = mean([reshape(mtmcc(j,2,p,1,:),1,18);reshape(mtmcc(j,3,p,1,:),1,18)]);           
            deltas(j,p,:) = (k2 - k1)./k1;                          
        end
        figure(4); hold on
        %scatter(deltas(j,2,:),deltas(j,3,:),'MarkerEdgeColor',color{j});% curveture against velocity
        scatter(deltas(j,2,:),expansionsobs(j,2,:),'MarkerEdgeColor',color{j});% expansion against velocity

        compensation(j,:) = sqrt(deltas(j,2,:).^2 + deltas(j,3,:).^2);
        angle = atan2d(deltas(j,3,:),deltas(j,2,:)); angle = angle(:)';
        compensation(j,(angle < 0 & angle > -90)) = ...
            -1 *compensation(j,(angle < 0 & angle > -90));
        compensation2(j,:) = sqrt(reshape(deltas(j,2,:),1,18).^2 + reshape(expansionsobs(j,2,:),1,18).^2);
        angle = atan2d(reshape(expansionsobs(j,2,:),1,18),reshape(deltas(j,2,:),1,18)); angle = angle(:)';
        compensation2(j,(angle < 0 & angle > -90)) = ...
            -1 *compensation2(j,(angle < 0 & angle > -90));   
     end 
    
     figure(5); hold on
    for ss = 1 : 18
        plot(1:6, compensation2(:,ss),':k','LineWidth',0.5)
    end
    errorbar(1:6, nanmean(compensation2,2),nanstd(compensation2,0,2),'-k','LineWidth',2.5)
    xlabel('Obstalce positions'); ylabel('Compensation')
    set(gca,'XTick',1:6)
    set(gca, 'XTickLabels',{'Most Rightward','Right','Center-L','Center-R','Left','Most Leftward'})
    xlim([0.5 6.5])   
    
end

    
% Run statistical analysis

% for with visual feedback
if vf == 1
    mtmc2 = zeros(3,3,3,2,18); 
    c = 0; % with/without visual feedback (0=> with, 1=> without)
    valid = 1; % check if there is enough data on both right and left direction
    obs = [ 5 1 ; 4 2; 3 3 ];
    dir = [0 1 ; 0 1 ; 0 1];   
    for s = 1 : 18
        filename = strcat('Comp3_s',num2str(s),'.mat');
        load(filename)
        for cond = 1:3 
            temp = params(1:Ns(cond+c*3),:,cond+c*3);
            for j = 1 : 3
                % no collision, no repeat, and either right to left -ward
                % directions, good trials, and no changed direction or predicted movement
                qnc1 = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == dir(j,1) & temp(:,8) == obs(j,1) & temp(:,24) == 0; 
                qnc2 = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == dir(j,2) & temp(:,8) == obs(j,2) & temp(:,24) == 0; 
                mtmc2(j,cond,1,1,s) = mean([temp(qnc1,13);temp(qnc2,13)]-[temp(qnc1,12);temp(qnc2,12)]);
                mtmc2(j,cond,1,2,s) = std([temp(qnc1,13);temp(qnc2,13)]-[temp(qnc1,12);temp(qnc2,12)]);
                mtmc2(j,cond,2,1,s) = mean([temp(qnc1,14);temp(qnc2,14)]);
                mtmc2(j,cond,2,2,s) = std([temp(qnc1,14);temp(qnc2,14)]);
                mtmc2(j,cond,3,1,s) = mean([temp(qnc1,19);temp(qnc2,19)]);
                mtmc2(j,cond,3,2,s) = std([temp(qnc1,19);temp(qnc2,19)]);
            end
        end    
    end
    save('corr2.mat','mtmc2')
end

if poolplotting2 == 1
    color = {'r', 'b', 'g'};
    load('corr2.mat')
    deltas = zeros(3,3,18); % obstacles, params (MT,Vel, MC), subjects
    
    for j = 1:3
        for p = 1:3
            k1 = reshape(mtmc2(j,1,p,1,:),1,18);
            k2 = mean([reshape(mtmc2(j,2,p,1,:),1,18);reshape(mtmc2(j,3,p,1,:),1,18)]);
            deltas(j,p,:) = (k2 - k1)./k1;            
        end
        figure(5); hold on
        scatter(deltas(j,2,:),deltas(j,3,:),'MarkerEdgeColor',color{j})
        for k = 1:length(deltas(j,2,:))
            text(deltas(j,2,k),deltas(j,3,k),num2str(k))
        end
    end
    
end
for s = 1:18
    figure(1);subplot(6,3,s); hold on % plot movement time
    errorbar(1:3,mtmc2(1:3,1,1,1,s),mtmc2(1:3,1,1,2,s),'-k') % HR = 0
    k1 = mean(reshape(mtmc2(1:3,2:3,1,1,s), 3 ,2),2); %HR !=0 , mean
    k2 = mean(reshape(mtmc2(1:3,2:3,1,2,s), 3 ,2),2); %HR !=0 , std
    errorbar(1:3,k1,k2,':b')
    ylim([0,0.8])
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'1&5','2&4','3'})
    
    figure(2);subplot(6,3,s); hold on % plot maximum curvature
    errorbar(1:3,mtmc2(1:3,1,2,1,s),mtmc2(1:3,1,2,2,s),'-k') % HR = 0
    k1 = mean(reshape(mtmc2(1:3,2:3,2,1,s), 3 ,2),2); %HR !=0 , mean
    k2 = mean(reshape(mtmc2(1:3,2:3,2,2,s), 3 ,2),2); %HR !=0 , std
    errorbar(1:3,k1,k2,':b')
    ylim([0,10])
    set(gca,'XTick',1:3)
    set(gca, 'XTickLabels',{'1&5','2&4','3'})
end
figure(1)
figname = 'Movement Time - visual feedback '; 
savefig(gcf,figname) 
saveas(gcf,figname,'tiff') 

figure(2)
figname = 'Maximum Curvature - visual feedback '; 
savefig(gcf,figname) 
saveas(gcf,figname,'tiff') 
close all
