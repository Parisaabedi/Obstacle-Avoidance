%% This function plot the individual data for Obstacle Avoidance
% Written by PA, June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
subs = 1:7; subs = [subs 9];
cond = {'F_';'NF_'};
colors = { 'g' ; 'k' ; 'r'};
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
markers = {'o';'+';'*';'>';'s';'d';'x';'^'};
% individual participants
% Initialization (HR conditions, mean/std, target, left/right,F/NF,subjects)
ima = zeros(3,2,5,2,2,8); % initilal movement angle    
mcurv = zeros(3,2,5,2,2,8); % movement curvature
DfromO = zeros(3,2,5,2,2,8); % Distance from Obstacle
NCs = zeros(3,5,2,2,8); % number of collisions
NR = zeros(3,5,2,2,8); % number of repeats
imas = { ima; ima; ima};
mcurvs = { mcurv; mcurv; mcurv};
DfromOs = { DfromO; DfromO; DfromO};
NRs = {NR; NR};
for s = 1:7 % for all subjects    
    for c = 1 : 2 % for both with and without feedback conditions
        sID = strcat('s',num2str(subs(s)));
        filename = strcat('Comp',cond{c,1},sID,'.mat');
        load(filename)
        for i = 1 : 3
            for j = 1 : 2
                for t = 1 : 5
                    % no collision, no repeat, and either right to left -ward
                    % directions, no weird mvs => do this for each specific obstacle
                    qnc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 0 & ...
                          Comp{i}(:,7) == (j-1) & Comp{i}(:,8) == t & Comp{i}(:,9) == 0; 
                    imas{1,1}(i,1,t,j,c,s) = mean(Comp{i}(qnc,3));
                    imas{1,1}(i,2,t,j,c,s) = std(Comp{i}(qnc,3));
                    mcurvs{1,1}(i,1,t,j,c,s) = mean(Comp{i}(qnc,4));
                    mcurvs{1,1}(i,2,t,j,c,s) = std(Comp{i}(qnc,4));
                    DfromOs{1,1}(i,1,t,j,c,s) = mean(Comp{i}(qnc,5));
                    DfromOs{1,1}(i,2,t,j,c,s) = std(Comp{i}(qnc,5));
                    qrc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 1 & ...
                          Comp{i}(:,7) == (j-1) & Comp{i}(:,8) == t & Comp{i}(:,9) == 0; 
                    NRs{1,1}(i,t,j,c,s) = sum(qrc);
                    % collision, no repeat, and either right or left -ward
                    % directions => do this for each specific obstacle
                    qc = Comp{i}(:,1) == 1 & Comp{i}(:,7) == (j-1) & ...
                         Comp{i}(:,8) == t  & Comp{i}(:,9) == 0; 
                    imas{2,1}(i,1,t,j,c,s) = mean(Comp{i}(qc,3));
                    imas{2,1}(i,2,t,j,c,s) = std(Comp{i}(qc,3));
                    mcurvs{2,1}(i,1,t,j,c,s) = mean(Comp{i}(qc,4));
                    mcurvs{2,1}(i,2,t,j,c,s) = std(Comp{i}(qc,4));
                    DfromOs{2,1}(i,1,t,j,c,s) = mean(Comp{i}(qc,5));
                    DfromOs{2,1}(i,2,t,j,c,s) = std(Comp{i}(qc,5));
                    NCs(i,t,j,c,s) = sum(qc);
                    % weird movements, no repeat, and either right or left -ward
                    % directions => do this for each specific obstacle
                    qcw = Comp{i}(:,1) == 0 & Comp{i}(:,7) == (j-1) & ...
                          Comp{i}(:,8) == t  & Comp{i}(:,9) == 1;
                    imas{3,1}(i,1,t,j,c,s) = mean(Comp{i}(qcw,3));
                    imas{3,1}(i,2,t,j,c,s) = std(Comp{i}(qcw,3));
                    mcurvs{3,1}(i,1,t,j,c,s) = mean(Comp{i}(qcw,4));
                    mcurvs{3,1}(i,2,t,j,c,s) = std(Comp{i}(qcw,4));
                    DfromOs{3,1}(i,1,t,j,c,s) = mean(Comp{i}(qcw,5));
                    DfromOs{3,1}(i,2,t,j,c,s) = std(Comp{i}(qcw,5));
                    qrc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 1 & ...
                          Comp{i}(:,7) == (j-1) & Comp{i}(:,8) == t & Comp{i}(:,9) == 1;
                    NRs{2,1}(i,t,j,c,s) = sum(qrc);
                end
            end
        end
    end
end

%%
for c = 1 : 2
    
    %%
    % Plot Initial movement angles   
    % for each obstacle
    figure(1); subplot(2,1,1);hold on
    xa = 5:-1:1;
    h = zeros(1,3);
    for i = 1:3 % for each HR condition
        t1 = reshape(imas{1,1}(i,1,:,1,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(imas{1,1}(i,2,:,1,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        q = t1(1:8,1) ~= 0;
        h(1,i)=errorbar(1:5,abs(nanmean(t1(q,:))-90),nanmean(t2(q,:)),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            if q(ss) == 1
                scatter(1:5, abs(90-t1(ss,:)),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss}); 
            end
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'},'Location','best')
    %legend('Location','southeast')    
    ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    %title('Initial movement angle modulation by varying head roll - Rightward MVs')
    title('Rightward MVs')

    figure(1); subplot(2,1,2);hold on
    for i = 1:3 % for each HR condition
        t1 = reshape(imas{1,1}(i,1,:,2,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(imas{1,1}(i,2,:,2,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        q = t1(1:8,1) ~= 0;
        errorbar(1:5,abs(nanmean(t1(q,:))-90),nanmean(t2(q,:)),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            if q(ss) == 1
                scatter(1:5, abs(90-t1(ss,:)),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});
            end
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    %title('Initial movement angle modulation by varying head roll - Leftward MVs')
    title('Leftward MVs')

    figname = strcat('Overall-', cond{c,1}, 'IMA-OBs');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')

    %%
    % Plot curvature
    % for each obstacle
    figure(2); subplot(2,1,1);hold on 
    xa = 5:-1:1;
    h = zeros(1,3);
    for i = 1:3 % for each HR condition
        t1 = reshape(mcurvs{1,1}(i,1,:,1,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(mcurvs{1,1}(i,2,:,1,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        h(1,i)=errorbar(1:5,nanmean(t1),nanmean(t2),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            scatter(1:5, t1(ss,:),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});       
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'},'Location','best')   
    ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    %title('Maximum movement curvature modulation by varying head roll - Rightward MVs')
    title('Rightward MVs')

    figure(2); subplot(2,1,2);hold on
    for i = 1:3 % for each HR condition
        t1 = reshape(mcurvs{1,1}(i,1,:,2,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(mcurvs{1,1}(i,2,:,2,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        errorbar(1:5,nanmean(t1),nanmean(t2),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            scatter(1:5, t1(ss,:),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});       
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    xlabel('Obstacle position'); ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    %title('Maximum movement curvature modulation by varying head roll - Leftward MVs')
    title('Leftward MVs')
    
    figname = strcat('Overall-',cond{c,1},'-MC-OBS');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    
    %%
    % Plot distance from obstacle 
    % for each obstacle
    figure(3); subplot(2,1,1);hold on 
    xa = 5:-1:1;
    h = zeros(1,3);
    for i = 1:3 % for each HR condition
        t1 = reshape(DfromOs{1,1}(i,1,:,1,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(DfromOs{1,1}(i,2,:,1,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        h(1,i)=errorbar(1:5,nanmean(abs(t1)),nanmean(t2),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            scatter(1:5, abs(t1(ss,:)),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});       
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'},'Location','best')   
    ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
    plot(xa,[nan nan nan nan nan])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'','','','',''}))
    xlim([0,6])
    %title('Distance form obstacle modulation by varying head roll - Rightward MVs')  
    title('Rightward MVs')            

    figure(3); subplot(2,1,2);hold on
    for i = 1:3 % for each HR condition
        t1 = reshape(DfromOs{1,1}(i,1,:,2,c,:),5,8);t1 = t1';t1 = fliplr(t1);
        t2 = reshape(DfromOs{1,1}(i,2,:,2,c,:),5,8);t2 = t2';t2 = fliplr(t2);
        errorbar(1:5,nanmean(abs(t1)),nanmean(t2),'Color',colors{i});
        for ss = 1 : 8 % for each obstacle
            scatter(1:5, abs(t1(ss,:)),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});       
            %plot(1:5, abs(90-t1(ss,:)),'Color',colorss(i,:),'Marker',markers{ss});       
        end
    end  
    xlabel('Obstacle position');ylabel('Distance from Obstacle'); 
    plot(xa,[nan nan nan nan nan])
    xlim([0,6])
    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    %title('Distance form obstacle modulation by varying head roll - Leftward MVs')            
    title('Leftward MVs')            
    
    figname = strcat('Overall-',cond{c,1},'-dfo-OBS');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    %%
    % Plot Number of collisition
    figure(4); subplot(1,2,1); hold on
    HR = [-30 0 30];
    for i = 1:3
        t1 = reshape(NCs(i,:,1,c,:),5,8);t1 = sum(t1);
        bar(HR(i), mean(t1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
        for ss = 1 : 8
            scatter(HR(i),t1(ss),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss})
        end
    end
    xlabel('Head Rotations (deg)'); ylabel('Obstacle collision - Rightward MVs')
    
    figure(4); subplot(1,2,2); hold on
    for i = 1:3
        t1 = reshape(NCs(i,:,2,c,:),5,8);t1 = sum(t1);
        bar(HR(i), mean(t1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
        for ss = 1 : 8
            scatter(HR(i),t1(ss),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss})
        end
    end
    xlabel('Head Rotations (deg)'); ylabel('Obstacle collision - Leftward MVs')
    %supertitle('Obstacle collision modulation by varying head roll')

    figname = strcat('Overall-',cond{c,1},'-NC');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
     %%   
    % Plot Number of right vs left -ward movements 
    figure(5); hold on               
    rl = zeros(3,5,2);
    h = zeros(1,3);
    pr = zeros(3,5);
    for i = 1:3
        r = reshape(NRs{1,1}(i,:,1,c,:),5,8);r = flipud(r);
        l = reshape(NRs{1,1}(i,:,2,c,:),5,8);l = flipud(l);
        pr = r.*100./(r+l); 
        h(1,i)=plot(1:5,nanmean(pr'),'Color',colors{i}); 
        for ss = 1 : 8
            scatter(1:5, pr(:,ss),'MarkerEdgeColor',colorss(i,:),'Marker',markers{ss});
        end
    end        
    legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'},'Location','best')
    %legend('Location','southeast') 
    xlabel('Obstacle position'); ylabel('Percentage of Rightward MVs')    

    set(gca,'XTick',fliplr(xa))
    set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
    %title('Movement Directions')
    xlim([0,6])
    figname = strcat('Overall-',cond{c,1},'LR-MDir');
    savefig(gcf,figname)
    saveas(gcf,figname,'tiff')
    close all
end
