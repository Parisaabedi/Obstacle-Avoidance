%% Obstacle Avoidacne Data analysis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HR = [-30 0 30];
cond1 = { 'FL'; 'F0'; 'FR'};
cond2 = {'NFL' ; 'NF0' ; 'NFR'};
colors = { 'g' ; 'k' ; 'r'};
init = 0;
pooled = 0;
poolploting = 0;
plotting = 0;
trajectories = 1;
if init == 1 
    % columns : 1-collision, 2-repeat, 3-initial movment angle, 4-movement curvature, 
    % 5-distance from obstacle, 6-movement duration, 7- movement direction
    % with regard to the obstacle(right = 0 , left = 1)
    % 8- Obstacle order
    Comp = cell(3,1);
    files = dir;

    for j = 1 : length(files)
        k = strfind(files(j).name,'_marked');
        for i = 1 : 3
            % load the data        
            if ~isempty(k)
                k1 = strfind(files(j).name,['5' cond1{i}]);
                if ~isempty(k1)   
                    load(files(j).name)
                    tk = size(param,1) - size(D.repeat,1);
                    if abs(tk) > 0
                        D.repeat = [D.repeat ; zeros(tk,1)];     
                    end
                    if strfind(filename,'_marked_marked')
                        m = matfile(filename(1:end-7), 'Writable', true);

                    else
                        m = matfile(filename, 'Writable', true);
                    end   
                    if size(param,2) < 27
                        param = [param D.repeat];
                        m.param = param;                                
                    end
                    q =  param(:,24) == 0; % Select valid trials 
                    os = D.obstacle_posi(:,1);
                    os = fliplr(os);                    
                    for ir = 1 : size(param,1)
                        if param(ir,25) == 0
                            Rarr = find(D.t{ir} >= D.handON{ir}, 1 ):find(D.t{ir} <= D.handOFF{ir}, 1, 'last' ); % movement period  
                            no = D.obstacleorder(ir);
                            temp = find(D.handscreen{ir}(Rarr,2) >= 30-D.obstacle_posi(no,2));                             
                            param(i,23) = D.handscreen{ir}(Rarr(1)+temp(1),1)-os(no,1); % Disntace from obstacle
                        end
                    end
                    Comp{i,1} = zeros(sum(q),8);
                    Comp{i,1}(:,1) = param(q,25); % collisions
                    Comp{i,1}(:,2) = param(q,27); % repeated trail (too long)
                    Comp{i,1}(:,3) = param(q,8);  % Initial movement angle
                    Comp{i,1}(:,4) = param(q,22); % Maximum curvature
                    Comp{i,1}(:,5) = param(q,23); % Distance from obstacle
                    Comp{i,1}(:,6) = param(q,17); % Movement duration
                    Comp{i,1}(:,7) = param(q,26); % Movement direction based on obstacle
                    Comp{i,1}(:,8) = D.obstacleorder(q); % obstacle orders

                end
            end

        end
    end
    save('CompF_s5.mat','Comp')
elseif init == 2
    load('Comp.mat')
end
%% Pooling
if pooled == 1
    files = dir;
    % columns : 1-collision, 2-repeat, 3-initial movment angle, 4-movement curvature, 
    % 5-distance from obstacle, 6-movement duration, 7-subject ID
    Comps = cell(3,1); 
    mtrajs = zeros(80,5*2,3);
    strajs = zeros(80,5*2,3);
    for j = 1 : length(files)
        k = strfind(files(j).name,'CompNF');
        if ~isempty(k)
            load(files(j).name)
            sidp = strfind(files(j).name,'_s');
            sid = str2double(files(j).name(sidp+2:end-4));
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
    save('CompNFs.mat','Comps')
    %save('trajs.mat','mtrajs','strajs')
elseif pooled == 2
    load('Comps.mat')
    %load('trajs.mat')
end
%%
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

%%
if plotting == 1
    subs = 1:7; subs = [subs 9];
    cond = {'F';'NF'};
    for s = 1 : 8
        for c = 1 : 2
            sID = strcat('s',num2str(subs(s)));
            filename = strcat('Comp',cond{c,1},'_',sID,'.mat');
            load(filename)
            % individual participants
            % Initialization (condition, mean/std, target, left/right)
            ima = zeros(3,2,5,2); % initilal movement angle    
            mcurv = zeros(3,2,5,2); % movement curvature
            DfromO = zeros(3,2,5,2); % Distance from Obstacle
            NC = zeros(3,5,2); % number of collisions
            NR = zeros(3,5,2); % number of repeats
            for i = 1 : 3
                for j = 1 : 2
                    for t = 1 : 5
                        % no collision, no repeat, and either right to left -ward
                        % directions, no weird mvs => do this for each specific obstacle
                        qnc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 0 & ...
                              Comp{i}(:,7) == (j-1) & Comp{i}(:,8) == t & Comp{i}(:,9) == 0; 
                        % collision, no repeat, and either right or left -ward
                        % directions => do this for each specific obstacle
                        qc = Comp{i}(:,1) == 1 & Comp{i}(:,7) == (j-1) & ...
                             Comp{i}(:,8) == t  & Comp{i}(:,9) == 0; 
                        % weird movements, no repeat, and either right or left -ward
                        % directions => do this for each specific obstacle
                        qcw = Comp{i}(:,1) == 0 & Comp{i}(:,7) == (j-1) & ...
                              Comp{i}(:,8) == t  & Comp{i}(:,9) == 1;
                        ima(i,1,t,j,c) = mean(Comp{i}(qnc,3));
                        ima(i,2,t,j,c) = std(Comp{i}(qnc,3));
                        mcurv(i,1,t,j,c) = mean(Comp{i}(qnc,4));
                        mcurv(i,2,t,j,c) = std(Comp{i}(qnc,4));
                        DfromO(i,1,t,j,c) = mean(Comp{i}(qnc,5));
                        DfromO(i,2,t,j,c) = std(Comp{i}(qnc,5));
                        NC(i,t,j,c) = sum(qc);
                    end
                end
            end
            %%
            % Plot Initial movement angles
            % For each head angle
            figure(1); subplot(2,2,1); hold on
            for i = 1 : 3
                bar(HR(i), abs(nanmean(ima(i,1,:,1,c),3)-90),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('initial movement angle (deg)')    
            figure(1); subplot(2,2,2);  hold on
            for i = 1: 3
                bar(HR(i), nanmean(ima(i,2,:,1,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (deg)')    
            title('Initial movement angle modulation by varying head roll - Rightward movements') 

            figure(1); subplot(2,2,3); hold on
            for i = 1 : 3
                bar(HR(i), abs(nanmean(ima(i,1,:,2,c),3)-90),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('initial movement angle (deg)')    
            figure(1); subplot(2,2,4);  hold on
            for i = 1: 3
                bar(HR(i), nanmean(ima(i,2,:,2,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (deg)')
            title('Initial movement angle modulation by varying head roll - Leftward movements') 

            figname = strcat(sID,cond{c,1},'-IMA');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            % for each obstacle
            figure(12); subplot(2,1,1);hold on
            xa = [250 150 50 -50 -150];
            xad = [5.5 15.5 25.5];
            colorss = [204 255 204; 204 204 204; 255 102 102]./255;
            h = zeros(2,3);
            for i = 1:3
                t1 = reshape(ima(i,1,:,1,c),5,1);
                t2 = reshape(ima(i,1,:,2,c),5,1);
                h(1,i)=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                h(2,i)=bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end  
            legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
                         'HR = 0-R','HR = 0-L',...
                         'HR = 30CW-R','HR = 30CW-L',...
                    })
            legend('Location','northeastoutside')    
            ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''}))
            title('Initial movement angle modulation by varying head roll')

            figure(12); subplot(2,1,2);hold on
            clear h;
            for i = 1:3
                t1 = reshape(ima(i,2,:,1,c),5,1);
                t2 = reshape(ima(i,2,:,2,c),5,1);
                h=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end 
            legend(h,{'HR = 30CCW-R'})
            legend('Location','northeastoutside');legend('hide')    
            xlabel('Obstacle position'); ylabel('Standard Deviations (deg)')
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))

            figname = strcat(sID,cond{c,1},'-IMA-OBs');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            %%
            % Plot curvature
            figure(2); subplot(2,2,1); hold on
            for i = 1:3
                bar(HR(i), nanmean(mcurv(i,1,:,1,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('movement curvature (cm)')
            figure(2); subplot(2,2,2); hold on
            for i = 1: 3
                bar(HR(i), nanmean(mcurv(i,2,:,1,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (cm)')    
            title('Movement curvature modulation by varying head roll - Rightward movements')

            figure(2); subplot(2,2,3); hold on
            for i = 1:3
                bar(HR(i), nanmean(mcurv(i,1,:,2,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('movement curvature (cm)')
            figure(2); subplot(2,2,4); hold on
            for i = 1: 3
                bar(HR(i), nanmean(mcurv(i,2,:,2,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (cm)')
            title('Movement curvature modulation by varying head roll - Leftward movements')

            figname = strcat(sID,cond{c,1},'-MC');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            % for each obstacle
            figure(22); subplot(2,1,1);hold on
            xa = [250 150 50 -50 -150];
            xad = [5.5 15.5 25.5];
            colorss = [204 255 204; 204 204 204; 255 102 102]./255;
            h = zeros(2,3);
            for i = 1:3
                t1 = reshape(mcurv(i,1,:,1,c),5,1);
                t2 = reshape(mcurv(i,1,:,2,c),5,1);
                h(1,i)=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                h(2,i)=bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end  
            legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
                         'HR = 0-R','HR = 0-L',...
                         'HR = 30CW-R','HR = 30CW-L',...
                    })
            legend('Location','northeastoutside')    
            ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''}))
            title('Maximum movement curvature modulation by varying head roll')

            figure(22); subplot(2,1,2);hold on
            clear h;
            for i = 1:3
                t1 = reshape(mcurv(i,2,:,1,c),5,1);
                t2 = reshape(mcurv(i,2,:,2,c),5,1);
                h=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end 
            legend(h,{'HR = 30CCW-R'})
            legend('Location','northeastoutside');legend('hide')    
            xlabel('Obstacle position'); ylabel('Standard Deviations (cm)')
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))

            figname = strcat(sID,cond{c,1},'-MC-OBS');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            %%
            % Plot distance from obstacle
            figure(3); subplot(2,2,1); hold on
            for i = 1: 3
                bar(HR(i), abs(nanmean(DfromO(i,1,:,1,c),3)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Distance from obstacle (cm)')
            figure(3); subplot(2,2,2); hold on
            for i = 1:3
                bar(HR(i), nanmean(DfromO(i,2,:,1,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (cm)')
            title('Distance from obstacle modulation by varying head roll - Rightward movements')

            figure(3); subplot(2,2,3); hold on
            for i = 1: 3
                bar(HR(i), abs(nanmean(DfromO(i,1,:,2,c),3)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Distance from obstacle (cm)')
            figure(3); subplot(2,2,4); hold on
            for i = 1:3
                bar(HR(i), nanmean(DfromO(i,2,:,2,c),3),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('Standard deviation (cm)')
            title('Distance from obstacle modulation by varying head roll - Leftward movements')

            figname = strcat(sID,cond{c,1},'-dfo');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')
        %     
            % for each obstacle
            figure(32); subplot(2,1,1);hold on
            xa = [250 150 50 -50 -150];
            xad = [5.5 15.5 25.5];
            colorss = [204 255 204; 204 204 204; 255 102 102]./255;
            h = zeros(2,3);
            for i = 1:3
                t1 = reshape(DfromO(i,1,:,1,c),5,1);
                t2 = reshape(DfromO(i,1,:,2,c),5,1);
                h(1,i)=bar(xa+xad(i), abs(t1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                h(2,i)=bar(xa-xad(i), abs(t2),'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end  
            legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
                         'HR = 0-R','HR = 0-L',...
                         'HR = 30CW-R','HR = 30CW-L',...
                    })
            legend('Location','northeastoutside')    
            ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''}))
            title('Distance form obstacle modulation by varying head roll')

            figure(32); subplot(2,1,2);hold on
            clear h;
            for i = 1:3
                t1 = reshape(DfromO(i,2,:,1,c),5,1);
                t2 = reshape(DfromO(i,2,:,2,c),5,1);
                h=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end 
            legend(h,{'HR = 30CCW-R'})
            legend('Location','northeastoutside');legend('hide')    
            xlabel('Obstacle position'); ylabel('Standard Deviations (cm)')
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))

            figname = strcat(sID,cond{c,1},'-dfo-OBS');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            %%
            % Plot Number of collisition
            figure(4); subplot(1,2,1); hold on
            for i = 1:3
                bar(HR(i), sum(NC(i,:,1,c)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('obstacle collision - rightward movements')
            figure(4); subplot(1,2,2); hold on
            for i = 1:3
                bar(HR(i), sum(NC(i,:,2,c)),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',8)
            end
            xlabel('Head Rotations (deg)'); ylabel('obstacle collision - leftward movements')
            supertitle('Obstacle collision modulation by varying head roll')

            figname = strcat(sID,cond{c,1},'-NC');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            % for each obstacle
            figure(42); hold on
            xa = [250 150 50 -50 -150];
            xad = [5.5 15.5 25.5];
            colorss = [204 255 204; 204 204 204; 255 102 102]./255;
            h = zeros(2,3);
            for i = 1:3
                t1 = reshape(NC(i,:,1,c),5,1);
                t2 = reshape(NC(i,:,2,c),5,1);
                h(1,i)=bar(xa+xad(i), t1,'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                h(2,i)=bar(xa-xad(i), t2,'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end  
            legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
                         'HR = 0-R','HR = 0-L',...
                         'HR = 30CW-R','HR = 30CW-L',...
                    })
            legend('Location','northeastoutside')    
            ylabel('Nomber of collisions');  xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
            title('Variation of collision number by varying head roll')    


            figname = strcat(sID,cond{c,1},'-NC-OBS');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')

            %%   
            % Plot Number of right vs left -ward movements 
            figure(5); hold on
            xa = [250 150 50 -50 -150];
            xad = [5.5 15.5 25.5];
            colorss = [204 255 204; 204 204 204; 255 102 102]./255;
            rl = zeros(3,5,2);
            h = zeros(2,3);
            for i = 1:3
                for j = 1 : 5
                    qnc = Comp{i}(:,1) == 0 & Comp{i}(:,2) == 0 & Comp{i}(:,8) == j ; % no collision, no repeat
                    rl(i,j,1) = sum(Comp{i}(qnc,7)==0);            
                    rl(i,j,2) = sum(Comp{i}(qnc,7)==1);            
                end 
                h(1,i)=bar(xa+xad(i), rl(i,:,1),'FaceColor',colors{i},'EdgeColor',colors{i},'BarWidth',0.1);
                h(2,i)=bar(xa-xad(i), rl(i,:,2),'FaceColor',colorss(i,:),'EdgeColor',colors{i},'BarWidth',0.1);
            end        
            legend([h(1,1) h(2,1) h(1,2) h(2,2) h(1,3) h(2,3)],{'HR = 30CCW-R','HR = 30CCW-L',...
                         'HR = 0-R','HR = 0-L',...
                         'HR = 30CW-R','HR = 30CW-L',...
                    })
            legend('Location','northeastoutside')
            xlabel('Head Rotations'); ylabel('number of movements')
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
            title('Modulation of movement direction by rolling the head')

            figname = strcat(sID,cond{c,1},'LR-MDir');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')
            close all
    
        end
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