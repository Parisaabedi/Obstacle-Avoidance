plotting = 1;
colors = { 'g' ; 'k' ; 'r'};
if plotting == 1
    subs = 1:7; subs = [subs 9];
    cond = {'F';'NF'};
    for s = 5:8
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
            for i = 1:3
                for j = 1:2 
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
            % for each obstacle
            figure(12); subplot(2,1,1);hold on  
            xa = 5:-1:1;            
            h = zeros(1,3);
            for i = 1:3
                t1 = flipud(reshape(ima(i,1,:,1,c),5,1));
                t2 = flipud(reshape(ima(i,2,:,1,c),5,1));
                h(1,i)=errorbar(1:5,abs(t1-90),t2,'Color',colors{i});
                
            end  
            legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
            legend('Location','southeast')    
            ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''}))
            xlim([0,6])
            title('Initial movement angle modulation by varying head roll - Rightward MVs')

            figure(12); subplot(2,1,2);hold on
            for i = 1:3
                t1 = flipud(reshape(ima(i,1,:,2,c),5,1));
                t2 = flipud(reshape(ima(i,2,:,2,c),5,1));
                errorbar(1:5,abs(t1-90),t2,'Color',colors{i});
                
            end              
            xlabel('Obstacle position'); ylabel('Initial movement angle (deg)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            xlim([0,6])
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))
            title('Initial movement angle modulation by varying head roll - Leftward MVs')

            figname = strcat(sID,cond{c,1},'-IMA-OBs'); 
            savefig(gcf,figname) 
            saveas(gcf,figname,'tiff')

            %%
            % Plot curvature
            % for each obstacle
            figure(22); subplot(2,1,1);hold on 
            xa = 5:-1:1;
            h = zeros(1,3);
            for i = 1:3
                t1 = flipud(reshape(mcurv(i,1,:,1,c),5,1));
                t2 = flipud(reshape(mcurv(i,2,:,1,c),5,1));
                h(1,i)=errorbar(1:5,t1,t2,'Color',colors{i});                
            end  
            legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
            legend('Location','southeast')           
            ylabel('Maximum curvature (cm)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            xlim([0,6]);
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''}))
            title('Maximum movement curvature modulation by varying head roll - Rightward MVs')

            figure(22); subplot(2,1,2);hold on
            for i = 1:3
                t1 = flipud(reshape(mcurv(i,1,:,2,c),5,1));
                t2 = flipud(reshape(mcurv(i,2,:,2,c),5,1));
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
            % Plot distance from obstacle            
            % for each obstacle
            figure(32); subplot(2,1,1);hold on       
            xa = 5:-1:1;
            h = zeros(1,3);
            for i = 1:3
                t1 = flipud(reshape(DfromO(i,1,:,1,c),5,1));
                t2 = flipud(reshape(DfromO(i,2,:,1,c),5,1));
                h(1,i)=errorbar(1:5,abs(t1),t2,'Color',colors{i});                
            end
            legend([h(1,1) h(1,2) h(1,3)],{'HR = 30CCW','HR = 0','HR = 30CW'})
            legend('Location','southeast')           
            ylabel('Distance from Obstacle (cm)'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            xlim([0,6]);
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'','','','',''})) 
            title('Distance form obstacle modulation by varying head roll - Rightawrd MVs')

            figure(32); subplot(2,1,2);hold on
            for i = 1:3
                t1 = flipud(reshape(DfromO(i,1,:,2,c),5,1));
                t2 = flipud(reshape(DfromO(i,2,:,2,c),5,1));
                errorbar(1:5,abs(t1),t2,'Color',colors{i});                
            end
            xlabel('Obstacle position');ylabel('Distance from Obstacle'); % xlabel('Obstacle position');
            plot(xa,[nan nan nan nan nan])
            xlim([1,6]);
            set(gca,'XTick',fliplr(xa))
            set(gca, 'XTickLabels',fliplr({'Most Rightward','Right','Center','Left','Most Leftward'}))            
            title('Distance form obstacle modulation by varying head roll - Leftward MVs')            

            figname = strcat(sID,cond{c,1},'-dfo-OBS');
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

            figname = strcat(sID,cond{c,1},'-NC');
            savefig(gcf,figname)
            saveas(gcf,figname,'tiff')
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
end