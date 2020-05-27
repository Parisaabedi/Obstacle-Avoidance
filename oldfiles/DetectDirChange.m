%% This function detects the changes in the trajectory direction
files = dir;
HR = [-30 0 30];
cond1 = { 'FL'; 'F0'; 'FR'};
cond2 = {'NFL' ; 'NF0' ; 'NFR'};
colors = { 'g' ; 'k' ; 'r'};
colorss = [204 255 204; 204 204 204; 255 102 102]./255;
trajectories = 1;
trajx = 1000*ones(1000,300,3);
trajy = 1000*ones(1000,300,3);
for j = 1 : length(files)
    k = strfind(files(j).name,'_marked');
    color = [ '-g'; '-k'; '-r'];
%         figure(10); hold on; % row trajectories        
    for i = 1 : 3
        %figure(50); hold on
        % load the data        
        if ~isempty(k)
            k1 = strfind(files(j).name,['9' cond2{i}]);
            if ~isempty(k1)   
                load(files(j).name)
                load('CompNF_s9.mat')
                tti = 1;
                for tt = 1 : N
                    if param(tt,24) == 0 && param(tt,25) == 0   % select good trials                            
                        on = find(D.t{tt,1}(:) > D.handON{1,tt},1);
                        off1 = find(D.t{tt,1}(:) < D.handOFF{1,tt});
                        off = off1(end);                       
                        trajx(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:off,1);
                        trajy(1:off-on+1,tti,i) = D.handscreen{1,tt}(on:off,2);
%                         figure(5*i);%hold on
%                         plot(D.handscreen{1,tt}(on:off,1),D.handscreen{1,tt}(on:off,2),'Color',colorss(i,:))
%                         figure(6*i); subplot(2,1,1)
%                         plot(1:off-on+1,D.handXv{tt}(on:off,1),1:off-on+1,D.handZv{tt}(on:off,1))   
%                         figure(6*i); subplot(2,1,2)
%                         plot(1:off-on+1,D.handXa{tt}(on:off,1),1:off-on+1,D.handZa{tt}(on:off,1))
                        dtrajn = diff(find(D.handXv{tt}(on:off,1)<0));
                        dtrajp = diff(find(D.handXv{tt}(on:off,1)>0));
                        if ~isempty(find(dtrajn>1, 1)) || ~isempty(find(dtrajp>1, 1))
                            if ~isempty(find(dtrajn>1, 1))
                                p1 = find(D.handXv{tt}(on:off,1)<0);
                                p2 = find(dtrajn>1, 1);
                                pointer = p1(p2);
                                dt = D.t{tt,1}(on+pointer) - D.handON{1,tt};
                                if dt > 0.1 && find(dtrajn>1,1) < length(dtrajn)-5
                                    Comp{i,1}(tt,9) = 1;
                                else
                                    D.handON{1,tt} = D.t{tt,1}(on+find(dtrajn>1, 1));% check later if you really want to change the onset time! 
                                end 
                            end
                            if ~isempty(find(dtrajp>1, 1))
                                p1 = find(D.handXv{tt}(on:off,1)>0);
                                p2 = find(dtrajp>1, 1);
                                pointer = p1(p2);
                                dt = D.t{tt,1}(on+pointer) - D.handON{1,tt};
                                if dt > 0.1 && find(dtrajp>1,1) < length(dtrajp)-5
                                    Comp{i,1}(tt,9) = 1;
                                else
                                    D.handON{1,tt} = D.t{tt,1}(on+find(dtrajp>1, 1));% check later if you really want to change the onset time! 
                                end
                            end
                            
                        end                        
                        c = 2;
                    end
                end  
                save('CompNF_s9.mat','Comp')
            end            
        end
        
    end        
end