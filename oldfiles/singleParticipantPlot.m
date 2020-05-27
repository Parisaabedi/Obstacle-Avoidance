%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by PAK, Jan. 2019
% Single participants plots for the result section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
% participant #2
color = { 'k'; 'g'; 'r'};
load('Comp3_s2.mat')
cs2 = zeros(3,2);
om2 = zeros(3,2);
for i = 4:6
    temp = params(1:Ns(i),:,i);
    ttraj = ntraj(:,1:Ns(i),:,i);  
    tvel = nvelacc(:,1:Ns(i),1,i);
    for j = 1:2
        t = 3;        
        qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                  temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0;
        pq = find(qnc == 1);
        om2(i-3,j) =  sum(temp(:,7) == (j-1) & temp(:,8) == t);
        for k = 1 : sum(qnc)
            figure(5); hold on; % trajectory plots
            plot(ttraj(:,pq(k),1)-18.8750,ttraj(:,pq(k),2),...
                'Linewidth',0.5, 'LineStyle',':','Color',color{i-3})
            figure(j); hold on; % Velocity Profile
            tsr = (temp(pq(k),13)-temp(pq(k),12))/2000;
            ts = 0 : tsr : (temp(pq(k),13)-temp(pq(k),12));
            plot(ts,tvel(:,pq(k)),'Linewidth',0.5, 'LineStyle',':','Color',color{i-3})            
            
        end
        figure(5); hold on; 
        plot(nanmean(ttraj(:,qnc,1),2)-18.8750,nanmean(ttraj(:,qnc,2),2),...
                'Linewidth',2.5, 'LineStyle','-','Color',color{i-3})
        scatter([18.8750-18.8750 0 0],[30-11 9 29],'ro')
        figure(j); hold on; % Velocity Profile
        tsr = mean((temp(qnc,13)-temp(qnc,12)))/2000;
        ts = 0 : tsr : mean((temp(qnc,13)-temp(qnc,12)));
        plot(ts,nanmean(tvel(:,qnc),2),'Linewidth',2.5, 'LineStyle','-','Color',color{i-3})
        cs2(i-3,j) = sum( temp(:,5) == 1 & temp(:,7) == (j-1) & temp(:,8) == t); % collision number
    end   
    
end

% participant # 16
load('Comp3_s16.mat')
cs16 = zeros(3,2);
om16 = zeros(3,2);
for i = 4:6
    temp = params(1:Ns(i),:,i);
    ttraj = ntraj(:,1:Ns(i),:,i);  
    tvel = nvelacc(:,1:Ns(i),1,i);
    for j = 1:2
        t = 3;        
        qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                  temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0;
        
        pq = find(qnc == 1);
        om16(i-3,j) = sum(temp(:,7) == (j-1) & temp(:,8) == t);
        for k = 1 : sum(qnc)
            figure(5); hold on; % trajectory plots
            plot(ttraj(:,pq(k),1)-18.8750,ttraj(:,pq(k),2),...
                'Linewidth',0.5, 'LineStyle',':','Color',color{i-3})
            figure(j); hold on; % Velocity Profile
            tsr = (temp(pq(k),13)-temp(pq(k),12))/2000;
            ts = 0 : tsr : (temp(pq(k),13)-temp(pq(k),12));
            plot(ts,tvel(:,pq(k)),'Linewidth',0.5, 'LineStyle',':','Color',color{i-3})            
            
        end
        figure(5); hold on; 
        plot(nanmean(ttraj(:,qnc,1),2)-18.8750,nanmean(ttraj(:,qnc,2),2),...
                'Linewidth',2.5, 'LineStyle','-','Color',color{i-3})
        scatter([18.8750-18.8750 0 0],[30-11 9 29],'ro')
        figure(j); hold on; % Velocity Profile
        tsr = mean((temp(qnc,13)-temp(qnc,12)))/2000;
        ts = 0 : tsr : mean((temp(qnc,13)-temp(qnc,12)));
        plot(ts,nanmean(tvel(:,qnc),2),'Linewidth',2.5, 'LineStyle','-','Color',color{i-3})
        cs16(i-3,j) = sum( temp(:,5) == 1 & temp(:,7) == (j-1) & temp(:,8) == t); % collision number
    end
    
end