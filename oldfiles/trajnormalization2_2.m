% This function calculate the averaged trajectory by first normalizing the
% tranjectory data
% first intropolating the missing points and then normalizing
% Written by PA, May 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mtraj,straj,scheck] = trajnormalization2_2(trajx,trajy,l,o,td,tn,s)
    os = [15.275 17.0750 18.875 20.6750 22.4750]; % Obstacle positions
    os = fliplr(os);
    scheck = 0;
    q = l ~= 0;
    tempx = zeros(s,sum(q));
    tempy = zeros(s,sum(q));
    for i = 1 : sum(q)
        yq = 9.5:(29-9.5)/(s-1):29;   
        tt = diff(trajy(1:l(i),i))./diff(trajx(1:l(i),i));
        q = isnan(tt)| tt == inf ;
        if sum(q)
            t = find(q==1);
            t2 = trajy(t,i) < 10;
            if sum(t2) 
                t3 = find(t2 == 1,1,'last');
                trajy(1:t(t3)-1,i) = nan;
                trajx(1:t(t3)-1,i) = nan;
            end
            trajy(q,i) = nan;
            trajx(q,i) = nan;
            
            %l(i) = l(i) - sum(q);
            t1 = trajy(1:l(i),i);
            t2 = trajx(1:l(i),i);
            temp1 = t1(~isnan(t1(:)));
            temp2 = t2(~isnan(t2(:)));
            %i
        else
            temp1 = trajy(1:l(i),i);
            temp2 = trajx(1:l(i),i);
        end
        if temp1(1)-yq(1) > .2
            temp1 = [yq(1); temp1];
            temp2 = [18.875 ; temp2];
        end
        if ~issorted(temp1)
            %pause
            temp = diff(temp1);
            q = temp <= 0;
            temp11 = temp1(~q);
            temp22 = temp2(~q);
            %clear temp1; clear temp2;
            temp1 = temp11; temp2 = temp22;
        end
        %temp2 = temp2;
        tempx(:,i) = interp1(temp1,temp2,yq','spline','extrap');            
        tempy(:,i) = yq';
        if max(tempx(:,i))> 40 || abs(min(tempx(:,i))) > 40
            qt = find((tempx(:,i))> 35 | abs((tempx(:,i))) > 35);
            if qt(1) == 1
                qt(1) = 2;
                tempx(1,i) = 18.875;
            end
            if length(qt)<5 && length(qt)>1
                tempx(qt,i) = ones(length(qt),1).* mean([tempx(qt(1)-1,i);tempx(qt(end)+1,i)]);
            elseif length(qt) == 1 && qt ~= 2
                tempx(qt,i) = tempx(qt-1,i)+(tempx(qt-1,i) - tempx(qt-2,i));
            else
                tempx(qt,i) = nan;
            end
            tn(i)
            pause
        end
        tt = find(tempy(:,i) > 19,1);
        if tempx(tt,i) > os(o(i))
            if td(i) == 1
                scheck = scheck + 1;
                tn(i)
                pause
            end
            td(i) = 0;
        else
            if td(i) == 0 
                scheck = scheck + 1;
                tn(i)
                pause
            end
            td(i) = 1;
        end
        
    end
    mtraj = zeros(length(yq),5,2); straj = mtraj;
    oprime = o(o~=0);
    for i = 1:5
        q = oprime' == i;
        %tt = tempx(:,q); 
        for j = 1 : 2
            tt = tempx(:,q); 
            q2 = td(q) == j-1; % right side or left side
            temp = nanmedian(tt(:,q2),2);
            temp = repmat(temp,[1,sum(q)])-tt(:,:);
            tt(abs(temp(:,1:sum(q)))>6) = nan;
            mtraj(:,i,j) = nanmean(tt(:,q2),2); %mtraj(:,2) = mean(tempy,2);
            straj(:,i,j) = nanstd(tt(:,q2),0,2); %straj(:,2) = std(tempy,0,2);
        end
    end
    
%     mtraj(:,1) = mean(tempx,2); mtraj(:,2) = mean(tempy,2);
%     straj(:,1) = std(tempx,0,2); straj(:,2) = std(tempy,0,2);
   
    
end