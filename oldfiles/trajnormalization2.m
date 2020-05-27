function [mtraj,straj] = trajnormalization2(trajx,trajy,l,s)
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
            i
        else
            temp1 = trajy(1:l(i),i);
            temp2 = trajx(1:l(i),i);
        end
        tempx(:,i) = interp1(temp1,temp2,yq','spline');            
        tempy(:,i) = yq';
        
    end
    mtraj(:,1) = mean(tempx,2); mtraj(:,2) = mean(tempy,2);
    straj(:,1) = std(tempx,0,2); straj(:,2) = std(tempy,0,2);
   
    
end