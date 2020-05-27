function [mtraj,straj] = trajnormalization(trajx,trajy,l,s)
    q = l ~= 0;
    tempx = zeros(s,sum(q));
    tempy = zeros(s,sum(q));
    for i = 1 : sum(q)
        handpath = sqrt((trajx(1:l(i)-1,i)-trajx(2:l(i),i)).^2) + ((trajy(1:l(i)-1,i)-trajy(2:l(i),i)).^2);
        shandpath = sum(handpath);
        normfactor = floor(shandpath / s);
        tempx(1,i) = trajx(1,i); tempy(1,i) = trajy(1,i);
        jp = 1;
        for j = 2 : s-1
            thp = sqrt((trajx(jp,i)-trajx(jp+1:l(i),i)).^2) + ((trajy(jp,i)-trajy(jp+1:l(i),i)).^2);
            tf = find(thp> normfactor | thp == normfactor,1);
            jp = jp + tf;
            tempx(j,i) = trajx(jp,i);  tempy(j,i) = trajy(jp,i);
        end
        tempx(end,i) = trajx(l(i),i); tempy(end,i) = trajy(l(i),i);
    end
    mtraj(:,1) = mean(tempx,2); mtraj(:,2) = mean(tempy,2);
    straj(:,1) = std(tempx,0,2); straj(:,2) = std(tempy,0,2);
   
    
end