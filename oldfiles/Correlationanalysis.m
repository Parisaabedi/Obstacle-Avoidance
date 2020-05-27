%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by PAK (Dec. 2018)
% plotting the correlation between several parameters:
% Max curvature and Reaction time
% Movement time and peak velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcrtcor = zeros(6,5,2,18); % correlation (cond, obstacle, direction,subs)
for s = 1:18
    filename = strcat('Comp3_s',num2str(s),'.mat');
    load(filename)    
    for cond = 1:6
        temp = params(1:Ns(cond),:,cond);       
        for j = 1 : 2
            for t = 1 : 5
                figure(s);subplot(6,5,t+(cond-1)*5); hold on
                qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0; % repeats are not allowed
                rts = temp(qnc,12);
                mcs = temp(qnc,19);
                scatter(rts,(-2*(j-1)+1)*mcs)
                covtemp = cov(rts,mcs);
                mcrtcor(cond,t,j,s) = covtemp(1,2);
            end
        end
    end
    figname = strcat(['RTMCcorr-s# ',num2str(s)]); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')  
end
close all

mcmtcor = zeros(6,5,2,18); % correlation (cond, obstacle, direction,subs)
for s = 1:18
    filename = strcat('Comp3_s',num2str(s),'.mat');
    load(filename)    
    for cond = 1:6
        temp = params(1:Ns(cond),:,cond);       
        for j = 1 : 2
            for t = 1 : 5
                figure(s);subplot(6,5,t+(cond-1)*5); hold on
                qnc = temp(:,4) == 0 & temp(:,5) == 0 & temp(:,6) == 0 &...
                      temp(:,7) == (j-1) & temp(:,8) == t & temp(:,24) == 0; % repeats are not allowed
                mts = temp(qnc,13)-temp(qnc,12);
                mcs = temp(qnc,19);
                scatter(mts,(-2*(j-1)+1)*mcs)
                covtemp = cov(mts,mcs);
                mcmtcor(cond,t,j,s) = covtemp(1,2);
            end
        end
    end
    figname = strcat(['MTMCcorr-s# ',num2str(s)]); 
    savefig(gcf,figname) 
    saveas(gcf,figname,'tiff')  
end
close all