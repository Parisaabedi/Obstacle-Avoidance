%% Analysis of rotation and expansion biases
% Written by PA, Nov. 2018

%% Match the relevant obstacles together
% mirrored obstacles are matched to have both rightward and leftward MVs (3, 1&5, 2&4)

% Initialization
alfabeta = zeros(3,2,18); %  3 pairs of obstacles, 2 parameters, #subjects
y = 19; % point of interest along y-axis
c = 1; % with/without visual feedback (0=> with, 1=> without)
valid = 1; % check if there is enough data on both right and left direction
obs = [3 3 ; 5 1 ; 4 2];
% Load the data
load('WholeComp3.mat')
j = 1;
for i = 1 : 18
    for j = 1:3 
        if j == 1
            ec = zeros(3,2);
            for t = 1 : 3
                rrd = 100*(reshape(RD(t+c*3,3,:,i),[],1)./sum(reshape(RD(t+c*3,3,:,i),[],1)));
                ec(t,:) = rrd > 10;            
            end
            if sum(ec(:)==0) > 0
                valid = 0; 
            else
                valid = 1;
            end
        else
            valid = 1;
        end
        if valid == 1
            % find pointers 
            p1 = find(trajc(:,2,1+c*3,1,obs(j,1),1,i) == y | trajc(:,2,1+c*3,1,obs(j,1),1,i) > y,1); 
            p1p = find(trajc(:,2,1+c*3,1,obs(j,2),2,i) == y | trajc(:,2,1+c*3,1,obs(j,2),2,i) > y,1); 
            p3 = find(trajc(:,2,3+c*3,1,obs(j,1),1,i) == y | trajc(:,2,3+c*3,1,obs(j,1),1,i) > y,1);
            p3p = find(trajc(:,2,3+c*3,1,obs(j,2),2,i) == y | trajc(:,2,3+c*3,1,obs(j,2),2,i) > y,1);
            p33 = find(trajc(:,2,2+c*3,1,obs(j,1),1,i) == y | trajc(:,2,2+c*3,1,obs(j,1),1,i) > y,1);
            p33p = find(trajc(:,2,2+c*3,1,obs(j,2),2,i) == y | trajc(:,2,2+c*3,1,obs(j,2),2,i) > y,1);

            % calculate deltaX for different directions and conditions
            dxrhr = trajc(p3,1,3+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
            dxlhr = trajc(p3p,1,3+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);
            dxrhl = trajc(p3,1,2+c*3,1,obs(j,1),1,i) - trajc(p1,1,1+c*3,1,obs(j,1),1,i);
            dxlhl = trajc(p3p,1,2+c*3,1,obs(j,2),2,i) - trajc(p1p,1,1+c*3,1,obs(j,2),2,i);

            % calculate alfa and beta
            alfabeta(j,1,i) = ((dxrhr-dxrhl)+(dxlhr-dxlhl))/4; % alfa
            alfabeta(j,2,i) = ((dxrhr+dxrhl)-(dxlhr+dxlhl))/4; % beta
        end
    end    
end

%% Plot the result
ddata = reshape(alfabeta(2,:,:),2,18);
figure; hold on; scatter(ddata(1,:),ddata(2,:))

%% stat analysis
% check if the expansion is consistently greater than zero for all participants
q = reshape(alfabeta(1,2,:),1,18);
qp = q == 0;
[h1,p1] = ttest(q(~qp)); % obstacle 3
[h2,p2] = ttest(reshape(alfabeta(2,2,:),1,18)); % obstacle 1&5
[h3,p3] = ttest(reshape(alfabeta(3,2,:),1,18)); % obstacle 2&4

