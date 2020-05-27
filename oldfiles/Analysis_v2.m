%% Data Analysis (Trajectory analysis)
% Trajectory data analysis for Obstacle Avoidance Project
% Written by PAK (August2018)

%% Load the file 
[filename,pathname] = uigetfile('*.mat','Please choose the DataSet you want to load...');
    cd(pathname);
    % trn = 1;
    FILES = filename;
    leng = length(filename);
    
    if exist([filename(1:leng-4) '.mat'],'file') == 2         
        if contains(filename,'marked')
            load([filename(1:leng-4) '.mat']);
        end
    end
%%  Intropolation (increase sampling rate to 120Hz)
figure(1); hold on
for i = 1 : D.trialnum
    tON = D.eventtime{i}(2); % time of target presentation
    %ton = find(D.t{i}(:) >= tON, 1 ); 
    tEM = D.handOFF{i}; % End of movement (targethit/collision)
    %tem = find(D.t{i}(:) >= tEM, 1 ); 
    Rarr = find(D.t{i} >= tON, 1 ):find(D.t{i} <= tEM, 1, 'last' );
    sr = length(Rarr) / (tEM - tON);
    x = nan(length(Rarr)*2-1,2);
    x(1:2:end,:)= D.handscreen{1,i}(Rarr(1):Rarr(end),:);
    % intropolate to find the missing data points
    xint = inpaint_nans(x,2);
    figure(1); plot(xint(:,1),xint(:,2))
end

%% Filtering, ON-OFF calculation, Velocity and Acceleration calculation


%% Standardization (removing a baseline from trajectories)

%% Normalization

%% Calculate different parameters
