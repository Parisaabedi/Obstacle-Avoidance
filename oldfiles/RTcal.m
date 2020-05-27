%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by PAK
% Recalculate the Reaction Times based on the method proposed 
% by Smeets and Brenner, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the Data
load('Comp_s1.mat')
for i = 1 : 6
    n = Ns(i);
    for j = 1 : n
        %% Calculate Velocity
        % read the trajectory data
        handx = ntraj(:,j,1,i);
        handy = ntraj(:,j,1,i);
        % calculate the velocity
        win = 0.005;
        sfr = 200/(params(j,13,i)-params(j,12,i)); % sampling frequency
        handXv = GRAIdiff(1/sfr,win,handx);
        handYv = GRAIdiff(1/sfr,win,handy);
        
        %% Calculate Acceleration    
        % check if you need to smooth your velocity calculation
        handXa = GRAIdiff(1/sfr,win,handXv);
        handYa = GRAIdiff(1/sfr,win,handYv);
        
        %% Calculate RT based on Threshold       
        vel = sqrt(handXv.^2 + handYv .^ 2);
        acc = sqrt(handXa.^2 + handYa .^ 2);

        %% Calculate RT based on Extrapolation
    end
end