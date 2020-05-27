% Correct the direction problem!

%% 
for i = 1 : 18
    % load the previous data
    filename = strcat('comp_s',num2str(i),'.mat');
    load(filename)
    
    % copy the direction
    temp = zeros(300,6);
    for j = 1 : 6
        temp(:,j) = params(:,7,j);
    end
    % load the new data
    filename = strcat('comp2_s',num2str(i),'.mat');
    load(filename)
    
    % paste the new data
    for j = 1 : 6
        params(:,7,j) = temp(:,j);
    end
    
    save(filename,'params','Ns','ntraj','nvelacc','straj')
end