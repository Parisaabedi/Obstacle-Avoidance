%% test the bias + expansion hypothesis
% written by PA, Nov. 2018

%% Generate a trajectory like path
y = -2: 0.002:2;
x1 = -y .^ 2 + 4;
y = y + 2;
x1p = -x1;
figure(1); subplot(1,2,1); plot(x1,y,'k',x1p,y,'k')

%% add the rotational biases
x2 = zeros(length(x1),1);y2 = x2;
x2p = x2;y2p = x2;
d = sqrt(x1.^2+y.^2);
alfa = 10;
for i = 1 : length(x1)
    a = atan2d(y(i),x1(i));
    a = a+alfa;
    y2(i) = sind(a) * d(i);
    x2(i) = cosd(a) * d(i);
    a = atan2d(y(i),x1p(i));
    a = a+alfa;
    y2p(i) = sind(a) * d(i);
    x2p(i) = cosd(a) * d(i);
end
figure(1);subplot(1,2,1); hold on; plot(x2,y2, 'r',x2p,y2p, 'r')

x22 = zeros(length(x1),1);y22 = x22;
x22p = x22;y22p = x22;
for i = 1 : length(x1)
    a = atan2d(y(i),x1(i));
    a = a-alfa;
    y22(i) = sind(a) * d(i);
    x22(i) = cosd(a) * d(i);
    a = atan2d(y(i),x1p(i));
    a = a-alfa;
    y22p(i) = sind(a) * d(i);
    x22p(i) = cosd(a) * d(i);
end
figure(1); subplot(1,2,1); hold on; plot(x22,y22, 'g',x22p,y22p, 'g')
%% add the expantion biases
beta = 1.05;
x3 = beta*x2; % R,HR
x3p = beta*x2p; % L,HR
figure(1); subplot(1,2,1); hold on; plot(x3,y2, ':r',x3p,y2p, ':r')
x33 = beta*x22; % R,HL
x33p = beta*x22p; % L,HL
figure(1); subplot(1,2,1); hold on; plot(x33,y22, ':g',x33p,y22p, ':g')

%% calculate dalta X
p1 = find(y == 1 | y > 1,1);
p2 = find(y2 == 1 | y2 > 1,1);
p2p = find(y2p == 1| y2p > 1,1);
p3 = find(y2 == 1 | y2 > 1,1);
p3p = find(y2p == 1|y2p > 1,1);
p33 = find(y22 == 1|y22 > 1,1);
p33p = find(y22p == 1|y22p > 1,1);
dxrhr = x3(p3) - x1(p1);
dxlhr = x3p(p3p) - x1p(p1);
dxrhl = x33(p33) - x1(p1);
dxlhl = x33p(p33p) - x1p(p1);
%% Calculate alfa
p1 = ceil(size(x1,2)/2);
alfaprime = ((dxrhr-dxrhl)+(dxlhr-dxlhl))/4;

% check if the calculation is correct
dxrhr2 = x2(p3) - x1(p1);
dxlhr2 = x2p(p3p) - x1p(p1);
dxrhl2 = x22(p33) - x1(p1);
dxlhl2 = x22p(p33p) - x1p(p1);
alfapp = ((dxrhr2-dxrhl2)+(dxlhr2-dxlhl2))/4;
%% calculate beta
betaprime = ((dxrhr+dxrhl)-(dxlhr+dxlhl))/4;

% check if the calculation is correct
xb = beta*x1; % R,HR
xbp = beta*x1p;
dxrhr22 = xb(p1) - x1(p1);
dxlhr22 = xbp(p1) - x1p(p1);
dxrhl22 = xb(p1) - x1(p1);
dxlhl22 = xbp(p1) - x1p(p1);
betapp = ((dxrhr22+dxrhl22)-(dxlhr22+dxlhl22))/4;

%% The other way, firt expansion and then rotation
figure(1); subplot(1,2,2); plot(x1,y,'k',x1p,y,'k')
%% add the expantion biases
%beta = 1.2;
xx2 = beta*x1; % R
xx2p = beta*x1p; % L
figure(1); subplot(1,2,2); hold on; plot(xx2,y, ':k',xx2p,y, ':k')

%% add the rotational biases
xx3 = zeros(length(x1),1);yy3 = xx3;
xx3p = xx3;yy3p = xx3;
d = sqrt(xx2.^2+y.^2);
%alfa = 3;
for i = 1 : length(xx2)
    a = atan2d(y(i),xx2(i));
    a = a-alfa;
    yy3(i) = sind(a) * d(i);
    xx3(i) = cosd(a) * d(i);
    a = atan2d(y(i),xx2p(i));
    a = a-alfa;
    yy3p(i) = sind(a) * d(i);
    xx3p(i) = cosd(a) * d(i);
end
figure(1);subplot(1,2,2); hold on; plot(xx3,yy3, 'r',x2p,yy3p, 'r')

xx33 = zeros(length(x1),1);yy33 = xx33;
xx33p = xx33;yy33p = xx33;
for i = 1 : length(x1)
    a = atan2d(y(i),xx2(i));
    a = a+alfa;
    yy33(i) = sind(a) * d(i);
    xx33(i) = cosd(a) * d(i);
    a = atan2d(y(i),xx2p(i));
    a = a+alfa;
    yy33p(i) = sind(a) * d(i);
    xx33p(i) = cosd(a) * d(i);
end
figure(1); subplot(1,2,2); hold on; plot(xx33,yy33, 'g',xx33p,yy33p, 'g')

%% calculate dalta X
ycross = 2;
p1 = find(y == ycross | y > ycross,1);
p2 = p1; p2p = p1; 
p3 = find(yy3 == ycross | yy3 > ycross,1);
p3p = find(yy3p == ycross|yy3p > ycross,1);
p33 = find(yy33 == ycross | yy33 > ycross,1);
p33p = find(yy33p == ycross|yy33p > ycross,1);
dxrhr = xx3(p3) - x1(p1);
dxlhr = xx3p(p3p) - x1p(p1);
dxrhl = xx33(p33) - x1(p1);
dxlhl = xx33p(p33p) - x1p(p1);
%% Calculate alfa
p1 = ceil(size(x1,2)/2);
alfaprime2 = ((dxrhr-dxrhl)+(dxlhr-dxlhl))/4;

% check if the calculation is correct (for this method I can't though)
% dxrhr2 = x2(p3) - x1(p1);
% dxlhr2 = x2p(p3p) - x1p(p1);
% dxrhl2 = x22(p33) - x1(p1);
% dxlhl2 = x22p(p33p) - x1p(p1);
% alfapp2 = ((dxrhr2-dxrhl2)+(dxlhr2-dxlhl2))/4;
%% calculate beta
betaprime2 = ((dxrhr+dxrhl)-(dxlhr+dxlhl))/4;

% check if the calculation is correct
% xb = beta*x1; % R,HR
% xbp = beta*x1p;
% dxrhr22 = xb(p1) - x1(p1);
% dxlhr22 = xbp(p1) - x1p(p1);
% dxrhl22 = xb(p1) - x1(p1);
% dxlhl22 = xbp(p1) - x1p(p1);
% betapp2 = ((dxrhr22+dxrhl22)-(dxlhr22+dxlhl22))/4;

%% Solve the system of non-linear equations
% Calculate the known parameters
acr = atan2d(y(p1),x1(p1));
acl = atan2d(y(p1),x1p(p1));

acrp = atan2d(y(p3),xx2(p3));
aclp = atan2d(y(p3p),xx2p(p3p));
% create the hyperparameter vector
hp = [y(p3) , xx3(p3), xx33(p33), xx3p(p3p),xx33p(p33p),acr,acl];
x0 = [acrp aclp -5 1.2];
k = alfabeta(x0,hp);
val = fsolve(@(x)alfabeta(x,hp),x0);