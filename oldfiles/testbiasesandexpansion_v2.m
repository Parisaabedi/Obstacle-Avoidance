%% test the bias + expansion hypothesis
% written by PA, Nov. 2018

%% Generate a trajectory like path
y = 0:0.005:20;
x1 = ((y-20).*y)./25;
x1 = -x1;
x1p = -x1;
figure(1); subplot(1,2,1); plot(x1,y,'k',x1p,y,'k')

%% add the rotational biases
x2 = zeros(length(x1),1);y2 = x2;
x2p = x2;y2p = x2;
% d = sqrt(x1.^2+y.^2);
alfa = 5;
R = [cosd(alfa) -sind(alfa);sind(alfa) cosd(alfa)];
R2 = [cosd(-alfa) -sind(-alfa);sind(-alfa) cosd(-alfa)];
for i = 1 : length(x1)
    v = [x1(i) ;y(i)];
    vr = R * v;
    y2(i) = vr(2);
    x2(i) = vr(1);
    v = [x1p(i) ;y(i)];
    vr = R * v;
    y2p(i) = vr(2); 
    x2p(i) = vr(1);
end
% x2 = ((y2-20).*y2)./25;
% x2p = ((y2p-20).*y2p)./25;
figure(1);subplot(1,2,1); hold on; plot(x2,y2, 'r',x2p,y2p, 'r')

x22 = zeros(length(x1),1);y22 = x22;
x22p = x22;y22p = x22;
for i = 1 : length(x1)
    v = [x1(i) ;y(i)];
    vr = R2 * v;
    y22(i) = vr(2);
    x22(i) = vr(1);
    v = [x1p(i) ;y(i)];
    vr = R2 * v;
    y22p(i) = vr(2); 
    x22p(i) = vr(1);
end
figure(1); subplot(1,2,1); hold on; plot(x22,y22, 'g',x22p,y22p, 'g')
%% add the expantion biases
beta = 1.5;
x3 = beta*x2; % R,HR
x3p = beta*x2p; % L,HR
figure(1); subplot(1,2,1); hold on; plot(x3,y2, ':r',x3p,y2p, ':r')
x33 = beta*x22; % R,HL
x33p = beta*x22p; % L,HL
figure(1); subplot(1,2,1); hold on; plot(x33,y22, ':g',x33p,y22p, ':g')

%% calculate dalta X
ybase = 10;
p1 = find(y == ybase | y > ybase,1);
p2 = find(y2 == ybase | y2 > ybase,1);
p2p = find(y2p == ybase| y2p > ybase,1);
p3 = find(y2 == ybase | y2 > ybase,1);
p3p = find(y2p == ybase|y2p > ybase,1);
p33 = find(y22 == ybase|y22 > ybase,1);
p33p = find(y22p == ybase|y22p > ybase,1);
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
dxrhr22 = xb(p3) - x1(p1);
dxlhr22 = xbp(p3p) - x1p(p1);
dxrhl22 = xb(p33) - x1(p1);
dxlhl22 = xbp(p33p) - x1p(p1);
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
    v = [xx2(i);y(i)];
    vr = R * v;
    yy3(i) = vr(2);
    xx3(i) = vr(1);
    v = [xx2p(i);y(i)];
    vr = R * v;
    yy3p(i) = vr(2);
    xx3p(i) = vr(1);
end
figure(1);subplot(1,2,2); hold on; plot(xx3,yy3, 'r',xx3p,yy3p, 'r')

xx33 = zeros(length(x1),1);yy33 = xx33;
xx33p = xx33;yy33p = xx33;
for i = 1 : length(x1)
    v = [xx2(i);y(i)];
    vr = R2 * v;
    yy33(i) = vr(2);
    xx33(i) = vr(1);
    v = [xx2p(i);y(i)];
    vr = R2 * v;
    yy33p(i) = vr(2);
    xx33p(i) = vr(1);
end
figure(1); subplot(1,2,2); hold on; plot(xx33,yy33, 'g',xx33p,yy33p, 'g')

%% calculate dalta X
ycross = 10;
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