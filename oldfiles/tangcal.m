%Example
t=0:0.01:10;
y=sin(t);
plot(t,y)
%-------------------------
dy=diff(y)./diff(t);
k=220; % point number 220
tang=(t-t(k))*dy(k)+y(k);
hold on
plot(t,tang)
scatter(t(k),y(k))
hold off