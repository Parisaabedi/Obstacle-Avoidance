%% Non-linear functions for seperating rotatoinal and expansion biases
% x => 1 alfacpR, 2 alfacpL, 3 alfa, 4 beta
% HP => 1 Yc, 2 dxrhr, 3 dxrhl, 4 dxlhr, 5 dxlhl, 6 alfacr, 7 alfaxl
% written by PA, Nov. 2018

%%
function F = alfabeta(x,hp)
    
    F(1) = tand(x(1)) - hp(1)/(x(4)*(hp(2)+hp(3))*cosd(hp(6)));
    F(2) = hp(1)*2*tand(x(1)) - (hp(2)+hp(3))*(tand(x(1))^2+tan(x(3))^2);
    F(3) = hp(1)*2*tand(x(2)) - (hp(4)+hp(5))*(tand(x(2))^2+tan(x(3))^2);
    F(4) = tand(x(2)) - hp(1)/(x(4)*(hp(4)+hp(5))*cosd(hp(7)));
end

