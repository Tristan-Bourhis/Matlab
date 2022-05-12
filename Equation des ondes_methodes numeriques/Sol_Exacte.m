function [V] = Sol_Exacte(x,t)

c=2;
L=1;
V = cos((pi/L).*c.*t).*sin((pi/L).*x)+(1/4).*cos(((10.*pi)/L).*c.*t).*sin(((10.*pi)/L).*x);

end