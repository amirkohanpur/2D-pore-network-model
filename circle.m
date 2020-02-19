% Plot a circle with (x,y) center and diameter d
function circle(x,y,d)
t=0:0.01:2*pi; 
xp=0.5*d*cos(t);
yp=0.5*d*sin(t);
plot(x+xp,y+yp,'black');
end