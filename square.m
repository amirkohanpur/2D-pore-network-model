% Plot a square with (x,y) center and side length l
function square(x,y,l)
line([x-(l/2),x+(l/2)],[y+(l/2),y+(l/2)]);
line([x-(l/2),x-(l/2)],[y+(l/2),y-(l/2)]);
line([x+(l/2),x+(l/2)],[y-(l/2),y+(l/2)]);
line([x+(l/2),x-(l/2)],[y-(l/2),y-(l/2)]);
end