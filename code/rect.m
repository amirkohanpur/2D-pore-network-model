% Plot a rectangle
% (x,y) is center
% lx is horizontal side and ly is vertical side
function rect(x,y,lx,ly)
line([x-(lx/2),x+(lx/2)],[y+(ly/2),y+(ly/2)]);
line([x-(lx/2),x-(lx/2)],[y+(ly/2),y-(ly/2)]);
line([x+(lx/2),x+(lx/2)],[y-(ly/2),y+(ly/2)]);
line([x+(lx/2),x-(lx/2)],[y-(ly/2),y-(ly/2)]);
end