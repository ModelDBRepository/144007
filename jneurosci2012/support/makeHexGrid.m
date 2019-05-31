function grid = makeHexGrid(width, height, spacing)
% function grid = makeHexGrid(width, height, spacing)
% This function returns the centers of a regularly spaced hexagonal grid
% spacing is the distance between points in the grid

nPoints = width*height;
grid = zeros(nPoints,2);
for i = 1:nPoints
    col = mod(i-1,width);
    row = floor((i-1)/width);
    if(mod(row,2))
        xoffset = spacing/2;
    else
        xoffset = 0;
    end
    
    grid(i, 1) = col*spacing + xoffset;
    grid(i, 2) = spacing/2*sqrt(3)*row;
end
%plot(grid(:,1), grid(:,2), 'o');