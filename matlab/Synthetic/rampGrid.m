function p = rampGrid(x,y,a,b,c,resolution)
p = x-x; % create a Zero vector with the same resolution than x
for j=0:resolution
    dy = -0.5 + j/resolution;
    for i=0:resolution
        dx = -0.5 + i/resolution;
        p = p + (a*(x+dx) + b*(y+dy) + c > 0);
    end
end
p = p /(resolution+1) / (resolution+1);
end