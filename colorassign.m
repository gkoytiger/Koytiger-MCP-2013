function [ colorvec ] = colorassign( aff )
%Greg Koytiger June 2012
aff = aff/1000 * 2;
%Creates Red -> Yellow gradient    
if aff <= 1
    colorvec = [1 aff 0];

%Creates Yellow -> Cyan gradient       

elseif aff < 2
    colorvec = [1-(aff-1) 1 aff-1];

    
end
    
    
end

