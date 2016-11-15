function theta=qtan(y1,x1,y2,x2)
    
theta=atan((y2-y1)/(x2-x1));
    
%90 degrees
if (x2==x1)&&(y2>y1);
    theta = pi/2;
   
%180 degrees
elseif (x2<x1)&&(y2==y1);
    theta = pi;
    
%270 degrees
elseif (x2==x1)&&(y2<y1);
    theta =-pi/2;  
    
%Zero
elseif (x2>x1)&&(y2==y1);
    theta =0;
    
    %First Quadrant
elseif (x2>x1)&&(y2>y1);
    theta =abs(theta);
    
%Second Quadrant
elseif (x2<x1)&&(y2>y1);
    theta = pi - abs(theta);
    
%Third Quadrant
elseif (x2<x1)&&(y2<y1);
    theta = -pi + abs(theta);
    
%Fourth Quadrant
elseif (x2>x1)&&(y2<y1);
    theta = -abs(theta); 
    
end

end