function [endPoint, borderCount, dirDegree] = createNewPoint(origin, dist, dirDegree, rectSize, Cfg)
% origin is a vector with two coordinates, dist is mesured in pixel per frame, dir is the direction to move in as a cos and sin
    dirRadiant = deg2rad(dirDegree);
    % if size(dist,2) ~= 1
    %     disp ('jj');
    % end 
    cs   = [cos(dirRadiant), sin(dirRadiant)];
    dxdy = cs * dist;  % increment of coordinate that needs to be add to the origin
    
    endPoint = origin + dxdy;
    borderCount = 0;
    
    %AH: note that the origin coord are in absolute, so we deal with coord
    %frame [0 0 14 14] lets suppose, which then at the end outside this
    %function is translated to the real coords.. now they compare with kind
    %of an inside margin which include the dot width by two [w/2 w/2 14-w/2
    %14-w/2].. and then recalls the function for no reason of recursion it
    %seems ?? because the actual intersection point is not returned.. they
    %return directly the next point after intersection, the dist is the
    %speed thay they multiply by the norm direction vestor which they get
    %from the angle.. if it does intersect with the borders, they find this
    %intersection and cosider it as origin and the speed is now the dist
    %remaining outside the rect and they recall again the function..

    %AH: if i create already a control for not having to bounce i should
    %not need this function, i can just use the code above for next point
    
    if endPoint(1) < round(Cfg.dot_w/2,6) || endPoint(2) < round(Cfg.dot_w/2,6) || ...
        endPoint(1) > rectSize(1)- round(Cfg.dot_w/2,6) || endPoint(2) > rectSize(2)-round(Cfg.dot_w/2,6) % if dot is touching rect border
        % if the end point is outside of the rect then 
        % - first will be taken the intersection pointbetween the segment origin-endPoint and the rect
        % - it will be evaluete the length of the segment out side the rect
        % - a new direction is computed depending on which wall it bounced
        % - at the end the function will recall this func again but as origin the intersection point found, and the new direction
        
        %figure();
        xlimit = [Cfg.dot_w/2 rectSize(1)-round(Cfg.dot_w/2,6)];
        ylimit = [Cfg.dot_w/2 rectSize(2)-round(Cfg.dot_w/2,6)];
        xbox = xlimit([1 1 2 2 1]);
        ybox = ylimit([1 2 2 1 1]);
        %mapshow(xbox,ybox,'DisplayType','polygon','LineStyle','none');

        xLine = [origin(1) endPoint(1)];
        yLine = [origin(2) endPoint(2)];
        %mapshow(xLine,yLine,'Marker','+');

        % this function is used to calculate the intersection between the  
        % 'box' (my rectangle by also considering the dimention of the ball)
        [xIntersection, yIntersection] = polyxpoly(xLine, yLine, xbox, ybox);

        % if the speed is really high, then we will have a high value of
        % dist, it is then traslated in the fact that it bounces more then
        % one time. in this case the second time it bounces the 'origin'
        % will be exactly on the edge of the box used for the intersection 
        % with the line, and this will then produced two points of intersection, but we are intrested only in the first one 
        
        if size(xIntersection,1) > 1
            indx = find(round(xIntersection,6) == round(origin(1),6));
            xIntersection(indx) = [];
            yIntersection(indx) = [];
        end

        %mapshow(xIntersection,yIntersection,'DisplayType','point','Marker','o');

        pair = [xIntersection, yIntersection; endPoint];
        
        dist_Intersection_EndPoint = pdist(pair,"euclidean");  % this is the length of the segment that is out side the rect
        %dddddd(zz) = dist_Intersection_EndPoint 
        
        %zz=zz+1;
        %         if size(dist_Intersection_EndPoint) == 0
        %             disp ('ahhhhhhhhhh');
        %         end

        if isempty(dist_Intersection_EndPoint)
            disp("here")
         %dist_Intersection_EndPoint = 0;  % Set a default value if empty
        end
        
        if size(dist,2) ~= 1
            disp ('jj');
        end
        % for all this condition i need to round the numebrs to have 6 digits, because there is a problem with the precesion of numbers
        if round(xIntersection,6) <= round(Cfg.dot_w/2,6)       % left border
            dirDegree = 180-dirDegree;
        elseif round(yIntersection,6) <= round(Cfg.dot_w/2,6)    % top
            dirDegree = 360-dirDegree;
        elseif round(xIntersection,6) >= round(rectSize(1)-Cfg.dot_w/2,6)    % right
            dirDegree = 180-dirDegree;
        elseif round(yIntersection,6) >= round(rectSize(2)-Cfg.dot_w/2,6)    % bottom
            dirDegree = 360-dirDegree;
        end

        
        [endPoint, borderCount, dirDegree] = createNewPoint([xIntersection, yIntersection], dist_Intersection_EndPoint, dirDegree, rectSize, Cfg);
        borderCount = 1;
    end
end
