function P = sliceDelaunay(d, varargin)
% SLICEDELAUNAY returns the cross-section points of a plane for a delaunay
% triangulation structure.
%
% SYNTAX
%   P = SLICEDELAUNAY(d, planeName, planePosition) returns points (P) of the
%   cross-sectional area that crosses the plane named planeName (possible values = 'x', 'y' or 'z')
%   at a specific position (double value)
%   P = SLICEDELAUNAY(d, normalOfPlane, pointOnPlane) returns points (P) of the
%   cross-sectional area that crosses the plane described by its normal
%   (3x1 or 1x3 vector) and any point (3x1 or 1x3 vector) that lies on this plane.
% 
% EXAMPLE
%   The following example cut a sphere using SLICEDELAUNAY function
%     % Create a sphere
%     [x,y,z] = sphere(31);
%     d = delaunayTriangulation(x(:), y(:), z(:));
%     % Lets plot it
%     trisurf(d.convexHull, d.Points(:,1), d.Points(:,2), d.Points(:,3)); 
%     axis equal; hold on;
%     
%     % Cut a slice on yz plane at x=0.25
%     px = sliceDelaunay(d, 'x', 0.25);
%     plot3(px(1,:), px(2,:), px(3,:), 'r.', 'markersize', 25)
% 
%     % Cut a slice on plane with a custom equation (normal = [1 1 1], V0 = [0 0 0.25])
%     px = sliceDelaunay(d, [1 1 1], [0 0 .25]);
%     plot3(px(1,:), px(2,:), px(3,:), 'b.', 'markersize', 25)
%
% NOTE
%   If multitoolbox (available online at 
%   http://www.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled) 
%   is found, SLICEDELAUNAY uses it. If it is not found, it uses for loop
%   which greatly reduces efficiency
%
% DISCLAIMER and COPYRIGHT
%   SLICEDELAUNAY uses an home-updated version of plane_line_intersect available
%   at http://www.mathworks.com/matlabcentral/fileexchange/17751-straight-line-and-plane-intersection 
%   
%   This function is written by Pariterre
%   Version 1 : April 13th, 2016
    
    % Make sure type is okay
    if ~isa(d, 'delaunayTriangulation')
        error('First argument must be of delaunayTriangulation type');
    end
    
    % Get the triangles of delaunay structure
    t = d.convexHull;
    p = d.Points;
    
    % Find all possible lines of the mesh (between each triangle points 1,
    % 2 et 3), starting points are on 1st, 3rd and 5th columns while end
    % points are on 2nd, 4th and 6th (This will be usefull later)
    t_tp = t(:,[1 2 1 3 2 3]); 
    % Find what is the cutting plane
    if ischar(varargin{1})
        useStandardizedPlane = true;
        if strcmp(varargin{1}, 'x')
            plane = [1 0 0; varargin{2} 0 0];
        elseif strcmp(varargin{1}, 'y')
            plane = [0 1 0; 0 varargin{2} 0];
        elseif strcmp(varargin{1}, 'z')
            plane = [0 0 1; 0 0 varargin{2}];
        else
            errordlg('Plan non reconnu');
            error('Plan non reconnu');
        end
    else % If normal/point formalism is used
        useStandardizedPlane = false;
        % Test for dimensions
        if ~(numel(varargin{1})==3 && length(varargin{1}) == 3 && numel(varargin{2})==3 && length(varargin{2}) == 3)
            error('Normal and point must be a 3x1 or 1x3 vector')
        end
        if size(varargin{1},1) == 3
            varargin{1} = varargin{1}';
        end
        if size(varargin{2},1) == 3
            varargin{2} = varargin{2}';
        end
        
        % set-up the plane
        plane = [varargin{1}; varargin{2}];
    end
        
    % Compute intersection between the plane and each lines
    [p_inter, check] = plane_line_intersect(plane(1,:)', plane(2,:)', p(t_tp(:,1:2:end),:)', p(t_tp(:,2:2:end),:)');
    % Remove every outside crossing 
    p_inter = p_inter(:,check==1);
    % Compute the convex hull in the desired plane
    if size(p_inter,2) <= 3
        % If there is 3 points or less, the convex hull is all of them
        % independetly of ther order
        p2 = p_inter; 
    else
        if useStandardizedPlane
            % If it's a standard plane (align with global reference frame, we
            % can just remove the axes of that plane to compute the convex hull
            p2 = p_inter(~plane(1,:),:);
        else
            % If its a plane not aligned on global axis, we must align point on a global reference frame in
            % order to compute the convex hull
            % Find a matrix of rotation
            Z = plane(1,:)'; % Z is the normal of the plane
            % Create any point which is not collinear to Z
            X = [1 0 0]';
            if dot(X,Z)/(norm(X)*norm(Z)) >0.95
                X = double(~X);
            end
            % Calculate, recalculate and normalize axes
            Y = cross(Z,X);
            X = cross(Y,Z);
            Z = Z/norm(Z); 
            X = X/norm(X);
            Y = Y/norm(Y);
            
            % Create the rotation matrix
            R = [X, Y, Z];
            
            % Rotate each element
            p2 = R' * p_inter;
            
            % Remove the z axis
            p2 = p2(1:2, :);
        end
    end
    
    % Compute convexe hull
    if isempty(p2)
        P = double.empty(3,0);
        return;
    else
        idx = convhull(p2(1,:), p2(2,:));
    end
    % Returns points
    P = p_inter(:,idx);
    
end