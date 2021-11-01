function [I,check]=plane_line_intersect(n,V0,P0,P1)
% PLANE_LINE_INTERSECT computes the intersection of a plane and a
% segment(or a straight line), using multitoolbox when available. If not
% available, a for loop is used instead which may greatly reduces the function
% efficiency.
% Multitoolbox is available online at : 
% http://www.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled
%
% INPUTS: 
%       n: normal vector of the Plane (3x1 vector) 
%       V0: any point that belongs to the Plane (3x1 vector) 
%       P0: end point 1 of the segment P0P1 (3xN matrix)
%       P1:  end point 2 of the segment P0P1 (3xN matrix)
%
% OUTPUTS:
%      I    is the points of intersection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3 => the intersection lies outside the segment P0P1
%
% EXAMPLE:
%   Determine the intersection of following the plane x+y+z+3=0 with the
%   segment P0P1: 
%   The plane is represented by the normal vector 
%     n=[1 1 1]';
%   and an arbitrary point that lies on the plane
%     V0=[1 1 -5]'
%   The segment is represented by the following 100 pairs of points
%     P0=rand(3, 1);
%     P1=rand(3, 1)*10;
%   Function returns all the points where pairs cross the plane
%     [I,check]=plane_line_intersect(n,V0,P0,P1);
%
% DISCLAIMER and COPYRIGHT
%   This function was originally written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%   Updated to accept multitoolbox by Pariterre
%     1st updated version : April 13th, 2016
%
%
% 
    % Check if multitoolbox is available
    persistent useMulti
    if isempty(useMulti)
        if exist('multiprod', 'file')
            useMulti = true;
        end
    end
    % Make sure n and V0 are column vectors
    if ~(size(n,1) == 3 && size(n,2) == 1)
        error('Normal must be a 3x1 vector')
    end
    if ~(size(V0,1) == 3 && size(V0,2) == 1)
        error('Point on the plane must be a 3x1 vector')
    end
    % Do some calculation
    nPoints = size(P0,2);
    u = P1-P0;
    w = P0 - repmat(V0, [1 nPoints]);
    if useMulti
        D = multiprod(n', u);
        N = -multiprod(n',w);
    else
        D = nan(1,nPoints);
        N = nan(1,nPoints);
        for i = 1:size(u,2)
            D(i) = dot(n,u(:,i));
            N(i) = -dot(n,w(:,i));
        end
    end
    
    % Check for parallel cases
    idxParallel = abs(D) < 10^-7;
    idxOnPlane = idxParallel & N==0;
    idxDisjoint = idxParallel & N~=0;
    % Ccompute the intersection parameter
    sI = nan(1,nPoints);
    sI(~idxParallel) = N(~idxParallel) ./ D(~idxParallel);
    I = P0+ repmat(sI, [3,1]).*u;
    % Check if intersection is between the points or outside them
    idxIntersectBetween = (sI(1,:) > 0 & sI(1,:) < 1);
    idxIntersectOutside = ~idxIntersectBetween & ~idxParallel;
    
    % Prepare check vector
    check = zeros(1, nPoints);
    check(idxOnPlane) = 2;
    check(idxDisjoint) = 0;
    check(idxIntersectBetween) = 1;
    check(idxIntersectOutside) = 3;
end