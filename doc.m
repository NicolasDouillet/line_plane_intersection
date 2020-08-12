%% line_plane_intersection
%
% Function to compute the intersection point
% between the (N,u) line and the (M,n) plane of the 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%% Syntax
%
% [I,rc] = line_plane_intersection(u, N, n, M);
%
% [I,rc] = line_plane_intersection(u, N, n, M, verbose);
%
%% Description
%
% [I,rc] = line_plane_intersection(u, N, n, M) computes the coordinates of I,
% the intersection point between the line (u,N) and the plane (n,M).
%
% In the most generic case, I is a point in the 3D space, but when
% the line is stricly parallel to the plane, I is the empty set, and when
% the line is included in the plane, I is a function handle corresponding
% to the system of parametric equations of the line.
%
% [I,rc] = line_plane_intersection(u, N, n, M, verbose) displays a message in
% console when verbose is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73853-planes-intersection planes_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73922-lines-intersection-3d-2d?s_tid=prof_contriblnk lines_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73490-point-to-plane-distance?s_tid=prof_contriblnk point_to_plane_distance> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73478-point-to-line-distance-3d-2d?s_tid=prof_contriblnk point_to_line_distance>
%
%% Principle
%
% Based on solving Descartes plane equation :
%
% ax + by + cz + d = 0, where n = [a, b, c] is a vector normal to the plane,
%
% combined with the parametric equations system of a 3D line :
%
% x(t) = x0 + at 
% y(t) = y0 + bt
% z(t) = z0 + ct
%
% where N0 = [x0, y0, z0] is a point belonging to the line, and u = [a, b, c], a vector directing this line.
%
%% Input arguments
%
% - u : real row or column vector double. numel(u) = 3. One director vector of the parametric line.
%
% - N : real row or column vector double. numel(N) = 3. One point belonging to the line.
%
% - n : real row or column vector double. numel(n) = 3. One normal vector to the plane.
%
% - M : real row or column vector double. numel(M) = 3. One point belonging to the plane.
%
% - verbose : logical *true(1)/false(0), to enable/disable the verbose mode.
%
%% Output arguments
%
% - I = [xI yI zI], real row or column vector double, the intersection point.
%
% - rc : return code, integer scalar in the set {1,2,3}.
%
%        0 : void / [] intersection
%
%        1 : point intersection (unique).
%
%        2 : line intersection
%
%        rc return code is necessary to distinguish between cases where
%        (N,u) line and the (M,n) plane intersection is a single point
%        and where it is the line itself.
%
%% Example #1
% One unique intersection point
n = [1 1 1];
M = n;
u = [1 0 0];
N = u; % (N,u) = (OX) axis
[I,rc] = line_plane_intersection(u, N, n, M) % one unique intersction point expected : I = [3 0 0], rc = 1

%% Example #2
% Line and plane are strictly // ; no intersection
n = [0 0 1];
M = [0 0 0]; % (M,n) = (XOY) plan
u = [1 2 0];
N = [0 0 6];
[I,rc] = line_plane_intersection(u, N, n, M, false) % line strictly // plane =>  I = [], rc = 0 expected

%% Example #3
% Line is included in the plane
n = [1 1 1];
M = (1/3)*[1 1 1];
u = [1 1 -2];
N = [0.5 0.5 0];
[I,rc] = line_plane_intersection(u, N, n, M) % line belongs to the plane, rc = 2 expected