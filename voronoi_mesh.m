% This DEMO calculates a Voronoi diagram with arbitrary points in arbitrary
% polytope/polyheron in 2D/3D

clear all;close all;clc

case_str = 'before';
% case_str = 'during';  %20,24,28,32,36    44,52,60 (4)  70,80,90(2)
% case_str = 'after';

% flag_boundary_included = 0; % if boundary is included, big mesh with 1 pos and 3 bndpts can be ran
flag_boundary_included = 1;
boxexpansion = 1;
n = [2000]; % number of voronoi cells
m = 1;  % number of trials for same number of voronoi cells
        % (for n=1, cell edge can vary from (xmax-xmin)/2 to (xmax-xmin)) 

        
%% generate random samples
num_bndpoints_inclusion = 3; % number of boundary points (8 in total for cubic boundary) included in delaundry triangulation, because at least 4 points are needed
d = 3;          % dimension of the space
tol = 1e-07;            % tolerance value used in "inhull.m" (larger value high precision, possible numerical error)
    
        
% box
case_box = 1;  % full simulation volume, for bigger mesh size
% case_box = 3;  % full simulation volume, for bigger mesh size
% case_box = 2;  % only volume for full loop, for smaller mesh size
% case_box = 4;  % only volume for partial loop, for smaller mesh size

flag_figure = 1;
% flag_figure = 0;

flag_outside = 0; % flag used to indicate if vorvx is outside of bounding box




for index_seed = 1:size(n,2)
ii = n(index_seed);

jj = 1;
while jj<=m
yourFolder = strcat(case_str,num2str(ii));
if not(isfolder(yourFolder))
    mkdir(yourFolder);
end


if case_box == 1
    xmax = 252;
    xmin = -250;
    ymax = 251;
    ymin = -251;
    zmax = 380;
    zmin = -380;
end




if case_box == 2
    % full loop, extrema calculated by main.py
    xmax = 77;
    xmin = -79;
    ymax = 68;
    ymin = -73;
    zmax = -57;
    zmin = -154;
end

if case_box == 3
    % expanded box
    xmax = ceil(252*boxexpansion);
    xmin = floor(-250*boxexpansion);
    ymax = ceil(251*boxexpansion);
    ymin = floor(-251*boxexpansion);
    zmax = ceil(380*boxexpansion);
    zmin = floor(-380*boxexpansion);
end


if case_box == 4
    % partial loop, extrema calculated by main.py
%     xmax = -4;
%     xmin = -108;
%     ymax = -48;
%     ymin = -151;
%     zmax = -137;
%     zmin = -190;
    xmax = 14;
    ymax = -28;
    zmax = -117;
    
    xmin = -250;
    ymin = -251;
    zmin = -380;
end


bnd0 = [xmin ymin zmax; xmax ymin zmin;  xmin ymax zmin; xmax ymax zmin; xmin ymax zmax; xmax ymin zmax; xmax ymax zmax; xmin ymin zmin];       % generate boundary point-candidates

K = convhull(bnd0);
bnd_pnts = bnd0(K,:);   % take boundary points from vertices of convex polytope formed with the boundary point-candidates
bnd_pnts = unique(bnd_pnts,'rows');


rng('shuffle')

pos0 = [];
for i=1:ii
    pos_x = randi([xmin xmax],1,1);
    pos_y = randi([ymin ymax],1,1);
    pos_z = randi([zmin zmax],1,1);
    pos01_temp = [pos_x,pos_y,pos_z];
    pos0 = [pos0;pos01_temp];
end


indice_bndpoints_inclusion = randperm(8,num_bndpoints_inclusion);


if flag_boundary_included == 1
    for i =1:num_bndpoints_inclusion
        pos0 = [pos0;bnd_pnts(indice_bndpoints_inclusion(i),:)];
    end
end

%% take points that are in the boundary convex polytope
in = inhull(pos0,bnd0,[],tol); 
% inhull.m is written by John D'Errico that efficiently check if points are
% inside a convex hull in n dimensions
% We use the function to choose points that are inside the defined boundary
u1 = 0;
for i = 1:size(pos0,1)
    if in(i) ==1
        u1 = u1 + 1;
        pos(u1,:) = pos0(i,:);
    end
end
%% 
% =========================================================================
% INPUTS:
% pos       points that are in the boundary      n x d matrix (n: number of points d: dimension) 
% bnd_pnts  points that defines the boundary     m x d matrix (m: number of vertices for the convex polytope
% boundary d: dimension)
% -------------------------------------------------------------------------
% OUTPUTS:
% vornb     Voronoi neighbors for each generator point:     n x 1 cells
% vorvx     Voronoi vertices for each generator point:      n x 1 cells
% =========================================================================

[vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);

for i = 1:size(pos,1)
    if max(vorvx{i}(:,1))>xmax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,1))<xmin-0.1
        flag_outside = 1;
    end
    if max(vorvx{i}(:,2))>ymax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,2))<ymin-0.1
        flag_outside = 1;
    end
    if max(vorvx{i}(:,3))>zmax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,3))<zmin-0.1
        flag_outside = 1;
    end
end

if flag_outside==1
    fprintf(strcat(num2str(jj),'outside \n'))
%         continue
end



%% PLOT
for i = 1:size(vorvx,2)
    col(i,:)= rand(1,3);
end



if flag_figure == 1

figure('position',[0 0 600 600],'Color',[1 1 1]);
hold on;
for i = 1:size(pos,1)
    K = convhulln(vorvx{i});
    flag_outside = 0;
%     scatter3(vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
    scatter3(pos0(i,1),pos0(i,2),pos0(i,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');

    if max(vorvx{i}(:,1))>xmax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,1))<xmin-0.1
        flag_outside = 1;
    end
    if max(vorvx{i}(:,2))>ymax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,2))<ymin-0.1
        flag_outside = 1;
    end
    if max(vorvx{i}(:,3))>zmax+0.1
        flag_outside = 1;
    end
    if min(vorvx{i}(:,3))<zmin-0.1
        flag_outside = 1;
    end

    if flag_outside == 1
        continue
    end

    trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',0)

end
hold off
axis('equal')
xlabel('X');ylabel('Y');zlabel('Z');
axis off

end



% Convert cell to a table and use first row as variable names
vorvx_namestring = ['n',num2str(ii),'_vorvx',num2str(jj),'.txt'];
vorvx_fileID = fopen(strcat(case_str,num2str(ii),'/',vorvx_namestring),'w');

for i = 1:size(vorvx,2)
    M = cell2mat(vorvx(i));
    A1 = M(:,1);
    A2 = M(:,2);
    A3 = M(:,3);

    fprintf(vorvx_fileID, '[');
    for j = 1:size(A1)
        formatSpec1 = '[%8.4f,';
        fprintf(vorvx_fileID,formatSpec1,A1(j));

        formatSpec2 = '%8.4f,';
        fprintf(vorvx_fileID,formatSpec2,A2(j));

        formatSpec3 = '%8.4f],\n';
        fprintf(vorvx_fileID,formatSpec3,A3(j));
    end
    fprintf(vorvx_fileID, ']');
    fprintf(vorvx_fileID, ',');
end
fclose(vorvx_fileID);



jj = jj+1;
end
end



