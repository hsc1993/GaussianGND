% This DEMO calculates a Voronoi diagram with arbitrary points in arbitrary
% polytope/polyheron in 2D/3D

clear all;close all;clc

epsilon = 1e-4;
case_str = 'cubic';


n = [4]; % number of middle cells (not cells on the edge)
m = 4;  % partition number for cell edge length from lowlim to highlim 
        % (for n=1, cell edge can vary from (xmax-xmin)/2 to (xmax-xmin)) 
% num_edgevariate_ramdomshift = 10;


flag_boundary_included = 0; % if boundary is included, big mesh with 1 pos and 3 bndpts can be ran
% flag_boundary_included = 1;


% box
case_box = 1;  % full simulation volume, for bigger mesh size
% case_box = 2;  % only volume for loop, for smaller mesh size
flag_figure = 1;
% flag_figure = 0;
flag_outside = 0; % flag used to indicate if vorvx is outside of bounding box

%% generate random samples
num_bndpoints_inclusion = 3; % number of boundary points (8 in total for cubic boundary) included in delaundry triangulation, because at least 4 points are needed
d = 3;          % dimension of the space
tol = 1e-07;            % tolerance value used in "inhull.m" (larger value high precision, possible numerical error)


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

bnd0 = [xmin ymin zmax; xmax ymin zmin;  xmin ymax zmin; xmax ymax zmin; xmin ymax zmax; xmax ymin zmax; xmax ymax zmax; xmin ymin zmin];       % generate boundary point-candidates

x_boxsize = xmax-xmin;
y_boxsize = ymax-ymin;
z_boxsize = zmax-zmin;


cellsize_matrix = [];
for ii = 1:size(n,2)

n_middle_cell = n(ii);  % number of discretizations for volume for regular(cubic) discretization 
nx_middle_cell = n_middle_cell;
ny_middle_cell = n_middle_cell;
nz_middle_cell = n_middle_cell;

yourFolder = strcat(case_str,num2str(n_middle_cell));
if not(isfolder(yourFolder))
    mkdir(yourFolder);
end


x_distance_lowlim = x_boxsize/(nx_middle_cell);
x_distance_highlim = x_boxsize/(nx_middle_cell-1);
y_distance_lowlim = y_boxsize/(ny_middle_cell);
y_distance_highlim = y_boxsize/(ny_middle_cell-1);
z_distance_lowlim = z_boxsize/(nz_middle_cell);
z_distance_highlim = z_boxsize/(nz_middle_cell-1);

x_distance_increment = (0.99*x_distance_highlim-x_distance_lowlim)/m;
y_distance_increment = (0.99*y_distance_highlim-y_distance_lowlim)/m;
z_distance_increment = (0.99*z_distance_highlim-z_distance_lowlim)/m;


for jj = 2:m+1
        x_distance = x_distance_lowlim+(jj-1)*x_distance_increment;
        y_distance = y_distance_lowlim+(jj-1)*y_distance_increment;
        z_distance = z_distance_lowlim+(jj-1)*z_distance_increment;    

        n_xcell_middle = floor(x_boxsize/x_distance);
        n_ycell_middle = floor(y_boxsize/y_distance);
        n_zcell_middle = floor(z_boxsize/z_distance);

        x_max_to_grid = (x_boxsize - n_xcell_middle*x_distance);
        y_max_to_grid = (y_boxsize - n_ycell_middle*y_distance);
        z_max_to_grid = (z_boxsize - n_zcell_middle*z_distance);

%         rng('shuffle')
%         x_shift = rand*(x_max_to_grid);
%         y_shift = rand*(y_max_to_grid);
%         z_shift = rand*(z_max_to_grid);       
        x_shift = 0;
        y_shift = 0;
        z_shift = 0;
        
        x_grid_min = xmin+x_shift;
        y_grid_min = ymin+y_shift;
        z_grid_min = zmin+z_shift;

        x_grid_max = xmax-x_max_to_grid+x_shift;
        y_grid_max = ymax-y_max_to_grid+y_shift;
        z_grid_max = zmax-z_max_to_grid+z_shift;


        n_xcellvertice = n_xcell_middle+1;
        n_ycellvertice = n_ycell_middle+1;
        n_zcellvertice = n_zcell_middle+1;

        x_edge_sites = [xmin,xmax];
        y_edge_sites = [ymin,ymax];
        z_edge_sites = [zmin,zmax];
        x_regular_sites = linspace(x_grid_min,x_grid_max,n_xcellvertice);
        y_regular_sites = linspace(y_grid_min,y_grid_max,n_ycellvertice);
        z_regular_sites = linspace(z_grid_min,z_grid_max,n_zcellvertice);
        cellsize = ((x_regular_sites(2)-x_regular_sites(1))*(y_regular_sites(2)-y_regular_sites(1))*(z_regular_sites(2)-z_regular_sites(1)))^(1/3);

        if ~any(x_regular_sites(:)==xmin)
            x_sites = [xmin,x_regular_sites];
        else 
            x_sites = x_regular_sites;
        end
        if ~any(y_regular_sites(:)==ymin)
            y_sites = [ymin,y_regular_sites];
        else 
            y_sites = y_regular_sites;
        end
        if ~any(z_regular_sites(:)==zmin)
            z_sites = [zmin,z_regular_sites];
        else 
            z_sites = z_regular_sites;
        end


        if ~any(x_regular_sites(:)==xmax)
            x_sites = [x_regular_sites,xmax];
        else 
            x_sites = x_regular_sites;
        end
        if ~any(y_regular_sites(:)==ymax)
            y_sites = [y_regular_sites,ymax];
        else 
            y_sites = y_regular_sites;
        end
        if ~any(z_regular_sites(:)==zmax)
            z_sites = [z_regular_sites,zmax];
        else 
            z_sites = z_regular_sites;
        end

        n_xcell = size(x_sites,2)-1;
        n_ycell = size(y_sites,2)-1;
        n_zcell = size(z_sites,2)-1;

        count = 0;
        for i = 1:n_xcell
            for j = 1:n_ycell
                for k = 1:n_zcell
                    count = count + 1;
                    vorvx_cube_cell{count} = [[x_sites(i),y_sites(j),z_sites(k)];
                        [x_sites(i),y_sites(j),z_sites(k+1)];
                        [x_sites(i),y_sites(j+1),z_sites(k)];
                        [x_sites(i),y_sites(j+1),z_sites(k+1)];
                        [x_sites(i+1),y_sites(j),z_sites(k)];
                        [x_sites(i+1),y_sites(j),z_sites(k+1)];
                        [x_sites(i+1),y_sites(j+1),z_sites(k)];
                        [x_sites(i+1),y_sites(j+1),z_sites(k+1)]];
                end
            end
        end
    


%% PLOT
for i = 1:size(vorvx_cube_cell,2)
    col(i,:)= rand(1,3);
end

if flag_figure == 1

figure('position',[0 0 600 600],'Color',[1 1 1]);
hold on;
for i = 1:size(vorvx_cube_cell,2)
    K = convhulln(vorvx_cube_cell{i});
    flag_outside = 0;
%   scatter3(vorvx_cube_cell{i}(:,1),vorvx_cube_cell{i}(:,2),vorvx_cube_cell{i}(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');

    if max(vorvx_cube_cell{i}(:,1))>xmax+0.1
        flag_outside = 1;
    end
    if min(vorvx_cube_cell{i}(:,1))<xmin-0.1
        flag_outside = 1;
    end
    if max(vorvx_cube_cell{i}(:,2))>ymax+0.1
        flag_outside = 1;
    end
    if min(vorvx_cube_cell{i}(:,2))<ymin-0.1
        flag_outside = 1;
    end
    if max(vorvx_cube_cell{i}(:,3))>zmax+0.1
        flag_outside = 1;
    end
    if min(vorvx_cube_cell{i}(:,3))<zmin-0.1
        flag_outside = 1;
    end

    if flag_outside == 1
        continue
    end

    trisurf(K,vorvx_cube_cell{i}(:,1),vorvx_cube_cell{i}(:,2),vorvx_cube_cell{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',0)

end
hold off
axis('equal')
xlabel('X');ylabel('Y');zlabel('Z');
axis off
end

vorvx_namestring = ['n',num2str(n_middle_cell),'_vorvx',num2str(jj-1),'.txt'];
vorvx_fileID = fopen(strcat(yourFolder,'/',vorvx_namestring),'w');

for i = 1:size(vorvx_cube_cell,2)
    M = cell2mat(vorvx_cube_cell(i));
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
cellsize_matrix = [cellsize_matrix,cellsize];
end

end

cellsize_matrix