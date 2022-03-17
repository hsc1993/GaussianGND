% This DEMO calculates a Voronoi diagram with arbitrary points in arbitrary
% polytope/polyheron in 2D/3D

clear all;close all;clc

epsilon = 1e-4;
case_str = 'before_cubic_random';
% case_str = 'during_cubic_random';
% case_str = 'after_cubic_random';

n = [5,6]; % number of middle cells (not cells on the edge)
m = 10;  % partition number for cell edge length from lowlim to highlim 
        % (for n=1, cell edge can vary from (xmax-xmin)/2 to (xmax-xmin)) 

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

    for jj = 1:m
        x_sites = xmin + (xmax-xmin).*rand(n_middle_cell,1);
        y_sites = ymin + (ymax-ymin).*rand(n_middle_cell,1);
        z_sites = zmin + (zmax-zmin).*rand(n_middle_cell,1);
        
        x_sites = [xmin; x_sites; xmax];
        y_sites = [ymin; y_sites; ymax];
        z_sites = [zmin; z_sites; zmax];


        count = 0;
        for i = 1:size(x_sites,1)-1
            for j = 1:size(y_sites,1)-1
                for k = 1:size(z_sites,1)-1
                    count = count + 1;
                    vorvx_cube_cell{count} = [[x_sites(i),y_sites(j),z_sites(k)];
                        [x_sites(i),y_sites(j),z_sites(k+1)];
                        [x_sites(i),y_sites(j+1),z_sites(k)];
                        [x_sites(i),y_sites(j+1),z_sites(k+1)];
                        [x_sites(i+1),y_sites(j),z_sites(k)];
                        [x_sites(i+1),y_sites(j),z_sites(k+1)];
                        [x_sites(i+1),y_sites(j+1),z_sites(k)];
                        [x_sites(i+1),y_sites(j+1),z_sites(k+1)]];
                    cellsize = (abs(x_sites(i+1)-x_sites(i))*abs(y_sites(i+1)-y_sites(i))*abs(z_sites(i+1)-z_sites(i)))^(1/3);
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

vorvx_namestring = ['n',num2str(n_middle_cell),'_vorvx',num2str(jj),'.txt'];
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