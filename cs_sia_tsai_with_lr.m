%Load climate data already pre-treated

precip = load('precip.m');
sat = load('sat.m');


%Melting component: Tsai et al., 2018 generalised to the entire grid
% Time t, in days
% Air temperature T_a(t), sampled at discrete time points corresponding with t, in degrees C
% Hp, percolation layer thickness, in meters
% k, effective conductivity of air, in W/m/degC
% h0, assumed height of surface temperature measurement, in meters
% L, latent heat of fusion of water, in J/kg
% Cp, specific heat of ice, in J/kg/degC
% rho, density of ice, in kg/m^3
% dt, output time resolution, in days
% Alternatively: Instead of k, DDF can be provided, in which case k=h0*rho*L*DDF
% Note: For faster performance, t and Ta can be defined functionally,
% rather than passed at each time step.

%Computes the accumulation    

accu = zeros(6144,360);
%for w = 1:12
    for j = 1:1:6144
        for k = 1:1:360
            if sat(j,k) < 0
                accu(j,k) = precip(j,k);
            end
        end
    end
%end

% dimension of a cs face
cx = 32;

% read the geometric arrays that give grid cell area, sides, 
% and centre spacings
subfold='GridMask_Laure';
land=read_cs_bin([subfold '/mask_32x6x32.bin'],1,'real*8',cx);
topo=read_cs_bin([subfold '/topo_malai_v2_ISOCORRECT.bin'],1,'real*8',cx);
RAC=read_cs_bin([subfold '/RAC.data'],1,'real*8',cx);
DXC=read_cs_bin([subfold '/DXC.data'],1,'real*8',cx);
DYC=read_cs_bin([subfold '/DYC.data'],1,'real*8',cx);
DXG=read_cs_bin([subfold '/DXG.data'],1,'real*8',cx);
DYG=read_cs_bin([subfold '/DYG.data'],1,'real*8',cx);
XC=read_cs_bin([subfold '/XC.data'],1,'real*8',cx);
YC=read_cs_bin([subfold '/YC.data'],1,'real*8',cx);


% physical parameters 

% Glen's law parameter: Pa^3 per year
% this corresponds to *very* soft ice
% +17 C according the Cuffey and Paterson relation!
% and should be made a bit smaller..
Aglen = 5e-15;

%Lapse rate
lapse_r = -0.0065 %K/m


% Exponent in Glens Law
nglen = 3;

% density x gravity
rhog = (920)*9.8;

Csia = 2*Aglen/(nglen+2)*(rhog)^nglen;

% model time step
dt = 1;

topol_file = 'topol.txt';

%%% GET CS TOPOLOGY INFORMATION %%%%%%

% Open the file for reading
fid = fopen(topol_file, 'r');

% Initialize an empty array to store the objects
objArray = [];

% N = 1
% S = 2
% E = 3
% W = 4

% classdef DataRow
%     properties
%         Field1 = Edge number (1 thru 6)
%         Field2 = Face
%         Field3 = Corresponding edge of neighbor facet
%         Field4 = Neighbor facet
%     end
% end

% Read each line from the file
while ~feof(fid)
    line = fgetl(fid); % Read a line
    
    if ischar(line)
        % Split the line into space-delimited tokens
        tokens = strsplit(line);
        
        % Check if the line has enough tokens
        if length(tokens) >= 7
            % Create a new DataRow object and assign fields
            %newRow = DataRow;
            newRow.Field1 = tokens{1};
            newRow.Field2 = tokens{3};
            newRow.Field3 = tokens{5};
            newRow.Field4 = tokens{7};
            
            % Append the new object to the array
            objArray = [objArray; newRow];
        end
    end
end

% Close the file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dict = containers.Map( ...
    {'N.Edge' 'S.Edge' 'E.Edge' 'W.Edge'},[1 2 3 4]);

% Now turn these objects into an array for easy access

%% (FACET, EDGE) --> NEIGHBOR FACET, NEIGHBOR EDGE

facets = zeros(6,4,2);
for i=1:24;
   
     facet = str2num(objArray(i).Field2);
     neighborFacet = str2num(objArray(i).Field4);
     edge = dict(objArray(i).Field1);
     neighborEdge = dict(objArray(i).Field3);
     facets(facet,edge,:) = [neighborFacet neighborEdge];
   
end

%% GENERATE LIST OF unique cell id's (Degrees Of Freedom)
dofs = reshape(1:(cx*cx*6),[32 6 32]);

%% DOF, NORTH/SOUTH/EAST/WEST NEIGHBOR
%% following is set of matlab commands to generate
%% a N x 5 matrix, where 
%%   first column is dof ID
%%   2nd   column is ID of north neigbor
%%   3rd   column is ID of south neigbor
%%   4th   column is ID of east  neigbor
%%   5th   column is ID of west  neigbor
%%
%% a python implementation may differ

dof_matrix = zeros(max(max(max(dofs))),5);
dof_matrix(:,1) = 1:(cx^2*6);


doffaces = zeros(32,32,6);
for i=1:6;
    doffaces(:,:,i) = squeeze(dofs(:,i,:));
end

for i=1:6;
    east_dofs = zeros(size(doffaces(:,:,1)));
    west_dofs = east_dofs;
    noth_dofs = east_dofs;
    soth_dofs = east_dofs;

    east_dofs(:,1:end-1) = doffaces(:,2:end,i);
    west_dofs(:,2:end) = doffaces(:,1:end-1,i);
    noth_dofs(1:end-1,:) = doffaces(2:end,:,i);
    soth_dofs(2:end,:) = doffaces(1:end-1,:,i);

    % north edge
    neighbor = facets(i,1,1);
    edge = facets(i,1,2);
    switch(edge)
        case(1);
            noth_dofs(end,:) = ...
                doffaces(end,end:-1:1,neighbor);
        case(2);
            noth_dofs(end,:) = ...
                doffaces(1,:,neighbor);
        case(3);
            noth_dofs(end,:) = ...
                doffaces(:,end,neighbor);
        case(4);
            noth_dofs(end,:) = ...
                doffaces(end:-1:1,1,neighbor);
    end

    % south edge
    neighbor = facets(i,2,1);
    edge = facets(i,2,2);
    switch(edge)
        case(1);
            soth_dofs(1,:) = ...
                doffaces(end,:,neighbor);
        case(2);
            soth_dofs(1,:) = ...
                doffaces(1,end:-1:1,neighbor);
        case(3);
            soth_dofs(1,:) = ...
                doffaces(end:-1:1,end,neighbor);
        case(4);
            soth_dofs(1,:) = ...
                doffaces(:,1,neighbor);
    end

    % east edge
    neighbor = facets(i,3,1);
    edge = facets(i,3,2);
    switch(edge)
        case(1);
            east_dofs(:,end) = ...
                doffaces(end,:,neighbor);
        case(2);
            east_dofs(:,end) = ...
                doffaces(1,end:-1:1,neighbor);
        case(3);
            east_dofs(:,end) = ...
                doffaces(end:-1:1,end,neighbor);
        case(4);
            east_dofs(:,end) = ...
                doffaces(:,1,neighbor);
    end

    % west edge
    neighbor = facets(i,4,1);
    edge = facets(i,4,2);
    switch(edge)
        case(1);
            west_dofs(:,1) = ...
                doffaces(end,end:-1:1,neighbor);
        case(2);
            west_dofs(:,1) = ...
                doffaces(1,:,neighbor);
        case(3);
            west_dofs(:,1) = ...
                doffaces(:,end,neighbor);
        case(4);
            west_dofs(:,1) = ...
                doffaces(end:-1:1,1,neighbor);
    end

    current_dofs = doffaces(:,:,i);
    dof_matrix(current_dofs(:),2) = noth_dofs(:);
    dof_matrix(current_dofs(:),3) = soth_dofs(:);
    dof_matrix(current_dofs(:),4) = east_dofs(:);
    dof_matrix(current_dofs(:),5) = west_dofs(:);


end
%% FINISH GENERATING DOF_MATRIX

%% CONVERGE GRID METRICS AND OTHER CS FIELDS TO VECTORS

dxc = cs_to_vec(dofs,DXC);
dyc = cs_to_vec(dofs,DYC);
dxg = cs_to_vec(dofs,DXG);
dyg = cs_to_vec(dofs,DYG);
rac = cs_to_vec(dofs,RAC);
land_vec = cs_to_vec(dofs,land);
xc = cs_to_vec(dofs,XC);
yc = cs_to_vec(dofs,YC);
topo_vec = cs_to_vec(dofs,topo);

flow2north = (land_vec == 1) & (land_vec(dof_matrix(:,2))==1);
flow2south = (land_vec == 1) & (land_vec(dof_matrix(:,3))==1);
flow2east  = (land_vec == 1) & (land_vec(dof_matrix(:,4))==1);
flow2west  = (land_vec == 1) & (land_vec(dof_matrix(:,5))==1);

%% INITIALISE THICKNESS (h_vec)
%% this way of initialisation would change in a coupled model
h_vec = 1e3*ones(length(dof_matrix),1);
h_vec(land_vec==0) = 0;

%% SAVE INITIAL THICKNESS
h_vec0 = h_vec;
sat_iter = sat;

for it=1:10000;

% begin building linear system
    
    surf_vec = topo_vec + h_vec;
    
    % east neighbor minus west
    dsdx_c = (surf_vec(dof_matrix(:,4)) - surf_vec(dof_matrix(:,5))) ./ ...
             (dxc(dof_matrix(:,1)) + dxc(dof_matrix(:,5)));
    % north neighbor minus south
    dsdy_c = (surf_vec(dof_matrix(:,2)) - surf_vec(dof_matrix(:,3))) ./ ...
             (dyc(dof_matrix(:,1)) + dyc(dof_matrix(:,3)));
    
    D_vec = Csia * (dsdx_c.^2 + dsdy_c.^2 + 1e-8) .^ ((nglen-1)/2) .* ...
        h_vec .^ (nglen+2);
    
    % diffusivity at south and west faces
    Dx_vec = .5 * (D_vec(dof_matrix(:,1)) + D_vec(dof_matrix(:,5)));
    Dy_vec = .5 * (D_vec(dof_matrix(:,1)) + D_vec(dof_matrix(:,3)));
    
    Dx_west = Dx_vec;
    Dx_east = Dx_vec(dof_matrix(:,4));
    Dy_south = Dy_vec;
    Dy_north = Dy_vec(dof_matrix(:,2));
    
    Brhs = h_vec;
    
    Iland = land_vec==1;
    I = Iland;
    Amat00 = ones(size(h_vec));
    Amat00(land_vec==1) = Amat00(Iland) + ...
        dt./rac(I) .* ( ...
            + dyg(I)./dxc(dof_matrix(I,5)) .* Dx_west(I)  ...
            + dyg(dof_matrix(I,4))./dxc(I) .* Dx_east(I)  ...
            + dxg(I)./dyc(dof_matrix(I,3)) .* Dy_south(I) ...
            + dxg(dof_matrix(I,2))./dyc(I) .* Dy_north(I) ...
            );
    Brhs(I) = Brhs(I) + ...
        dt./rac(I) .* ( ...
            - dyg(I)./dxc(dof_matrix(I,5)) .* Dx_west(I) .* ...
                (topo_vec(I) - topo_vec(dof_matrix(I,5)))  ...
            + dyg(dof_matrix(I,4))./dxc(I) .* Dx_east(I) .* ...
                (topo_vec(dof_matrix(I,4)) - topo_vec(I))  ...
            - dxg(I)./dyc(dof_matrix(I,3)) .* Dy_south(I) .* ...
                (topo_vec(I) - topo_vec(dof_matrix(I,3)))  ...
            + dxg(dof_matrix(I,2))./dyc(I) .* Dy_north(I) .* ...
                (topo_vec(dof_matrix(I,2)) - topo_vec(I)) ...
            );
    
    
    %Modification of the TSAI code for the melting
    %Lapse rate correction

    for i = 1:size(sat,2)
        for j = 1:size(sat,1)
            sat_iter(j,i) = sat_iter(j,i) + lapse_r * h_vec(j);
        end
    end

    melt_map = zeros(6144,360);%initialize the melting grid for each cell and for each day

    for j = 1:1:6144
        t = [1:1:360];
        Ta = sat(j,:);
        Hp=10; k=40; h0=2; L=334*10^3; Cp=2100; rho=920; dt = 1.0;
        
        % Calculates rate factor for ODE, C_rate, and DDF, degree day factor (with unit conversions)
        C_rate = k/(h0*rho*Cp*Hp);
        DDF = (1000*3600*24*k/(h0*rho*L)); %*360 pour faire a l'année
        
        % Alternatively, uses DDF as input and calculates k. If desired, change 0 to 1
        %if 0
            % DDF can be defined instead of k, in units of mm/day/degC
         %   DDF = 10; %= daily DDF 
            % Calculates k, including unit conversion to SI units
          %  k = DDF/(1000*3600*24)*h0*rho*L;
        %end
        
        % Defines start and end times in seconds, ODE tolerances
        t0=3600*24*t(1); tf=3600*24*t(end);
        options = odeset('RelTol',1e-5,'AbsTol',1e-6);
        
        % Runs ODE solver to solve for percolation temperature Tp
        [t_out,Tp] = ode45(@(t_out,y) percolation_ode(t_out,y,C_rate,t,Ta),[t0:3600*24*dt:tf],Ta(1),options);
        t_out_days = t_out/(3600*24);
        %%sat
        % Calculates melt rate, given Tp(t) and Ta(t)
        melt_rate = 0*t_out;
        for i=1:length(t_out)
            if Tp(i)>0
                melt_rate(i) = DDF*max(interp1(t,Ta,t_out_days(i)),0);
            end
        end
    
        %melt = mean(melt_rate);
        melt_map(j,:) = melt_rate;
        
    end

    %The SMB is calculated for every grid cell, every day
    smb_1 = zeros(6144,360);
    smb_1 = accu - melt_map;
    
    %Compute the accumulative SMB
    smb_accumulative = zeros(6144,1);
    smb_accumulative = sum(smb_1,2);
    
    %Remove values on ocean
    path(path,'/Users/laure/Desktop/Modele_icesheet');
    land=read_cs_bin(['GridMask_Laure/mask_32x6x32.bin'],1,'real*8',32);
    land = reshape(land,[6144 1]);
    iii = find(land==0);
    smb_accumulative(iii) = NaN;
    smb_accumulative_c = smb_accumulative(~isnan(smb_accumulative));
    
    [a_smb ela] =  smb(surf_vec(I),yc(I), smb_accumulative(I));
    
    Brhs(I) = Brhs(I) + dt*a_smb;
    
    
    AmatNorth = zeros(size(h_vec));
    I = flow2north;
    AmatNorth(I) = AmatNorth(I) ...
        - dt./rac(I) .* ...
            dxg(dof_matrix(I,2))./dyc(I) .* Dy_north(I);
    
    AmatSouth = zeros(size(h_vec));
    I = flow2south;
    AmatSouth(I) = AmatSouth(I) ...
        - dt./rac(I) .* ...
            dxg(I)./dyc(dof_matrix(I,3)) .* Dy_south(I);
    
    AmatEast = zeros(size(h_vec));
    I = flow2east;
    AmatEast(I) = AmatEast(I) ...
        - dt./rac(I) .* ...
            dyg(dof_matrix(I,4))./dxc(I) .* Dx_east(I);
    
    AmatWest = zeros(size(h_vec));
    I = flow2west;
    AmatWest(I) = AmatWest(I) ...
        - dt./rac(I) .* ...
            dyg(I)./dxc(dof_matrix(I,5)) .* Dx_west(I);
    

%% Acols give coefficients in a row for weights of corresponding cell and 
%% its neighbors 
%% .. although these do not correspond to matrix diagonals
    Acols = [Amat00 AmatNorth AmatSouth AmatEast AmatWest];
    func_handle = @(x) cgfunc(x,dof_matrix,Acols);

%% A DIRECT SOLVE WOULD REQUIRE FORMING THE MATRIX
%% WHICH DOES NOT HAVE A REGULAR SPARSITY PATTERN IN CS
%% SO matrix *action* is implemented as a function, and an indirect solver is used
    [h_vecnew,cgflag,cgrelres,cgiter] = cgs(func_handle,Brhs);
    h_vecnew(h_vecnew<0) = 0;
    
    disp(['timestep ' num2str(it)]);
    disp('finished cg, [flag relerr iter]: ')
    disp(['[' num2str(cgflag) ',' num2str(cgrelres) ',' num2str(cgiter) ']']);

    h_vec = h_vecnew;
end

H = vec_to_cs(h_vec);
H0 = vec_to_cs(h_vec0);


figure(1); crossmap(H); %icesheet height only

%Isostatic correction
rho_ice = 920;
rho_lito = 3300;
ration = rho_ice / rho_lito;
wss = ration * H;
topo_new = topo - wss + H;

figure(2); crossmap(topo_new);%total topography
%%------------------------------------------------------%%%


function vec = cs_to_vec (cs_dof, vec_cs)
    vec_unordered = vec_cs(:);
    dof_unordered = cs_dof(:);
    [temp I] = sort(dof_unordered);
    vec = vec_unordered(I);
end

function cs = vec_to_cs (vec)
    cs=reshape(vec,[32 6 32]);
end


%% FUNCTION WHICH IMPLEMENTS MATRIX ACTION
function y = cgfunc (x, varargin)
    
    dofs = varargin{1};
    A_columns = varargin{2};
    x_vals = x(dofs);
    y = sum(A_columns.*x_vals,2);
end


function [a, ela] = smb(surf, lat, smb_accu)
    ela = zeros(size(lat)); 
    ela(lat < -30) = 5500 + 100 * lat(lat < -30);    
    ela(lat > 30) = 5500 - 100 * lat(lat > 30);     
    ela(lat >= -30 & lat <= 30) = 2500;             

    
    elev_diff = surf - ela;

    
    if isequal(size(smb_accu), size(surf))
        a = smb_accu;  
    else
        a(elev_diff > 0) = 1 * elev_diff(elev_diff > 0);   
        a(elev_diff < 0) = 2 * elev_diff(elev_diff < 0);   
    end

    
    %a(a > 1000) = 1000;  %???????  
    %a(a < -20) = -20;  %???????
end

    
    
   
