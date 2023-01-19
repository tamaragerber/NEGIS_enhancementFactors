% estimate eigenvalues from CP

clear all
close all

load CP.mat;

sm = shaperead('shear_margin_outline.shp');
wgs84 = almanac('earth','wgs84','meters');

figure()
plot(sm(1).X,sm(1).Y);hold on
plot(sm(2).X,sm(2).Y)
plot(sm(3).X,sm(3).Y)

[smlat1, smlon1] = polarstereo_inv(sm(1).X,sm(1).Y, 6378137.0,0.08181919,70,-45);
[smlat2, smlon2] = polarstereo_inv(sm(2).X,sm(2).Y, 6378137.0,0.08181919,70,-45);
[smlat3, smlon3] = polarstereo_inv(sm(3).X,sm(3).Y, 6378137.0,0.08181919,70,-45);

d = distance(smlat1(1:end-1),smlon1(1:end-1),smlat1(2:end),smlon1(2:end),wgs84);
smdist1 = [0,cumsum(d)]*1e-3;
smdist1_int = [smdist1(1):0.5:smdist1(end-1)];
smlat1_int = interp1(smdist1(1:end-1),smlat1(1:end-1),smdist1_int);
smlon1_int = interp1(smdist1(1:end-1),smlon1(1:end-1),smdist1_int);

d = distance(smlat2(1:end-1),smlon2(1:end-1),smlat2(2:end),smlon2(2:end),wgs84);
smdist2 = [0,cumsum(d)]*1e-3;
smdist2_int = [smdist2(1):0.5:smdist2(end-1)];
smlat2_int = interp1(smdist2(1:end-1),smlat2(1:end-1),smdist2_int);
smlon2_int = interp1(smdist2(1:end-1),smlon2(1:end-1),smdist2_int);

d = distance(smlat3(1:end-1),smlon3(1:end-1),smlat3(2:end),smlon3(2:end),wgs84);
smdist3 = [0,cumsum(d)]*1e-3;
smdist3_int = [smdist3(1):0.5:smdist3(end-1)];
smlat3_int = interp1(smdist3(1:end-1),smlat3(1:end-1),smdist3_int);
smlon3_int = interp1(smdist3(1:end-1),smlon3(1:end-1),smdist3_int);

polygx = [sm(1).X,sm(2).X,flip(sm(3).X),sm(1).X(1)];
polygy = [sm(1).Y,sm(2).Y,flip(sm(3).Y),sm(1).Y(1)];

r = find(isnan(polygx)==1);
polygx(r) =[];
polygy(r) = [];
plot(polygx,polygy);

save('icestream_in.mat','polygx','polygy', 'smlat1_int','smlat2_int','smlat3_int','smlon1_int','smlon2_int','smlon3_int')

coord = load('cp_coordinates.mat');

lat = coord.coord(:,1);
lon = coord.coord(:,2);
[cpx cpy] = polarstereo_fwd(lat,lon, 6378137.0,0.08181919,70,-45);


% find cp inside shear margins
% -----------------------------------------------------------------------
in = inpolygon(cpx,cpy,polygx,polygy);

cpicestream = find(in==1);
cpoutside = find(in==0);

% find delta lambda
% -------------------------------------------------------------------------
for j=1:length(in)
    try
    dellam(j) = abs(CP{j}.slope_dlambda);
    catch
    dellam(j) = NaN;    
    end
    if dellam(j) >1
        dellam(j)= NaN;
    end
    
    if in(j) == 1 %inside ice stream
        lambda_x(j) = 0;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
        
    else %outside ice stream 
        lambda_x(j) = 0.1;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
    end
    
    d = distance(lat(j),lon(j),[smlat1_int,smlat2_int,smlat3_int],[smlon1_int, smlon2_int, smlon3_int],wgs84);
    if min(d)<2000 % vicinity of sm
        lambda_x(j) = (1-dellam(j))/3;
        lambda_z(j) = (1-dellam(j))/3;
        lambda_y(j) = 1-lambda_z(j)-lambda_x(j);
    end
        
        
    if lambda_x(j)>1 ||lambda_y(j)>1||lambda_z(j)>1 ||lambda_x(j)<0 ||lambda_y(j)<0||lambda_z(j)<0
        lambda_x(j) = NaN;
        lambda_y(j) = NaN;
        lambda_z(j) = NaN;
    end

end

%% plot

for j=1:length(in)
    i=in(j);
    [x,y] = polarstereo_fwd(CP{j}.lat,CP{j}.lon, 6378137.0,0.08181919,70,-45);
    scatter(x,y,25,lambda_x(j),'filled');
end
colorbar()


csv_table(:,1) = lat;
csv_table(:,2) = lon;
csv_table(:,3) = cpx;
csv_table(:,4) = cpy;
csv_table(:,5) = lambda_x';
csv_table(:,6) = lambda_y';
csv_table(:,7) = lambda_z';

nn = find(isnan(lambda_y)==1);
csv_table(nn,:) = [];
writematrix(csv_table,'a2_cp.csv')


%% Eigenvalues from birefrincence data
clear all
close all
wgs84 = almanac('earth','wgs84','meters');

sm = shaperead('shear_margin_outline.shp');

load('icestream_in.mat')

D = readmatrix('combined_D.csv');
L = readmatrix('combined_L.csv');
Q = readmatrix('combined_Q.csv');

% remove 0 points 
o = find(D(:,3) == 0);
D(o,:) = [];
o = find(L(:,3) == 0);
L(o,:) = [];
o = find(Q(:,3) == 0);
Q(o,:) = [];

% Q:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = Q(:,1);
lon = Q(:,2);
[x,y] = polarstereo_fwd(lat,lon, 6378137.0,0.08181919,70,-45);
dellam = Q(:,3);

figure()
scatter(x,y,25,dellam,'filled');hold on
plot(sm(1).X,sm(1).Y,'.k');hold on
plot(sm(2).X,sm(2).Y,'.k')
plot(sm(3).X,sm(3).Y,'.k')

% find data inside shear margins
% -----------------------------------------------------------------------
in = inpolygon(x,y,polygx,polygy);

cpicestream = find(in==1);
cpoutside = find(in==0);

% find delta lambda
% -------------------------------------------------------------------------
for j=1:length(in)
   
    if in(j) == 1 %inside ice stream
        lambda_x(j) = 0;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
        
    else %outside ice stream 
        lambda_x(j) = 0.1;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
    end
    
    d = distance(lat(j),lon(j),[smlat1_int,smlat2_int,smlat3_int],[smlon1_int, smlon2_int, smlon3_int],wgs84);
    if min(d)<2000 % vicinity of sm
        lambda_x(j) = (1-dellam(j))/3;
        lambda_z(j) = (1-dellam(j))/3;
        lambda_y(j) = 1-lambda_z(j)-lambda_x(j);
    end
        
    if lambda_x(j)>1 ||lambda_y(j)>1||lambda_z(j)>1 ||lambda_x(j)<0 ||lambda_y(j)<0||lambda_z(j)<0
        lambda_x(j) = NaN;
        lambda_y(j) = NaN;
        lambda_z(j) = NaN;
    end

end

csv_tableQ(:,1) = lat;
csv_tableQ(:,2) = lon;
csv_tableQ(:,3) = x;
csv_tableQ(:,4) = y;
csv_tableQ(:,5) = lambda_x';
csv_tableQ(:,6) = lambda_y';
csv_tableQ(:,7) = lambda_z';

nn = find(isnan(lambda_y)==1);
csv_tableQ(nn,:) = [];

writematrix(csv_tableQ,'a2_birQ.csv')

clear lat lon x y lambda_x lambda_y lambda_z

%%
% L:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = L(:,1);
lon = L(:,2);
[x,y] = polarstereo_fwd(lat,lon, 6378137.0,0.08181919,70,-45);
dellam = L(:,3);

figure()
scatter(x,y,25,dellam,'filled');hold on
plot(sm(1).X,sm(1).Y,'.k');hold on
plot(sm(2).X,sm(2).Y,'.k')
plot(sm(3).X,sm(3).Y,'.k')

% find data inside shear margins
% -----------------------------------------------------------------------
in = inpolygon(x,y,polygx,polygy);

cpicestream = find(in==1);
cpoutside = find(in==0);

% find delta lambda
% -------------------------------------------------------------------------
for j=1:length(in)
   
    if in(j) == 1 %inside ice stream
        lambda_x(j) = 0;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
        
    else %outside ice stream 
        lambda_x(j) = 0.1;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
    end
    
    d = distance(lat(j),lon(j),[smlat1_int,smlat2_int,smlat3_int],[smlon1_int, smlon2_int, smlon3_int],wgs84);
    if min(d)<2000 % vicinity of sm
        lambda_x(j) = (1-dellam(j))/3;
        lambda_z(j) = (1-dellam(j))/3;
        lambda_y(j) = 1-lambda_z(j)-lambda_x(j);
    end
        
    if lambda_x(j)>1 ||lambda_y(j)>1||lambda_z(j)>1 ||lambda_x(j)<0 ||lambda_y(j)<0||lambda_z(j)<0
        lambda_x(j) = NaN;
        lambda_y(j) = NaN;
        lambda_z(j) = NaN;
    end

end

csv_tableL(:,1) = lat;
csv_tableL(:,2) = lon;
csv_tableL(:,3) = x;
csv_tableL(:,4) = y;
csv_tableL(:,5) = lambda_x';
csv_tableL(:,6) = lambda_y';
csv_tableL(:,7) = lambda_z';

nn = find(isnan(lambda_y)==1);
csv_tableL(nn,:) = [];

writematrix(csv_tableL,'a2_birL.csv')

clear lat lon x y lambda_x lambda_y lambda_z 

%% D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = D(:,1);
lon = D(:,2);
[x,y] = polarstereo_fwd(lat,lon, 6378137.0,0.08181919,70,-45);
dellam = D(:,3);
figure()
scatter(x,y,25,dellam,'filled');hold on
plot(sm(1).X,sm(1).Y,'.k');hold on
plot(sm(2).X,sm(2).Y,'.k')
plot(sm(3).X,sm(3).Y,'.k')

% find data inside shear margins
% -----------------------------------------------------------------------
in = inpolygon(x,y,polygx,polygy);

cpicestream = find(in==1);
cpoutside = find(in==0);

% find delta lambda
% -------------------------------------------------------------------------
for j=1:length(in)
   
    if in(j) == 1 %inside ice stream
        lambda_x(j) = 0;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
        
    else %outside ice stream 
        lambda_x(j) = 0.1;
        lambda_y(j) = lambda_x(j) + dellam(j);
        lambda_z(j) = 1-lambda_y(j)-lambda_x(j);
    end
    
    d = distance(lat(j),lon(j),[smlat1_int,smlat2_int,smlat3_int],[smlon1_int, smlon2_int, smlon3_int],wgs84);
    if min(d)<2000 % vicinity of sm
        lambda_x(j) = (1-dellam(j))/3;
        lambda_z(j) = (1-dellam(j))/3;
        lambda_y(j) = 1-lambda_z(j)-lambda_x(j);
    end
        
    if lambda_x(j)>1 ||lambda_y(j)>1||lambda_z(j)>1 ||lambda_x(j)<0 ||lambda_y(j)<0||lambda_z(j)<0
        lambda_x(j) = NaN;
        lambda_y(j) = NaN;
        lambda_z(j) = NaN;
    end

end

csv_tableD(:,1) = lat;
csv_tableD(:,2) = lon;
csv_tableD(:,3) = x;
csv_tableD(:,4) = y;
csv_tableD(:,5) = lambda_x';
csv_tableD(:,6) = lambda_y';
csv_tableD(:,7) = lambda_z';

nn = find(isnan(lambda_y)==1);
csv_tableD(nn,:) = [];

writematrix(csv_tableD,'a2_birD.csv')

clear lat lon x y lambda_x lambda_y lambda_z 

%% %% Eigenvalues from ice flow model
clear all
close all


mod = shaperead('theta_flow_model.shp');
for i = 1:length(mod)
    x(i) = mod(i).x_nps;
    y(i) = mod(i).y_nps;
    [lat(i), lon(i)] = polarstereo_inv(x(i),y(i), 6378137.0,0.08181919,70,-45);
    
    T = [mod(i).E_hor_max,0,0;
        0,mod(i).E_hor_min,0
        0, 0, mod(i).E_vert];
    R = [cosd(mod(i).theta_flow), -sind(mod(i).theta_flow), 0
        sind(mod(i).theta_flow), cosd(mod(i).theta_flow), 0
        0, 0, 1];
    T_rot = R*T*R';
    
    lambda_hmin(i) = mod(i).E_hor_min;
    lambda_hmax(i) = mod(i).E_hor_max;
    lambda_z(i) = mod(i).E_vert;
    t= mod(i).theta_flow;
    if t<90
        theta(i) = t+90;
    elseif t>90
        theta(i) = t-90;
    end
end

csv_table_mod(:,1) =lat;
csv_table_mod(:,2) = lon;
csv_table_mod(:,3) = x;
csv_table_mod(:,4) = y;
csv_table_mod(:,5) = lambda_hmin;
csv_table_mod(:,6) = lambda_hmax;
csv_table_mod(:,7) = lambda_z;
csv_table_mod(:,8) = theta;


writematrix(csv_table_mod,'a2_model.csv') 
