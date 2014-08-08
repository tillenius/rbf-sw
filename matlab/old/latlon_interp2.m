function he = latlon_interp2(la,th,h,lae,the)
% Interpolate DWD data from lat-lon grid to arbitrary points
% NOTE! Longitude now in [0,2pi)
% INPUT:
% la - longitude vector (1st column of DWD data) in [0,2pi)
% th - latitude vector (2nd column of DWD data) in (-pi/2,pi/2)
% h - height field vector (3rd column of DWD data)
% lae - longitude vector of evaluation points in [0,2pi]
% the - latitude vector of evaluation points in [-pi/2,pi/2]
% OUTPUT:
% he - height field vector at evaluation points

ne = length(la);
res = sqrt(ne*2)/2;
L = reshape(la,res*2,res);
T = reshape(th,res*2,res);
H = reshape(h,res*2,res);
h_la = L(2,1)-L(1,1);
h_th = T(1,1)-T(1,2);

L_pad = [L(:,1:2) L L(:,1:2)]; %pad top and bottom
L_pad = [L_pad(1,:) - 2*h_la ; L_pad(1,:) - h_la ; L_pad ; ...
    L_pad(end,:) + h_la ; L_pad(end,:) + 2*h_la ];  %pad sides
% T_pad = [T(:,1) + 2*h_th , T(:,1) + h_th , T , T(:,end) - h_th , ...
%     T(:,end)- 2*h_th ]; %pad top and bottom
T_pad = [pi-T(:,2:-1:1), T, -pi-T(:,end:-1:end-1)];
T_pad = [T_pad(1:2,:) ; T_pad ; T_pad(1:2,:)]; %pad sides

% Pad top and bottom with shifted columns
H_pad = [[H(res+1:end,2:-1:1);H(1:res,2:-1:1)] ...
         H [H(res+1:end,end:-1:end-1);H(1:res,end:-1:end-1)]];
% Pad left and right with two rows
H_pad = [H_pad(end-1:end,:); H_pad; H_pad(1:2,:)];

he = interp2(L_pad.',T_pad.',H_pad.',lae,the,'cubic');

% % OLD JUNK
% % % [la,is] = sort(la); h = h(is); th = th(is);
% % % i1 = 1:res; i2 = res+1:res+2; i3 = res; i4 = res+1:2*res;
% % % Lp = [L(i1,:); L(i2,:)+2*pi; L(i3,:)-2*pi; L(i4,:)];
% % % Tp = [T(i1,:); T(i2,:); T(i3,:); T(i4,:)];
% % % Hp = [H(i1,:); H(i2,:); H(i3,:); H(i4,:)];
% % % 
% % % i1 = 2:-1:1; i2 = 1:res; i3 = res:-1:res-1; i4 = [res+4:2*res+3 1:res+3];
% % % Tp = [pi-Tp(:,i1) Tp(:,i2) -pi-Tp(:,i3)];
% % % Lp = [Lp(:,i1) Lp(:,i2) Lp(:,i3)];
% % % Hp = [Hp(i4,i1) Hp(:,i2) Hp(i4,i3)];
% % % 
% % % surf(L,T,Hp), shading interp, pause;
% % 
% % % la_u = unique(la); th_u = unique(th); %find unique coords
% % 
% % %%% PAD WITH GHOST DATA, NEEDED FOR INTERPOLATION
% % % ia=find(la==la_u(end)|la==la_u(end-1)); %find indexes for points on right edge
% % % add shifted to h, la and th
% % % h = [h(ia);h];
% % % la = [la(ia)-2*pi;la];
% % % th = [th(ia);th];
% % 
% % % ia=find(la==la_u(1)); %find indexes for points on left edge
% % % ia=find(la==la_u(1)|la==la_u(2)); %find indexes for points on left edge
% % % add shifted to h, la and th
% % % h = [h(ia);h];
% % % la = [la(ia)+2*pi;la];
% % % th = [th(ia);th];
% % 
% % % ia = find(th==th_u(2));
% % % ia = find(th==th_u(1)|th==th_u(2)); %find indexes for points on bottom edge
% % % ia = find(th==th_u(2)|th==th_u(3)); %find indexes for points on bottom edge
% % % add shifted to h, la and th
% % % h = [h(ia);h];
% % % h = [h(ia);h];
% % % la = [la(ia);la];
% % % th = [-pi-th(ia);th];
% % 
% % % ia = find(th==th_u(end-1));
% % % ia = find(th==th_u(end)|th==th_u(end-1)); %find indexes for points on top edge
% % % ia = find(th==th_u(end-1)|th==th_u(end-2)); %find indexes for points on top edge
% % % add shifted to h, la and th
% % % h = [h(ia(end:-1:1));h];
% % % h = [h(ia);h];
% % % la = [la(ia);la];
% % % th = [pi-th(ia);th];
% % 
% % %%% SETUP THE NEW GRID, INTERSECT AND RESHAPE SOLUTION TO MATCH
% % % la_u = unique(la); th_u = unique(th); %find unique coords
% % % [L,T] = meshgrid(la_u,th_u);
% % % [c,ia,ib] = intersect([la th],[L(:) T(:)],'rows');
% % % h = h(ia);
% % % H = reshape(h,size(L));
% % % whos H
% % % surf(L,T,H); shading interp; pause;
% % 
% % %%% INTERPOLATE TO EVALUATION POINTS
% % % he = interp2(L,T,H,lae,the,'cubic');