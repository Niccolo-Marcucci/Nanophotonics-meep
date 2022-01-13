% clc
clear all
close all

addpath('Near-to-Far-Field-Transformation')

name = 'spiral_outcoupler_787f164995_design_TM_gd3_buriedDBR_onSiO2_positive_N5_charge1_D5000nm_simend01.0e-04-res66_220112-014717_nearfield_t30.0.mat';
name = 'polSplitter_ffd522835c_design_TM_gd3_buriedDBR_onSiO2_positive_N9_Dphi60_sigma1-res66_220112-104136_nearfield_t30.0.mat';

fields = load(['data/',name]);

volume = [fields.Lx, fields.Ly, 0]*1e-6;

Ex = fields.Ex;
Ey = fields.Ey;

Ex = Ex(2:end,:);
Ey = Ey(:,2:end);
% Ex = Ex(:,2:end);
% Ey = Ey(2:end,:);

N = size(Ex);
x = linspace(-1,1,N(1))*volume(1);
y = linspace(-1,1,N(2))*volume(2);

[X,Y] = meshgrid(y,x);

% imagesc(abs(Ex).^2+abs(Ey).^2)

c0 = physconst('lightspeed');

wavelength = 0.57*1e-6;

freq = c0/wavelength;


obj_near = Field2D(X, Y, Ex, Ey, freq, 'm');
obj_near.plotNearField()
% drawnow

farGridNumber = 1e2;

kx = linspace(-1, 1, farGridNumber)*1;
ky = linspace(-1, 1, farGridNumber)*1;

tic
obj_far = obj_near.getFarFieldCPU({kx, ky});
toc
% 
% figure
% obj_far.plotFarField(2,5,.2)
%%
Ex_far = obj_far.fx;
Ey_far = obj_far.fy;
% save([name(1:end-4),'_MAT.mat'],"Ex_far","Ey_far");
% 
% error
Ux = obj_far.kx;
Uy = obj_far.ky;
ux = linspace(-1,1,farGridNumber);
uy = linspace(-1,1,farGridNumber);

ER = sqrt(2)/2*Ex_far + sqrt(2)/2*Ey_far*exp(-1i*pi/2);
EL = sqrt(2)/2*Ex_far + sqrt(2)/2*Ey_far*exp(+1i*pi/2);

S3 = 1i*(Ex_far.*conj(Ey_far)-Ey_far.*conj(Ex_far));
S0 = (abs(Ex_far).^2+abs(Ey_far).^2);
chi = 0.5*asin( real(S3)./S0);

E = sqrt(real(Ex_far).^2+real(Ey_far).^2)+1i*sqrt(imag(Ex_far).^2+imag(Ey_far).^2);

% color maps
col_vec=linspace(0,1,256);
map_wave=[col_vec' col_vec' ones(256,1);
          ones(255,1) col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'];
map_wave_dark=[zeros(255,1) zeros(255,1) col_vec(end-1:-1:1)';
               col_vec' zeros(256,1) zeros(256,1)];
map_yell_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)' zeros(255,1) ;
               col_vec' zeros(256,1) zeros(256,1)];
map_oran_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'.^2 zeros(255,1) ;
               col_vec' col_vec([1:end/2,end/2:-1:1])' zeros(256,1)];
map_intensity=[ones(256,1) ones(256,1) col_vec(end:-1:1)';
               ones(255,1) col_vec(end-1:-1:1)' zeros(255,1)];
%%
fig=figure('units','centimeters','outerposition',[0 0 25 25]);
Rho = sqrt(Ux.^2+Uy.^2);
Theta = atan2(Uy,Ux);
rho = Rho(1:end-1,1:end-1);
theta = Theta(1:end-1,1:end-1);

PsiR = angle(ER);
dxPsiR = (diff(PsiR')./diff(Ux'))';
dxPsiR = dxPsiR(1:end-1,:);
dyPsiR = (diff(PsiR)./diff(Uy));
dyPsiR = dyPsiR(:,1:end-1);
mR = (dxPsiR .* rho .* sin(theta) - dyPsiR .* rho .* cos(theta)) ;
mR(mR>8) = NaN;
PsiL = angle(EL);
dxPsiL = (diff(PsiL')./diff(Ux'))';
dxPsiL = dxPsiL(1:end-1,:);
dyPsiL = (diff(PsiL)./diff(Uy));
dyPsiL = dyPsiL(:,1:end-1);
mL = (dxPsiL .* rho .* sin(theta) - dyPsiL .* rho .* cos(theta)) ;
mL(abs(mL)>8) = NaN;
%     plot_surf(rho,theta,(mL),'hsv',"Left circular polarization OAM","symmetric",5,0);
%     plot_surf(ux,uy,round(m),'hsv',"Right circular polarization phase","symmetric",5);
subplot(2,3,1)
plot_surf(ux,uy,abs(ER).^2,'hot',"Right circular polarization intensity");
subplot(2,3,4)
plot_surf(ux,uy,abs(EL).^2,'hot',"Left circular polarization intensity");
subplot(2,3,2)
plot_surf(ux,uy,angle(ER),'jet',"Right circular polarization OAM","symmetric",5);
subplot(2,3,5)
plot_surf(uy,ux,angle(EL),'jet',"Left circular polarization OAM","symmetric",5);


% rho_grid = linspace(min(min(rho)),max(max(rho)),1000);
% theta_grid = linspace(min(min(theta)),max(max(theta)),1000);
% [X,Y] = meshgrid(rho_grid,theta_grid);
% mR_fit = scatteredInterpolant(rho(:),theta(:),mR(:),'nearest');
% mR_grid = mR_fit(X,Y);
% mR_grid(isnan(mR_grid)) = 0;
% 
% mL_fit = scatteredInterpolant(rho(:),theta(:),mL(:),'nearest');
% mL_grid = mL_fit(X,Y);
% mL_grid(isnan(mL_grid)) = 0;
% 
% subplot(2,3,3)
% plot(rho_grid,sum(mR_grid)/length(theta_grid))
% %     plot_surf(rho,theta,round(mR),'hsv',"Right circular polarization OAM","symmetric",5,0);
% xlabel('\rho')
% ylabel('OAM')
% ylim([-5 0]);
% xlim([0, 0.2]);
% nicePlot
% subplot(2,3,6)
% plot(rho_grid,sum(mL_grid)/length(theta_grid))
% %     plot_surf(rho,theta,round(mL),'hsv',"Left circular polarization OAM","symmetric",5,0);
% xlabel('\rho')
% ylabel('OAM')
% ylim([-5 0]);
% xlim([0, 0.2]);
% nicePlot
function plot_surf(ux,uy,quantity,map,picture_title,symmetry,massimo,image)
    if nargin < 8 || image 
        imagesc(ux,uy,quantity);
    else
        x = linspace(min(min(ux)),max(max(ux)),1000);
        y = linspace(min(min(uy)),max(max(uy)),1000);
        [X,Y] = meshgrid(x,y);
        fit_IM = scatteredInterpolant(ux(:),uy(:),quantity(:),'nearest');
%         s=surface(ux,uy,quantity);
%         s.EdgeColor='none';
        imagesc(x,y,fit_IM(X,Y));

    end
    set(gca,'YDir','normal') 
    ax=gca;
    if nargin > 3
        colormap(ax,map)
    if nargin > 4
        title(picture_title)
    if nargin > 5 
        if symmetry == "symmetric"
            c = max(abs([min(min(quantity)),max(max(quantity))]));
            caxis([-c c])
        elseif symmetry == "unitary"
            caxis([0 1])
        else
            error('Wrong "symmetry" signment')
        end
    if nargin > 6
        if symmetry == "symmetric"
            caxis([-1 1]*massimo)
        elseif symmetry == "unitary"
            caxis([0 1]*massimo)
        else
            error('Wrong "symmetry" signment')
        end
    end;end;end;end
    colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
end
