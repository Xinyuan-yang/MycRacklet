% Author : Fabian Barras, 2013
% This function draws CT-Diagram from input file CT_Diagram.dat
% No variable to edit
% Modified by Fatima Fekak

close all;
clear all;
clc;

%ff: noeuds indices in time
CT_Diagram = load('ST_Diagram_id.cra');
%fic = fopen('ST_Diagram_id.cra');
%CT_Diagram = textscan(fic, '%f');
%fclose(fic);
%CT_Diagram = load('ST_Diagram_shear_velo_jump.cra');

% ff: simulation parameters
fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);
nele = parameters{2}(1);
dom_size = parameters{2}(5);

nut = parameters{2}(9);
nub = parameters{2}(10);

cst = parameters{2}(11);
csb = parameters{2}(12);
%csb = sqrt(2*(1-nut)/(1-2*nut))*cst; %cpt

psi = parameters{2}(21);
phi = parameters{2}(22);

crt = (0.862+1.14*nut)/(1+nut)*cst;
crb = (0.862+1.14*nub)/(1+nub)*csb;

str = sprintf('ST-Diagram | psi = %d & phi = %d | Domain size : %d m',psi,phi,dom_size);

x = linspace(0,1,nele);
y = load('Timer_ST_Diagram_id.cra');
maxy = max(y);

figure(1)
imagesc (x,y,CT_Diagram);
set(gca,'LineWidth',2,'Fontsize',23);
xlabel('x/X','Fontsize',23)
ylabel('c_st/X','Fontsize',23)

title(str,'Fontsize',23)
axis xy;
save2pdf('st_diag.pdf',gcf,600); 
saveas(gcf,'st_diag.fig');

%displacement jumps

% normal_displacement_jumps = load('normal_displacement_jumps.cra');
% 
% figure(2)
% imagesc (x,y,normal_displacement_jumps);
% set(gca,'LineWidth',2,'Fontsize',23);
% xlabel('x/X','Fontsize',23)
% ylabel('c_st/X','Fontsize',23)
% 
% title(str,'Fontsize',23)
% axis xy;
% save2pdf('disp_jumps.pdf',gcf,600); 
% saveas(gcf,'disp_jumps.fig');

%velocity jumps

shear_velocity_jumps = load('ST_shear_velocity_jumps.cra');

figure(2)
imagesc (x,y,shear_velocity_jumps);

%To change de colormap limits
% colorbar
% lim = caxis;
% caxis([0 0.01])

set(gca,'LineWidth',2,'Fontsize',23);
xlabel('x/X','Fontsize',23)
ylabel('c_st/X','Fontsize',23)
title(str,'Fontsize',23)
axis xy;
save2pdf('shear_velocity_jumps.pdf',gcf,600); 
saveas(gcf,'shear_velocity_jumps.fig');




% hold on;
% 
% csx = [0.2 0.1];
% csy = [0.1*maxy 0.1*maxy+0.1];
% 
% plot(csx,csy,'r','Linewidth',4);
% 
% csx = [0.2 0.1];
% csy = [0.1*maxy 0.1*maxy+0.1*cst/csb];
% 
% plot(csx,csy,'g','Linewidth',4);
% 
% csx = [0.5 0.6];
% csy = [0.3*maxy 0.3*maxy+0.1*cst/crt];
% 
% plot(csx,csy,'k','Linewidth',4);
% 
% csx = [0.5 0.6];
% csy = [0.3*maxy 0.3*maxy+0.1*cst/crb];
% 
% plot(csx,csy,'c','Linewidth',4);
% 
% csx = [0.8 0.9];
% csy = [0.1*maxy 0.1*maxy+0.1];
% 
% %legend('c_{s}^{+}','c_{s}^{-}','c_{R}^{+}','c_{R}^{-}');
% 
% plot(csx,csy,'r','Linewidth',4);
% 
% csx = [0.8 0.9];
% csy = [0.1*maxy 0.1*maxy+0.1*cst/csb];
% 
% plot(csx,csy,'g','Linewidth',4);
% 
% csx = [0.5 0.4];
% csy = [0.3*maxy 0.3*maxy+0.1*cst/crt];
% 
% plot(csx,csy,'k','Linewidth',4);
% 
% csx = [0.5 0.4];
% csy = [0.3*maxy 0.3*maxy+0.1*cst/crb];
% 
% plot(csx,csy,'c','Linewidth',4);
