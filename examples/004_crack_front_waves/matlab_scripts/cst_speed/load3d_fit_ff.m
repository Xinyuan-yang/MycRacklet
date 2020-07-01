% Author : Damien Spielmann 2015
% This function plot the loading in time and the iteration from the files "loads.dat" and "literation.dat"
% No variable to edit
% Modified by Fatima Fekak

close all;
clear all;
clc;

%A supprimer
% x1 = linspace(0,1);
% x2 = linspace(3/4,1);
% y1 = sin(2*pi*x1);
% y2 = sin(2*pi*x2);
% 
% figure(1)
% 
% % plot on large axes
% plot(x1,y1)
% 
% % create smaller axes in top right, and plot on it
% axes('Position',[.6 .5 .3 .3])
% box on
% plot(x2,y2)


%ff: x, y and z loadings in time / or in space
loading = load('loading.cra');

% ff: simulation parameters
fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);
nele = parameters{2}(1);
dom_size = parameters{2}(5);
beta = parameters{2}(4);

E = parameters{2}(7);

nut = parameters{2}(9);
nub = parameters{2}(10);

cst = parameters{2}(11);
csb = parameters{2}(12);

delta_c = parameters{2}(13);

%csb = sqrt(2*(1-nut)/(1-2*nut))*cst; %cpt

%ff: equal to 1 for homogeneous media
csbcst = csb/cst;

cdt=sqrt(2*(1-nut)/(1-2*nut));%already divided by cs+
cdb=sqrt(2*(1-nub)/(1-2*nub))*csbcst;%""

crt = (0.862+1.14*nut)/(1+nut);
crb = (0.862+1.14*nub)/(1+nub)*csbcst;

psi = parameters{2}(21);
phi = parameters{2}(22);

a0 = parameters{2}(17);

n_str = parameters{2}(15);
s_str = parameters{2}(16);

%ff: number of time outputs / or number of elements
nts=loading(1);
%ff: crack velocity => toutes les vitesses de Cracklett sont par rapport à
%Cs
Vcrack = loading(2);

x_timer = load('Timer_ST_Diagram_id.cra'); %time=Cst/X 
x_space = 1:1:nts ; % number of elements

%Read the file
for h=0:nts-1      
   taux(h+1)=loading(3*h+1+2)/s_str;
   tauy(h+1)=loading(3*h+2+2)/n_str; 
   tauz(h+1)=loading(3*h+3+2)/s_str; 
end

% Morrissey and Rice formula
% % %loading mat files
% % load r_cra_zone_list.mat;
% % load r_coh_zone_list.mat;
% % % Theoritical loading
% % M = E/(2*(1-nut*nut));
% % v0=0.9; 
% % %const = 1/(1-(v0*(1/crt))); %already divided by cs+
% % const = 2.25;
% % %for i=1:size(x_tip,1)
% % for i=1:size(r_coh_zone,1)
% %     tau(i) = const * sqrt(M*delta_c*s_str/(r_coh_zone(i,1) * pi)) / s_str;    
% % end
% % 
% % %
% % figure(2);
% % plot(x_timer(1:end),taux(1:end),'o r','LineWidth',1.5);
% % hold on ;
% % %plot(x_timer(1:end),tau(1:end),'* b','LineWidth',1.5);
% % plot(r_coh_zone(1:end,2),tau(1:end),'* b','LineWidth',1.5);
% % xlabel('Cst/X','Fontsize',15)
% % ylabel('loading','Fontsize',15)
% % str1 = sprintf('psi = %d & phi = %d',psi,phi);
% % str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
% % title({str1; str2},'Fontsize',15)
% % legend('\tau_{x}/\tau^{str}_{s}','Location','northeast');
% % save2pdf('loading_VS_time.pdf',gcf,600);
% % saveas(gcf,'loading_VS_time.fig');


% 
%using Polynomial curve fitting
find = false;
ind = 1;
i=2;
while (find==false)
    if (psi== 0)
        if (tauy(i)==tauy(1))
            i=i+1;
        else
            find = true;
            ind = i;
        end
    else
        if (taux(i)==taux(1))
            i=i+1;
        else
            find = true;
            ind = i;
        end
    end        
end

%time control

% p = polyfit(x_timer(ind:end)',taux(ind:end),14);
% taux_fit = taux ;
% taux_fit(ind:end) = polyval(p,x_timer(ind:end));
% 
% figure(2);
% plot(x_timer,taux_fit,'--b','LineWidth',1.5);
% hold on;
% plot(x_timer,taux,'-r','LineWidth',1.5);
% xlabel('Cst/X','Fontsize',15)
% ylabel('loading','Fontsize',15)
% str1 = sprintf('psi = %d & phi = %d',psi,phi);
% str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
% title({str1; str2},'Fontsize',15)
% legend('\tau^{fit}_{x}/\tau^{str}_{s}','\tau_{x}/\tau^{str}_{s}','Location','northeast');
% save2pdf('loadingXFIT_VS_time.pdf',gcf,600);
% saveas(gcf,'loadingXFIT_VS_time.fig');
% 
% fileID = fopen('loading_fit.cra','w');
% fprintf(fileID,'%i\n',nts);
% fprintf(fileID,'%f\n',Vcrack);
% for i=1:nts
%     fprintf(fileID,'%f\n',taux_fit(i)*s_str);
%     fprintf(fileID,'%f\n',tauy(i)*n_str);
%     fprintf(fileID,'%f\n',tauz(i)*s_str);
% end
% fclose(fileID);

% space_control
propagation_domain1 = int16(0.7 * nts) ;

% fit de taux
% just for mode II and mixed mode 
px1 = polyfit(x_space(1:propagation_domain1),taux(1:propagation_domain1),5);
taux_fit = taux ;
%
propagation_domain2 = int16(0.7 * nts) ;
taux_fit(ind:propagation_domain2) = polyval(px1,x_space(ind:propagation_domain2));
%
figure(2);
plot(x_space(1:end)./nts,taux_fit(1:end),'* b','LineWidth',1.5);
hold on;
plot(x_space(1:end)./nts,taux(1:end),'o r','LineWidth',1.5);
xlabel('x/X','Fontsize',15)
ylabel('loading','Fontsize',15)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title({str1; str2},'Fontsize',15)
legend('\tau^{fit}_{x}/\tau^{str}_{s}','\tau_{x}/\tau^{str}_{s}','Location','northeast');
save2pdf('loadingXFIT_VS_space.pdf',gcf,600);
saveas(gcf,'loadingXFIT_VS_space.fig');

% fit de tauy
% just for mode I and mixed mode 
py1 = polyfit(x_space(ind:propagation_domain1),tauy(ind:propagation_domain1),5);
tauy_fit = tauy ;
%
tauy_fit(ind:propagation_domain2) = polyval(py1,x_space(ind:propagation_domain2));
%
figure(22);
plot(x_space(1:end)./nts,tauy_fit(1:end),'* b','LineWidth',1.5);
hold on;
plot(x_space(1:end)./nts,tauy(1:end),'o r','LineWidth',1.5);
xlabel('x/X','Fontsize',15)
ylabel('loading','Fontsize',15)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title({str1; str2},'Fontsize',15)
legend('\tau^{fit}_{y}/\tau^{str}_{n}','\tau_{y}/\tau^{str}_{n}','Location','northeast');
save2pdf('loadingYFIT_VS_space.pdf',gcf,600);
saveas(gcf,'loadingYFIT_VS_space.fig');

%
fileID = fopen('loading_fit.cra','w');
fprintf(fileID,'%i\n',nts);
fprintf(fileID,'%f\n',Vcrack);
for i=1:nts
    if (psi == 90)
    fprintf(fileID,'%f\n',taux_fit(i)*s_str);
    fprintf(fileID,'%f\n',tauy(i)*n_str);
    fprintf(fileID,'%f\n',tauz(i)*s_str);
    end
    if (psi == 0)
    fprintf(fileID,'%f\n',taux(i)*s_str);
    fprintf(fileID,'%f\n',tauy_fit(i)*n_str);
    fprintf(fileID,'%f\n',tauz(i)*s_str);       
    end
    if (psi ~= 0) && (psi ~= 90)
    fprintf(fileID,'%f\n',taux_fit(i)*s_str);
    fprintf(fileID,'%f\n',tauy(i)*n_str);
    fprintf(fileID,'%f\n',tauz(i)*s_str);
    end
    
end
fclose(fileID);

% % Mixed-mode
% fileID = fopen('loading_fit.cra','w');
% fprintf(fileID,'%i\n',nts);
% fprintf(fileID,'%f\n',Vcrack);
% for i=1:nts
%     fprintf(fileID,'%f\n',taux_fit(i)*s_str);
%     fprintf(fileID,'%f\n',tauy_fit(i)*n_str);
%     fprintf(fileID,'%f\n',tauz(i)*s_str);
% end
% fclose(fileID);


%x_tip as a function of time
x_tip = load('x_tip.cra');

figure(5);
plot(x_timer,x_tip/nele,'--b','LineWidth',1.5);
xlabel('Cst/X','Fontsize',15)
ylabel('x_{tip}/X','Fontsize',15)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title({str1; str2},'Fontsize',15)
save2pdf('xtip_VS_time.pdf',gcf,600);
saveas(gcf,'xtip_VS_time.fig');


% %loading as a function of x_tip
% figure(6);
% plot(x_tip/nele,taux_fit,'--b','LineWidth',1.5);
% hold on;
% plot(x_tip/nele,taux,'-r','LineWidth',1.5);
% xlabel('x/X','Fontsize',15)
% ylabel('loading','Fontsize',15)
% str1 = sprintf('psi = %d & phi = %d',psi,phi);
% str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
% title({str1; str2},'Fontsize',15)
% legend('\tau^{fit}_{x}/\tau^{str}_{s}','\tau_{x}/\tau^{str}_{s}','Location','northeast');
% save2pdf('loadingXFIT_VS_space.pdf',gcf,600);
% saveas(gcf,'loadingXFIT_VS_space.fig');
% %number of iteration as a function of x_tip
 iterations = load('iterations.cra');
%
find = false;
ind = 0;
i=1;
while (find==false)
    if (iterations(i,2)==0)
        i=i+1;        
    else         
        find = true;
        ind = i;
        
    end
end
%
figure(7);
plot(iterations(ind:end,1)/nele,iterations(ind:end,2),'--b','LineWidth',1.5);
xlabel('x/X','Fontsize',15)
ylabel('iterations','Fontsize',15)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title({str1; str2},'Fontsize',15)
save2pdf('iterations_VS_space.pdf',gcf,600);
saveas(gcf,'iterations_VS_space.fig');


%axis([0 maxy 0 cdt])

%title('load','Fontsize',23);

%Dilatational waves
% cdt=sqrt(2*(1-nut)/(1-2*nut));%c'est deja divisé par cs+
% cdb=sqrt(2*(1-nub)/(1-2*nub))*csbcst;%gleichfalls
% 
% xx = [0 maxy];
% csy = [csbcst csbcst];
% crby = [crb crb];
% sqrcsy = [sqrt(2) sqrt(2)]; 
% crty = [crt crt];

% plot(xx,[csbcst,csbcst],'--r')
% %plot(xx,[cdt,cdt],'--k');
% %plot(xx,[cdb,cdb],'--r');
% plot(xx,[crt,crt],':k');
% plot(xx,[crb,crb],':r');

% it = load('literation.dat');
% for h=0:nts-1      %nts-1 si c est pour le temps 
% 
%    itt(h+1)=it(h+1,1);
% end
% hold off
% hold on
% %plot(y,lo,'--b');%y si c'est pour le temps
% plot(x,itt,'-r','LineWidth',2);%y si c'est pour le temps
% 
% ylim([0 1.5]);
% axis([0 max(x)*1.2 0 1.5]);
% xlabel('c_{s}^{+}t/Z','Fontsize',38)
% ylabel('\tau /\tau_{c}','Fontsize',38);
% set(gca,'Color', 'none','LineWidth',2,'Fontsize',32);
% h_legend=legend('\tau_{x}','\tau_{y}','\tau_{z}','it/n','Location','northeast');
% %set(h_legend,'FontSize',40);
% legend('boxoff');

