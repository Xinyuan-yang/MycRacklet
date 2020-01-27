% Fatima Fekak

close all;
clear all;
clc;

% simulation parameters
fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);
nele = parameters{2}(1);
dom_size = parameters{2}(5);
beta = parameters{2}(4);

nut = parameters{2}(9);
nub = parameters{2}(10);

cst = parameters{2}(11);
csb = parameters{2}(12);
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

%
x_timer = load('Timer_ST_Diagram_id.cra'); %time=Cst/X 

%number of iteration as a function of x_tip
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
plot(iterations(ind:end,1)/nele,iterations(ind:end,2),'--b','LineWidth',2);
xlabel('x/X','Fontsize',23)
ylabel('iterations','Fontsize',23)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('iterations_VS_space.pdf',gcf,600); 
saveas(gcf,'iterations_VS_space.fig');

%computing v0 = vcrack/cs = 1/(it*beta) as a function of 

v0 = 1./(iterations(:,2)*beta);

figure(4);
plot(iterations(ind:end,1)/nele,v0(ind:end),'--b','LineWidth',2);
xlabel('x/X','Fontsize',23)
ylabel('v_{crack}/Cs','Fontsize',23)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('vcrack_VS_space.pdf',gcf,600);
saveas(gcf,'vcrack_VS_space.fig');

%x_tip as a function of time
x_tip = load('x_tip.cra');

figure(1);
plot(x_timer,x_tip/nele,'*','LineWidth',2);
xlabel('Cst/X','Fontsize',23)
ylabel('x_{tip}/X','Fontsize',23)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('x_tip_VS_time.pdf',gcf,600);
saveas(gcf,'x_tip_VS_time.fig');

%Cst/x = f(x_tip/X) => eauivqlent to ct_diag
figure(2);
plot(x_tip/nele,x_timer,'*','LineWidth',2);
ylabel('Cst/X','Fontsize',23)
xlabel('x_{tip}/X','Fontsize',23)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('time_VS_x_tip.pdf',gcf,600);
saveas(gcf,'time_VS_x_tip.fig');

