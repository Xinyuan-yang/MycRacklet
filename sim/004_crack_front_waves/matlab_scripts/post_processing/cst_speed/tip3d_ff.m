% Author : Damien Spielmann, 2015
% This function first captures different tips from cracking index ST_Diagram type.
% At the break point, check if the capture is OK and press F5 to continue.
% Then it prints the evolution of propagation speed and energy dissipation.
% Modified by fekak fatima for super-shear mode II

close all;
clear all;
clc;

%% Variables to edit

sampling = 1; % Spatial jump between samples (i = i + sampling) 
l_contact = 0; % Is there contact in the left area ?
r_contact = 0; % Is there contact in the right area ?
rate = 0.2; %0.2 % Part of samples used in the polynomial approximation of CT_Diagram lines

%Start of the left and right observation zone. By default = 0.5
r_start = 0.0;
%% 
%CT_Diagram = load('ST_Diagram_id3d.dat');
CT_Diagram = load('ST_Diagram_id.cra');
fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);

nex = parameters{2}(1);
nez = parameters{2}(2);
nts = parameters{2}(3);
X = parameters{2}(5);
Z = parameters{2}(6);
%Z=X/4;
nut = parameters{2}(9);
nub = parameters{2}(10);

cst=parameters{2}(11);
csb=parameters{2}(12);

csbcst = csb/cst;

tau_n_cr = parameters{2}(15);
max_n_open = parameters{2}(16);

psi = parameters{2}(21);
phi = parameters{2}(22);

Gc = 0.5*tau_n_cr*max_n_open;


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
% %choose direction and line                              %
% direction=0;                                            %
% line=nez/2;                                             %
%                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cdt=sqrt(2*(1-nut)/(1-2*nut));%c'est deja divisé par cs+
% cdb=sqrt(2*(1-nub)/(1-2*nub))*csbcst;%gleichfalls
% 
% crt = (0.862+1.14*nut)/(1+nut);
% crb = (0.862+1.14*nub)/(1+nub)*csbcst;
% 
 ct_size = size(CT_Diagram);
% 
 dx = X/(nex-1);
% 
% space = linspace(0+dx,1-dx, nex);
 space = 0:dx:X;
 test = linspace(0,X, nex);
 time = load('Timer_ST_Diagram_id.cra'); %Cst/X
% 
% %Read the file and put it into 2d
% 
%     for h=0:ct_size(1)-1
%         for i=0:nex-1
%             ind(h+1,i+1)=CT_Diagram(h+1,line*nex+i+1);%+1 matlab indice
%         end
%     end
% 
% % Right side
% 
coh = 1;
cra = 1;
conts = 1;
conte = 1;
% CT_Diagram=ind;
%%

lim1 = 0.;
lim2 = 0.75 ;
i=1;
while (i <= nex)
    
  ind_coh = 0;
  ind_cra = 0;
  ind_conts = 0;
  ind_conte = 0;

    for t = 1:ct_size(1)-1
       
        if (((CT_Diagram(t,i)==0)&&(CT_Diagram(t+1,i)~=0))||((CT_Diagram(t,i)==5)&&(CT_Diagram(t+1,i)~=5)))
           if(space(i)>= lim1*X)&&(space(i)<= lim2*X)
            r_coh_zone(coh,1) = space(i);
            r_coh_zone(coh,2) = time(t);
            ind_coh = 1;
           end
                        
        end
        
        if (((CT_Diagram(t,i)==1)&&(CT_Diagram(t+1,i)~=1))||((CT_Diagram(t,i)==6)&&(CT_Diagram(t+1,i)~=6)))
           if(space(i)>= lim1*X)&&(space(i)<= lim2*X)
            r_cra_zone(cra,1) = space(i);
            r_cra_zone(cra,2) = time(t);
            ind_cra = 1;
           end             
        end
           
    end

    coh = coh + ind_coh;
    cra = cra + ind_cra;
     
    i = i+sampling;
    
end
%
save('r_coh_zone_list.mat','r_coh_zone');
save('r_cra_zone_list.mat','r_cra_zone');

%
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
plot(iterations(ind:end,1)/nex,iterations(ind:end,2),'--b','LineWidth',1.5);
xlabel('x/X','Fontsize',15)
ylabel('iterations','Fontsize',15)
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('iterations_VS_space.pdf',gcf,600);
saveas(gcf,'iterations_VS_space.fig');

%
figure(1);
xlabel('x/X','Fontsize',38);
ylabel('c_{s}t/X','Fontsize',38);

%title('If tips are correctly captured, press F5, else Shift-F5','Fontsize',23)

%hold on;
%axis([0 X]);

%plot(r_coh_zone(:,1)./X,r_coh_zone(:,2)./X,'*');
plot(r_cra_zone(:,1)./X,r_cra_zone(:,2),'*c');

%x = [0 X];
%cdtt = [+X/cdt -X/cdt+X/cdt];
%cdbb = [+X/cdb -X/cdb+X/cdb];

%%plot(x./Z,cdtt./Z,'--k');
%plot(x,cdbb,'--r');


%%
%influence
% yn=0;
% for i=1:length(r_coh_zone(:,1))
%     if r_coh_zone(i,2)>=-r_coh_zone(i,1)/cdt+X/cdt %%if the top material is decisive
%         xinfl=r_coh_zone(i,1);
%         tinfl=r_coh_zone(i,2);
%         yn=1;
%         break;
%     end
%     if r_coh_zone(i,2)>=-r_coh_zone(i,1)/cdb+X/cdb %%if the top bottom is decisive
%         xinfl=r_coh_zone(i,1);
%         tinfl=r_coh_zone(i,2);
%         yn=1;
%         break;
%     end
% end
% 
% if yn==1
%     tinflplot = [0 tinfl];
%     xinflplot = [0 xinfl];
%     plot(xinflplot./Z,[tinfl,tinfl]./Z,'--b');
%     plot([xinfl,xinfl]./Z,tinflplot./Z,'--b');
%     legend('Cohesive Tips','Crack Tips','c_{D}^{+} infl','D inf.','location','NorthWest');
% else
%     legend('Cohesive Tips','Crack Tips','location','Northeast');
% end
% set(gca,'Color', 'none','LineWidth',2,'Fontsize',32);
% legend('boxoff');
%keyboard;
%%


% Velocity computation
%dust c est ce qu on utilise pas car bool=0 on veut que les velocities

%
[r_coh]  = velociraptor(r_coh_zone,rate);
[r_cra]  = velociraptor(r_cra_zone,rate);


figure(2);
% x = [0 max(max(time))];
% csy = [csbcst csbcst];
% crby = [crb crb];
% sqrcsy = [sqrt(2) sqrt(2)]; 
% crty = [crt crt];

%axis([0,max(max(time)) 0 1]);

%Right velocities

%title('Tip velocities','Fontsize',23)

%hold on;
%plot(r_coh(:,2),r_coh(:,1),'b','Linewidth',2);
plot(r_cra(:,2),r_cra(:,1)./X,'* r','Linewidth',2);
xlabel('C_st/X','Fontsize',38);
ylabel('v/C_s','Fontsize',38);
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('velocity_VS_time.pdf',gcf,600);
saveas(gcf,'velocity_VS_time.fig');
%
figure(3);
% x = [0 max(max(time))];
% csy = [csbcst csbcst];
% crby = [crb crb];
% sqrcsy = [sqrt(2) sqrt(2)]; 
% crty = [crt crt];

%axis([0,max(max(time)) 0 1]);

%Right velocities

%title('Tip velocities','Fontsize',23)

%hold on;
%plot(r_coh(:,2),r_coh(:,1),'b','Linewidth',2);
plot(r_coh(:,2),r_coh(:,1)./X,'* r','Linewidth',2);
xlabel('C_st/X','Fontsize',38);
ylabel('v_{coh}/C_s','Fontsize',38);
str1 = sprintf('psi = %d & phi = %d',psi,phi);
%str2 = sprintf('Crack velocity : %d *Cs m/s',Vcrack);
title(str1,'Fontsize',15)
save2pdf('velocity_coh_VS_time.pdf',gcf,600);
saveas(gcf,'velocity_coh_VS_time.fig');



%plot(x./Z,[cst/cst cst/cst],'--r','Linewidth',2);
%plot(x./Z,csy,':r','Linewidth',2);
%plot(x./Z,crty,'--g','Linewidth',2);
%plot(x,crby,':g','Linewidth',2);

%% 
% %dilatational influence plot
% if yn==1
%     for i=1:length(r_coh(:,2))
%      if r_coh(i,2)==tinfl
%            indexi=i;
%            break;
%      end
%     end
%     vinflplot = [0 r_coh(indexi,1)];
%     %plot(tinflplot,[r_coh(indexi,1),r_coh(indexi,1)],'--b');
%     plot([tinfl,tinfl]./Z,vinflplot,'-.b');
%     legend('Crack front','c_{s}^{+}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}','c_{D} infl.','Location','SouthEast');%'southoutside','Orientation','horizontal');
% else
%     legend('Crack front','c_{s}^{+}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}','Location','SouthEast');%'southoutside','Orientation','horizontal');
% end
%     %legend('Cohesive front', 'Crack front','c_{s}^{-}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}', 'c_{R}^{-}/c_{s}^{+}',23,'Location','SouthEast');%'southoutside','Orientation','horizontal');
% xlabel('c_{s}^{+}t/Z','Fontsize',38);
% ylabel('v_{0}/c_{s}^{+}','Fontsize',38);
% %axis([0,max(max(time))./Z 0 max(mean(r_coh(:,1))*1.05,csbcst)])
% set(gca,'Color', 'none','LineWidth',2,'Fontsize',32);
% legend('boxoff');
% hold off;
% saveas(gcf,'vfront.fig');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3);
% x = [0 max(max(time))];
% csy = [csbcst csbcst];
% crby = [crb crb];
% sqrcsy = [sqrt(2) sqrt(2)]; 
% crty = [crt crt];
% 
% %axis([0,max(max(time)) 0 1]);
% 
% %Right velocities
% 
% %title('Tip velocities','Fontsize',23)
% 
% hold on;
% plot(r_coh(:,2)./Z,r_coh(:,1),'b','Linewidth',2);
% 
% plot(x./Z,[cst/cst cst/cst],'--r','Linewidth',2);
% %plot(x./Z,csy,':r','Linewidth',2);
% plot(x./Z,crty,'--g','Linewidth',2);
% %plot(x,crby,':g','Linewidth',2);
% 
% 
% %dilatational influence plot
% if yn==1
%     for i=1:length(r_coh(:,2))
%      if r_coh(i,2)==tinfl
%            indexi=i;
%            break;
%      end
%     end
%     vinflplot = [0 r_coh(indexi,1)];
%     %plot(tinflplot,[r_coh(indexi,1),r_coh(indexi,1)],'--b');
%     plot([tinfl,tinfl]./Z,vinflplot,'-.b');
%     legend('Cohesive front','c_{s}^{+}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}','c_{D} infl.','Location','SouthEast');%'southoutside','Orientation','horizontal');
% else
%     legend('Cohesive front','c_{s}^{+}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}','Location','SouthEast');%'southoutside','Orientation','horizontal');
% end
%     %legend('Cohesive front', 'Crack front','c_{s}^{-}/c_{s}^{+}','c_{R}^{+}/c_{s}^{+}', 'c_{R}^{-}/c_{s}^{+}',23,'Location','SouthEast');%'southoutside','Orientation','horizontal');
% xlabel('c_{s}^{+}t/Z','Fontsize',38);
% ylabel('v_{0}/c_{s}^{+}','Fontsize',38);
% %axis([0,max(max(time))./Z 0 max(mean(r_coh(:,1))*1.05,csbcst)])
% set(gca,'Color', 'none','LineWidth',2,'Fontsize',32);
% legend('boxoff');
% saveas(gcf,'vcohesive.fig');
%% 
