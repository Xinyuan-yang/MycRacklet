% Author : Damien Spielmann, 2015
% This function computes the crack front velocity fluctuations between
% the simulation without the asperity and the one with asperity

close all;
clear all;
%clearvars -except CT_Diagram CT_Diagram_model
clc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                       %
sampling = 1;       % Spatial jump between samples (i = i + sampling)                                                  %
samplingz = 1;      % Smooth the output pit on out of samplingz of the z lines                                         %
rate = 0.2;        % Part of samples used in the polynomial approximation of CT_Diagram lines. Initially: 0.2         %
display = 10;       % To display the number of FW (if display=1->all computed FW. If 2-> 1 all 2 FW will be plotted)   %
amplitude = 1.;    % Amplification factor of the pulse's Amplitude                                                    %
CO_CR = 0;          % 0 for cohesive or 1 for crack front                                                              %
spanx = 15;         %smoothing variable must be odd (default is 5)                                                     %
spanz = 15;         %smoothing variable must be odd (default is 5)                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%loading mat files
load r_cra_list.mat;
%
load r_cra_model_list.mat;
%
load drawAsp_list.mat;
%
load r_cra_zone_list.mat

%parameters = load('Parameters.dat');
 
fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);

nex = parameters{2}(1);
nez = parameters{2}(2);
nts = parameters{2}(3);0.


X = parameters{2}(5);
Z = parameters{2}(6);

%a0=parameters(12);

nut = parameters{2}(9);
nub = parameters{2}(10);

cst=parameters{2}(11);
csb=parameters{2}(12);

csbcst = csb/cst;

tau_n_cr = parameters{2}(15);
max_n_open = parameters{2}(16);

Gc = 0.5*tau_n_cr*max_n_open;
dGcrit=0.1*Gc;

cdt=sqrt(2*(1-nut)/(1-2*nut));%already divided by cs+
cdb=sqrt(2*(1-nub)/(1-2*nub))*csbcst;%""

crt = (0.862+1.14*nut)/(1+nut);%already divided by cs+
crb = (0.862+1.14*nub)/(1+nub)*csbcst;

dx = X/nex;
dz = Z/nez;

space = linspace(0,X, nex);
spacez = linspace(0,Z, nez/samplingz);
%
timer = load('Timer_ST_Diagram_id.cra');
nb_steps = size(timer,1);
start_step = round(nb_steps / 3) ; %time step at which we want to start the post-processing
end_step = nb_steps;   %time step at which we want to end the post-processing
time = timer(start_step:end_step,1);
ntss=size(time,1);
%
% % fic = fopen('ST_Diagram_id.cra', 'rb');
% % nb_bytes=4;
% % position = nex*nez*(start_step-1)*nb_bytes;
% % status = fseek(fic,position,'bof');
% % nb = nex *nez ;
% % ST_Diagram = fread(fic, [nb ntss], 'uint');
% % %position = ftell(fic);
% % fclose(fic);
% % CT_Diagram = ST_Diagram';
% % %ntss=size(CT_Diagram,1);%already multiplied by cs 
% % clear ST_Diagram;
% % %
% % 
% % if CO_CR==0
% %     a=5;b=0;
% % else
% %     a=6;b=1;
% % end
% % 
% % coh = 1;
% % cra1 = 2000000;
% % %Compute the crack front propagation velocity of the simulation with the asperity
% % for line=(nez/2-1):samplingz:nez-1
% %    %draw asperity
% %    for ii=1:nex
% %        if CT_Diagram(1,line*nex+ii)==5 %(CT_Diagram(1,line*nex+ii)==5)&&(CT_Diagram(1,line*nex+ii)~=5)%
% %             
% %           drawAsp(coh,1) = space(ii); 
% %           drawAsp(coh,2) = spacez(line/samplingz+1);
% %           coh=coh+1;
% %        end
% %    end
% % 
% % cra = 1;
% % i=1;
% % 
% % while (i <= nex)
% % 
% %   
% % %  if CT_Diagram(1,line*nex+i)==5 %(CT_Diagram(1,line*nex+ii)==5)&&(CT_Diagram(1,line*nex+ii)~=5)%
% % %             
% % %           drawAsp(coh,1) = space(i)/X; 
% % %           drawAsp(coh,2) = spacez(line+1)/X;
% % %           coh=coh+1;
% % %         end
% %  
% % 
% %     for t = 1:ntss-1
% %         
% % 
% %         if ((((CT_Diagram(t,line*nex+i)==a)&&(CT_Diagram(t+1,line*nex+i)~=a))||((CT_Diagram(t,line*nex+i)==b)&&(CT_Diagram(t+1,line*nex+i)~=b)))&&(cra<cra1))%(CT_Diagram(t,i)==1)||(CT_Diagram(t,i)==6)%
% %            
% %             r_cra_zone(cra,1) = space(i);
% %             r_cra_zone(cra,2) = time(t);
% %             cra = cra + 1; 
% %         end
% %            
% %     end
% % 
% % 
% %      
% %     i = i+sampling;
% %     
% % end
% % 
% % cra1=cra;
% %  
% % [r_cra((line-(nez/2-1))/samplingz+1,:,:)]  = velociraptor_xz(r_cra_zone,rate);
% % end
% % save('r_cra_list.mat','r_cra');
% % save('drawAsp_list.mat','drawAsp');
% % 
% % clear CT_Diagram;
% % clear r_cra_zone;
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Compute the crack front propagation velocity of the simulation without the asperity (model)
% % CT_Diagram_model_all = load('ST_Diagram_id_model.cra');
% % CT_Diagram_model = CT_Diagram_model_all(start_step:end_step,:);
% % 
% % %
% % cra2 = 1;
% % i=1;
% % while (i <= nex)
% %     
% % 
% %     for t = 1:ntss-1
% %         
% %         %%needed to draw the influence line on cohesive front
% % %          if (((CT_Diagram_model(t,i)==0)&&(CT_Diagram_model(t+1,i)~=0))||((CT_Diagram_model(t,i)==5)&&(CT_Diagram_model(t+1,i)~=5)))
% % %            
% % %             r_coh_zone(coh,1) = space(i);
% % %             r_coh_zone(coh,2) = time(t);X = parameters{2}(5);
% % 
% % %             ind_coh = 1;
% % %         end
% %         
% %         if (((CT_Diagram_model(t,i)==a)&&(CT_Diagram_model(t+1,i)~=a))||((CT_Diagram_model(t,i)==b)&&(CT_Diagram_model(t+1,i)~=b))&&(cra2<cra))%(CT_Diagram_model(t,i)==1)||(CT_Diagram_model(t,i)==6)%
% %            
% %             r_cra_zone(cra2,1) = space(i);
% %             r_cra_zone(cra2,2) = time(t);
% %             cra2 = cra2 + 1;
% %         end
% %            
% %     end
% % 
% %     i = i+sampling;
% %     
% % end
% % clear CT_Diagram_model;
% % 
% % [r_cra_model]  = velociraptor_xz(r_cra_zone,rate);
% % cracracksize=size(r_cra_model);
% % save('r_cra_model_list.mat','r_cra_model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
x0Asp=min(drawAsp(:,1));%plot(drawAsp(:,1)./Z,drawAsp(:,2)./Z,'*k','Linewidth',1);

x1Asp=max(drawAsp(:,1));

%cracracksize=size(r_cra_model);
cracracksize=size(r_cra);

%positions in x of the begining and the end of FW post-treatment
%x_start = 0.4995 ; %x/X
%x_start = 0.3187 ;
x_start = 0.2 ;
ind_start = 1;
for i=1:size(r_cra_zone,1)
    if (abs((r_cra_zone(i,1) - (x_start * X))) < 1.e-4)
        ind_start = i;
        break ;
    end
end
%
%x_end = 0.4076 ;
x_end =  0.4457;
ind_end = size(r_cra_zone,1);
for i=1:size(r_cra_zone,1)
    if (abs((r_cra_zone(i,1) - (x_end * X))) < 1.e-4)
        ind_end = i;
        break ;
    end
end
%
figure(1);
hold on;
spacez2 = linspace(Z/2,Z, nez/2/samplingz);
%compute the pulse G (velocity fluctuation between the crack front propagation velocity of the simulation without the asperity and with the asperity)
for line=1:nez/samplingz/2
      %G(line,:)= smooth((amplitude).*(squeeze(r_cra(line,1:size(r_cra_model,1),1))'-squeeze(r_cra_model(:,1)))./ X + squeeze(r_cra(1,1:size(r_cra_model,1),2))'./X,spanx);%-(amplitude).*((-1/(crt))./((1-squeeze(r_cra_model(:,1))/(crt)).^2)).*(Gc*(squeeze(r_cra(line,:,1))'-squeeze(r_cra_model(:,1)))).*(1/dGcrit)+squeeze(r_cra_model(:,2))/X;
      G(line,:)= smooth((amplitude).*(squeeze(r_cra(line,ind_start:ind_end,1))'-squeeze(r_cra_model(ind_start:ind_end,1)))./ X + squeeze(r_cra(1,ind_start:ind_end,2))'./X,spanx);       
      G2(line,:)=squeeze(r_cra_model(ind_start:ind_end,2))./X; %lines representing the crack front without the asperity (no fluctuations)
end

if CO_CR==0
    for i=1:0.95*size(G,2)/display
        GG(:,i*display) = smooth(G(:,i*display),spanz);
        %plot(G(:,i*display),spacez./Z,'b');
            plot(GG(:,i*display),spacez2./Z,'b');
            plot(G2(:,i*display),spacez2./Z,':k');
    end
else
    for i=1:0.95*size(G,2)/display
        GG(:,i*display) = smooth(G(:,i*display));
        plot(GG(:,i*display),spacez2./Z,'c');
        plot(G2(:,i*display),spacez2./Z,':k');
    end
end

for i=1:cracracksize(1)
    diff(i)=-min(r_cra(:,i,1)-r_cra_model(i,1)) ./ X;
end

%plot asperity
plot(drawAsp(:,1)./X,drawAsp(:,2)./Z,'.k','Linewidth',1);

%Rayleigh speed

% xn=max(r_cra_zone(:,1));
% t0=min(drawAsp(:,2));
% tn=max(r_cra_zone(:,2));

x0=min(r_cra_zone(:,1));
xn=max(r_cra_zone(:,1));
t0=min(r_cra_zone(:,2));%le temps est deja multiplier par cs - ff: est divisé par X
tn=max(r_cra_zone(:,2));

%damien
% v_long=(xn-x0)/(tn-t0); %attention it's actually not true since at the beginning (before the asperity) the propagation is faster than after the propagation. So this is a mean value
% %Pulse's velocity along the Crack Front to propagate at CR relative to the
% %asperityr_cra
% v_transt=sqrt(crt*crt-v_long*v_long);
% v_transb=sqrt(crb*crb-v_long*v_long);
% %time needed to reach the asperity%plot(drawAsp(:,1)./Z,drawAsp(:,2)./Z,'*k','Linewidth',1);
% t0Asp=t0+(x0Asp-x0)/v_long;blanchisserie epfl
% %position on the crack front where it should be if travelling at CR
% zncft=v_transt*(tn-t0Asp);
% zncfb=v_transb*(tn-t0Asp);
% zR = [Z/2 Z/2+zncft];
% zR2 = [Z/2 Z/2-zncfb];
% xx = [0+x0Asp xn];

%fatima
%v_long=0.3 ;
v_long=(xn-x0)./(tn-t0)./X;
v_transt=sqrt(crt*crt-v_long*v_long);
zn = (xn-x0Asp) * (v_transt ./ v_long) + Z/2;
zR = [Z/2 zn];
xx = [x0Asp xn];

crplott=plot(xx./X,zR./Z,'--g','Linewidth',1);
%crplotb=plot(xx./Z,zR2./Z,':g','Linewidth',2);

%Shear wave speed

%damien
% vs_transt=sqrt(1-v_long*v_long);3284
% v_transb=sqrt(csbcst*csbcst-v_long*v_long);
% zncft=vs_transt*(tn-t0Asp);
% zncfb=v_transb*(tn-t0Asp);
% zS = [Z/2 Z/2+zncft];
% zS2 = [Z/2 Z/2-zncfb];
% xx = [0+x0Asp xn];

%fatima
vs_transt=sqrt(1-v_long*v_long);
zn = (xn-x0Asp) * (vs_transt ./ v_long) + Z/2;
zS= [Z/2 zn];
xx = [x0Asp xn];

csplott=plot(xx./X,zS./Z,'--r','Linewidth',1);
%csplotb=plot(xx./Z,zS2./Z,':r','Linewidth',2);


%Dilatational wave speed

%damien
% v_transt=sqrt(cdt*cdt-v_long*v_long);
% v_transb=sqrt(cdb*cdb-v_long*v_long);
% zncft=v_transt*(tn-t0Asp);
% zncfb=v_transb*(tn-t0Asp);
% zD = [Z/2 Z/2+zncft];
% zD2 = [Z/2 Z/2-zncfb];
% xx = [0+x0Asp xn];

%fatima
vd_transt=sqrt(cdt*cdt-v_long*v_long);
zn = (xn-x0Asp) * (vd_transt ./ v_long) + Z/2;
zD= [Z/2 zn];
xx = [x0Asp xn];

cdplott=plot(xx./X,zD./Z,'--k','Linewidth',1);
%cdplotb=plot(xx./Z,zD2./Z,':k','Linewidth',2);

% %Dilatational influence plot
% dil_infl_index=0;
% for i=1:length(r_cra_zone(:,1))
%     if r_cra_zone(i,2)>=-r_cra_zone(i,1)/cdb+X/cdb %%if the bottom is decisive
%         dxinfl=r_cra_zone(i,1);
%         dil_infl_index=i;
%         break;
%     end
%     if r_cra_zone(i,2)>=-r_cra_zone(i,1)/cdt+X/cdt %%if the top material is decisive but the top is always set as the stiffener so it will be the most decisive. Darum at the second position in the loop
%         dxinfl=r_cra_zone(i,1);
%         dil_infl_index=i;%plot(drawAsp(:,1)./Z,drawAsp(:,2)./Z,'*k','Linewidth',1);
% 
%         break;
%     end
% 
% end
% %plot(r_cra_zone(:,1),diff,':b');



% %Shear wave influence plot
% in_domain=0;
% for i=1:length(r_cra_zone(:,1))
%     if r_cra_zone(i,2)>=-r_cra_zone(i,1)/(csb/cst)+X/(csb/cst) 
%         sxinfl=r_cra_zone(i,1);
%         break;
%     end
%     if r_cra_zone(i,2)>=-r_cra_zone(i,1)/(cst/cst)+X/(cst/cst) 
%         sxinfl=r_cra_zone(i,1);
%         break;
%         in_domain=1;
%     end
% 
% end
% 
% if (in_domain==1)&&(dil_infl_index~=0)
%     Dil_inf=plot([dxinfl/Z,dxinfl/Z],[0,Z/Z],'-.k','Linewidth',2);
%     Shr_inf=plot([sxinfl/Z,sxinfl/Z],[0,Z/Z],'-.r','Linewidth',2);
%     h_legend=legend([crplott,csplott,cdplott,Dil_inf,Shr_inf],'c_{R}^{+}','c_{S}^{+}','c_{D}^{+}','c_{D}^{+} influence','c_{s}^{+} influence','location','northeast');
% elseif (in_domain==0)&&(dil_infl_index~=0)
%     Dil_inf=plot([dxinfl/Z,dxinfl/Z],[0,Z/Z],'-.k','Linewidth',2);
%     h_legend=legend([crplott,csplott,cdplott,Dil_inf],'c_{R}^{+}','c_{S}^{+}','c_{D}^{+}','c_{D}^{+} influence','location','northeast');
% 
% else
%     h_legend=legend([crplott,csplott,cdplott],'c_{R}^{+}','c_{S}^{+}','c_{D}^{+}','location','northeast');
% end


%h_legend=legend([crplott,csplott,cdplott],'$c_{R}$','$c_{S}$','$c_{D}$','location','best');
h_legend=legend([cdplott],'$c_{D}$','location','best');
set(h_legend,'Interpreter','Latex');
set(h_legend,'FontSize',15);
%axis equal;
%axis([0.25 X/2/X Z/Z/2 Z/Z]);
axis([x0Asp/X x_end Z/Z/2 Z/Z]);
%axis([0 x_end Z/Z/2 Z/Z]);
ylabel('$\frac{z}{L_z}$','Fontsize',15,'interpreter','Latex')
xlabel('$\frac{(v-v_{0})}{c_{S}} + const. \frac{v_{0}t}{L_x}$','Fontsize',15,'interpreter','Latex');%\Delta G_{rest}/\Delta G_{crit}+v_{0}*t'\
legend('boxoff');
%set(gca, 'Color', 'none','Fontsize',24);
saveas(gcf,'CFW_positions.fig');
save2pdf('CFW_positions.pdf',gcf,600);
%%
%export_fig test -tiff -png -m2.5 -transparent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%plot the amplitude normalized with cs
figure(2)
hold on
plot(r_cra_zone(ind_start:ind_end,1)./X,smooth(diff(ind_start:ind_end),spanx),'-b','Linewidth',1.5);
%plot(r_cra_zone(1:end_post,1)./X,diff(1:end_post),':b','Linewidth',2);
%plot(r_cra_zone(:,1)./Z,diff,':b','Linewidth',2);
% if (dil_infl_index~=0)
%     Dil_inf=plot([dxinfl/Z,dxinfl/Z],[0,max(diff)],'-.k','Linewidth',2);
% end
ylabel('$\frac{\max(\Delta v)}{c_{S}}$','Fontsize',15,'interpreter','Latex');%V=(v-vo)
xlabel('$\frac{v_{0}t}{L_x}$','Fontsize',15,'interpreter','Latex');
%axis([0 (xn+0.25)/Z 0 max(diff)]);
% if (dil_infl_index~=0)
%     legend([Dil_inf],'c_{D}^{+} influence','location','northeast');
%     legend('boxoff');
% end
%set(gca, 'Color', 'none','Fontsize',24); % Sets axes background
saveas(gcf,'CFW_amplitude.fig');
save2pdf('CFW_amplitude.pdf',gcf,600);
%%
% figure(3)
% %loglog(r_cra_zone(ind_start:ind_end,1)./X,diff(ind_start:ind_end),':b','Linewidth',2);
% X_TEST = r_cra_zone(ind_start:ind_end,1)./X;
% Y_TEST = smooth(diff(ind_start:ind_end),spanx);
% plot(log(X_TEST),log(Y_TEST)) ;
% % hold on;
% % loglog(X_TEST, X_TEST .^(-0.5), 'r');
% %loglog(r_cra_zone(1:end_post,1)./X,smooth(diff(1:end_post),spanx),':b','Linewidth',2);
% grid on;
% saveas(gcf,'CFW_log_amplitude.fig');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %plot the amplitude normalized with v0
% figure(3)
% hold on
% plot(r_cra_zone(:,1)./Z,smooth(diff./squeeze(r_cra_model(:,1))',spanx),':b','Linewidth',2);
% %plot(r_cra_zone(:,1)./Z,diff,':b','Linewidth',2);
% if (dil_infl_index~=0)
%     Dil_inf=plot([dxinfl/Z,dxinfl/Z],[0,max(diff./squeeze(r_cra_model(:,1))')],'-.k','Linewidth',2);
% end
% ylabel('\Delta v^{max}/v_{0}','Fontsize',26);%V=(v-vo)
% xlabel('v_{0}t/Z','Fontsize',26);
% axis([0 (xn+0.25)/Z 0 max(diff./squeeze(r_cra_model(:,1))')]);
% if (dil_infl_index~=0)
%     legend([Dil_inf],'c_{D}^{+} influence','location','northeast');
%     legend('boxoff');r_cra_zone
% end
% set(gca, 'Color', 'none','Fontsize',24); % Sets axes background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %plot the amplitude for a fix z
% figure(5)
% hold on
% chosen_line=0.4;
% lini=1;%=int64((0.5-chosen_line)*nez);
% if (dil_infl_index~=0)
%     amp=(squeeze(r_cra(lini,1:dil_infl_index))'-squeeze(r_cra_model(1:dil_infl_index,1)))./Z;
%     plot(r_cra_zone(1:dil_infl_index,1)./Z,amp,':b','Linewidth',2);
% else
%     amp=(squeeze(r_cra(lini,1:size(r_cra_model,1)))'-squeeze(r_cra_model(:,1)))./Z;
%     plot(r_cra_zone(:,1)./Z,amp,':b','Linewidth',2);
% end
% ylabel('\Delta v/c_{S}^{+}','Fontsize',26);
% xlabel('v_{0}t/Z','Fontsize',26);
% if (dil_infl_index~=0)
%     axis([x1Asp/Z (dxinfl*1.12)/Z 0 max(amp)]);
% else
%     axis([x1Asp/Z (xn)/Z 0 max(amp)]);
% end
% legend(sprintf('Amplitude at z=%0.2fZ',chosen_line),'location','northeast');
% legend('boxoff');
% set(gca, 'Color', 'none','Fontsize',24); % Sets axes background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%compute and plot front waves velocity data and interpolation until the Cd
%influence
%compute Front Wave velocity
figure(6)
hold on
dil_infl_index = 0;
t1Asp=t0+(max(drawAsp(:,1))-x0)/v_long;%time when the CF has reached the end of the asperity
t_ind_start = r_cra_zone(ind_start,2);
t_start = max(t1Asp,t_ind_start);

cc=1;
maxline=0;
if (dil_infl_index~=0)
    for i=ind_start:ind_end
      Gmin=min(G(:,i-ind_start+1));
      for line=1:nez/samplingz/2
             GL=G(line,i-ind_start+1);
             if ((GL == Gmin)&&(r_cra_zone(i,1)>=x1Asp)&&(r_cra_zone(i,1)<dxinfl)&&(line>maxline)) %(Gmin<Z/2/vs_transt*v_long/Z+min(drawAsp(:,1))/Z))%%the point for the computation of the front wave velocity are taken between the asperity and the influence of the dilatational line
             GGmin(cc,1) = r_cra_zone(i,1)-x1Asp;
             GGmin(cc,2) = spacez(line);%Z/2-
             GGmin(cc,3) = r_cra_zone(i,2)-t_start; 
             %if line>0.85*nez/2
                    maxline=line;
             %end
             cc=cc+1;
             break;%ATTENTION avec ce break seulement le pic le plus loin de l axe x=Z/2 est pris
             end
      end
    end
else
     for i=ind_start:ind_end
      Gmin=min(G(:,i-ind_start+1));      
      for line=1:(nez/samplingz/2)
             GL=G(line,i-ind_start+1);
             if ((GL == Gmin)&&(r_cra_zone(i,1)>=x1Asp)&&(line>maxline)) %(Gmin<Z/2/vs_transt*v_long/Z+min(drawAsp(:,1))/Z))%%             
             GGmin(cc,1) = r_cra_zone(i,1);
             GGmin(cc,2) = spacez(line);%Z/2-
             %GGmin(cc,3) = (r_cra_zone(i,2)-t_start) * X; 
             GGmin(cc,3) = r_cra_zone(i,2) * X; 
             %if line>0.85*nez/2
                    maxline=line;
             %end
             cc=cc+1;
             break;
             end
      end
    end
end
%Based on position
% popo = polyfit((GGmin(:,1)),GGmin(:,2),1);
% vfwx=popo(1)*v_long; %tiangle semblable
% vfwx_asp=sqrt(vfwx*vfwx+v_long*v_long);
%
% xnfwx=vfwx*(tn-t0Asp);
% z = [Z/2 Z/2+xnfwx];
% xx = [0+x0Asp xn];
% vfwxplot=plot(xx./Z,z./Z,':y','Linewidth',2);


%Based on mean values
% %v_fw is the front wave velocity in the z direction
% vfwm_all=(GGmin(:,2))./(GGmin(:,3)); 
% vfwm = mean(vfwm_all);
% %v_fw_asp is the front wave velocity / asperity
% vfwm_asp=sqrt(v_long*v_long+vfwm*vfwm);
% 
% vfwm_asp/crt
%
% xnfwm=vfwm*(tn-t0Asp);
% z = [Z/2 Z/2+xnfwm];
% xx = [0+x0Asp xn];
% vfwmplot=plot(xx./Z,z./Z,':g','Linewidth',2);

% %Based on time
% popov = polyfit((GGmin(:,3)),GGmin(:,2),1);
% v_long2=polyfit(r_cra_zone(:,2),r_cra_zone(:,1),1);
% v_long2=v_long2(1);
% vfwv = popov(1);
% vfwv_asp=sqrt(vfwv*vfwv+v_long*v_long);
% %vfwcrt=vfwv_asp/crt;
% vfwcrt=vfwv/(sqrt(crt*crt-v_long*v_long));
% 
% xnfwv=vfwv*(tn-t0Asp);
% z = [Z/2 Z/2+xnfwv];
% xx = [0+x0Asp xn];
% vfwvplot=plot(xx./Z,z./Z,'m','Linewidth',1);

%linear
[popov1,S1] = polyfit((GGmin(:,3)),GGmin(:,2),1);
[y1,delta1]=polyval(popov1,(GGmin(:,3)),S1);
appro1=plot(GGmin(:,1)./X,y1./Z+Z/Z/2,'--r','Linewidth',1);
[popovmax] = polyfit((GGmin(:,3)),y1+delta1,1);
[ymax]=polyval(popovmax,GGmin(:,3));
[popovmin] = polyfit((GGmin(:,3)),y1-delta1,1);
[ymin]=polyval(popovmin,GGmin(:,3));
firstcoeffmax=(ymax(size(GGmin,1))-y1(1))/(GGmin(size(GGmin,1),3)-GGmin(1,3));
firstcoeffmin=(ymin(size(GGmin,1))-y1(1))/(GGmin(size(GGmin,1),3)-GGmin(1,3));
% ymax=y1+delta1;
% ymin=y1-delta1;

% [y2]=polyval([firstcoeffmax popov1(2)],(GGmin(:,3)));
% appro2=plot(GGmin(:,1)./Z+x1Asp/Z,y2./Z+Z/Z/2,':k','Linewidth',1);
% [y3]=polyval([firstcoeffmin popov1(2)],(GGmin(:,3)));
% appro3=plot(GGmin(:,1)./Z+x1Asp/Z,y3./Z+Z/Z/2,':k','Linewidth',1);

delta1=(firstcoeffmax-firstcoeffmin)/4;

%par rapport a l axe x
%vfwcrt1=popov1(1)/sqrt(crt^2-v_long^2);
%par rapport a l asperity
vfwAcrt1=sqrt(v_long^2+popov1(1)^2)/crt;
%deltaA=sqrt(v_long^2+(popov1(1)+delta1)^2)/crt-vfwAcrt1;
vfwAcdt1=sqrt(v_long^2+popov1(1)^2)/cdt;
deltaA=sqrt(v_long^2+(popov1(1)+delta1)^2)/cdt-vfwAcdt1;
%Z component
vfwzcrt1=popov1(1)/crt;
%deltaz=delta1/crt;
vfwzcdt1=popov1(1)/cdt;
deltaz=delta1/cdt;

%Asperity
plot(drawAsp(:,1)./X,drawAsp(:,2)./Z,'.k','Linewidth',1);

%cubic
% [popov3,S3] = polyfit((GGmin(:,3)),GGmin(:,2),3);
% [y3,delta3]=polyval(popov3,(GGmin(:,3)),S3);
% appro3=plot(GGmin(:,1)./Z+x1Asp/Z,y3./Z+Z/Z/2,'-g','Linewidth',1);
% %par rapport a l asperity
% %vfwcrt3=(GGmin(0.15*(cc-1):(cc-1),3).*GGmin(0.15*(cc-1):(cc-1),3).*3.*popov3(1)+GGmin(0.15*(cc-1):(cc-1),3).*2.*popov3(2)+popov3(3))/(sqrt(crt*crt-v_long*v_long));
% %par rapport a l asperity
% vfwAcrt3=sqrt(v_long^2+(GGmin(:,3).*GGmin(:,3).*3.*popov3(1)+GGmin(:,3).*2.*popov3(2)+popov3(3)).^2)/crt;
% vfwAcrt3max=max(vfwAcrt3);
% vfwAcrt3min=min(vfwAcrt3);
% vfwAcrt3mean=mean(vfwAcrt3);
% 
% vfwAcdt3=sqrt(v_long^2+(GGmin(:,3).*GGmin(:,3).*3.*popov3(1)+GGmin(:,3).*2.*popov3(2)+popov3(3)).^2)/cdt;
% vfwAcdt3max=max(vfwAcdt3);
% vfwAcdt3min=min(vfwAcdt3);
% vfwAcdt3mean=mean(vfwAcdt3);

data=plot(GGmin(:,1)./X,GGmin(:,2)./Z+Z/Z/2,'xb','Linewidth',1);  

%Rayleigh wave speed
crplott=plot(xx./X,zR./Z,'--g','Linewidth',2);

%Shear wave speed
csplott=plot(xx./X,zS./Z,'--r','Linewidth',2);

%Dilatational wave speed
cdplott=plot(xx./X,zD./Z,'--k','Linewidth',2);

%Dilatational influence plot
% if (dil_infl_index~=0)
%     Dil_inf=plot([dxinfl/Z,dxinfl/Z],[Z/Z/2,Z/Z],'-.k','Linewidth',1);
% end

%axis equal;
ylabel('$\frac{z}{L_z}$','Fontsize',15,'interpreter','Latex');
xlabel('$\frac{v_{0}t}{L_x}$','Fontsize',15,'interpreter','Latex');
%axis([x0Asp/X X/2/X Z/Z/2 Z/Z]);
axis([x0Asp/X x_end Z/Z/2 Z/Z]);
if (dil_infl_index~=0)
    h_legend=legend([crplott,csplott,cdplott,Dil_inf,data,appro1],'c_{R}^{+}','c_{S}^{+}','c_{D}^{+}','c_{D}^{+} influence','front wave position',sprintf('v_{FWz}/c_{R}^{+}=%0.3f ± %0.3f \nv_{FW}/c_{R}^{+}=%0.3f ± %0.3f',vfwzcrt1,deltaz,vfwAcrt1,deltaA),'location','northeast');
else
    h_legend=legend([crplott,csplott,cdplott,data,appro1],'$c_{R}$','$c_{S}$','$c_{D}$','front wave position',sprintf('$c/c_{R}=%0.3f \\pm %0.3f$ \n$c_f/c_{R}=%0.3f \\pm %0.3f$',vfwzcrt1,deltaz,vfwAcrt1,deltaA),'location','best');
    %h_legend=legend([cdplott,data,appro1],'$c_{D}$','front wave position',sprintf('$c/c_{D}=%0.3f \\pm %0.3f$ \n$c_f/c_{D}=%0.3f \\pm %0.3f$',vfwzcdt1,deltaz,vfwAcdt1,deltaA),'location','best');
end
legend('boxoff');
set(h_legend,'Interpreter','Latex');
set(h_legend,'FontSize',15);
%set(gca, 'Color', 'none','Fontsize',24); % Sets axes background
saveas(gcf,'CFW_velocity.fig');
save2pdf('CFW_velocity.pdf',gcf,600);
%,sprintf('v_{FW} v_{FW}/c_{R}^{+}=%0.2f',vfwcrt)
%linearSD3=100*mean(delta1)
%cubicSD3=100*mean(delta3)
%,sprintf('v_{FW} v_{FW}/c_{R}^{+}=%0.2f',vfwcrt)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute and plot front waves velocity data and interpolation until it
%reaches the duplicated domain
%compute Front Wave velocity
% figure(4)
% hold on
% 
% t1Asp=t0+(max(drawAsp(:,1))-x0)/v_long;%time when the CF has reached the end of the asperity
% 
% cc=1;
% maxline=0;
% for i=1:cracracksize(1)
%     Gmin=min(G(:,i));
%     for line=1:nez/samplingz/2
%         GL=G(line,i);
%         if ((GL == Gmin)&&(r_cra_zone(i,1)>=x1Asp)&&(line>maxline)) %(Gmin<Z/2/vs_transt*v_long/Z+min(drawAsp(:,1))/Z))%%the point for the computation of the front wave velocity are taken between the asperity and the influence of the dilatational line
%             GGmin(cc,1) = r_cra_zone(i,1)-x1Asp;
%             GGmin(cc,2) = spacez(line);
%             GGmin(cc,3) = r_cra_zone(i,2)-t1Asp; 
%             %if line>0.85*nez/2
%                 maxline=line;
%             %end
%             cc=cc+1;
%             break;
%         end
%     end
% end
% 
% [popov1,S1] = polyfit((GGmin(:,3)),GGmin(:,2),1);
% [y1,delta1]=polyval(popov1,(GGmin(:,3)),S1);
% appro1=plot(GGmin(:,1)./Z+x1Asp/Z,y1./Z+Z/Z/2,'-c','Linewidth',1);
% %par rapport a l asperity
% %vfwcrt1=popov1(1)/sqrt(crt^2-v_long^2)
% %par rapport a l asperity
% vfwAcrt1=sqrt(v_long^2+popov1(1)^2)/crt;
% vfwAcdt1=sqrt(v_long^2+popov1(1)^2)/cdt;
% 
% [popov3,S3] = polyfit((GGmin(:,3)),GGmin(:,2),3);
% [y3,delta3]=polyval(popov3,(GGmin(:,3)),S3);
% appro3=plot(GGmin(:,1)./Z+x1Asp/Z,y3./Z+Z/Z/2,'-g','Linewidth',1);
% %par rapport a l asperity
% %vfwcrt3=(GGmin(0.15*(cc-1):(cc-1),3).*GGmin(0.15*(cc-1):(cc-1),3).*3.*popov3(1)+GGmin(0.15*(cc-1):(cc-1),3).*2.*popov3(2)+popov3(3))/(sqrt(crt*crt-v_long*v_long));
% %par rapport a l asperity
% vfwAcrt3=sqrt(v_long^2+(GGmin(:,3).*GGmin(:,3).*3.*popov3(1)+GGmin(:,3).*2.*popov3(2)+popov3(3)).^2)/crt;
% vfwAcrt3max=max(vfwAcrt3);
% vfwAcrt3min=min(vfwAcrt3);
% vfwAcrt3mean=mean(vfwAcrt3);
% 
% vfwAcdt3=sqrt(v_long^2+(GGmin(:,3).*GGmin(:,3).*3.*popov3(1)+GGmin(:,3).*2.*popov3(2)+popov3(3)).^2)/cdt;
% vfwAcdt3max=max(vfwAcdt3);
% vfwAcdt3min=min(vfwAcdt3);
% vfwAcdt3mean=mean(vfwAcdt3);
% 
% data=plot(GGmin(:,1)./Z+max(drawAsp(:,1))/Z,GGmin(:,2)./Z+Z/Z/2,'xb','Linewidth',1);  
% 
% %Rayleigh wave speed
% crplott=plot(xx./Z,zR./Z,'--g','Linewidth',2);
% 
% %Shear wave speed
% csplott=plot(xx./Z,zS./Z,'--r','Linewidth',2);
% 
% %Dilatational wave speed
% cdplott=plot(xx./Z,zD./Z,'--k','Linewidth',2);
% 
% %Dilatational influence plot
% Dil_inf=plot([dxinfl/Z,dxinfl/Z],[Z/Z/2,Z/Z],'-.k','Linewidth',1);
% 
% axis equal;
% ylabel('z/Z','Fontsize',23);
% xlabel('x/Z','Fontsize',23);
% axis([0 1.5*(xn)/Z Z/Z/2 Z/Z]);
% h_legend=legend([crplott,csplott,cdplott,Dil_inf,data,appro1,appro3],'c_{R}^{+}','c_{S}^{+}','c_{D}^{+}','c_{D}^{+} influence','front wave position',sprintf('linear inter. v_{FW}/c_{R}^{+}=%0.2f',vfwAcrt1),sprintf('cubic inter. v_{FW}/c_{R}^{+}=%0.2f',vfwAcrt3mean),'location','northeast');
% legend('boxoff');
% set(h_legend,'FontSize',11);
% set(gca, 'Color', 'none','Fontsize',22); % Sets axes background
% %,sprintf('v_{FW} v_{FW}/c_{R}^{+}=%0.2f',vfwcrt)
% linearSD4=100*mean(delta1)
% cubicSD4=100*mean(delta3)
toc/60

