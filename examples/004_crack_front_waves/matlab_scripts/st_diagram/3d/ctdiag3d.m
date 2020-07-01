% Author : Damien Spielmann 2015
% This function draws CT-Diagram from input file CT_Diagram.dat
% No variable to edit
% modified by fatima fekak

close all;
clear all;
clc;

fic = fopen('Parameters.cra');
parameters = textscan(fic, '%s %f');
fclose(fic);

nex = parameters{2}(1);
nez = parameters{2}(2);
X = parameters{2}(5);
Z = parameters{2}(6);

nut = parameters{2}(9);
nub = parameters{2}(10);

cst = parameters{2}(11);
csb = parameters{2}(12);

csbcst = csb/cst;

cdt=sqrt(2*(1-nut)/(1-2*nut));%c'est deja divisé par cs+
cdb=sqrt(2*(1-nub)/(1-2*nub))*csbcst;%gleichfalls

crt = (0.862+1.14*nut)/(1+nut);%c'est deja divisé par cs+
crb = (0.862+1.14*nub)/(1+nub)*csbcst;%gleichfalls

psi = parameters{2}(21);
phi = parameters{2}(22);


timer = load('Timer_ST_Diagram_id.cra');
%timer = load('Timer_ST_Diagram_shear_velo_jump.cra');
nb_steps = size(timer,1);
start_step = round(nb_steps / 3) ; %time step at which we want to start the post-processing
%end_step = nb_steps;   %time step at which we want to end the post-processing
end_step = round(2*nb_steps / 3) ;
time = timer(start_step:end_step,1);
ntss=size(time,1);
%
%For St_diagram_id
fic = fopen('ST_Diagram_id.cra', 'rb');
nb_bytes=4;
position = nex*nez*(start_step-1)*nb_bytes;
status = fseek(fic,position,'bof');
nb = nex * nez ;
ST_Diagram = fread(fic, [nb ntss], 'uint');
%position = ftell(fic);
fclose(fic);


% %For ST_Diagram_shear_velo_jump
% fic = fopen('ST_Diagram_shear_velo_jump.cra','rb');
% nb_bytes=8;
% position = nex*nez*(start_step-1)*nb_bytes;
% status = fseek(fic,position,'bof');
% nb = nex * nez ;
% ST_Diagram = fread(fic, [nb ntss], 'double');
% %position = ftell(fic);
% fclose(fic);
 
% ff
x = linspace(0,1,nex);
z = linspace(0,1,nez);

sampling = 20 ;
n=1;
k=10000;

% create the video writer with 1 fps
 %writerObj = VideoWriter('myVideo.avi');
 %writerObj.FrameRate = 1;

 % set the seconds per image
 %secsPerImage = 2;

 % open the video writer
 %open(writerObj);

while n<=size(ST_Diagram,2)
    for i=0:nez-1
        ind(:,i+1) = ST_Diagram(nex*i+1:nex*i+nex,n);
    end   
    imagesc(x,z,ind');
    axis xy;
    %To change de colormap limits
    %colorbar
    %lim = caxis;
    %caxis([0 0.0001])
    %caxis([0 0.001])
    
    nomFichier = sprintf('./png/st_diag_%d.png',k);
    saveas(gcf, nomFichier);
    %frame = getframe(gca);
    %size(frame.cdata)

     %for v=1:secsPerImage
     %    writeVideo(writerObj, frame);
     %end
    n = n + sampling ;
    k = k+1 ;
end

%k=1;

% while k<=nts2
%     for i=0:nez-1
%         ind(:,i+1) = ST_Diagram_suite(nex*i+1:nex*i+nex,k);
%     end   
%     imagesc(x,z,ind');
%     axis xy;
%     frame = getframe;
% 
%      %for v=1:secsPerImage
%          writeVideo(writerObj, frame);
%      %end
%     k = k + sampling ;
% end

 % close the writer object
% close(writerObj);

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
% %choose direction and line                              %
% direction=0;                                            %
% line=nez/2;                                             %
%                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %Read the file
% if direction==0
%     for h=0:nts-1
%         for i=0:nex-1
%             ind(h+1,i+1)=CT_Diagram(h+1,line*nex+i+1);%+1 to get uy and +1 matlab indice
%         end
%     end
%     x = linspace(0,X,nex);
% else
%     for h=0:nts-1
% 	    for j=0:nez-1
%             ind(h+1,j+1)=CT_Diagram(h+1,line+j*nex+1);
%         end
%     end
%     x = linspace(0,Z,nez);
%     X=Z;
% end
% 
% 
% %str = sprintf('c-t Diagram | psi = %d & phi = %d | Domain size : %d m',psi,phi,X);
% 
% 
% 
% y = load('Timer_ST_Diagram_id.cra');
% maxy = max(y);
% 
% imagesc (x./Z,y./Z,ind);
% 
% %title('ct-diagram cracking index','Fontsize',23);
% %title(str,'Fontsize',23)
% axis xy;
% 
% hold on;
% 
% %Dilatational waves
% x = [0 X];
% cdtt = [+X/cdt -X/cdt+X/cdt];
% cdbb = [+X/cdb -X/cdb+X/cdb];
% 
% plot(x./Z,cdtt./Z,'--k');
% %plot(x,cdbb,':k');
% 
% 
% %Shear waves
% csx = [0.2 0.3];
% csy = [0.7*maxy 0.7*maxy+0.1];
% plot(csx./Z,csy./Z,'--r','Linewidth',2);
% 
% % csx = [0.5 0.7];
% % csy = [0.7*maxy 0.7*maxy+0.2*cst/csb];
% % plot(csx,csy,':r','Linewidth',2);
% 
% 
% %Rayleigh wavex
% crx = [0.1 0.2];
% cry = [0.2*maxy 0.2*maxy+0.1/crt];
% plot(crx./Z,cry./Z,'--g','Linewidth',2);
% 
% % crx = [0.2 0.4];
% % cry = [0.2*maxy 0.2*maxy+0.2/crb];
% % plot(crx,cry,':g','Linewidth',2);
% 
% 
% xlabel('x/Z','Fontsize',38)
% ylabel('c_s^{+}t/Z','Fontsize',38)
% set(gca,'Color', 'none','LineWidth',2,'Fontsize',32);
% lg=legend('c_{d}^{+} infl.','c_{s}^{+}','c_{R}^{+}','Location','northeast');
% %set(lg,'FontSize',32);
% legend('boxoff');
