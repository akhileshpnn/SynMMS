%% Prepare for normalization
radiusGUV=318;
sizeSquare=320;
sizeProt=3;
testIm=zeros(sizeSquare,'uint8');
protInitPos=ceil(sqrt(radiusGUV^2/2));
for u=1:sizeProt+80; testIm(protInitPos+u-1,protInitPos+u-1)=1; end

protInitPos=ceil(sqrt(radiusGUV^2/5));
for u=1:sizeProt; testIm(1*protInitPos+floor(1*u/2)-1,2*protInitPos+u-1)=1; end
protInitPos=ceil(sqrt(radiusGUV^2/13));
%for u=1:sizeProt; testIm(2*protInitPos+floor(2*u/3)-1,3*protInitPos+u-1)=1; end
for u=1:sizeProt+40; testIm(floor(3*radiusGUV/5+3*u/4)-1,floor(4*radiusGUV/5)+u-1)=1; end

SE=strel('disk',sizeProt,0);
testIm=imdilate(testIm,SE);
SE=strel('square',3);
testIm=testIm+254*imerode(testIm,SE);
for u=1:sizeSquare
    for v=1:sizeSquare
        if (sqrt((u-1)^2+(v-1)^2)<radiusGUV-1) 
            testIm(u,v)=255;
        elseif (sqrt((u-1)^2+(v-1)^2)<radiusGUV && testIm(u,v)==0)
            testIm(u,v)=1;
        end
    end
end
clear u v SE ans
PM_Pos=(testIm==1);
cm=colormap(gray(256));cm(2,:)=[0 0 1];colormap(cm);
figure(1); imagesc(testIm); axis image

%%

if ~exist('SS_In','var')
    if exist('Cell_Init_New.png','file')
        SS_In=imread('Cell_Init_New.png');
    else
        SS_In=testIm;
    end
end
figure(1); imagesc(SS_In); axis image

SS_Init=zeros(size(SS_In,1),size(SS_In,2),9);
Msk_Int=SS_In>0;
Msk_PM=SS_In==1;
Msk_Cyt=SS_In==255;

PM_Init=1;
Cyt_Init=1;
SS_Init(:,:,1)=Msk_Int *PM_Init;
SS_Init(:,:,2)=Msk_Int*Cyt_Init;

SS_P=maltesCA;
SS_P.StochVal=0.0;
%Diff-Rads: 70 um^2/s ~ 6, 20 um^2/s ~ 3, but: 4 ~ 34 um^2/s
%=> 1px = 10 um & 1 time step = 1 s
%=> 0.01 um / px & 1e-6 s / time step
SS_P.DiffRad=[4,4,4,3,3,3,3,2,2];
SS_P.Reps=99;

% moderate kinase/pptase activity, no cytosolic kinase, immediate
% dissociation of p-stat-tub, no association of p-stat and free tub
%SS_P.IntPar=[1e-2,0e-5,1e-1,1,1e-3,0e-3,1e-2];
% 10x Reaction for all Diffusion!!!
% sim params: # x 1e6; e.g. kinase 1e-6 => 100%/s (max 2/s)
% pptase max@4uM: 0.1/s => 1e-7
% conc. Stathmin = 10uM; high-affinity site: k_on 16/s; k_off 4.8e-2/s 
% phosphorylated high-affinity site: k_on 1.8/s; k_off 2.3/s 
SS_P.IntPar=[1e-3,2e-5,1e-1,  0.08,0.00048, 0.000128,0.016,...
                             0.023,0.0186,    0.0017,0.009];
% steady-state distribution:
complex_SS=0.5*(PM_Init+Cyt_Init+SS_P.IntPar(5)/SS_P.IntPar(7)-...
          sqrt((PM_Init+Cyt_Init+SS_P.IntPar(5)/SS_P.IntPar(7))^2-...
              4*PM_Init*Cyt_Init));
SS_Init(:,:,4)=complex_SS*Msk_Int;
SS_Init(:,:,1)=( PM_Init-complex_SS)*Msk_Int;
SS_Init(:,:,2)=(Cyt_Init-complex_SS)*Msk_Int;

SS_Init_Zero=SS_Init;

%%
fehlerGrenze=1e-8; 
fehler=1; Reps=0;
clear maltesCA O V; 
set(groot,'defaultAxesColorOrder','remove');

tic;
while fehler>fehlerGrenze
  [O, ~]=maltesCA(SS_P,SS_In,SS_Init);
  SS_Init=squeeze(O(:,:,:,:,end));
  fehler=(SS_Init-squeeze(O(:,:,:,:,end-1)))./SS_Init;
  fehler(SS_Init==0)=0;
  fehler=sqrt(sum(fehler(:).*fehler(:)))/numel(fehler);
  disp(fehler);
  Reps=Reps+SS_P.Reps;
end
toc;

plotData=squeeze(O(:,:,:,:,end));

%%
width=500; height=300;
fig=figure(3); clf;
fig.Position=[50 50 width height];
fig.PaperPosition = [0 0 width/72 height/72];
fig.PaperUnits = 'inches';


freeTub=squeeze(plotData(:,:,1));
for u=sizeSquare:-1:1; gradPlotXVal(u,4)= sqrt((u-1)^2+(floor(3*u/4))^2); end
gradPlotXVal(:,1)=sqrt(2)*(0:sizeSquare-1);
for u=sizeSquare:-1:1; gradPlotXVal(u,3)= sqrt((u-1)^2+(floor(1*u/2))^2); end
for u=sizeSquare:-1:1; gradPlotXVal(u,2)= sqrt((u-1)^2+(floor(3*u/4))^2); end
gradPlotXVal=0.05*gradPlotXVal;

for u=sizeSquare:-1:1; gradPlotData(u,4)= freeTub(u,floor(3*u/4)+1); end
for u=sizeSquare:-1:1; gradPlotData(u,1)= freeTub(floor(u),u); end
for u=sizeSquare:-1:1; gradPlotData(u,3)= freeTub(floor(1*u/2)+1,u); end
for u=sizeSquare:-1:1; gradPlotData(u,2)= freeTub(floor(3*u/4)+1,u); end

title('free tubulin gradient in presence of protrusions')
pH=plot(gradPlotXVal,gradPlotData, 'Linewidth',3);
ylabel('concentration of free tubulin');
xlabel('distance from membrane [µm]');

axes('Position',[0.1 0.4 0.5 0.5]);
colormap(gray(256));
imagesc(flip(freeTub)), axis image, colorbar('Location','westoutside')
set(gca,'xtick',[])
set(gca,'ytick',[])
annotation('arrow',[0.244,0.45],[.403,.743],'Color',get(pH(1),'Color'),'LineWidth',2);
annotation('arrow',[0.244,0.480],[.403,.692],'Color',get(pH(2),'Color'),'LineWidth',2);
annotation('arrow',[0.244,0.506],[.403,.617],'Color',get(pH(3),'Color'),'LineWidth',2);
annotation('arrow',[0.244,0.526],[.403,.550],'Color',get(pH(4),'Color'),'LineWidth',2);
print('Grad2D','-dpng','-r0')

clear width height fig pH
