%% Prepare for normalization
%clear
if ~exist('SS_In','var')
    if exist('Cell_New_Init.png','file')
        SS_In=imread('Cell_New_Init.png');
    else
        SS_In=uint8(255*ones(301,1));
        SS_In(151)=1;
    end
end
cm=colormap(gray(256));cm(2,:)=[0 0 1];colormap(cm);
nSpecies=9;
nParams=9;

SS_Init=zeros(size(SS_In,1),size(SS_In,2),nSpecies);
Msk_Int=SS_In>0;
Msk_PM=SS_In==1;
Msk_Cyt=SS_In==255;

SS_Init(:,:,1)=2.0*Msk_Int;
SS_Init(:,:,2)=Msk_Int;

SS_P=maltesCA;
SS_P.StochVal=0.0;
SS_P.IntPar=ones(nParams,1);
SS_P.DiffRad=4.0*ones(nSpecies,1);
SS_P.Reps=999;

% moderate kinase/pptase activity, no cytosolic kinase, immediate
% dissociation of p-stat-tub, no association of p-stat and free tub
SS_P.IntPar=[2e-4,2e-5,1e-2,1,1e-3,0e-3,1e-2,1,1e-1];

% steady-state distribution:
%complex_SS=0.5*(SS_Init(6,1,1)+SS_Init(6,1,2)+SS_P.IntPar(5)/SS_P.IntPar(7)-...
%    sqrt((SS_Init(6,1,1)+SS_Init(6,1,2)+SS_P.IntPar(5)/SS_P.IntPar(7))^2-...
%    4*SS_Init(6,1,1)*SS_Init(6,1,2)));
%SS_Init(:,:,4)=complex_SS*Msk_Int;
%SS_Init(:,:,1)=(SS_Init(6,1,1)-complex_SS)*Msk_Int;
%SS_Init(:,:,2)=(SS_Init(6,1,2)-complex_SS)*Msk_Int;

%SS_Init_Zero=SS_Init;
figure(1); plot(SS_In)

%Species: tubulin, stathmin, pStathmin,
% stat-tub01, pStat-tub01, stat-tub10, pStat-tub10, stat-tub11, pStat-tub11


%%
fehlerGrenze=1e-9;
fehler=1; Reps=0;
clear maltesCA O V; 

%Sample parameter, only 1 tub binding site
%SS_P.IntPar=[1e-2,0e-5,1e-1,1,1e-3,0e-3,1e-2,1,1];
SS_P.IntPar=[2e-4,2e-5,1e-2,1,1e-3,0e-3,1e-2,1,1e-2];

%from paper in s-1, on-rates in µM-1s-1:
%k_on01    1.6, k_on10   0.9, k_on01_p  0.0128, k_on10_p  0.17
%k_off01 0.048, k_off10 1.86, k_off01_p      8, k_off10_p 2.3
%pptase: 5µM @ 0.02 µM-1s-1 => 0.1
%kinase: 2µM @ 0.001µM-1s-1 => 0.002
SS_P.DiffRad=[4,4,4,3,3,3,3,2,2];
%9Species
SS_P.IntPar=[1e-3,2e-5,1e-2,  0.08,0.00048, 0.000128,0.016,...
                             0.023,0.0186,    0.0017,0.009];
%5Species
SS_P.IntPar=[1e-3,2e-5,1e-1,  0.08,0.00048, 0.000128,0.016,...
                             0.00,0.00,    0,0];

SS_P.Reps=100;
SS_Init=zeros(size(SS_In,1),size(SS_In,2),nSpecies);
%Conc. Tubulin for 5 Spec 13.6µM
SS_Init(:,1)=10.0*Msk_Int;
SS_Init(:,2)=10.0*Msk_Int;
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

%%
%set(groot,'defaultAxesColorOrder','remove');
%    colorOrder=colormap(hsv(11)); 
colorOrder=zeros(11,3);
colorOrder(1,:)=[1 0 0];colorOrder(10,:)=[0.5 0 0];
colorOrder(2,:)=[0 0 1];colorOrder(11,:)=[0 0 0.5];
colorOrder(3,:)=[0 0.5 1];
colorOrder(4,:)=[0 1 0];colorOrder( 5,:)=[1 0 1];
colorOrder(6,:)=[0 1 0.5];colorOrder( 7,:)=[1 0.3 0.6];
colorOrder(8,:)=[0 1 1];colorOrder( 9,:)=[1 0.6 0.3];
set(groot,'defaultAxesColorOrder',colorOrder);

Data=squeeze(O(:,:,:,:,end));
Data(:,10)=sum(Data(:,[1 4 5 6 7 8 9 8 9]),2);
Data(:,11)=sum(Data(:,[2 3 4 5 6 7 8 9]),2);

xTickData=0.04*((0:size(Data,1)-1)-size(Data,1)/2);

width=300; height=250;
fig=figure(1); clf;
fig.Position=[50 50 width height];
fig.PaperPosition = [0 0 width/72 height/72];
fig.PaperUnits = 'inches';

plot(xTickData,Data,'--','LineWidth',2);
legend({'tubulin','stathmin','p-stat','stat-tub01','p-stat-tub01',...
        'stat-tub10','p-stat-tub10','stat-tub11','p-stat-tub11','tub_{tot}','stat_{tot}'});
xlim([0,4])
ylim([0,25])
xlabel('distance [µm]');
ylabel('concentration [µM]');
title('5 Species, 10.0µM tubulin')
Data_5Spec_10_0microM=Data;
Data_pStat(:,1)=sum(Data(:,[3 5 7 9]),2);
print('SupFig_6_5Species_10_0microM','-dsvg','-r0')