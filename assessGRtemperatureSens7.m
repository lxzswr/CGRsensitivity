% this script will...

% relate variations in the growth rate
% to tropical temperatures
% and precipitation and/or alpha

% NOTES FOR REMI:
% could update this to GCP 2019 version
% alpha is being read in (alpha = AET/PET) but not used
% Soil moisture could also be read in.
warning off;

inputdirClim= '/Volumes/RemiLBNL/project7/Data/CRU4.01/';

% load the PFT map
filename=strcat('/Volumes/RemiLBNL/project7/Data/LCMODIS.mat');
    tmp=load(filename);
    env.PFTs=tmp;
    env.PFTs=env.PFTs(:,3:end);
    clear tmp filename
% Mig2011:  ENF DBF GRA CRO     SAV     SHB     EBF MF WET
% MODIS:    1   4   10  (12,14) (8,9)   (6,7)   (2) (5) (11)

%%
if ~exist('GCPdata','var')
    
    GCPdata=importdata('/Volumes/RemiLBNL/project7/Data/Global_Carbon_Budget_2017v1.1.xlsx');
    %GCPdata=importdata('/Volumes/RemiLBNL/project7/Data/Global_Carbon_Budget_2019v1.0.xlsx');

    % not the indices (1,2,3,4,5 below) likely need to be changed for updated version
    
    % extract the data
    GCPdata.yearsGCP=GCPdata.data.GlobalCarbonBudget(:,1);
    GCPdata.ffEmissions=(GCPdata.data.GlobalCarbonBudget(:,2));
    GCPdata.lucEmissions=(GCPdata.data.GlobalCarbonBudget(:,3));
    GCPdata.totalEmissions=(GCPdata.ffEmissions+GCPdata.lucEmissions);
    GCPdata.growthRateGCP=(GCPdata.data.GlobalCarbonBudget(:,4));
    GCPdata.oceanUptake=(GCPdata.data.GlobalCarbonBudget(:,5));
    GCPdata.landUptake=GCPdata.totalEmissions-GCPdata.growthRateGCP-GCPdata.oceanUptake;
    
    GCPdata.GCPmodels=GCPdata.data.TerrestrialSink;
    
    GCPdata.airbornF=GCPdata.growthRateGCP./GCPdata.totalEmissions;
    
    % the model data is 4:17
    GCPdata.modelz=GCPdata.GCPmodels(:,4:end);
    GCPdata.GCPresSink=GCPdata.GCPmodels(:,2);
    
    for iii=1:17
        GCPdata.ACO2modelz(:,iii)=GCPdata.totalEmissions-GCPdata.modelz(:,iii)-GCPdata.oceanUptake;
    end
    %GCPdata.ACO2modelz=GCPdata.totalEmissions-GCPdata.modelz(:,:)-GCPdata.oceanUptake;
    
end

%%
% % load the monthly CO2 data from Mauna Loa
if ~exist('CO2','var')

monthlyCO2data=importdata('/Volumes/RemiLBNL/project7/Data/monthly_in_situ_co2_mlo 2018 no header.csv');
% Interpolate missing data (from some month before summer of 1964):
ind = monthlyCO2data(:, 2) == -99.99;
monthlyCO2data(ind, 2) = interp1(monthlyCO2data(~ind, 1), monthlyCO2data(~ind, 2), monthlyCO2data(ind, 1));

mData.t = monthlyCO2data(:, 1);
mData.CO2 = monthlyCO2data(:, 2);
clear ind monthlyCO2data
end

%%
% from 1959 to 2016, load CRU climate data 
% and extract tropical temps, alpha
if ~exist('env.alpha','var')

years=1959:2016;
countx=1;
mData.tropicalT=[];
mData.tropicalAlpha=[];
mData.tropicalPAR=[];
mData.tropicalPP=[];
mData.tropicalSWC=[];
mData.tropicalVPD=[];

aData.CO2gr=GCPdata.growthRateGCP;
aData.tropicalT=[];
aData.tropicalAlpha=[];
aData.tropicalPAR=[];
aData.tropicalPP=[];
aData.tropicalSWC=[];
aData.tropicalVPD=[];

aData.globalT=[];
aData.globalAlpha=[];
aData.globalPAR=[];
aData.globalPP=[];
aData.globalSWC=[];
aData.globalVPD=[];

for year=years
    
    yearstr=num2str(year);
       
    %%%%%%% START LOAD ENVIRONMENTAL VARIABLES
    % load ALPHA
    filename=strcat(inputdirClim,'alpha/cru_alpha_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.alpha=tmp.(names{:}); %csvread(strcat(inputdirClim,'alpha/cru_alpha_',num2str(year),'.csv'),R,C);
    env.alpha=env.alpha(:,3:end);
    
    % load temp (is in degrees C)
    filename=strcat(inputdirClim,'tmp2/cru_tmp_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.temp=tmp.(names{:}); 

    % load PAR 
    filename=strcat(inputdirClim,'par/cru_par_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.par=tmp.(names{:}); 
   
    % load PP
    filename=strcat(inputdirClim,'pre2/cru_pre_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.pre=tmp.(names{:}); 
    
    % load SWC
    filename=strcat(inputdirClim,'SWC/cru_SWC_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.swc=tmp.(names{:}); 
    
    % load VPD
    filename=strcat(inputdirClim,'vpd2/cru_vpd_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.vpd=tmp.(names{:}); 
    
    % topical forest index is '2'
%   indX=env.PFTs==2;

    % here taking tropical by latitude (from the temperature matrix)
    indX=env.temp(:,2)<23 & env.temp(:,2)>-23;
    
    mData.tropicalT=vertcat(mData.tropicalT,nanmean(env.temp(indX,3:end))');
    mData.tropicalAlpha=vertcat(mData.tropicalAlpha,nanmean(env.alpha(indX,1:end))');
    mData.tropicalPAR=vertcat(mData.tropicalPAR,nanmean(env.par(indX,3:end))');
    mData.tropicalPP=vertcat(mData.tropicalPP,nanmean(env.pre(indX,3:end))');
    mData.tropicalSWC=vertcat(mData.tropicalSWC,nanmean(env.swc(indX,3:end))');
    mData.tropicalVPD=vertcat(mData.tropicalVPD,nanmean(env.vpd(indX,3:end))');
    
    aData.tropicalT=vertcat(aData.tropicalT,nanmean(nanmean(env.temp(indX,3:end))));
    aData.tropicalAlpha=vertcat(aData.tropicalAlpha,nanmean(nanmean(env.alpha(indX,1:end))));
    aData.tropicalPAR=vertcat(aData.tropicalPAR,nanmean(nanmean(env.par(indX,3:end))));
    aData.tropicalPP=vertcat(aData.tropicalPP,nanmean(nansum(env.pre(indX,3:end),2)));%use sum for pp
    aData.tropicalSWC=vertcat(aData.tropicalSWC,nanmean(nanmean(env.swc(indX,3:end))));
    aData.tropicalVPD=vertcat(aData.tropicalVPD,nanmean(nanmean(env.vpd(indX,3:end))));
    
    
    
    aData.globalT=vertcat(aData.globalT,nanmean(nanmean(env.temp(:,3:end))));
    aData.globalAlpha=vertcat(aData.globalAlpha,nanmean(nanmean(env.alpha(:,1:end))));
    aData.globalPAR=vertcat(aData.globalPAR,nanmean(nanmean(env.par(:,3:end))));
    aData.globalPP=vertcat(aData.globalPP,nanmean(nansum(env.pre(:,3:end),2)));%use sum for pp
    aData.globalSWC=vertcat(aData.globalSWC,nanmean(nanmean(env.swc(:,3:end))));
    aData.globalVPD=vertcat(aData.globalVPD,nanmean(nanmean(env.vpd(:,3:end))));
    
    
    countx=countx+1;
end 
clear countx tmp indX names year
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finished loading data
%   Now detrend annual data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aData.TTDT=detrend(aData.tropicalT);
aData.TADT=detrend(aData.tropicalAlpha);
aData.TPDT=detrend(aData.tropicalPAR);
aData.TPPDT=detrend(aData.tropicalPP);
aData.TSWCDT=detrend(aData.tropicalSWC);
aData.TVPDDT=detrend(aData.tropicalVPD);
aData.CO2DT=detrend(aData.CO2gr);

%consider the effect of volcaneo
% cd('/Volumes/RemiLBNL/project7/code');
% vol_years=[1963; 1964; 1982; 1985; 1992; 1993]; %according to P.Cox
% for i=1:length(vol_years)
%     index=years==vol_years(i);
%     aData.CO2gr(index)=NaN;
% end
% aData.CO2DT=detrend(aData.CO2gr);


cd('/Volumes/RemiLBNL/project7/emd');
aData.TTEMD=emd(aData.tropicalT);
aData.TAEMD=emd(aData.tropicalAlpha);
aData.TPEMD=emd(aData.tropicalPAR);
aData.TPPEMD=emd(aData.tropicalPP);
aData.TSWCEMD=emd(aData.tropicalSWC);
aData.TVPDEMD=emd(aData.tropicalVPD);
aData.CO2EMD=emd(aData.CO2gr);
aData.TDEMD=emd(aData.tropicalVPD);



for iii=1:17
    aData.modelzDT(:,iii)=detrend(GCPdata.modelz(:,iii));
end

for iii=1:17
    aData.ACO2modelzDT(:,iii)=detrend(GCPdata.ACO2modelz(:,iii));
end




% get 4 month-lagged precipitation, annual values.
for yy = 1:57

 aData.tropicalPPlag(yy,1) = nansum(mData.tropicalPP((yy*12 - 4): (yy*12 + 8 - 1)));

end
 
aData.TPPlagDT=detrend(aData.tropicalPPlag);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finished loading data
%   Now deseasonalize and detrend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the growth rate
mData.CO2gr=vertcat(0.8,diff(mData.CO2));

% remove the seasonal cycle of the growth rate (12 month, 58 years)
mData.CO2grSeasonal=repmat(mean(reshape(mData.CO2gr,[12,58]),2),58,1);
mData.CO2grDeseason=mData.CO2gr-(mData.CO2grSeasonal-mean(mData.CO2grSeasonal));

% detrend the deseasonalized growth rate data
mData.CO2grDeseasonDetrend = detrend(mData.CO2grDeseason);

% convert to annual contribution
for ii=1:length(mData.CO2grDeseason)-12
    mData.CO2grAnn(ii)=sum(mData.CO2gr(ii:ii+12));
    mData.CO2grDSAnn(ii)=sum(mData.CO2grDeseason(ii:ii+12));
    mData.CO2grDSAnnDT(ii)=sum(mData.CO2grDeseasonDetrend(ii:ii+12));
end




% deseasonalize and detrend Tropical temperature
mData.TTseasonal=repmat(mean(reshape(mData.tropicalT,[12,58]),2),58,1); %the first T stands for Tropic, the second T stands for Temperature
mData.TTdeseason=mData.tropicalT-(mData.TTseasonal-mean(mData.TTseasonal));

mData.TTdeseasonDetrend = detrend(mData.TTdeseason);
for ii=1:length(mData.CO2grDeseason)-12
    mData.TTAnn(ii)=mean(mData.tropicalT(ii:ii+12)); %temperature is mean
    mData.TTDSAnn(ii)=mean(mData.TTdeseason(ii:ii+12));
    mData.TTDSAnnDT(ii)=mean(mData.TTdeseasonDetrend(ii:ii+12));
end
% mData.TTdsSmooth=smooth(mData.TTDSAnn,12);
% mData.TTdsDTSmooth=smooth(mData.TTDSAnnDT,12);



% now similarly for Alpha
% remove seasonal cycle from tropical alpha
mData.TAseasonal=repmat(mean(reshape(mData.tropicalAlpha,[12,58]),2),58,1);
mData.TAdeseason=mData.tropicalAlpha-(mData.TAseasonal-mean(mData.TAseasonal));

mData.TAdeseasonDetrend = detrend(mData.TAdeseason);
for ii=1:length(mData.CO2grDeseason)-12
    mData.TAAnn(ii)=mean(mData.tropicalAlpha(ii:ii+12));
    mData.TADSAnn(ii)=mean(mData.TAdeseason(ii:ii+12));
    mData.TADSAnnDT(ii)=mean(mData.TAdeseasonDetrend(ii:ii+12));
end
% mData.TAdsSmooth=smooth(mData.TADSAnnDT,12);



% now similarly for PAR
% remove seasonal cycle from tropical PAR
mData.TPseasonal=repmat(mean(reshape(mData.tropicalPAR,[12,58]),2),58,1);
mData.TPdeseason=mData.tropicalPAR-(mData.TPseasonal-mean(mData.TPseasonal));

mData.TPdeseasonDetrend = detrend(mData.TPdeseason);
for ii=1:length(mData.CO2grDeseason)-12
    mData.TPAnn(ii)=mean(mData.tropicalPAR(ii:ii+12));
    mData.TPDSAnn(ii)=mean(mData.TPdeseason(ii:ii+12));
    mData.TPDSAnnDT(ii)=mean(mData.TPdeseasonDetrend(ii:ii+12));
end
% mData.TPdsSmooth=smooth(mData.TPDSAnn,12);
% mData.TPdsDTSmooth=smooth(mData.TPDSAnnDT,12);


% now similarly for PRE
mData.TPPseasonal=repmat(mean(reshape(mData.tropicalPP,[12,58]),2),58,1); %the first T stands for Tropic, the second T stands for Temperature
mData.TPPdeseason=mData.tropicalPP-(mData.TPPseasonal-mean(mData.TPPseasonal));

mData.TPPdeseasonDetrend = detrend(mData.TPPdeseason);
for ii=1:length(mData.CO2grDeseason)-12
    mData.TPPAnn(ii)=sum(mData.tropicalPP(ii:ii+12));
    mData.TPPDSAnn(ii)=sum(mData.TPPdeseason(ii:ii+12));
    mData.TPPDSAnnDT(ii)=sum(mData.TPPdeseasonDetrend(ii:ii+12));
end


% now similarly for SWC, however SWC only has annual value
% mData.TSWCseasonal=repmat(mean(reshape(mData.tropicalSWC,[12,58]),2),58,1); %the first T stands for Tropic, the second T stands for Temperature
% mData.TSWCdeseason=mData.tropicalSWC-(mData.TSWCseasonal-mean(mData.TSWCseasonal));
% 
% mData.TSWCdeseasonDetrend = detrend(mData.TSWCdeseason);
% for ii=1:length(mData.CO2grDeseason)-12
%     mData.TSWCAnn(ii)=sum(mData.tropicalSWC(ii:ii+12));
%     mData.TSWCDSAnn(ii)=sum(mData.TSWCdeseason(ii:ii+12));
%     mData.TSWCDSAnnDT(ii)=sum(mData.TSWCdeseasonDetrend(ii:ii+12));
% end

mData.TSWCAnn=repmat(mean(reshape(mData.tropicalSWC,[1,58]),2),58,1);

figure;
plot(mData.CO2grDSAnn)
title('CO2 growth rate deseasonalized, annualized')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finished data processing
%   Now analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% Try repeat the running mean regression analysis of Wang et al. 2015
Y1=mData.CO2grDeseasonDetrend;
Y2=mData.TTdeseasonDetrend;
%y=Y1';
y=Y1;

%X = [ones(size(Y1)) Y2 Y3]; % Temperature (Y2), Alpha (Y3), PAR (Y4)
X = [ones(size(Y1)) Y2]; % Temperature (Y2), Alpha (Y3), PAR (Y4)

% DEFINE THE WINDOW LENGTH (15 years * 12 Months here) it is actually a 30
% years window.
window=15*12;
b=nan(size(X,2),length(X));
countx=1;
for ii=(window+1):length(X)-window-1
[b(:,countx),~,~,~,~] = regress(y(ii-window:ii+window,1),X(ii-window:ii+window,:)); %later on probably add bint
countx=countx+1;
end

% plot the intercept
figure;plot(12*b(1,:),'.')
title('Intercept of temperature sensitivity')

% plot the slope
figure;plot(12*b(2,:),'.')
title('Temperature sensitivity (slope)')









%%%%%%%%%%%% Import NEE, GPP and Reco from TRENDY and CMIP5, and FluxCOM %%%%%%%%%%%%%%%%%

%for TRENDY, NEE has been added through the spreadsheet 
%provide annual NEE, GPP and Reco, NEEDT, GPPDT and RecoDT
%provide monthly NEE, GPP and Reco
dataHome='/Volumes/RemiLBNL/project7/TRENDY/'; %9 models, NC format
cd(dataHome);

clear Trendy
%unit is KgC/m2/h
%resolution: 1 degree, annual

%load in area 0.5 D
filename_area=strcat('/Volumes/RemiLBNL/project6/code4Remi/areaGrid.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename_area);
areaGrid=rot90(areaGrid,3);

areaGrid1D=[];

%convert 0.5 degree to 1 degree
for iii=1:360
    for jjj=1:180
        areaGrid1D(iii,jjj)=nansum(nansum(areaGrid((iii-1)*2+1:iii*2,(jjj-1)*2+1:jjj*2),1),2);
    end
end



Trendy.Models={'CLM4CN','HYLAND','LPJ_GUESS','LPJ','OCN','ORCHIDEE','SDGVM','TRIFFID','VEGAS'};

for iii=1:9
    files=dir(strcat(dataHome,Trendy.Models{iii},'*'));
    filename=files(1).name;
    Trendy.GPP2(:,:,:,iii)=ncread(filename,'gpp');
    Trendy.RA2(:,:,:,iii)=ncread(filename,'ra');
    Trendy.RH2(:,:,:,iii)=ncread(filename,'rh');
    Trendy.NBP2(:,:,:,iii)=ncread(filename,'nbp');
end

%convert -99999 to NaN
Trendy.GPP2(Trendy.GPP2<-1000)=NaN;
Trendy.RA2(Trendy.RA2<-1000)=NaN;
Trendy.RH2(Trendy.RH2<-1000)=NaN;
Trendy.NBP2(Trendy.NBP2<-1000)=NaN;

%convert unit to Kg/yr
for iii=1:110
    for jjj=1:9
        Trendy.GPP2(:,:,iii,jjj)=Trendy.GPP2(:,:,iii,jjj).*areaGrid1D.*24.*365;
        Trendy.RA2(:,:,iii,jjj)=Trendy.RA2(:,:,iii,jjj).*areaGrid1D.*24.*365;
        Trendy.RH2(:,:,iii,jjj)=Trendy.RH2(:,:,iii,jjj).*areaGrid1D.*24.*365;
        Trendy.NBP2(:,:,iii,jjj)=Trendy.NBP2(:,:,iii,jjj).*areaGrid1D.*24.*365;
    end
end

%calculate annual GPP for each year and each model, unit it PgC
Trendy.annualGPP=squeeze(nansum(nansum(Trendy.GPP2,1),2))./(10^12);
Trendy.annualRA=squeeze(nansum(nansum(Trendy.RA2,1),2))./(10^12);
Trendy.annualRH=squeeze(nansum(nansum(Trendy.RH2,1),2))./(10^12);
Trendy.annualNBP=squeeze(nansum(nansum(Trendy.NBP2,1),2))./(10^12);

Trendy.annualNEP=Trendy.annualGPP-Trendy.annualRA-Trendy.annualRH;
Trendy.annualReco=Trendy.annualRA+Trendy.annualRH;

%the trend changed since 1959, so only consider 1959-2010, otherwise have
%to use EMD method.
for iii=1:9
    Trendy.GPPDT(:,iii)=detrend(Trendy.annualGPP(58:110,iii));
    Trendy.NEPDT(:,iii)=detrend(Trendy.annualNEP(58:110,iii));
    Trendy.RecoDT(:,iii)=detrend(Trendy.annualReco(58:110,iii));
    Trendy.NBPDT(:,iii)=detrend(Trendy.annualNBP(58:110,iii));
    
    if iii==6 %one model only has value from 1981
        Trendy.GPPDT(:,iii)=NaN;
        Trendy.NEPDT(:,iii)=NaN;
        Trendy.RecoDT(:,iii)=NaN;
        Trendy.NBPDT(:,iii)=NaN;
        
        Trendy.GPPDT(end-29:end,iii)=detrend(Trendy.annualGPP(81:110,iii));
        Trendy.NEPDT(end-29:end,iii)=detrend(Trendy.annualNEP(81:110,iii));
        Trendy.RecoDT(end-29:end,iii)=detrend(Trendy.annualReco(81:110,iii));
        Trendy.NBPDT(end-29:end,iii)=detrend(Trendy.annualNBP(81:110,iii));
    end
end






%for CMIP5, need to add NEE, GPP and Reco
%since all dataset have different coordinates, first interpolate them to 1
%degree. 
%provide annual NEE, GPP and Reco, NEEDT, GPPDT and RecoDT
%provide monthly NEE, GPP and Reco


dataHome='/Volumes/RemiLBNL/project7/CMIP5/';
cd(dataHome);

trig=0;

if trig==1
    load('CMIP52D.mat');
else
    clear CMIP5

    %need to use the models that have both GPP and NBP results
    %only consider historical run, 1850-2005

    CMIP5.Models={'BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    
    %CMIP5.Models={'bcc-csm1','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    
    
    %CMIP5.Models={'BNU-ESM','CCSM4','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};

    num=length(CMIP5.Models);

dataHome1='/Volumes/RemiLBNL/project7/CMIP5/Fluxes/';
dataHome2='/Volumes/RemiLBNL/project7/CMIP5/area/';

 for iii=1:num
        files=dir(strcat(dataHome1,'gpp_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;
        %read in first, and then interpolate to 1 D. this is monthly data 

        cd(dataHome1);
        tmpGPP2=ncread(filename,'gpp'); %if cannot read, check if cd to the directory
        tmpGPP2(tmpGPP2>1000)=NaN;

        files=dir(strcat(dataHome1,'npp_*',CMIP5.Models{iii},'*')); %BNU-ESM from bcc-csm...
        filename=files(1).name;

        tmpNPP2=ncread(filename,'npp');
        tmpNPP2=tmpNPP2.*100; %?
        tmpNPP2(tmpNPP2>1000)=NaN;
        
        
        files=dir(strcat(dataHome1,'rh_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;

        tmpRH2=ncread(filename,'rh');
        tmpRH2=tmpRH2.*100; %?
        tmpRH2(tmpRH2>1000)=NaN;
        
        
        files=dir(strcat(dataHome1,'nbp_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;

        tmpNBP2=ncread(filename,'nbp');
        %tmpNBP2=tmpNBP2.*100; %?
        tmpNBP2(tmpNBP2>1000)=NaN;
        
        

        if iii==8 %inmcm model GPP and NBP values are negative
            tmpGPP2=tmpGPP2.*(-1);
        end
        
        
       if iii==11 % MIROC-ESM fluxes are two orders of magnitude smaller than others
            tmpGPP2=tmpGPP2.*100;
            tmpNPP2=tmpNPP2.*100;
            tmpRH2=tmpRH2.*100;
            tmpNBP2=tmpNBP2.*100;
       end
        
       
        if iii==12 % MIROC-ESM-CHEM fluxes are two orders of magnitude smaller than others
            tmpGPP2=tmpGPP2.*100;
            tmpNPP2=tmpNPP2.*100;
            tmpRH2=tmpRH2.*100;
            tmpNBP2=tmpNBP2.*100;
        end


        tmplat=ncread(filename,'lat');
        tmplon=ncread(filename,'lon');
        
                
        %load the area and conversion factor of each model
        cd(dataHome2);
        files=dir(strcat(dataHome2,'areacella_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;
        
        tmparea=ncread(filename,'areacella');
        
        files=dir(strcat(dataHome2,'sftlf_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;
        
        tmpfl=ncread(filename,'sftlf');
        
    if size(tmpGPP2,3)<156*12
        %jstart=(156*12-size(tmpGPP2,3))./12+1;
        bridge_tmp=tmpGPP2;
        tmpGPP2(:,:,1:(156*12-size(tmpGPP2,3)))=NaN;
        tmpGPP2(:,:,(156*12-size(tmpGPP2,3))+1:156*12)=bridge_tmp;
        
        bridge_tmp=tmpNPP2;
        tmpNPP2(:,:,1:(156*12-size(tmpNPP2,3)))=NaN;
        tmpNPP2(:,:,(156*12-size(tmpNPP2,3))+1:156*12)=bridge_tmp;
        
        bridge_tmp=tmpRH2;
        tmpRH2(:,:,1:(156*12-size(tmpRH2,3)))=NaN;
        tmpRH2(:,:,(156*12-size(tmpRH2,3))+1:156*12)=bridge_tmp;
        
        bridge_tmp=tmpNBP2;
        tmpNBP2(:,:,1:(156*12-size(tmpNBP2,3)))=NaN;
        tmpNBP2(:,:,(156*12-size(tmpNBP2,3))+1:156*12)=bridge_tmp;
    end
        
        
    %convert unit from Kg/m2/s to PgC/yr
    
    for jjj=1:156 %in total 156 years 
        
        aGPP2=nanmean(tmpGPP2(:,:,(jjj-1)*12+1:jjj*12),3).*tmparea.*tmpfl.*24.*365.*3600./100;
        aNPP2=nanmean(tmpNPP2(:,:,(jjj-1)*12+1:jjj*12),3).*tmparea.*tmpfl.*24.*365.*3600./100;
        aRH2=nanmean(tmpRH2(:,:,(jjj-1)*12+1:jjj*12),3).*tmparea.*tmpfl.*24.*365.*3600./100;
        aNBP2=nanmean(tmpNBP2(:,:,(jjj-1)*12+1:jjj*12),3).*tmparea.*tmpfl.*24.*365.*3600./100;


        CMIP5.annualGPP(jjj,iii)=squeeze(nansum(nansum(aGPP2,2),1))./(10^12);
        CMIP5.annualNPP(jjj,iii)=squeeze(nanmean(nansum(aNPP2,2),1))./(10^12);
        CMIP5.annualRH(jjj,iii)=squeeze(nanmean(nansum(aRH2,2),1))./(10^12);
        CMIP5.annualNBP(jjj,iii)=squeeze(nanmean(nansum(aNBP2,2),1))./(10^12);
    end

 end


    %some years do not have NBP and GPP values, assign them to NaN;
    CMIP5.annualGPP(CMIP5.annualGPP==0)=NaN;
    CMIP5.annualNPP(CMIP5.annualNPP==0)=NaN;
    CMIP5.annualRH(CMIP5.annualRH==0)=NaN;
    CMIP5.annualNBP(CMIP5.annualNBP==0)=NaN;
    
    CMIP5.annualNEP=CMIP5.annualNPP-CMIP5.annualRH;
    CMIP5.annualReco=CMIP5.annualGPP-CMIP5.annualNPP+CMIP5.annualRH;


    %detrend, only consider 1959-2005
    for iii=1:num
        tmp1=CMIP5.annualGPP(end-46:end,iii);
        inD=isnan(tmp1);
        CMIP5.GPPDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];

        tmp1=CMIP5.annualNEP(end-46:end,iii);
        inD=isnan(tmp1);
        CMIP5.NEPDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];


        tmp1=CMIP5.annualReco(end-46:end,iii);
        inD=isnan(tmp1);
        CMIP5.RecoDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
        
        tmp1=CMIP5.annualNBP(end-46:end,iii);
        inD=isnan(tmp1);
        CMIP5.NBPDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
    end

    save('CMIP52D','CMIP5','-v7.3'); %takes a while
end






%load FLUXCOM results
dataHome='/Volumes/RemiLBNL/project7/FluxCom/annual/';
cd(dataHome);
clear FC;

filename_area=strcat('/Volumes/RemiLBNL/project6/code4Remi/areaGrid.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename_area);
areaGrid=flip(rot90(areaGrid,1)); %make sure the area grid is consistent with the GPP grid
%provide annual NEE, GPP and Reco, NEEDT, GPPDT and RecoDT
%provide monthly NEE, GPP and Reco
%for now, only load annual dataset
%from 1980 to 2013
%FC.Models={'_HB.ANN','_HB.MARS','_HB.RF','ANN','MARS','RF'};
FC.Models={'ANN','MARS','RF'}; %load these which applied to NEE and GPP simulation

num=length(FC.Models);


for iii=1:num
    
    
    for jjj=1:34
        files=dir(strcat(dataHome,'GPP.',FC.Models{iii},'*'));
        filename=files(jjj).name;
        %read in first, and then interpolate to 1 D. this is monthly data 

        tmpGPP2=ncread(filename,'GPP'); %if cannot read, check if cd to the directory
        tmpGPP2(tmpGPP2<-1000)=NaN;

        files=dir(strcat(dataHome,'NEE.',FC.Models{iii},'*'));
        filename=files(jjj).name;

        tmpNEE2=ncread(filename,'NEE');
        tmpNEE2(tmpNEE2<-1000)=NaN;     


        files=dir(strcat(dataHome,'TER.',FC.Models{iii},'*'));
        filename=files(jjj).name;

        tmpTER2=ncread(filename,'TER');
        tmpTER2(tmpTER2<-1000)=NaN;   

        FC.GPP2(:,:,jjj,iii)=tmpGPP2;
        FC.NEE2(:,:,jjj,iii)=tmpNEE2;
        FC.TER2(:,:,jjj,iii)=tmpTER2;
    end
   
    
end

%unit is g/m2/day, 
for iii=1:34
    for jjj=1:3
        FC.GPP2(:,:,iii,jjj)=FC.GPP2(:,:,iii,jjj).*areaGrid.*365; %is areaGrid the right direction?
        FC.NEE2(:,:,iii,jjj)=FC.NEE2(:,:,iii,jjj).*areaGrid.*365;
        FC.Reco2(:,:,iii,jjj)=FC.TER2(:,:,iii,jjj).*areaGrid.*365;
    end
end

%to PgC
FC.annualGPP=squeeze(nansum(nansum(FC.GPP2,1),2))./(10^15);
FC.annualNEE=squeeze(nansum(nansum(FC.NEE2,1),2))./(10^15);
FC.annualReco=squeeze(nansum(nansum(FC.Reco2,1),2))./(10^15);


for iii=1:3
    FC.GPPDT(:,iii)=detrend(FC.annualGPP(:,iii));
    FC.NEEDT(:,iii)=detrend(FC.annualNEE(:,iii));
    FC.RecoDT(:,iii)=detrend(FC.annualReco(:,iii));
end






%%for CMIP5, need modelled TAS and PP
%provide annual tropical temperature for 14 models, from 1959-2005

%dataHome='/Volumes/RemiLBNL/project7/CMIP5/Fluxes';

dataHome1='/Volumes/RemiLBNL/project7/CMIP5/Fluxes/';
dataHome2='/Volumes/RemiLBNL/project7/CMIP5/area/';
dataHome3='/Volumes/RemiLBNL/project7/CMIP5/prtas/';
dataHome='/Volumes/RemiLBNL/project7/CMIP5';
cd(dataHome);

trig=0;

if trig==1
    load('CMIP52D.mat');
else
    

    %need to use the models that have both GPP and NBP results
    %only consider historical run, 1850-2005
    %CMIP5.Models={'bcc-csm1','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    CMIP5.Models={'BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    CMIP5.Clstart=[1850,1850,1850,1956,0,1934,1934,1850,1850,0,1850,1850,1850,1850]; %1934.12
    %CMIP5.Models={'BNU-ESM','CCSM4','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};

    num=length(CMIP5.Models);

%initialize the output, all assigned as NaN
num_years=2005-1959+1;

bData.tropicalT=nan(num_years,num);
bData.TTDT=nan(num_years,num);
bData.tropicalPP=nan(num_years,num);
bdata.TPPDT=nan(num_years,num);

 for iii=1:num %for each model
     
        %read in the land mask . NOTE: the NETCDF deminsion is different
        %than Matlab's
        cd(dataHome2);
        files=dir(strcat(dataHome2,'sftlf_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;
        lcmask=ncread(filename,'sftlf'); %if cannot read, check if cd to the directory
        %select out the tropical land mask
        num_lat=size(lcmask,2);
        low_b=round(67./(180./num_lat));
        high_b=round(113./(180./num_lat));
        
%         lcmask2=lcmask;
%         lcmask2=NaN;
        lcmask(:,[1:low_b,high_b:end])=NaN;
        lcmask_final=lcmask;
     
        
        
        %read temperature first
        cd(dataHome3);
        files=dir(strcat(dataHome3,'tas_*',CMIP5.Models{iii},'*'));
        
        tmpTT=[];tmpTT2=[];
        
        if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'tas'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 1959 to 2005    
            switch CMIP5.Clstart(iii)

                case 1850
                    tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
                case 1934
                    tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
                    tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
                case 1956
                    tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);

            end
        
            %covert to annual tmp
            for iy=1:47 %1959-2005
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                bData.tropicalT(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1);
            end
        
        
        end
        
        
        
        
        
        
        %read precipiation
                
        cd(dataHome3);
        files=dir(strcat(dataHome3,'pr_*',CMIP5.Models{iii},'*'));
        
        tmpTT=[];tmpTT2=[];
        
        if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'pr'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 1959 to 2005    
            switch CMIP5.Clstart(iii)

                case 1850
                    tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
                case 1934
                    tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
                    tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
                case 1956
                    tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);

            end
        
            %covert to annual tmp
            for iy=1:47 %1959-2005
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                bData.tropicalPP(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1).*365.*24.*3600; %change kg/m2/s to mm/yr
            end
        
        
        end
  
  

 end





    %detrend
    for iii=1:num
        tmp1=bData.tropicalT(:,iii);
        inD=isnan(tmp1);
        bData.TTDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
        
        tmp1=bData.tropicalPP(:,iii);
        inD=isnan(tmp1);
        bData.TPPDT(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
    end

    save('CMIP5Climate','CMIP5','-v7.3'); %takes a while
end







%for FLUXCOM data, need CRU-NCEP data, from 1980 to 2016
% filename=strcat('/Volumes/RemiLBNL/project6/datasets_pro/pp_tair_cruncep.mat');    
% load(filename);
% cData.tropicalT=[];
% cData.tropicalPP=[];
% cData.TTDT=[];
% 
% pp_cruncep_year=yearPP;
% tair_cruncep_year=yearTair;
% 
% for i=1:37
%     
%     tmp_T=squeeze(yearTair(i,:,:)).*(lc>0);
%     tmp_PP=squeeze(yearPP(i,:,:)).*(lc>0);
%     
%     cData.tropicalT(i,1)=nanmean(reshape(tmp_T(154:226,:),[],1));
%     cData.tropicalPP(i,1)=nanmean(reshape(tmp_PP(154:226,:),[],1));
%      
% end
% 
% cData.TTDT=detrend(cData.tropicalT);
filename=strcat('/Volumes/RemiLBNL/project6/data4remi/ELM_data/ElNino2015.c20171016.nc');    

tmpTT=ncread(filename,'TSA'); %if cannot read, check if cd to the directory
tmpPP=ncread(filename,'RAIN');

%we know it is from 1979 to 2016 monthly data
for i=1:38
    tmp_y_T=squeeze(nanmean(tmpTT(:,:,(i-1)*12+1:i*12),3))-273.15;
    tmp_y_P=squeeze(nanmean(tmpPP(:,:,(i-1)*12+1:i*12),3)).*3600.*24.*365;
    
    
    cData.tropicalT(i,1)=nanmean(reshape(tmp_y_T(:,33:64),[],1));
    cData.tropicalPP(i,1)=nanmean(reshape(tmp_y_P(:,33:64),[],1));
end

cData.TTDT=detrend(cData.tropicalT);
cData.TPPDT=detrend(cData.tropicalPP);









dataHome='/Volumes/RemiLBNL/project6/data4Remi/Jena_globalGPP'; %add this for a later study




%%% add GRACE data, 1979-2016, monthly, though it is provided in 0.5
%%% degree, the effective resolution is 3 degree. the data only includes
%%% land
cd('/Volumes/RemiLBNL/project7');

filename=strcat('/Volumes/RemiLBNL/project7/GRACE_REC_v02.nc');    

rec_GRACE=ncread(filename,'rec_ensemble_mean'); %if cannot read, check if cd to the directory

for i=1:38
    TWS.y_tGRACE(i,1)=squeeze(nanmean(reshape(rec_GRACE(:,134:226,(i-1)*12+1:i*12),[],1))); %
end

TWS.dt_tGRACE=detrend(TWS.y_tGRACE);


for i=1:38
    TWS.y_GRACE(i,1)=squeeze(nanmean(reshape(rec_GRACE(:,:,(i-1)*12+1:i*12),[],1))); %
end

TWS.dt_GRACE=detrend(TWS.y_GRACE);












% add other precipitation datasets
% all observation-based: GPCC, PRECL, CPC, UDEL
% maybe not use CPC because it is short.


apData.GPCCtropicP = [];
apData.PRECLtropicP = [];
apData.CPCtropicP = [];
apData.UDELtropicP = [];
apData.GHCNtropicP = [];

% read in GPCC, monthly, 0.5 degree
% from 1891 to 2016
% only use 1959 to 2016
filename = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/precip.mon.total.v2018.nc';
rawdata = ncread(filename,'precip'); % mm
lat_range = ncread(filename,'lat');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1891).*12 + 1;
    high = (year-1891).*12 + 12;
    
    tmp_raw = rawdata(:,lat_ind,low:high);
    
    tmp_raw(tmp_raw<=0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    tmp_sum=nansum(tmp_raw,3);
    tmp_sum(tmp_sum==0)=NaN; % remove those ocean values
    
    apData.GPCCtropicP(i) = nanmean(reshape(tmp_sum,[],1));
    
end



% read in PRECL, monthly, 1 degree
filename = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/precip.mon.mean.1x1.nc';

rawdata = ncread(filename,'precip'); % mm/day
lat_range = ncread(filename,'lat');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1948).*12 + 1;
    high = (year-1948).*12 + 12;
    
    tmp_raw = rawdata(:,lat_ind,low:high).*30;
    
    tmp_raw(tmp_raw<=0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    tmp_sum=nansum(tmp_raw,3);
    tmp_sum(tmp_sum==0)=NaN; % remove those ocean values
    
    apData.PRECLtropicP(i) = nanmean(reshape(tmp_sum,[],1));
    
    %apdata.PRECLtropic(i) = nanmean(reshape(nansum(tmp_raw,3),[],1));
    
end


% read in CPC, 1 file per year, daily, 0.5 degree
dir_name = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/CPC/';
cd(dir_name);
files = dir(strcat(dir_name,'*nc')); % from 1979 to 2018

apData.CPCtropicP(1:20)=NaN; % 1959-1978 no values

for i = 1:length(files)-2 % 2017 and 2018 no need here
    
    filename = files(i).name;
    rawdata = ncread(filename,'precip');
    lat_range = ncread(filename,'lat');
    lat_ind = lat_range < 23 & lat_range > -23;
    
    tmp_raw = rawdata(:,lat_ind,:);
    tmp_raw(tmp_raw<0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    tmp_sum=nansum(tmp_raw,3);
    tmp_sum(tmp_sum==0)=NaN; % remove those ocean values
    
    apData.CPCtropicP(i+20) = nanmean(reshape(tmp_sum,[],1));
    
    
end


% read in UDEL

dir_name = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/UDEL/';
cd(dir_name);
files = dir(strcat(dir_name,'*')); % from 1900 to 2017
%load in land mask
% filename=strcat('/Volumes/RemiLBNL/project7/Data/CRU4.01/','alpha/cru_alpha_2000.mat');
% tmp=load(filename);
% names=fieldnames(tmp);
% tmp2=tmp.(names{:}); 
% land_latlon=tmp2(:,1:2);


for i = 62:length(files)-1 % only use 1959 to 2016 here, only use 1959 to 2016 here, also two empty files at the beginning
    
    filename = files(i).name;
    rawdata = importdata(filename);
    
    
    lat_ind = rawdata(:,2) < 23 & rawdata(:,2) > -23; % the tropics
    
    tmp_raw = rawdata(lat_ind,3:14);
    tmp_raw(tmp_raw<0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    apData.UDELtropicP(i-61) = nanmean(reshape(nansum(tmp_raw,2),[],1));
    
end




% read in GHCN-M, monthly, 5 degree, 1900-2015
filename = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/GHCN_precip.mon.total.nc';

%https://www.esrl.noaa.gov/psd/data/gridded/data.ghcngridded.html

rawdata = ncread(filename,'precip'); % mm/day
lat_range = ncread(filename,'lat');

lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1900).*12 + 1;
    high = (year-1900).*12 + 12;
    
    if year == 2016 | year == 2015
        apData.GHCNtropicP(i)=NaN;
    else
    
        tmp_raw = rawdata(:,lat_ind,low:high).*30;

        tmp_raw(tmp_raw<=0)=NaN;
        tmp_raw(tmp_raw>1000)=NaN;

        tmp_sum=nansum(tmp_raw,3);
        tmp_sum(tmp_sum==0)=NaN; % remove those ocean values

        apData.GHCNtropicP(i) = nanmean(reshape(tmp_sum,[],1));
     end
    
    %apdata.PRECLtropic(i) = nanmean(reshape(nansum(tmp_raw,3),[],1));
    
end



% detrend these datasests
apData.GPCCDT_P=detrend(apData.GPCCtropicP);
apData.PRECLDT_P=detrend(apData.PRECLtropicP);

ind = isnan(apData.CPCtropicP);
apData.CPCDT_P=detrend(apData.CPCtropicP,'linear');
apData.CPCDT_P(~ind)=detrend(apData.CPCtropicP(~ind),'linear');

apData.UDELDT_P=detrend(apData.UDELtropicP);
apData.GHCNDT_P=detrend(apData.GHCNtropicP);










% add other temperature datasets
% https://climatedataguide.ucar.edu/climate-data/global-land-precipitation-and-temperature-willmott-matsuura-university-delaware
% all observation-based: GISS, BEST, UDEL. CRUTEMP4 not the same to CRU? 
apData.GISStropicT = [];
apData.BESTtropicT = [];
apData.UDELtropicT = [];
apData.CRUTEMtropicT = [];

% read in GISS TEMP
% from 1880 to 2019.5, 2 degree
% only use 1959 to 2016
filename = '/Volumes/RemiLBNL/project7/more_temperature_dataset/gistemp250_GHCNv4.nc';
rawdata = ncread(filename,'tempanomaly'); % anomaly in K, same as anomaly in C
lat_range = ncread(filename,'lat');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1880).*12 + 1;
    high = (year-1880).*12 + 12;
    
    tmp_raw = rawdata(:,lat_ind,low:high);
    
%     tmp_raw(tmp_raw<-100)=NaN;
%     tmp_raw(tmp_raw>100)=NaN;
    
    tmp_mean=nanmean(tmp_raw,3);
    tmp_mean(tmp_mean==0)=NaN; % remove those ocean values
    
    apData.GISStropicT(i) = nanmean(reshape(tmp_mean,[],1)); % anomaly, not absolute T
    
end





% read in Berkeley TEMP (BEST)
% from 1750 to 2019.4, 1 degree
% only use 1959 to 2016
filename = '/Volumes/RemiLBNL/project7/more_temperature_dataset/BEST_Complete_TAVG_LatLong1.nc';
rawdata = ncread(filename,'temperature'); % anomaly in C
lat_range = ncread(filename,'latitude');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1750).*12 + 1;
    high = (year-1750).*12 + 12;
    
    tmp_raw = rawdata(:,lat_ind,low:high);
    
    
    tmp_mean=nanmean(tmp_raw,3);
    tmp_mean(tmp_mean==0)=NaN; % remove those ocean values
    
    apData.BESTtropicT(i) = nanmean(reshape(tmp_mean,[],1)); % anomaly, not absolute T
    
end



% read in UDEL

dir_name = '/Volumes/RemiLBNL/project7/more_temperature_dataset/UDEL_air_temp_2017/';
cd(dir_name);
files = dir(strcat(dir_name,'*')); % from 1900 to 2017, 0.5 degree

for i = 62:length(files)-1 % only use 1959 to 2016 here, also two empty files at the beginning
    
    filename = files(i).name;
    rawdata = importdata(filename);
    
    
    lat_ind = rawdata(:,2) < 23 & rawdata(:,2) > -23; % the tropics
    
    tmp_raw = rawdata(lat_ind,3:14);
    tmp_raw(tmp_raw<0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    apData.UDELtropicT(i-61) = nanmean(reshape(nanmean(tmp_raw,2),[],1));
    
end



% read in CRUTEM v4
% from 1850 to 2019.5, 5 degree
% only use 1959 to 2016
filename = '/Volumes/RemiLBNL/project7/more_temperature_dataset/CRUTEM.4.6.0.0.anomalies.nc';
rawdata = ncread(filename,'temperature_anomaly'); % anomaly in K, same as anomaly in C
lat_range = ncread(filename,'latitude');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1850).*12 + 1;
    high = (year-1850).*12 + 12;
    
    tmp_raw = rawdata(:,lat_ind,low:high);
    
%     tmp_raw(tmp_raw<-100)=NaN;
%     tmp_raw(tmp_raw>100)=NaN;
    
    tmp_mean=nanmean(tmp_raw,3);
    tmp_mean(tmp_mean==0)=NaN; % remove those ocean values
    
    apData.CRUTEMtropicT(i) = nanmean(reshape(tmp_mean,[],1)); % anomaly, not absolute T
    
end




% detrend these datasests
apData.GISSDT_T=detrend(apData.GISStropicT);
apData.BESTDT_T=detrend(apData.BESTtropicT);
apData.UDELDT_T=detrend(apData.UDELtropicT);
apData.CRUTEMDT_T=detrend(apData.CRUTEMtropicT);









% add future climate scenarios from CMIP5
%%for CMIP5, need modelled TAS and PP
%provide annual tropical temperature for 14 models, from 2006 to 2100

%dataHome='/Volumes/RemiLBNL/project7/CMIP5/Fluxes';

dataHome1='/Volumes/RemiLBNL/project7/CMIP5/Fluxes/';
dataHome2='/Volumes/RemiLBNL/project7/CMIP5/area/';
dataHome3='/Volumes/RemiLBNL/project7/CMIP5/RCPprtas/';
dataHome='/Volumes/RemiLBNL/project7/CMIP5';
cd(dataHome);

trig=0;

if trig==1
    load('CMIP5future2D.mat');
else
    

    %need to use the models that have both GPP and NBP results
    %CMIP5.Models={'bcc-csm1','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    CMIP5.Models={'BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    CMIP5.Clstart=[2006,2006,2006,2006,0,2006,2006,2006,2006,0,2006,2006,2006,2006];
    
    num=length(CMIP5.Models);

    %initialize the output, all assigned as NaN
    num_years=2100-2006+1;

    rcpData.tropicalT45=nan(num_years,num);
    rcpData.TTDT45=nan(num_years,num);
    rcpData.tropicalPP45=nan(num_years,num);
    rcpData.TPPDT45=nan(num_years,num);

    rcpData.tropicalT85=nan(num_years,num);
    rcpData.TTDT85=nan(num_years,num);
    rcpData.tropicalPP85=nan(num_years,num);
    rcpData.TPPDT85=nan(num_years,num);

 for iii=1:num %for each model
     
        %read in the land mask . NOTE: the NETCDF deminsion is different
        %than Matlab's
        cd(dataHome2);
        files=dir(strcat(dataHome2,'sftlf_*',CMIP5.Models{iii},'*'));
        filename=files(1).name;
        lcmask=ncread(filename,'sftlf'); %if cannot read, check if cd to the directory
        %select out the tropical land mask
        num_lat=size(lcmask,2);
        low_b=round(67./(180./num_lat));
        high_b=round(113./(180./num_lat));
        
%         lcmask2=lcmask;
%         lcmask2=NaN;
        lcmask(:,[1:low_b,high_b:end])=NaN;
        lcmask_final=lcmask;
     
        
        
        %read rcp45 temperature first
        cd(dataHome3);
        files=dir(strcat(dataHome3,'tas_*',CMIP5.Models{iii},'*_rcp45*'));
        
        tmpTT=[];tmpTT2=[];
        
        if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'tas'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 2006 to 2100   
            
            if strcmp(CMIP5.Models{iii},'HadGEM2-ES') % this model stops at 2099
                tmpTT_month=tmpTT;
                tmpTT_month(:,:,end+1:end+12)=tmpTT(:,:,end-12+1:end);
            else
                tmpTT_month=tmpTT;
            end
%             switch CMIP5.Clstart(iii)
% 
%                 case 1850
%                    tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
%                 case 1934
%                     tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
%                     tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
%                 case 1956
%                     tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);
% 
%             end
        
            %covert to annual tmp
            for iy=1:95 %2100-2006
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                rcpData.tropicalT45(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1)-273.15;
            end
        
        
        end
        
        
        
        
        
        
        %read precipiation
                
        cd(dataHome3);
        files=dir(strcat(dataHome3,'pr_*',CMIP5.Models{iii},'*_rcp45*'));
        
        tmpTT=[];tmpTT2=[];
        
        if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'pr'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 2006 to 2100    
            if strcmp(CMIP5.Models{iii},'HadGEM2-ES') % this model stops at 2099
                tmpTT_month=tmpTT;
                tmpTT_month(:,:,end+1:end+12)=tmpTT(:,:,end-12+1:end);
            else
                tmpTT_month=tmpTT;
            end
%             switch CMIP5.Clstart(iii)
% 
%                 case 1850
%                     tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
%                 case 1934
%                     tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
%                     tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
%                 case 1956
%                     tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);
% 
%             end
        
            %covert to annual tmp
            for iy=1:95 %1959-2005
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                rcpData.tropicalPP45(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1).*365.*24.*3600; %change kg/m2/s to mm/yr
            end
        
        
        end
  
  
        
        
         %read rcp85 temperature first
        cd(dataHome3);
        files=dir(strcat(dataHome3,'tas_*',CMIP5.Models{iii},'*_rcp85*'));
        
        tmpTT=[];tmpTT2=[];
        
    if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'tas'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 2006 to 2100   
            
            if strcmp(CMIP5.Models{iii},'HadGEM2-ES') % this model stops at 2099
                tmpTT_month=tmpTT;
                tmpTT_month(:,:,end+1:end+12)=tmpTT(:,:,end-12+1:end);
            else
                tmpTT_month=tmpTT;
            end
%             switch CMIP5.Clstart(iii)
% 
%                 case 1850
%                    tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
%                 case 1934
%                     tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
%                     tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
%                 case 1956
%                     tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);
% 
%             end
        
            %covert to annual tmp
            for iy=1:95 %2100-2006
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                rcpData.tropicalT85(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1)-273.15;
            end
        
        
    end
        
        
        
        
        
        
        %read rcp 85 precipiation
                
        cd(dataHome3);
        files=dir(strcat(dataHome3,'pr_*',CMIP5.Models{iii},'*_rcp85*'));
        
        tmpTT=[];tmpTT2=[];
        
        if length(files)>0 && CMIP5.Clstart(iii)>0
        
            for jjj=1:length(files)

                filename=files(jjj).name;
            %read in first, and then interpolate to 1 D. this is monthly data 

                tmpTT2=ncread(filename,'pr'); %if cannot read, check if cd to the directory
                tmpTT=cat(3,tmpTT,tmpTT2); % if there are 

            end
            
            
            %cut the year from 2006 to 2100    
            if strcmp(CMIP5.Models{iii},'HadGEM2-ES') % this model stops at 2099
                tmpTT_month=tmpTT;
                tmpTT_month(:,:,end+1:end+12)=tmpTT(:,:,end-12+1:end);
            else
                tmpTT_month=tmpTT;
            end
%             switch CMIP5.Clstart(iii)
% 
%                 case 1850
%                     tmpTT_month=tmpTT(:,:,(1958-1850+1)*12:end);
%                 case 1934
%                     tmpTT_month=tmpTT(:,:,(1958-1934+1)*12-11:end);
%                     tmpTT_month(:,:,end+1)=tmpTT_month(:,:,end);
%                 case 1956
%                     tmpTT_month=tmpTT(:,:,(1958-1956+1)*12:end);
% 
%             end
        
            %covert to annual tmp
            for iy=1:95 %1959-2005
                tmpTT_annual=squeeze(nanmean(tmpTT_month(:,:,(iy-1)*12+1:iy*12),3));
                tmpTT_annual_L=tmpTT_annual.*(lcmask_final>0);
                tmpTT_annual_L(tmpTT_annual_L==0)=NaN;
                rcpData.tropicalPP85(iy,iii)=nanmean(nanmean(tmpTT_annual_L,2),1).*365.*24.*3600; %change kg/m2/s to mm/yr
            end
        
        
        end
        
 end





    %detrend rcp 45
    for iii=1:num
        tmp1=rcpData.tropicalT45(:,iii);
        inD=isnan(tmp1);
        rcpData.TTDT45(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
        
        tmp1=rcpData.tropicalPP45(:,iii);
        inD=isnan(tmp1);
        rcpData.TPPDT45(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
    end
    
    
    
    %detrend rcp 85
    for iii=1:num
        tmp1=rcpData.tropicalT85(:,iii);
        inD=isnan(tmp1);
        rcpData.TTDT85(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
        
        tmp1=rcpData.tropicalPP85(:,iii);
        inD=isnan(tmp1);
        rcpData.TPPDT85(:,iii)=[detrend(tmp1(~inD));tmp1(inD)];
    end
    
    
     

    save('CMIP5futureClimate','rcpData','-v7.3'); %takes a while
end









%%%%%%%% LOAD CRU-JRA v2.0 data (For TRENDY) %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% END: LOAD CRU-JRA v2.0 data (For TRENDY) %%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%% LOAD MODIS LST %%%%%%%%%%%%%%%%%%%%%%%%%%


load('/Volumes/RemiLBNL/project7/Data/MODIS_LST.mat');
% % below is the code to process the raw data

% % ref : http://hdfeos.org/zoo
% 
% clear all;
% 
% import matlab.io.hdfeos.*
% import matlab.io.hdf4.*
% 
% 
% % open directory
% datadir = '/Volumes/RemiLBNL/MODIS/LSTraw/';
% cd(datadir);
% 
% file_list = dir(strcat(datadir,'*hdf'));
% 
% interested_field = {'LST_Day_CMG', 'QC_Day','LST_Night_CMG','QC_Night'};
% 
% MOD11_LST = [];
% 
% for i = 1: length(file_list)
%     
%     FILE_NAME=file_list(i).name;
%     GRID_NAME='MODIS_MONTHLY_0.05DEG_CMG_LST';
% 
%     % get file info
%     field_info = hdfinfo(FILE_NAME, 'eos');
% 
%     file_id = gd.open(FILE_NAME, 'rdonly');
% 
%     % Open Grid
%     grid_id = gd.attach(file_id, GRID_NAME);
% 
%     % Define the Data Field
%     DATAFIELD_NAME='LST_Day_CMG';
% 
%     % Read the dataset.
%     data = gd.readField(grid_id, DATAFIELD_NAME, [], [], []);
%     gd.detach(grid_id);
%     gd.close(file_id);
%     
%     % regrid it to 0.5 degree
%     data2 = imresize(data,0.1);
%     
% 
%     % reference: https://lpdaac.usgs.gov/documents/118/MOD11_User_Guide_V6.pdf
%     tmp = rot90(flip(data2),3).*0.02 - 273.15;
%     
%     MOD11_LST(:,:,i) = tmp;
%     
%     disp(i);
%     
% end
% 
% 
% 
% % convert to annual values
% start_m = 1;
% 
% 
% MOD11_LST_ann = [];
% for yy = 2000:2019
%     
%     if yy == 2000
%         n_month = 11;
%     else
%         n_month = 12;
%     end
%     
%     
%     MOD11_LST_ann(:,:,yy-1999) = nanmean(MOD11_LST(:,:,start_m: start_m + n_month -1),3);
%     start_m = start_m + n_month;
%     
% end


%%%%%%%% END: MODIS LST %%%%%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%% load in TRENDY v6 (1901-2016) %%%%%
dataHome = '/Volumes/RemiLBNL/project13_ModelGPPsen/data/TRENDYv6/';

Trendy_v6.Models={'CLM4.5','DLEM','ISAM','JULES','LPJ-wsl','ORCHIDEE','ORCHIDEE-MICT',...
    'SDGVM','VEGAS','VISIT'};

scenarios = {'S0','S1','S2','S3'};


areaGrid2 = rot90(flip(areaGrid),3);


%%%%% load in NBP;
for ss = 1:4 
    
    for iii = 1: length(Trendy_v6.Models)
        
        filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_nbp.mat');
        
        %filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_gpp.mat');
        
        
        % gridded monthly NBP
        gridded_NBP = load(filename);
        tmp1 = gridded_NBP.nbp_05;
        
        %tmp1 = gridded_NBP.gpp_05;
        
        tmp1(tmp1<-1000) = NaN;
        
        tmp2 = [];
        % gridded annual NBP
        for yyy = 1:116 % for each year, only store annual values.
        
            %convert unit from kg/m2/month to g/m2/yr, from 1901-2016, 116 years in total

            tmp2(:,:,yyy) = nansum(tmp1(:,:,(yyy-1)*12+1:yyy*12),3).*1000;
        end
        
        % annual NBP
        Trendy_v6.NBP(:,ss,iii) = squeeze(nansum(nansum(tmp2.*areaGrid2,1),2))./(10^15);       
        Trendy_v6.NBPDT(:,ss,iii) = detrend(Trendy_v6.NBP(end - 57:end,ss,iii)); % only from 1959 - 2016
        
        
        % annual tropical NBP
        Trendy_v6.TNBP(:,ss,iii) = squeeze(nansum(nansum(tmp2(134:226,:,:).*areaGrid2(134:226,:),1),2))./(10^15);       
        Trendy_v6.TNBPDT(:,ss,iii) = detrend(Trendy_v6.TNBP(end - 57:end,ss,iii));
        
    end
    
end


%%%%% load in GPP;
for ss = 1:4 
    
    for iii = 1: length(Trendy_v6.Models)
        
        %filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_nbp.mat');
        
        filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_gpp.mat');
        
        
        % gridded monthly NBP
        gridded_NBP = load(filename);
        %tmp1 = gridded_NBP.nbp_05;
        
        tmp1 = gridded_NBP.gpp_05;
        
        tmp1(tmp1<-1000) = NaN;
        
        tmp2 = [];
        % gridded annual NBP
        for yyy = 1:116 % for each year, only store annual values.
        
            %convert unit from kg/m2/month to g/m2/yr, from 1901-2016, 116 years in total

            tmp2(:,:,yyy) = nansum(tmp1(:,:,(yyy-1)*12+1:yyy*12),3).*1000;
        end
        
        % annual GPP
        Trendy_v6.GPP(:,ss,iii) = squeeze(nansum(nansum(tmp2.*areaGrid2,1),2))./(10^15);       
        Trendy_v6.GPPDT(:,ss,iii) = detrend(Trendy_v6.GPP(end - 57:end,ss,iii)); % only from 1959 - 2016
        
        
        % annual tropical NBP
        Trendy_v6.TGPP(:,ss,iii) = squeeze(nansum(nansum(tmp2(134:226,:,:).*areaGrid2(134:226,:),1),2))./(10^15);       
        Trendy_v6.TGPPDT(:,ss,iii) = detrend(Trendy_v6.TGPP(end - 57:end,ss,iii));
        
    end
    
end


%%%%% load in LAI;
for ss = 1:4 
    
    for iii = 1: length(Trendy_v6.Models)
        
        %filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_nbp.mat');
        
        filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{ss},'_lai.mat');
        
        if isfile(filename)
        
            % gridded monthly NBP
            gridded_NBP = load(filename);

            tmp1 = gridded_NBP.lai_05;
            
            tmp1(tmp1<-1000) = NaN;

            tmp2 = [];
            % gridded annual NBP
            for yyy = 1:116 % for each year, only store annual values.

               %from 1901-2016, 116 years in total

                tmp2(:,:,yyy) = nansum(tmp1(:,:,(yyy-1)*12+1:yyy*12),3);
            end

            % annual LAI (area weighted)
            Trendy_v6.LAI(:,ss,iii) = squeeze(nansum(nansum(tmp2.*areaGrid2,1),2))./(nansum(nansum(areaGrid2,1),2));       
            Trendy_v6.LAIDT(:,ss,iii) = detrend(Trendy_v6.LAI(end - 57:end,ss,iii)); % only from 1959 - 2016


            % annual tropical LAI
            Trendy_v6.TLAI(:,ss,iii) = squeeze(nansum(nansum(tmp2(134:226,:,:).*areaGrid2(134:226,:),1),2))./(nansum(nansum(areaGrid2(134:226,:),1),2));       
            Trendy_v6.TLAIDT(:,ss,iii) = detrend(Trendy_v6.TLAI(end - 57:end,ss,iii));
        end
        
    end
    
end



%%%%%% load LAI3g

load('/Users/xzluo/Dropbox/Data/lai3g_05m'); %lai3g_05, from 1981.5 to 2016

tmp2 = [];
for yyy = 1:35 %for each year, only store annual values, from 1982 to 2016

    tmp2(:,:,yyy) = nanmean(lai3g_05m(:,:,(yyy-1)*12 + 7 : yyy*12 + 6),3); % or use nansum?
    
end

aData.lai3g = squeeze(nansum(nansum(tmp2.*areaGrid2,1),2))./(nansum(nansum(areaGrid2,1),2)); 
aData.Tlai3g = squeeze(nansum(nansum(tmp2(134:226,:,:).*areaGrid2(134:226,:),1),2))./(nansum(nansum(areaGrid2(134:226,:),1),2));       


save('/Volumes/RemiLBNL/project7/Data/project7data2.mat','aData','apData','areaGrid2','CMIP5','env','FC','GCPdata','Trendy_v6','TWS','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% specify a group of variables that link 2D and 1D dataset

% read in those index that change 2D to 1D
lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;
m=length(lat);
n=length(lon);
    
lonlat2=[];
for ii = 1:m
    latz=lat(ii)*ones(n,1);
    lonz=lon';
    tmp=[lonz,latz];
    lonlat2=vertcat(lonlat2,tmp);
end
    
% find location of land in whole grid
filename=strcat('/Volumes/RemiLBNL/project6/data4Remi/Pmodel_data/LUE_out_m.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename);
lonlatland=lonlat;
indXClimate=ismember(lonlat2,lonlatland,'rows');

%load the 0.5 degree LC
filename=strcat('/Volumes/RemiLBNL/project6/data4remi/LCMODIS.mat');
load(filename)

ppp_names={'ENF','MF','DF','EBF','SAV','GRA','SH','CRO'};
ppp_nums=[1,5,4,2,9,10,7,12];

LC05D=landCoverPadded;
LC05D(LC05D==6)=7; %merge CSH and OSH
LC05D(LC05D==8)=9; %merge WSA and SAV
LC05D(LC05D==3)=4; %merge DNF and DBF

% load in area grid
filename_area=strcat('/Volumes/RemiLBNL/project6/code4Remi/areaGrid.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename_area);
areaGrid2 = areaGrid;

% tmp = reshape(areaGrid2(areaGrid2>0),[],1);
% area1D = tmp;

filename_area=strcat('/Volumes/RemiLBNL/project6/code4Remi/area.mat'); % the file 'LUE_out_m2' only has 31 years

load(filename_area);
area1D = area;

LC1D = reshape(LC05D(areaGrid2>0),[],1); % change LC map to 1 dimensional map

% only for the tropics
LC_Tropical = nan(size(LC1D));
LC_Tropical(LC1D == 2) = 1; % EBF
LC_Tropical(LC1D == 7 | LC1D == 9) = 2; % SH and SAV, (arid)
LC_Tropical(LC1D == 12 | LC1D == 14) = 3; % CRO
% 

%total_area = nansum(squeeze(area1D(indX,:)),1);
% %%%%% end: specify a group of variables that link 2D and 1D dataset














%%%%%% get monthly PP data to study drought.
% concantenate monthly data
inputdirClim= '/Volumes/RemiLBNL/project7/Data/CRU4.01/';
env.temp_month = [];
env.pre_month = [];
years = 1959:2016;

for year=years
    
    yearstr=num2str(year);
    
    % load temp (is in degrees C)
    filename=strcat(inputdirClim,'tmp2/cru_tmp_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.temp=tmp.(names{:}); 
    
    env.temp_month = horzcat(env.temp_month,env.temp(:,3:end));
    
    % load PP
    filename=strcat(inputdirClim,'pre2/cru_pre_',yearstr,'.mat');
    tmp=load(filename);
    names=fieldnames(tmp);
    env.pre=tmp.(names{:}); 
    env.pre_month = horzcat(env.pre_month,env.pre(:,3:end));

end


% deseason and detrend TT
env.dd_temp = [];
for i = 1: size(env.temp_month,1) % for every vegetated pixel
    
    raw_var = env.temp_month(i,1:end)';
    
    seasonal_mean = repmat(mean(reshape(raw_var,[12,58]),2),58,1); % get the monthly mean over the period, assign it to every year
    deseason_val = raw_var-(seasonal_mean-mean(seasonal_mean)); % remove the annual baseline and seasonal baseline, get anomalies
    
    dd_val = detrend(deseason_val); % 
    
    env.dd_temp(i,:) = dd_val;

end


% deseason and detrend precipitation
env.dd_pre = [];
for i = 1: size(env.pre_month,1) % for every vegetated pixel
    
    raw_var = env.pre_month(i,1:end)';
    
    seasonal_mean = repmat(mean(reshape(raw_var,[12,58]),2),58,1); % get the monthly mean over the period
    deseason_val = raw_var-(seasonal_mean-mean(seasonal_mean));
    
    dd_val = detrend(deseason_val); % 
    
    env.dd_pre(i,:) = dd_val;

end

% use three month three 
env.pre_month3 = movmean(env.pre_month,3,2);


%%%%%% quantify drought area, duration (frequency) over the area, intensity over the area for each year. (multi-year just sum it up)
% use precip only

env.drought_duration = [];
env.drought_area = [];
env.drought_intensity = [];

%cutoff_v = [1,5,10,20,30,40,50];
cutoff_v = [1,5,10,15,20,25,30,40,50];

for i = 1: size(env.pre_month,1) % for every pixel;
    
    
    dd_val = reshape(env.pre_month(i,:),[12,58]);
    
    for j = 1:size(cutoff_v,2)
        
        cutoff = prctile(dd_val,cutoff_v(j),2);
    
        ind_drought = dd_val < cutoff; % water deficit period;

        % drought duration (0-12 months)
        drought_duration = nansum(ind_drought,1);

        % drought area (0/1)
        drought_area = drought_duration > 0; % if there is a drought that happens

        % drought intensity (distance from the mean for every drought month)
        % tmp2 = dd_val.* ind_drought;

        tmp2 = (dd_val - nanmean(dd_val,1)).* ind_drought; % if use mean, the anomaly of drought could be positive...
        tmp2(~ind_drought) = NaN;
        tmp3 = drought_duration;
        tmp3(tmp3 == 0) = NaN;

        drought_intensity = nansum(tmp2,1)./tmp3;


        % store the results
        env.drought_duration(i,:,j) = drought_duration;
        env.drought_area(i,:,j) = drought_area.*area1D(i,1); % double check area1D to make sure its value makes sense
        env.drought_intensity(i,:,j) = drought_intensity;
        
        
    end
    
    
end


%commented on Aug 21, 2021. save('/Volumes/RemiLBNL/project7/Data/project7data3.mat','aData','apData','areaGrid2','CMIP5','env','FC','GCPdata','Trendy_v6','TWS','-v7.3');
save('/Volumes/RemiLBNL/project7/Data/env3.mat','env','-v7.3');











% %%%%%%%%%% reorganize other precip data, to get drought from other sources
% 
% 
% target_coord = env.temp(:,1:2); %lon and lat
% org_PRECL = nan(size(env.pre_month));
% org_GPCC = nan(size(env.pre_month));
% org_UDEL = nan(size(env.pre_month));
% 
% %%%% PRECL
% lat_1degree = -89.5:1:89.5;
% lon_1degree = -179.5:1:179.5;
% [tmp1, tmp2] = meshgrid(lon_1degree,lat_1degree);
% 
% temp = horzcat(reshape(tmp1,[],1),reshape(flip(tmp2),[],1)); % because it is one degree resolution
% 
% for i = 1:67420
%     
%     tmp = temp - target_coord(i,:);
%     [a,b] = min(abs(tmp(:,1)) + abs(tmp(:,2)));
%     
%     org_PRECL(i,:) = apData.PRECL_m(b(1),:);
%     
% end
% 
% 
% 
% %%%% GPCC
% lat_05degree = -89.75:0.5:89.75;
% lon_05degree = -179.75:0.5:179.75;
% [tmp1, tmp2] = meshgrid(lon_05degree,lat_05degree);
% 
% temp = horzcat(reshape(tmp1,[],1),reshape(flip(tmp2),[],1)); % because it is one degree resolution
% 
% for i = 1:67420
%     
%     tmp = temp - target_coord(i,:);
%     [a,b] = min(abs(tmp(:,1)) + abs(tmp(:,2)));
%     
%     org_GPCC(i,:) = apData.GPCC_m(b(1),:);
%     
% end
% 
% %%%% UDEL
% lat_05degree = -89.75:0.5:89.75;
% lon_05degree = -179.75:0.5:179.75;
% [tmp1, tmp2] = meshgrid(lon_05degree,lat_05degree);
% 
% temp2 = horzcat(reshape(tmp1,[],1),reshape(flip(tmp2),[],1));
% 
% %%%
% rawdata = importdata('/Volumes/RemiLBNL/project7/more_precipitation_dataset/UDEL/precip.2016');
% lonlatland = rawdata(:,1:2); % the lon and lat of only land in UDEL
% 
% lat=-89.75:0.5:89.75;
% lon=-179.75:0.5:179.75;
% m=length(lat);
% n=length(lon);
%     
% lonlat2=[];
% for ii = 1:m
%     latz=lat(ii)*ones(n,1);
%     lonz=lon';
%     tmp=[lonz,latz];
%     lonlat2=vertcat(lonlat2,tmp);
% end
% 
% indXLAND=ismember(lonlat2,lonlatland,'rows');
% %%%
% 
% temp = temp2(indXLAND,:);
% 
% for i = 1:67420
%     
%     tmp = temp - target_coord(i,:);
%     [a,b] = min(abs(tmp(:,1)) + abs(tmp(:,2)));
%     
%     org_UDEL(i,:) = apData.UDEL_m(b(1),:);
%     
% end
% 
%%%%%%%%%% End: reorganize other precip data, to get drought from other sources








%%%%%%%%% Calculate the precentile using other datasources %%%%%%
cutoff_v = [1,5,10,15,20,25,30,40,50];


% PRECL
apData.drought_duration_PRECL = [];
apData.drought_area_PRECL = [];
apData.drought_intensity_PRECL = [];

for i = 1: size(org_PRECL,1) % for every pixel;
    
    
    dd_val = reshape(org_PRECL(i,:),[12,58]);
    
    for j = 1:size(cutoff_v,2)
        
        cutoff = prctile(dd_val,cutoff_v(j),2);
    
        ind_drought = dd_val < cutoff; % water deficit period;

        % drought duration (0-12 months)
        drought_duration = nansum(ind_drought,1);

        % drought area (0/1)
        drought_area = drought_duration > 0; % if there is a drought that happens

        % drought intensity (distance from the mean for every drought month)
        % tmp2 = dd_val.* ind_drought;

        tmp2 = (dd_val - nanmean(dd_val,1)).* ind_drought; % if use mean, the anomaly of drought could be positive...
        tmp2(~ind_drought) = NaN;
        tmp3 = drought_duration;
        tmp3(tmp3 == 0) = NaN;

        drought_intensity = nansum(tmp2,1)./tmp3;


        % store the results
        apData.drought_duration_PRECL(i,:,j) = drought_duration;
        apData.drought_area_PRECL(i,:,j) = drought_area.*area1D(i,1); % double check area1D to make sure its value makes sense
        apData.drought_intensity_PRECL(i,:,j) = drought_intensity;
        
        
    end
    
    
end





% GPCC
apData.drought_duration_GPCC = [];
apData.drought_area_GPCC = [];
apData.drought_intensity_GPCC = [];

for i = 1: size(org_GPCC,1) % for every pixel;
    
    
    dd_val = reshape(org_GPCC(i,:),[12,58]);
    
    for j = 1:size(cutoff_v,2)
        
        cutoff = prctile(dd_val,cutoff_v(j),2);
    
        ind_drought = dd_val < cutoff; % water deficit period;

        % drought duration (0-12 months)
        drought_duration = nansum(ind_drought,1);

        % drought area (0/1)
        drought_area = drought_duration > 0; % if there is a drought that happens

        % drought intensity (distance from the mean for every drought month)
        % tmp2 = dd_val.* ind_drought;

        tmp2 = (dd_val - nanmean(dd_val,1)).* ind_drought; % if use mean, the anomaly of drought could be positive...
        tmp2(~ind_drought) = NaN;
        tmp3 = drought_duration;
        tmp3(tmp3 == 0) = NaN;

        drought_intensity = nansum(tmp2,1)./tmp3;


        % store the results
        apData.drought_duration_GPCC(i,:,j) = drought_duration;
        apData.drought_area_GPCC(i,:,j) = drought_area.*area1D(i,1); % double check area1D to make sure its value makes sense
        apData.drought_intensity_GPCC(i,:,j) = drought_intensity;
        
        
    end
    
    
end




% UDEL
apData.drought_duration_UDEL = [];
apData.drought_area_UDEL = [];
apData.drought_intensity_UDEL = [];

for i = 1: size(org_UDEL,1) % for every pixel;
    
    
    dd_val = reshape(org_UDEL(i,:),[12,58]);
    
    for j = 1:size(cutoff_v,2)
        
        cutoff = prctile(dd_val,cutoff_v(j),2);
    
        ind_drought = dd_val < cutoff; % water deficit period;

        % drought duration (0-12 months)
        drought_duration = nansum(ind_drought,1);

        % drought area (0/1)
        drought_area = drought_duration > 0; % if there is a drought that happens

        % drought intensity (distance from the mean for every drought month)
        % tmp2 = dd_val.* ind_drought;

        tmp2 = (dd_val - nanmean(dd_val,1)).* ind_drought; % if use mean, the anomaly of drought could be positive...
        tmp2(~ind_drought) = NaN;
        tmp3 = drought_duration;
        tmp3(tmp3 == 0) = NaN;

        drought_intensity = nansum(tmp2,1)./tmp3;


        % store the results
        apData.drought_duration_UDEL(i,:,j) = drought_duration;
        apData.drought_area_UDEL(i,:,j) = drought_area.*area1D(i,1); % double check area1D to make sure its value makes sense
        apData.drought_intensity_UDEL(i,:,j) = drought_intensity;
        
        
    end
    
    
end

















%%%%%%%%%%%%%%%% Get lagged precipitation %%%%%%%%%%%%%%%%%%%
% read in PRECL, monthly, 1 degree
filename = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/precip.mon.mean.1x1.nc';

rawdata = ncread(filename,'precip'); % mm/day
lat_range = ncread(filename,'lat');
lat_ind = lat_range < 23 & lat_range > -23;

for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1948).*12 + 1 - 4;
    high = (year-1948).*12 + 12 - 4;
    
    low = max(low,1);
    
    tmp_raw = rawdata(:,lat_ind,low:high).*30;
    
    tmp_raw(tmp_raw<=0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    tmp_sum=nansum(tmp_raw,3);
    tmp_sum(tmp_sum==0)=NaN; % remove those ocean values
    
    apData.PRECLlag(i) = nanmean(reshape(tmp_sum,[],1));
    
    
end





%%%%%%%%% read in GPCC, monthly, 0.5 degree
% from 1891 to 2016
% only use 1959 to 2016
filename = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/precip.mon.total.v2018.nc';
rawdata = ncread(filename,'precip'); % mm/month
lat_range = ncread(filename,'lat');
lat_ind = lat_range < 23 & lat_range > -23;


for year = 1959:2016
    
    i=year-1958;
    
    low = (year-1891).*12 + 1 - 4;
    high = (year-1891).*12 + 12 - 4;
    
    low = max(low,1);
    
    tmp_raw = rawdata(:,lat_ind,low:high);
    
    tmp_raw(tmp_raw<=0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    tmp_sum=nansum(tmp_raw,3);
    tmp_sum(tmp_sum==0)=NaN; % remove those ocean values
    
    apData.GPCClag(i) = nanmean(reshape(tmp_sum,[],1));
    
    
end



%%%%%%%%% read in UDEL

dir_name = '/Volumes/RemiLBNL/project7/more_precipitation_dataset/UDEL/';
cd(dir_name);
files = dir(strcat(dir_name,'*')); % from 1900 to 2017
%load in land mask
% filename=strcat('/Volumes/RemiLBNL/project7/Data/CRU4.01/','alpha/cru_alpha_2000.mat');
% tmp=load(filename);
% names=fieldnames(tmp);
% tmp2=tmp.(names{:}); 
% land_latlon=tmp2(:,1:2);

for i = 62:length(files)-1 %  only use 1959 to 2016 here, also two empty files at the beginning
    
    filename = files(i).name;
    rawdata = importdata(filename);
    
    % previous year
    filename = files(i-1).name;
    rawdata_pre = importdata(filename);
    
    lat_ind = rawdata(:,2) < 23 & rawdata(:,2) > -23; % the tropics
    
    tmp_raw = horzcat(rawdata_pre(lat_ind,11:14),rawdata(lat_ind,3:10));
    tmp_raw(tmp_raw<0)=NaN;
    tmp_raw(tmp_raw>1000)=NaN;
    
    apData.UDELlag(i-61) = nanmean(reshape(nansum(tmp_raw,2),[],1));
    
    
end


apData.GPCClagDT_P=detrend(apData.GPCClag);
apData.PRECLlagDT_P=detrend(apData.PRECLlag);
apData.UDELlagDT_P=detrend(apData.UDELlag);

%%save('/Volumes/RemiLBNL/project7/Data/project7data2.mat','aData','apData','areaGrid2','CMIP5','env','FC','GCPdata','Trendy_v6','TWS','-v7.3');
%save('/Volumes/RemiLBNL/project7/Data/project7data3.mat','aData','apData','areaGrid2','CMIP5','env','FC','GCPdata','Trendy_v6','TWS','-v7.3');

