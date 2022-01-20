% master_figures_v2_r1

% the figures used for revision



%%%%%%%%%%%% Figure R2. + climate variable + LULC from Yue 2020 %%%%%

filename = '/Users/xiangzhongluo/Documents/Work_folder/project7/ORC_LUC/FigshareData/Baseline/dft_eluc_orc_1850_2015.csv';
raw_table = readtable(filename);
tmp = [raw_table.Var1,raw_table.ELUC];

ind = tmp(:,1)>1958;
ELUC_CY = tmp(ind,:);


%%%%%%%%%%%%% Extended Data Fig. 2, STD of many fluxes %%%%%%%%%%%%%%%%%%%%%%%%%

fE2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 12], ...
    'OuterPosition', [2, 2, 25, 12]);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0,'SpacingHoriz',0.15,'MR',0.05,'ML',0.1,'MarginTop',0.05,'MarginBottom',0.2); 


%color_choice=[0,0,0; 141,160,203; 252,141,98;228,26,28 ;102,194,165]./255;

color_choice=[0,0,0; 141,160,203; 252,141,98;178,24,43;102,194,165]./255;

annual_fluxes2 = [];

tmp1 = GCPdata.data.Land0x2DUseChangeEmissions(:,4:5); % LUC emissions from two book keeping model 
tmp2 = GCPdata.data.OceanSink(:,4:11); % refill fc data to 1959-2016
tmp3 = vertcat(nan(23,1),aData.Tlai3g); % 1982-2016, LAI;
tmp4 = vertcat(nan(38,1),fire_em(1:20)'); % 1997-2020 (update to 2016)
tmp5 = vertcat(ELUC_CY(:,2),nan(1,1)); % 1959-2015, extended to 2016

annual_fluxes2 = horzcat(GCPdata.growthRateGCP, GCPdata.oceanUptake, GCPdata.lucEmissions,tmp1,tmp5,tmp2,tmp4,tmp3);
%1,1,1,2,1,

IAV_fluxes2 = [];
nIAV_fluxes2 = [];


lag = 20;

for i = 1:  size(annual_fluxes2,2) % get IAV (per 20 yrs) of each flux
    
    clear tmp;
    
    ydata = annual_fluxes2(:,i);
    
    if i == 1 % only if the CGR that consider NaN
    
        for j = 1:length(ydata)-lag

            tmp(j) = nanstd(ydata(j:j+lag)); 

        end
    
    else % for models are nan values
        for j = 1:length(ydata)-lag
            
             tmp(j) = std(ydata(j:j+lag));

        end
    end
     
      
    IAV_fluxes2(:,i) = tmp;
    
    % normalise IAV by the first available records
    
    if isnan(IAV_fluxes2(1,i))
        nIAV_fluxes2(:,i) = IAV_fluxes2(:,i)./IAV_fluxes2(24,i); % LAI data is shorter
    else
        nIAV_fluxes2(:,i) = IAV_fluxes2(:,i)./IAV_fluxes2(1,i);
    end
%     
end

hold on;

years = 1970:2007;


% line 1: IAV of CGR
xdata = years;
ydata = IAV_fluxes2(:,1);
    
faceColor = color_choice(1,:);
ppx1 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 2: IAV of Ocean uptake
xdata = years;
ydata = nanmean(IAV_fluxes2(:,7:14),2);
ydata_sd = nanstd(IAV_fluxes2(:,7:14),1,2);

faceColor = color_choice(2,:);
ppx2 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(8),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% line 3: IAV of land use emission

xdata = years;
ydata = IAV_fluxes2(:,3);
    
faceColor = color_choice(3,:);
ppx3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 4: IAV of CY land use emission

xdata = years;
ydata = IAV_fluxes2(:,6);
    
faceColor = color_choice(3,:);
ppx3b = plot(xdata,ydata,'--','LineWidth',2.5,'color',faceColor);


set(gca,'box','off');
ylabel('STD_{flux} (Pg C yr^-^1)','FontSize',12);
xlabel('Year','FontSize',12);


% line 4: IAV of fire emissions

lag2 = 10;
tmp = [];

for j = 1:length(fire_em)-lag2

    tmp(j) = nanstd(fire_em(j:j+lag2)); 

end


xdata = 2001:2014;
ydata = tmp;
    
faceColor = color_choice(4,:);
ppx4 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
%ylabel('STD_{Fire}','FontSize',12);
text(1971,1.35,'(a)');

%%%%%%%%% add another axis
yyaxis right
set(gca,'YColor','k');

% line 5: IAV of LAI
xdata = years;
ydata = IAV_fluxes2(:,16);
    
faceColor = color_choice(5,:);
ppx5 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
ylim([0 0.1]);

ylabel('STD_{LAI}','FontSize',12);


%legend([ppx1,ppx2,ppx3,ppx3b,ppx4,ppx5],{'CGR','Ocean','LULC','CY-LULC','Fire','LAI'},'FontSize',9,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);


%%%%% panel 2, add climate variability
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0,'SpacingHoriz',0.15,'MR',0.05,'ML',0.1,'MarginTop',0.05,'MarginBottom',0.2); 


color_choice=[166,97,26;123,50,148;128,205,193;1,133,113]./255;

annual_clim2 = [];

tmp1 = aData.TTDT; % tropical detrended tair 1959-2016
tmp2 = aData.TPPDT; % tropical detrended pp 1959-2016
tmp3 = aData.TVPDDT; % tropical detrended vpd 1959-2016
tmp4 = aData.TSWCDT; % tropical detrended swc 1959-2016


annual_clim2 = horzcat(tmp1,tmp2,tmp3,tmp4);
%1,1,1,2,1,

IAV_clim2 = [];
nIAV_clim2 = [];


lag = 20;

for i = 1:  size(annual_clim2,2) % get IAV (per 20 yrs) of each flux
    
    clear tmp;
    
    ydata = annual_clim2(:,i);
    
    if i == 1 % only if the CGR that consider NaN
    
        for j = 1:length(ydata)-lag

            tmp(j) = nanstd(ydata(j:j+lag)); 

        end
    
    else % for models are nan values
        for j = 1:length(ydata)-lag
            
             tmp(j) = std(ydata(j:j+lag));

        end
    end
     
      
    IAV_clim2(:,i) = tmp;
    
    % normalise IAV by the first available records
    
%     if isnan(IAV_clim2(1,i))
%         nIAV_clim2(:,i) = IAV_clim2(:,i)./IAV_clim2(24,i); % LAI data is shorter
%     else
        nIAV_clim2(:,i) = IAV_clim2(:,i)./IAV_clim2(1,i);
%      end
%     
end

hold on;

years = 1970:2007;


% line 1: IAV of Tair
xdata = years;
ydata = nIAV_clim2(:,1);
    
faceColor = color_choice(1,:);
ppc1 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 2: IAV of pp
xdata = years;
ydata = nIAV_clim2(:,2);
    
faceColor = color_choice(2,:);
ppc2 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 3: IAV of vpd
xdata = years;
ydata = nIAV_clim2(:,3);
    
faceColor = color_choice(3,:);
ppc3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 4: IAV of swc
xdata = years;
ydata = nIAV_clim2(:,4);
    
faceColor = color_choice(4,:);
ppc4 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);


set(gca,'box','off');
ylabel('normalised STD_{climate}','FontSize',12);
xlabel('Year','FontSize',12);

text(1971,1.75,'(b)');


legend([ppx1,ppx2,ppx3,ppx3b,ppx4,ppx5,ppc1,ppc2,ppc3,ppc4],{'CGR','Ocean','LULC','CY-LULC','Fire','LAI','Tair','Precip','VPD','SWC'},'FontSize',8,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);


set(fE2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureE2','-djpeg','-r600');     
     

%%%%%%%%%%%%% End: Extended Data Fig. R2, add in CY2020 paper %%%%%%%%%%%%%%%%%%%%%%%%% 









%load('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/project7data3.mat','env');
%%%%%%%%%%%%%%%%%%% Figure 2. Analyse sensitivity, leads to drought %%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) r2 for IAV_CGR and water indexes (TWS, SWC, precip...)
% (b) Annual drough affected area (month-weighted), cascade bar
% (c) Drought-affected area and IAV_CGR (scatter)
% (d) r2 between drought-affected area and IAV using different threshold

f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 20], ...
    'OuterPosition', [2, 2, 25, 20]);

% sub_panel figure
SpacingVert = 0.1;
SpacingHoriz = 0.11;
MR = 0.05;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.11;


color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
facecolor = color_choice(1,:);

%%%%%%% panel a
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

r_bar = []; % initiate results, store r2, p

%%% prepare the time series of climate data
xtmp = horzcat(vertcat(nan(20,1),TWS.y_tGRACE),...
    aData.tropicalSWC,...
    aData.tropicalPP,apData.GPCCtropicP', apData.PRECLtropicP', apData.UDELtropicP',...
    aData.tropicalVPD,...
    aData.tropicalT,apData.GISStropicT',apData.BESTtropicT',apData.UDELtropicT',apData.CRUTEMtropicT');


%non_auto_tmp = tmp(1:end-1,:) - tmp(2:end,:);

%tmp2 = movmean(tmp,20,'omitnan'); % omit nan values.
tmp2 = movmean(xtmp,20);
xdata = tmp2(11:48,:);
%ydata = gamma_T;
ydata = IAV_fluxes(:,1);

tmp_r1 = [];
tmp_r2 = [];

for  i = 1: size(xtmp,2)
    
    tmp_r1 = corrcoef(xdata(1:end-1,i),xdata(2:end,i),'rows','pairwise');
    tmp_r2(i,1) = tmp_r1(1,2);
    
end

non_auto_xdata = xdata(1:end-1,:) - xdata(2:end,:).*tmp_r2';
non_auto_ydata = ydata(1:end-1,:) - 0.895*ydata(2:end,:);


% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions'); % need to double
% check the location of dw function
% for i = 1:size(tmp,2)
%     mdl = fitlm(1:37,non_auto_xdata(:,i)); 
%     [p1,tmp_DW1] = dwtest(mdl,'exact','both');
%     DW1(i,1)=round(tmp_DW1.*100)./100;
% end
% 
% 
% mdl = fitlm(1:37,non_auto_ydata(:,1)); 
% [p1,tmp_DW1] = dwtest(mdl,'exact','both');
% disp(tmp_DW1);

% remove volcanoe years, newly added
[vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');
ydata(vol_ind) = NaN;



for i = 1: size(xdata,2)
    
    %p=polyfit(xdata(i,:),ydata,1);
    
    [r1,r2]=corrcoef(non_auto_xdata(:,i),non_auto_ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;    
    rsq=round(r1(2)^2*100)./100; 
    
    r_bar(i,1) = rsq;
    r_bar(i,2) = sig;
    
end

% initiate reorgnized bar data
reog_bar = [];
reog_bar(1,2,1) = r_bar(1,1); % r2 for TWS
reog_bar(1,2,2) = NaN; 

reog_bar(2,2,1) = r_bar(2,1); % r2 for SWC
reog_bar(2,2,2) = NaN; 

reog_bar(3,2,1) = nanmean(r_bar(3:6,1)); % r2 for Precip
reog_bar(3,2,2) = nanstd(r_bar(3:6,1))./sqrt(4); 

reog_bar(4,2,1) = r_bar(7,1); % r2 for VPD
reog_bar(4,2,2) = NaN; 

reog_bar(5,2,1) = nanmean(r_bar(8:12,1)); % r2 for T
reog_bar(5,2,2) = nanstd(r_bar(8:12,1))./sqrt(5); 


%%%%%%%% non auto-correlated
for i = 1: size(xdata,2)
    
    %p=polyfit(xdata(i,:),ydata,1);
    
    [r1,r2]=corrcoef(xdata(:,i),ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;    
    rsq=round(r1(2)^2*100)./100; 
    
    r_bar(i,1) = rsq;
    r_bar(i,2) = sig;
    
end

% initiate reorgnized bar data
%reog_bar = [];
reog_bar(1,1,1) = r_bar(1,1); % r2 for TWS
reog_bar(1,1,2) = NaN; 

reog_bar(2,1,1) = r_bar(2,1); % r2 for SWC
reog_bar(2,1,2) = NaN; 

reog_bar(3,1,1) = nanmean(r_bar(3:6,1)); % r2 for Precip
reog_bar(3,1,2) = nanstd(r_bar(3:6,1))./sqrt(4); 

reog_bar(4,1,1) = r_bar(7,1); % r2 for VPD
reog_bar(4,1,2) = NaN; 

reog_bar(5,1,1) = nanmean(r_bar(8:12,1)); % r2 for T
reog_bar(5,1,2) = nanstd(r_bar(8:12,1))./sqrt(5); 


%%%%% plot the bar figure
b = bar(reog_bar(:,:,1),'FaceColor',[166 206 227]/255,'EdgeColor','none');

b(1).FaceColor = 'flat';
for i = 1:3 % water
    b(1).CData(i,:) = [166 206 227]/255;
end

for i = 4
    b(1).CData(i,:) = [206 227 166]/255;
end

for i = 5
    b(1).CData(i,:) = [227 166 206]/255;
end

b(2).EdgeColor = 'flat';
b(2).LineWidth = 1.5;
for i = 1:3 % water
    b(2).CData(i,:) = [166 206 227]/255;
end

for i = 4
    b(2).CData(i,:) = [206 227 166]/255;
end

for i = 5
    b(2).CData(i,:) = [227 166 206]/255;
end

b(2).FaceColor = [255 255 255]/255;


hold on;

%%% error bar for non - auto-corrected
er = errorbar([1:5]-0.15, reog_bar(:,1,1),reog_bar(:,1,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold on;

%%% error bar for auto-corrected
er = errorbar([1:5]+0.15, reog_bar(:,2,1),reog_bar(:,2,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 



%%%% demonstrate significance
text([1:5]-0.2,reog_bar(1:5,1,1)+0.1,'**','FontSize',12);
text(1+0.12,reog_bar(1,2,1)+0.1,'*','FontSize',12);
text(3+0.12,reog_bar(3,2,1)+0.1,'**','FontSize',12);

set(gca, 'box', 'off');

factor_names = {'TWS';...
    'SWC';'MAP';'VPD';...
    'MAT'};

ylim([0 0.9]);

ylabel(['r^2 with ','STD_{CGR}'],'FontSize',12);
h = gca;
h.XTick = 1:length(factor_names);
h.XTickLabel = factor_names;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
%xlabel('Variables','FontSize',12);

text(0.7,0.9,'(a)','fontsize',12);

 
 
%%%%%%% panel b
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

cutoff_v = [1,5,10,15,20,25,30,40,50];
area_drought_20 = []; % iniate the matrix to store area affected by drought
area_drought_1 = [];
r_IAV = [];

for j = 1:size(cutoff_v,2)
    
    indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0; % index of the tropical land

    aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,j)),1);
    aData.drought_area = nansum(squeeze(env.drought_area(indX,:,j)),1);
    aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,j)),1);
    
    %M = movmean(aData.drought_area,20,'omitnan');
    total_area = nansum(squeeze(area1D(indX,:)),1);
    
    
    tmp = nansum(squeeze(env.drought_duration(indX,:,j)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
    aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
    aData.drought_area2 = tmp./total_area/12; % the drought affected area per month
    M = movmean(aData.drought_area2,20,'omitnan');
    
    area_drought_20(:,j) = M(11:48).*100;
    area_drought_1(:,j) = aData.drought_area2.*100;
    
    % get the correlation coefficient from CGR, TRENDY NBP and Fluxcom NBP
    % add TRENDY GPP, fluxcom GPP, LAI3g
    % add TRENDY Reco, fluxcom Reco.
    for mm = 1:size(IAV_fluxes,2)
       
        ydata = IAV_fluxes(:,mm);

        [r1,r2] = corrcoef(M(11:48),ydata,'rows','pairwise');
        sig=round(r2(2)*1000)./1000;    
        rsq=round(r1(2)^2*100)./100; 
        
        %r_IAV(mm,j) = rsq;
        r_IAV(mm,j) = round(r1(2)*100)./100;
    end
    

end

% plot out the cascade data.
color_choice = [165,15,21;203,24,29;251,106,74;252,174,145;254,229,217]./255;


% xdata = [1959:2016,flip(1959:2016)];
% baseline_ydata = zeros(58,1);

xdata = [1969:2006,flip(1969:2006)];
baseline_ydata = zeros(38,1);


for i = 1:4 % 1% 10%, 25%, 50%
    
    switch i
        case 1
            ydata = [log(area_drought_20(:,1));baseline_ydata]';
        case 2
            ydata = [log(area_drought_20(:,3));flip(log(area_drought_20(:,1)))]';
        case 3
            ydata = [log(area_drought_20(:,6));flip(log(area_drought_20(:,3)))]';
        case 4
            ydata = [log(area_drought_20(:,9));flip(log(area_drought_20(:,6)))]';
    end
    
%     if i == 1
%         ydata = [log(area_drought_20(:,i));baseline_ydata]';
%     elseif i == 5
%         ydata = [log(area_drought_20(:,7));flip(log(area_drought_20(:,4)))]';
%     else
%         ydata = [log(area_drought_20(:,i));flip(log(area_drought_20(:,i-1)))]';
%     end
    
    patch(xdata,ydata,color_choice(i+1,:),'EdgeColor','none');
    
end

hold on;
plot(1969:2006, log(area_drought_20(:,[1,3,6,9])),'-','LineWidth',2,'color','w'); % do not plot out 1%, as the area is small

xlim([1965,2010]);
ylim([0,log(50)]);

% set(gca, 'YScale', 'log');
set(gca, 'box', 'off');
h = gca;
h.YTick = [log(1),log(10),log(50)];
h.YTickLabel = {'1','10','50'};

ylabel('Area affected by drought (%)','FontSize',12);
xlabel('Years','FontSize',12);

text(1967,log(49),'(b)','fontsize',12);

% h2=colorbar;%
% set(h2,'YTick',0:2:10);
% set(h2,'YTickLabels',{'1','5','10','20','50'});
% ylabel(h2,{'Drought intensity'},'FontSize',8);



%%%%%%% panel c
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,3,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

xdata = area_drought_20(:,3);
ydata = IAV_fluxes(:,1);
xdata_all=[ones(size(xdata)) xdata];

%scatter(area_drought_20(:,3),IAV_fluxes(:,1));
[b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    
mdl = fitlm(xdata,ydata);

[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = 'k';

sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
hold on;

plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

text(7.5,0.9,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
text(7.5,0.84,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);


% text(-0.1,-0.1,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
% text(-0.1,0.1,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);


%%intercept 0 case

mdl = fitlm(xdata,ydata,'Intercept',false);
[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = [255,50,50]./255;

plot(xdata, xdata.*mdl.Coefficients.Estimate,'--','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

hold on;    

text(6.8,1.2,strcat('y =',num2str(round(mdl.Coefficients.Estimate,2)),'x'),'Color',faceColor,'FontSize',10);
text(6.8,1.14,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color',faceColor,'FontSize',10);

% text(-0.1,-0.1,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
% text(-0.1,0.1,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);


text(6.6,1.39,'(c)','fontsize',12);

xlabel('Area affected by extreme droughts (%)','FontSize',12);
ylabel('STD_{CGR} (PgC yr^{-1})','FontSize',12);


%%%%%%% panel d
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,4,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

color_choice=[0,0,0; 194,165,207;146,197,222]./255;


for mm = 1:3
    
    switch mm % CGR, TRENDY, FLUXCOM
        case 1
            ydata = r_IAV(1,:);
        case 2
%             ydata = nanmean(r_IAV(4:13,:),1);
%             ydata_sd = nanstd(r_IAV(4:13,:),1,1)./sqrt(10);
            
            ydata = nanmean(r_IAV(4:18,:),1);
            ydata_sd = nanstd(r_IAV(4:18,:),1,1)./sqrt(15);
        case 3
%             ydata = nanmean(r_IAV(14:16,:),1);
%             ydata_sd = nanstd(r_IAV(14:16,:),1,1)./sqrt(3);
            
            ydata = nanmean(r_IAV(19:21,:),1);
            ydata_sd = nanstd(r_IAV(19:21,:),1,1)./sqrt(3);
    end
    
    faceColor = color_choice(mm,:);
    
    xdata = cutoff_v;
    ppx(mm) = plot(xdata, ydata ,'-','LineWidth',2.5,'color',faceColor);
    
    hold on;
    
    if mm == 2 | mm ==3
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);
        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
    end
    
end

% ppx(4) = plot(cutoff_v, r_IAV(14,:) ,'--','LineWidth',1.5,'color',[194,165,207]./255);
% ppx(5) = plot(cutoff_v, r_IAV(17,:) ,'--','LineWidth',1.5,'color',[255,165,207]./255);

factor_names = {'50%';'40%';'30%';'';'20%';'';'10%';'5%';'1%'};

%ylabel('STD_{CGR} or STD_{NEE} explained (r^2)','FontSize',12);

ylabel({'r between area affected by extreme'; ' droughts and STD_{CGR} or STD_{NEE}'},'FontSize',11);

h = gca;
%h.XTick = 1:length(factor_names);
h.XTick = xdata;
h.XTickLabel = flip(factor_names);
%h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
xlabel('Drought Intensity','FontSize',12);

%ylim([0 1]);
ylim([-0.5 1]);
text(1,0.95,'(d)','fontsize',12);

set(gca,'box','off');

legend([ppx(1),ppx(2),ppx(3)],{'CGR','DGVMs','FLUXCOM'},'Orientation','vertical','Location','best','Box','off');
%legend([ppx(1),ppx(2),ppx(3),ppx(4),ppx(5)],{'CGR','TRENDY','FLUXCOM','OCN','SDGVM'},'Orientation','vertical','Location','best','Box','off');


set(f2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure2','-djpeg','-r600'); 


%%%%%%%%%%%%%%%%%% End: Figure 2. Analyse sensitivity %%%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%%%%%% Figure S3: autocorrelated correlations %%%%%%%%%%%%%%%%%%%

fs3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 12], ...
    'OuterPosition', [2, 2, 25, 12]);

SpacingHoriz = 0.15;
MarginBottom = 0.2;

% use 8 year interval to get gamma_T

clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
lag=10;

yearGCB = 1960:2016;
ydata_sen=aData.CO2DT(end-56:end); 
xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPDT(end-56:end) aData.TPDT(end-56:end)];

xlim_value=nan(3,1);

% remove volcanoe years
[vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');

ydata_sen(vol_ind) = NaN;
xdata_sen(vol_ind,:) = NaN;

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Users/xiangzhongluo/Dropbox/Code/project7/code');
        r1datalag(i,2)=bootstrap_slope2(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,100);

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
ydatalag=r1datalag(:,1);         
gamma_T = ydatalag;



% use 8 year interval to get gamma_W

clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;

yearGCB = 1960:2016;
ydata_sen=aData.CO2DT(end-56:end); 
xdata_sen=[ones(size(ydata_sen)) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100

xlim_value=nan(3,1);

% remove volcanoe years
[vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');

ydata_sen(vol_ind) = NaN;
xdata_sen(vol_ind,:) = NaN;

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Users/xiangzhongluo/Dropbox/Code/project7/code');
        r1datalag(i,2)=bootstrap_slope2(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,100);

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
ydatalag=r1datalag(:,1);         
gamma_W = ydatalag;


ydata = GCPdata.growthRateGCP; 

IAV_8CGR = [];
lag = 10;
for j = 1:length(ydata)-lag

    IAV_8CGR(j,1) = nanstd(ydata(j:j+lag)); 

end





%%%%%%% panel a: STDCGR and drought affected area

%%% show correlation between r_T and IAV_CGR
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

clear sxx;
%color_choice=[0.8,0,0;0,0,0.8];
xdata = IAV_8CGR(2:end-1,1) - 0.895*IAV_8CGR(3:end,1) - 0.1;
xdata_all=[ones(size(xdata)) xdata];

faceColor = color_choice2(1,:);
for i = 1:1 % plot against r_T 
    
    switch i
        case 1
            ydata = gamma_T(1:end-1,1) - 0.9363*gamma_T(2:end,1);
        case 2
            ydata = gamma_W(1:end-1,1) - 0.8582*gamma_W(2:end,1);
    end
    
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % not use 0.32
    
    mdl = fitlm(xdata,ydata);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);

    
    sxx(i) = scatter(xdata,ydata,10,faceColor,'filled');
    hold on;


    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);

    hold on;

    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');

    hold on;    
    
end


% ylim([0 6.5]);
% xlim([0.8 1.4]);
text(0,-1.25,strcat('r^2=',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'FontSize',11,'color',faceColor);


ylabel(strcat('adj.',{' '}, char(947),'_{CGR}','^T (PgC yr^{-1} K^{-1})'),'FontSize',12);
xlabel(strcat('adj.',{' '}, 'STD_{CGR} (PgC yr^-^1)'),'FontSize',12); 

xlim([-0.25 0.15]);
ylim([-2 2]);
text(-0.23,1.95,'(a)','fontsize',12);




 %%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 faceColor = color_choice2(2,:);
 for i = 2:2 % plot against  r_W
    
    switch i
        case 1
            ydata = gamma_T(1:end-1,1) - 0.9363*gamma_T(2:end,1);
        case 2
            ydata = gamma_W(1:end-1,1) - 0.8582*gamma_W(2:end,1);
    end
    
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % not use 0.32
    
    mdl = fitlm(xdata,ydata);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);

    
    sxx(i) = scatter(xdata,ydata,10,faceColor,'filled');
    hold on;

    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);

    hold on;

    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');

    hold on;    
    
 end

 
set(gca, 'YDir','reverse');
%ylim([-2 1.5]); 
text(0.0,0.8,strcat('r^2=',num2str(round(mdl.Rsquared.Ordinary,2)),', p>0.05'),'FontSize',11,'color',faceColor);

ylabel(strcat('adj.',{' '}, char(947),'_{CGR}','^W (PgC yr^{-1}100 mm^{-1})'),'FontSize',12);
xlabel(strcat('adj.',{' '}, 'STD_{CGR} (PgC yr^-^1)'),'FontSize',12); 
 


legend(sxx,{strcat('adj.',char(947),'_{CGR}','^T'); strcat('adj.', char(947),'_{CGR}','^W')},'location', 'southwest' ,'box','off','FontSize',8);

 


%%%%%%% panel b: STDCGR and drought affected area
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

xdata = area_drought_20(1:end-1,3) - area_drought_20(2:end,3);
ydata = IAV_fluxes(1:end-1,1) - IAV_fluxes(2:end,1);
xdata_all=[ones(size(xdata)) xdata];

%scatter(area_drought_20(:,3),IAV_fluxes(:,1));
[b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    
mdl = fitlm(xdata,ydata);

[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = 'k';

sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
hold on;

plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

text(0.1,-0.12,strcat('y =',num2str(round(b(2),2)),'x + ',num2str(round(b(1),2))),'Color','k','FontSize',10);
text(0.1,-0.15,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);


% text(-0.1,-0.1,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
% text(-0.1,0.1,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);


%%intercept 0 case

mdl = fitlm(xdata,ydata,'Intercept',false);
[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = [255,50,50]./255;

plot(xdata, xdata.*mdl.Coefficients.Estimate,'--','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

hold on;    

text(0.1,-0.18,strcat('y =',num2str(round(mdl.Coefficients.Estimate,2)),'x'),'Color',faceColor,'FontSize',10);
text(0.1,-0.21,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color',faceColor,'FontSize',10);
%text(-0.3,0.13,'autocorrelation corrected','Color',faceColor,'FontSize',10);


text(-0.38,0.14,'(b)','fontsize',12);

xlabel('adj. Area affected by extreme droughts (%)','FontSize',11);
ylabel(strcat('adj.',{' '}, 'STD_{CGR} (PgC yr^{-1})'),'FontSize',11);





set(fs3,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS3','-djpeg','-r600'); 

%%%%%%%%%%%%%%%% Figure S3: autocorrelated correlations %%%%%%%%%%%%%%%%%%%





%%%%%% Figure. S1. For reviewer 3 comments. %%%%%%%%%
%%%%%% examining long-term average ENSO MEI and CGR and STDCGR

% load in ENSO MEI
% source "https://psl.noaa.gov/enso/mei.old/table.html"
load('/Users/xiangzhongluo/Dropbox/Data/MEI_v1.mat');
% convert MEI to annual values, only need 1959 to 2016.
MEI_a = nanmean(MEI(10:67,2:end),2);

tmp = movmean(MEI_a,20);
MEI_20 = tmp(11:48,:);


% load in ENSO NINO34
% source "https://psl.noaa.gov/enso/mei.old/table.html"
load('/Users/xiangzhongluo/Dropbox/Data/nino34.mat');
% convert MEI to annual values, only need 1959 to 2016.
nino34_a = nanmean(nino34(12:69,2:end),2);

tmp = movmean(nino34_a,20);
nino34_20 = tmp(11:48,:);

f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 12], ...
    'OuterPosition', [2, 2, 18, 12]);

% sub_panel figure
SpacingVert = 0.1;
SpacingHoriz = 0.15;
MR = 0.05;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.13;

for i = 1:3
    
    switch i
        case 1
            xdata = MEI_a;
            %xdata = nino34_a;
            ydata = aData.CO2DT;
        case 2
            xdata = MEI_20;
            ydata = IAV_fluxes(:,1);
         case 3
            xdata = MEI_20;
            ydata = area_drought_20(:,3);
    end
    
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    subaxis(1,3,i,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    % linear regression
    xdata_all=[ones(size(xdata)) xdata];
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 

    mdl = fitlm(xdata,ydata);

    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);
    
    [r1,r2] = corrcoef(xdata,ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;   

    % scatter
    faceColor = 'k';
    sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
    hold on;

    % fitted line with uncertainty
    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
    hold on;

    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');
    
    xmin = min(xdata); xmax = max(xdata); xrange = xmax - xmin;
    ymin = min(ydata); ymax = max(ydata); yrange = ymax - ymin;
    
if b(1) < 0
    text(xmin + 0.5*xrange,ymin + 0.1*yrange,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
else
    text(xmin + 0.5*xrange,ymin + 0.1*yrange,strcat('y =',num2str(round(b(2),2)),'x+',num2str(round(b(1),2))),'Color','k','FontSize',10);
end

    text(xmin + 0.5*xrange,ymin + 0.05*yrange,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01 '),'Color','k','FontSize',10);
    text(xmin + 0.05*xrange,ymax,strcat('(',char(96+i),')'),'fontsize',12);
    
    ylim([ymin - 0.1*yrange,ymax + 0.05*yrange]);

    switch i
        case 1
            xlabel('MEI (annual value)','FontSize',12);
            ylabel('\DeltaCGR (PgC yr^{-1})','FontSize',12);
        case 2
            xlabel('MEI (20 yrs average)','FontSize',12);
            ylabel('STD_{CGR} (PgC yr^{-1})','FontSize',12);
        case 3
            xlabel('MEI (20 yrs average)','FontSize',12);
            ylabel('Extreme drought affected area (%)','FontSize',12);
    end
    
    
end

set(f2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS4','-djpeg','-r600'); 


%%%%%%%%% El Nino %%%%%%%%%%%%%%%%%%%%%%%%%%



    
  %%%% Figure S5. check the fluxcom and dgvm on the drought area - STDCGR 

%%%% prepare data
% Global STDNEE
globe_fluxes = [];
tmp1 = nanmean(squeeze(GCPdata.GCPmodels(end-57:end,4:18)),2); % DGVM
tmp2 = nanmean(vertcat(nan(21,3),FC.annualNEE,nan(3,3)),2); % FLUXCOM, refill fc data to 1959-2016

globe_fluxes = horzcat(GCPdata.growthRateGCP,tmp1,tmp2); %CGR, DGVMmean, FLUXCOMmean



% Tropic and Tropical drought area STDNEE

% convert drought duration to 2D

drough_du_2D = nan(360,720,58);

for i = 1: 58
    
    tmp = env.drought_duration(:,i,3);
    tmp2 = squeeze(tmp);    
    tmp3 = nan(360*720,1);
    tmp3(indXClimate) = tmp2;

    duration_dis = reshape(tmp3,720,360); 
    duration_dis = duration_dis./12; % weighted by time 

    % change desert to NaN;
    tmp_LC = rot90(LC05D,3);
    duration_dis(tmp_LC == 16) = NaN;

    drough_du_2D(:,:,i) = rot90(duration_dis,1);
end


% Tropical STDNEE
%tropic_fluxes = [];

tmp1 = squeeze(nanmean(Trendy_v6.NBP2,4)); % DGVM mean NBP;
tmp2 = tmp1(:,:,end-57:end); % DGVM from 1959 to 2016;
tmp3 = tmp2.*areaGrid; % times area
tmp4 = nansum(reshape(tmp3,[],58)); %global NBP
tmp5 = nansum(reshape(tmp3(134:226,:,:),[],58)); %tropical NBP

tmp6 = tmp3.*drough_du_2D;
tmp7 = nansum(reshape(tmp6(134:226,:,:),[],58)); %tropical NBP

globe_fluxes = horzcat(globe_fluxes, tmp4', tmp5', tmp7'); % DGVMs global + tropic + T_droughts

%FC
tmp1 = squeeze(nanmean(FC.NEE2,4)); % DFC mean NBP;
tmp1 = permute(tmp1,[2,1,3]);
tmp2 = tmp1(:,:,:); % FC from 1982 to 2013;
%tmp3 = tmp2.*areaGrid; % times area
tmp3 = tmp2.*1; % in dataprocessing already times area times area
tmp4 = nansum(reshape(tmp3,[],34)); %global NBP
tmp4b = vertcat(nan(21,1),tmp4',nan(3,1));% 1980:2013

tmp5 = nansum(reshape(tmp3(134:226,:,:),[],34)); %tropical NBP
tmp5b = vertcat(nan(21,1),tmp5',nan(3,1));% 1980:2013

tmp6 = tmp3.*drough_du_2D(:,:,22:55);% 1980:2013
tmp7 = nansum(reshape(tmp6(134:226,:,:),[],34)); %tropical NBP
tmp7b = vertcat(nan(21,1),tmp7',nan(3,1));% 1980:2013

globe_fluxes = horzcat(globe_fluxes, tmp4b, tmp5b, tmp7b); % DGVMs global + tropic + FC global + tropic

IAV_fluxes = [];
nIAV_fluxes = [];


lag = 20;

for i = 1:  size(globe_fluxes,2) % get IAV (per 20 yrs) of each flux
    
    clear tmp;
    
    ydata = globe_fluxes(:,i);
    
    if i == 1 % only if the CGR that consider NaN
    
        for j = 1:length(ydata)-lag

            tmp(j) = nanstd(ydata(j:j+lag)); 

        end
    
    else % for models are nan values
        for j = 1:length(ydata)-lag

            tmp(j) = std(ydata(j:j+lag));

        end
    end
        
      
    IAV_fluxes(:,i) = tmp;
    
    % normalise IAV by the first available records
    
    if isnan(IAV_fluxes(1,i))
        nIAV_fluxes(:,i) = IAV_fluxes(:,i)./IAV_fluxes(22,i); % fluxcom data is shorter
    else
        nIAV_fluxes(:,i) = IAV_fluxes(:,i)./IAV_fluxes(1,i);
    end

end


fs5=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 12], ...
    'OuterPosition', [2, 2, 20, 12]);

% sub_panel figure
SpacingVert = 0.1;
SpacingHoriz = 0.15;
MR = 0.05;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.2;

color_choice=[0 ,0 ,0; 253,174,97 ; 215,48,39]./255;

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,3,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

xdata = area_drought_20(:,3);
ydata = IAV_fluxes(:,1);

    xdata_all=[ones(size(xdata)) xdata];
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    mdl = fitlm(xdata,ydata);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);
    [r1,r2] = corrcoef(xdata,ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;   
    % scatter
    faceColor = 'k';
    sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
    hold on;
    % fitted line with uncertainty
    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
    hold on;
    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');
    xmin = min(xdata); xmax = max(xdata); xrange = xmax - xmin;
    ymin = min(ydata); ymax = max(ydata); yrange = ymax - ymin;
    if b(1) < 0
        text(6.6,0.5,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
    else
        text(6.6,0.5,strcat('y =',num2str(round(b(2),2)),'x+',num2str(round(b(1),2))),'Color','k','FontSize',10);
    end
    if sig < 0.01
       text(6.6,0.45,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);
    else
        text(6.6,0.45,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p=',num2str(sig)),'Color','k','FontSize',10);
    end    %text(xmin + 0.05*xrange,ymax,strcat('(',char(96+i),')'),'fontsize',12);
    
    %ylim([ymin - 0.1*yrange,ymax + 0.05*yrange]);
text(6.6,1.28,'(a)');
xlim([6.5 9]);    
ylim([0 1.3]);

ylabel('STD_{CGR} (PgC yr^{-1})','FontSize',10);
xlabel('Extreme drought affected area (%)','FontSize',10);

% DGVMs    
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,3,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


for i = 1:3
    
xdata = area_drought_20(:,3);
ydata = IAV_fluxes(:,3+i)./10^15;

    xdata_all=[ones(size(xdata)) xdata];
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    mdl = fitlm(xdata,ydata);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);
    [r1,r2] = corrcoef(xdata,ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;   
    % scatter
    faceColor = color_choice(i,:);
    sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
    hold on;
    % fitted line with uncertainty
    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
    hold on;
    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');
    xmin = min(xdata); xmax = max(xdata); xrange = xmax - xmin;
    ymin = min(ydata); ymax = max(ydata); yrange = ymax - ymin;
    if b(1) < 0
        text(6.6,0.5-0.14*(i-1),strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color',faceColor,'FontSize',10);
    else
        text(6.6,0.5-0.14*(i-1),strcat('y =',num2str(round(b(2),2)),'x+',num2str(round(b(1),2))),'Color',faceColor,'FontSize',10);
    end
    if sig < 0.01
       text(6.6,0.45-0.14*(i-1),strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color',faceColor,'FontSize',10);
    else
        text(6.6,0.45-0.14*(i-1),strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p=',num2str(sig)),'Color',faceColor,'FontSize',10);
    end
    %text(xmin + 0.05*xrange,ymax,strcat('(',char(96+i),')'),'fontsize',12);
    
    %ylim([ymin - 0.1*yrange,ymax + 0.05*yrange]);
end
text(6.6,1.28,'(b)');
xlim([6.5 9]);
ylim([0 1.3]);

ylabel('DGVMs STD_{NEE} (PgC yr^{-1})','FontSize',10);
xlabel('Extreme drought affected area (%)','FontSize',10);


% FLUXCOM
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,3,3,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


for i = 1:3
    
xdata = area_drought_20(22:35,3); %nan value not considered
ydata = IAV_fluxes(22:35,6+i)./10^15;

    xdata_all=[ones(size(xdata)) xdata];
    [b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    mdl = fitlm(xdata,ydata);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);
    [r1,r2] = corrcoef(xdata,ydata,'rows','pairwise');
    sig=round(r2(2)*1000)./1000;   
    % scatter
    faceColor = color_choice(i,:);
    sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
    hold on;
    % fitted line with uncertainty
    plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
    hold on;
    [xdata_s,I] = sort(xdata);
    h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
    set(h,'facealpha',0.3);
    set(h,'LineStyle','none');
    xmin = min(xdata); xmax = max(xdata); xrange = xmax - xmin;
    ymin = min(ydata); ymax = max(ydata); yrange = ymax - ymin;
    if b(1) < 0
        text(6.6,0.6+0.5-0.14*(i-1),strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color',faceColor,'FontSize',10);
    else
        text(6.6,0.6+0.5-0.14*(i-1),strcat('y =',num2str(round(b(2),2)),'x+',num2str(round(b(1),2))),'Color',faceColor,'FontSize',10);
    end
    if sig < 0.01
       text(6.6,0.6+0.45-0.14*(i-1),strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color',faceColor,'FontSize',10);
    else
        text(6.6,0.6+0.45-0.14*(i-1),strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p=',num2str(sig)),'Color',faceColor,'FontSize',10);
    end
    %text(xmin + 0.05*xrange,ymax,strcat('(',char(96+i),')'),'fontsize',12);
    
    %ylim([ymin - 0.1*yrange,ymax + 0.05*yrange]);
end
text(6.6,1.28,'(c)');
xlim([6.5 9]);
ylim([0 1.3]);

ylabel('FLUXCOM STD_{NEE} (PgC yr^{-1})','FontSize',10);
xlabel('Extreme drought affected area (%)','FontSize',10);

legend([sxx(1),sxx(2),sxx(3)],{'global','tropical','tropical drought affected region'},'FontSize',11,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);


set(fs5,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS5','-djpeg','-r600'); 
%%%% End Figure S5 %%%%%%%





%%%%% Figure R2. add CMIP5 data to plot STD
globe2_fluxes = [];
tmp1 = squeeze(GCPdata.GCPmodels(end-57:end,4:18)); % DGVM
tmp2 = vertcat(nan(21,3),FC.annualNEE,nan(3,3)); % FLUXCOM, refill fc data to 1959-2016
tmp3 = vertcat(CMIP5.annualNBP(end-46:end,:),nan(11,14)); % CMIP5, 1959-2005

globe2_fluxes = horzcat(GCPdata.growthRateGCP,tmp1,tmp2,tmp3); %CGR, DGVM, FLUXCOM, CMIP5

IAV_fluxes = [];
nIAV_fluxes = [];
lag = 20;

for i = 1:  size(globe2_fluxes,2) % get IAV (per 20 yrs) of each flux
    
    clear tmp;
    
    ydata = globe2_fluxes(:,i);
    
    if i == 1 % only if the CGR that consider NaN
    
        for j = 1:length(ydata)-lag

            tmp(j) = nanstd(ydata(j:j+lag)); 

        end
    
    else % for models are nan values
        for j = 1:length(ydata)-lag

            tmp(j) = std(ydata(j:j+lag));

        end
    end
        
      
    IAV_fluxes(:,i) = tmp;
    
    % normalise IAV by the first available records
    
    if isnan(IAV_fluxes(1,i))
        nIAV_fluxes(:,i) = IAV_fluxes(:,i)./IAV_fluxes(22,i); % fluxcom data is shorter
    else
        nIAV_fluxes(:,i) = IAV_fluxes(:,i)./IAV_fluxes(1,i);
    end

end


years = 1970:2007;
%color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
%color_choice=[1 ,0.2 ,0.2; 0,0,0; 0.2,0.2,1; 0,0.8,0];

fr1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 12, 10], ...
    'OuterPosition', [2, 2, 12, 10]);

% sub_panel figure
SpacingVert = 0.1;
SpacingHoriz = 0.15;
MR = 0.05;
ML = 0.15;
MarginTop = 0.03;
MarginBottom = 0.22;

color_choice=[0,0,0;200,100,100;150,150,150; 194,165,207;146,197,222]./255;

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,1,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


hold on;
% line 1: IAV of CGR
xdata = years;
ydata = nIAV_fluxes(:,1);
    
faceColor = color_choice(1,:);
ppx1 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 2: IAV of TRENDY
xdata = years;

ydata = nanmean(nIAV_fluxes(:,2:16),2);
ydata_sd = nanstd(nIAV_fluxes(:,2:16),1,2);

faceColor = color_choice(4,:);
ppx2 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(10),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% line 3: IAV of FLUXCOM
xdata = years;

ydata = nanmean(nIAV_fluxes(:,17:19),2);
ydata_sd = nanstd(nIAV_fluxes(:,17:19),1,2);

faceColor = color_choice(5,:);
ppx3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(3),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% line 4: IAV of CMIP5
xdata = years;

ydata = nanmean(nIAV_fluxes(:,20:end),2);
ydata_sd = nanstd(nIAV_fluxes(:,20:end),1,2);

faceColor = color_choice(2,:);
ppx4 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(14),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


set(gca,'box','off');
ylabel('normalized STD_{CGR} or STD_{NEE}','FontSize',11);
xlabel('Year','FontSize',11);
%text(1971,1.48,'(c)','fontsize',12);

legend([ppx1,ppx2,ppx3,ppx4],{'CGR','DGVMs','FLUXCOM','CMIP5'},'FontSize',10,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);


set(fr1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureR2','-djpeg','-r600'); 
%%%%%%%%%%%%%% Fig R2. Add CMIP5 %%%%%%%%%%%%%






%%%%%% Additional data process, use three month precipitation
% % use three month three 
% env.pre_month3 = movmean(env.pre_month,3,2);
% 
% 
% %%%%%% quantify drought area, duration (frequency) over the area, intensity over the area for each year. (multi-year just sum it up)
% % use precip only
% 
% env.drought_duration = [];
% env.drought_area = [];
% env.drought_intensity = [];
% 
% %cutoff_v = [1,5,10,20,30,40,50];
% cutoff_v = [1,5,10,15,20,25,30,40,50];
% 
% for i = 1: size(env.pre_month3,1) % for every pixel;
%     
%     
%     dd_val = reshape(env.pre_month3(i,:),[12,58]);
%     
%     for j = 1:size(cutoff_v,2)
%         
%         cutoff = prctile(dd_val,cutoff_v(j),2);
%     
%         ind_drought = dd_val < cutoff; % water deficit period;
% 
%         % drought duration (0-12 months)
%         drought_duration = nansum(ind_drought,1);
% 
%         % drought area (0/1)
%         drought_area = drought_duration > 0; % if there is a drought that happens
% 
%         % drought intensity (distance from the mean for every drought month)
%         % tmp2 = dd_val.* ind_drought;
% 
%         tmp2 = (dd_val - nanmean(dd_val,1)).* ind_drought; % if use mean, the anomaly of drought could be positive...
%         tmp2(~ind_drought) = NaN;
%         tmp3 = drought_duration;
%         tmp3(tmp3 == 0) = NaN;
% 
%         drought_intensity = nansum(tmp2,1)./tmp3;
% 
% 
%         % store the results
%         env.drought_duration(i,:,j) = drought_duration;
%         env.drought_area(i,:,j) = drought_area.*area1D(i,1); % double check area1D to make sure its value makes sense
%         env.drought_intensity(i,:,j) = drought_intensity;
%         
%         
%     end
%     
%     
% end


%commented on Aug 21, 2021. save('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/project7data3.mat','aData','apData','areaGrid2','CMIP5','env','FC','GCPdata','Trendy_v6','TWS','-v7.3');
%save('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/env3.mat','env','-v7.3');








%%%%%%%%%%%%%%%%%%% Figure 2. use 3 month average Analyse sensitivity, leads to drought %%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) r2 for IAV_CGR and water indexes (TWS, SWC, precip...)
% (b) Annual drough affected area (month-weighted), cascade bar
% (c) Drought-affected area and IAV_CGR (scatter)
% (d) r2 between drought-affected area and IAV using different threshold

load('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/env3.mat','env');

f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 12], ...
    'OuterPosition', [2, 2, 25, 12]);

% sub_panel figure
SpacingVert = 0.1;
SpacingHoriz = 0.11;
MR = 0.05;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.15;


color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
facecolor = color_choice(1,:);


%%% reassgine variables
cutoff_v = [1,5,10,15,20,25,30,40,50];
area_drought_20 = []; % iniate the matrix to store area affected by drought
area_drought_1 = [];
r_IAV = [];

for j = 1:size(cutoff_v,2)
    
    indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0; % index of the tropical land

    aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,j)),1);
    aData.drought_area = nansum(squeeze(env.drought_area(indX,:,j)),1);
    aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,j)),1);
    
    %M = movmean(aData.drought_area,20,'omitnan');
    total_area = nansum(squeeze(area1D(indX,:)),1);
    
    
    tmp = nansum(squeeze(env.drought_duration(indX,:,j)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
    aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
    aData.drought_area2 = tmp./total_area/12; % the drought affected area per month
    M = movmean(aData.drought_area2,20,'omitnan');
    
    area_drought_20(:,j) = M(11:48).*100;
    area_drought_1(:,j) = aData.drought_area2.*100;
    
    % get the correlation coefficient from CGR, TRENDY NBP and Fluxcom NBP
    % add TRENDY GPP, fluxcom GPP, LAI3g
    % add TRENDY Reco, fluxcom Reco.
    for mm = 1:size(IAV_fluxes,2)
       
        ydata = IAV_fluxes(:,mm);

        [r1,r2] = corrcoef(M(11:48),ydata,'rows','pairwise');
        sig=round(r2(2)*1000)./1000;    
        rsq=round(r1(2)^2*100)./100; 
        
        %r_IAV(mm,j) = rsq;
        r_IAV(mm,j) = round(r1(2)*100)./100;
    end
    

end
%%% reassgine variables



%%%%%%% panel a
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

xdata = area_drought_20(:,3);
ydata = IAV_fluxes(:,1);
xdata_all=[ones(size(xdata)) xdata];

%scatter(area_drought_20(:,3),IAV_fluxes(:,1));
[b,bint,~,~,~] = regress(ydata,xdata_all,0.05); % 
    
mdl = fitlm(xdata,ydata);

[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = 'k';

sxx(i) = scatter(xdata,ydata,15,faceColor,'filled');
hold on;

plot(xdata, xdata.*b(2)+b(1),'-','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

text(9.5,0.9,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
text(9.5,0.84,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);




%%intercept 0 case

mdl = fitlm(xdata,ydata,'Intercept',false);
[ypred,yci] = predict(mdl);
ydata_low=yci(:,1);
ydata_up=yci(:,2);

faceColor = [255,50,50]./255;

plot(xdata, xdata.*mdl.Coefficients.Estimate,'--','color',faceColor);
hold on;

[xdata_s,I] = sort(xdata);
h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], faceColor);
set(h,'facealpha',0.3);
set(h,'LineStyle','none');

hold on;    

text(7.8,1.2,strcat('y =',num2str(round(mdl.Coefficients.Estimate,2)),'x'),'Color',faceColor,'FontSize',10);
text(7.8,1.14,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color',faceColor,'FontSize',10);

% text(-0.1,-0.1,strcat('y =',num2str(round(b(2),2)),'x',num2str(round(b(1),2))),'Color','k','FontSize',10);
% text(-0.1,0.1,strcat('r^2 =',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'Color','k','FontSize',10);

ylim([0.8 1.4]);
text(7.6,1.37,'(a)','fontsize',12);

xlabel('Area affected by extreme droughts (%)','FontSize',12);
ylabel('STD_{CGR} (PgC yr^{-1})','FontSize',12);


%%%%%%% panel d
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

color_choice=[0,0,0; 194,165,207;146,197,222]./255;


for mm = 1:3
    
    switch mm % CGR, TRENDY, FLUXCOM
        case 1
            ydata = r_IAV(1,:);
        case 2
%             ydata = nanmean(r_IAV(4:13,:),1);
%             ydata_sd = nanstd(r_IAV(4:13,:),1,1)./sqrt(10);
            
            ydata = nanmean(r_IAV(4:18,:),1);
            ydata_sd = nanstd(r_IAV(4:18,:),1,1)./sqrt(15);
        case 3
%             ydata = nanmean(r_IAV(14:16,:),1);
%             ydata_sd = nanstd(r_IAV(14:16,:),1,1)./sqrt(3);
            
            ydata = nanmean(r_IAV(19:21,:),1);
            ydata_sd = nanstd(r_IAV(19:21,:),1,1)./sqrt(3);
    end
    
    faceColor = color_choice(mm,:);
    
    xdata = cutoff_v;
    ppx(mm) = plot(xdata, ydata ,'-','LineWidth',2.5,'color',faceColor);
    
    hold on;
    
    if mm == 2 | mm ==3
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);
        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
    end
    
end

% ppx(4) = plot(cutoff_v, r_IAV(14,:) ,'--','LineWidth',1.5,'color',[194,165,207]./255);
% ppx(5) = plot(cutoff_v, r_IAV(17,:) ,'--','LineWidth',1.5,'color',[255,165,207]./255);

factor_names = {'50%';'40%';'30%';'';'20%';'';'10%';'5%';'1%'};

%ylabel('STD_{CGR} or STD_{NEE} explained (r^2)','FontSize',12);

ylabel({'r between area affected by extreme'; ' droughts and STD_{CGR} or STD_{NEE}'},'FontSize',11);

h = gca;
%h.XTick = 1:length(factor_names);
h.XTick = xdata;
h.XTickLabel = flip(factor_names);
%h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
xlabel('Drought Intensity','FontSize',12);

%ylim([0 1]);
ylim([-0.5 1]);
text(1,0.95,'(b)','fontsize',12);

set(gca,'box','off');

legend([ppx(1),ppx(2),ppx(3)],{'CGR','DGVMs','FLUXCOM'},'Orientation','vertical','Location','best','Box','off');
%legend([ppx(1),ppx(2),ppx(3),ppx(4),ppx(5)],{'CGR','TRENDY','FLUXCOM','OCN','SDGVM'},'Orientation','vertical','Location','best','Box','off');


set(f2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure2_3month','-djpeg','-r600'); 


%%%%%%%%%%%%%%%%%% End: use 3 month average Figure 2. Analyse sensitivity %%%%%%%%%%%%%%%%%%%%%%

