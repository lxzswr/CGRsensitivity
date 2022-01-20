% new master figure to revise the sensitivity paper.
% by X. Luo. Sep 6th


%%% Load data %%%%%%%%%
% load('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/project7data.mat'); % acquired by running assessGRtemperatureSens6.
% load('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/project7data_sup.mat');

%/Users/xiangzhongluo/Documents/Work_folder
%/Volumes/RemiLBNL

%/Users/xiangzhongluo/Dropbox
%/Users/xiangzhongluo/Dropbox

load('/Users/xiangzhongluo/Documents/Work_folder/project7/Data/project7data3.mat');

%%%%%%%%%%%%%% Figure S1. T and W sensitivity, but different methods %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Get the rT and rW from many methods %%%%%%%
warning off;

sensCGR = [];

lag = 20; % 20-yr moving window.
% get Tsen
for m = 1:9
    
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1 vif1 rsquare1 aic1;
    
    switch m % swtich between methods
        case 1 % M1 CGR = MAT + MAP + RAD
            yearGCB = 1959:2016;
            ydata_sen=aData.CO2DT;
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT./100 aData.TPDT];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP',' + ', char(947),'_{CGR}','^R',char(916),'RAD');
            
        case 2 % M2 CGR = MAT + MAP + RAD + interactions
            yearGCB = 1959:2016;
            ydata_sen=aData.CO2DT;
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT./100 aData.TPDT normalize(aData.TTDT.*aData.TPPDT./100)];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP',' + ', char(947),'_{CGR}','^R',char(916),'RAD',...
                ' + ', char(947),'_inorm(MATxMAP)');
          
        case 3 % M3 CGR = MAT + GRACE
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end) aData.TPDT(end-36:end)];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS',' + ', char(947),'_{CGR}','^R',char(916),'RAD');

        case 4 % M4 CGR = MAT + MAPlag
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100 aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100            
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}',' + ', char(947),'_{CGR}','^R',char(916),'RAD');
            
        case 5 % M5 CGR = MAT + GRACE + interaction
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end) aData.TPDT(end-36:end) normalize(aData.TTDT(end-36:end)).*normalize(TWS.dt_tGRACE(end-36:end))];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS',' + ', char(947),'_{CGR}','^R',char(916),'RAD',' + ', char(947),'_inorm(MATxTWS)' );

        case 6 % M6 CGR = MAT + MAPlag + interaction
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100 aData.TPDT(end-56:end) normalize(aData.TPPlagDT(end-56:end)./100).*normalize(aData.TTDT(end-56:end))];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}',' + ', char(947),'_{CGR}','^R',char(916),'RAD',' + ', char(947),'_inorm(MATxMAP_{lag})');
          
        case 7 % M7 CGR = MAT; 
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end)];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT');
            
        case 8 % M8 CGR = GRACE; 
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end)];
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^{TWS}',char(916),'TWS');
            
        case 9 % M9 CGR = MAPlag;
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100]; %for pp-based sensitivity, times 100            
            ldtext = strcat(char(916),'CGR =', char(947),'_{CGR}','^W',char(916),'MAP_{lag}');

    end


    % remove volcanoe years
    vol_years=[1991; 1992; 1993];
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

        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable
        tmp = vif(xdata_sen(i:i+lag,2:end));
        if length(tmp)>1
            vif1(i,1) = nanmean(tmp(1:2)); % only check the VIF of the first two
        else
            vif1(i,1) = nanmean(tmp);
        end
        rsquare1(i,1) = mdl.Rsquared.Adjusted;
        aic1(i,1) = mdl.ModelCriterion.AIC;
        
        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
    len_rec = length(r1datalag); % length of available data, as TWS is short
    
    sensCGR.Tsen(1:len_rec,m) = r1datalag(:,1);
    sensCGR.Tsen_std(1:len_rec,m) = r1datalag(:,2);
    sensCGR.Tsig(1:len_rec,m) = sig1(:,1);
    sensCGR.Tvif(1:len_rec,m) = vif1(:,1);
    sensCGR.Trsq(1:len_rec,m) = rsquare1(:,1);
    sensCGR.Taic(1:len_rec,m) = aic1(:,1);
    
    sensCGR.yearlag(1:len_rec,m) = yearlag';
    
    %not use j = 8,9, as it is only Wsens
    
end



ldtext = cell(1,9);

% get Wsen
for m = 1:9
    
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    
    switch m % swtich between methods
        case 1 % M1 CGR = MAT + MAP + RAD
            yearGCB = 1959:2016;
            ydata_sen=aData.CO2DT;
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT./100 aData.TTDT aData.TPDT];
            ldtext{1} = strcat('M1:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP',' + ', char(947),'_{CGR}','^R',char(916),'RAD');
            
        case 2 % M2 CGR = MAT + MAP + RAD + interactions
            yearGCB = 1959:2016;
            ydata_sen=aData.CO2DT;
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT./100 aData.TTDT aData.TPDT normalize(aData.TTDT.*aData.TPPDT./100)];
            ldtext{2} = strcat('M2:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP',' + ', char(947),'_{CGR}','^R',char(916),'RAD',...
                ' + ', char(947),'_inorm(',char(916),'MATx',char(916),'MAP)');
          
        case 3 % M3 CGR = MAT + GRACE
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) aData.TPDT(end-36:end)];
            ldtext{3} = strcat('M3:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS',' + ', char(947),'_{CGR}','^R',char(916),'RAD');

        case 4 % M4 CGR = MAT + MAPlag
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen))  aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100            
            ldtext{4} = strcat('M4:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}',' + ', char(947),'_{CGR}','^R',char(916),'RAD');
            
        case 5 % M5 CGR = MAT + GRACE + interaction
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) aData.TPDT(end-36:end) normalize(aData.TTDT(end-36:end)).*normalize(TWS.dt_tGRACE(end-36:end))];
            ldtext{5} = strcat('M5:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS',' + ', char(947),'_{CGR}','^R',char(916),'RAD',' + ', char(947),'_inorm(',char(916),'MATx',char(916),'TWS)' );

        case 6 % M6 CGR = MAT + MAPlag + interaction
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end) normalize(aData.TPPlagDT(end-56:end)./100).*normalize(aData.TTDT(end-56:end))];
            ldtext{6} = strcat('M6:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}',' + ', char(947),'_{CGR}','^R',char(916),'RAD',' + ', char(947),'_inorm(',char(916),'MATx',char(916),'MAP_{lag})');
          
        case 7 % M7 CGR = MAT; 
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end)];
            ldtext{7} = strcat('M7:',char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT');
            
        case 8 % M8 CGR = GRACE; 
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end)];
            ldtext{8} = strcat('M8:',char(916),'CGR =', char(947),'_{CGR}','^{TWS}',char(916),'TWS');
            
        case 9 % M9 CGR = MAPlag;
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100]; %for pp-based sensitivity, times 100            
            ldtext{9} = strcat('M9:',char(916),'CGR =', char(947),'_{CGR}','^W',char(916),'MAP_{lag}');

    end


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

        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
    len_rec = length(r1datalag);
    
    sensCGR.Wsen(1:len_rec,m) = r1datalag(:,1);
    sensCGR.Wsen_std(1:len_rec,m) = r1datalag(:,2);
    sensCGR.Wsig(1:len_rec,m) = sig1(:,1);
    sensCGR.yearlag(1:len_rec,m) = yearlag';
    
    % not use j = 7, as it is only Tsen
    
end


% start plotting the figures
f1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 40], ...
    'OuterPosition', [2, 2, 18, 40]);

%consider the effect of volcaneo
%vol_years=[1963; 1964; 1982; 1985; 1991; 1992; 1993]; %according to P.Cox
vol_years=[1991; 1992; 1993];

% sub_panel figure
SpacingVert = 0.18;
SpacingHoriz = 0.11;
MR = 0.10;
ML = 0.08;
MarginTop = 0.035;
MarginBottom = 0.15;

color_choice=[255,100,100;100,100,100;100,100,255;5,113,176; 150,150,150;0,0,0]./255;


for p = 1:4
    
    % for senT panel
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    subaxis(4,2,(p-1)*2 + 1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    for i = 1:2 % each panel has two lines
        
        if 2*(p-1) + i == 8
            continue;
        end
        
        ydatalag=sensCGR.Tsen(:,2*(p-1) + i);
        ydatalag_sd=sensCGR.Tsen_std(:,2*(p-1) + i);
        sig1 = sensCGR.Tsig(:,2*(p-1) + i);
        xdatalag=sensCGR.yearlag(:,2*(p-1) + i);
        
        % remove 0 values
        ydatalag(ydatalag == 0) = NaN;
        ydatalag_sd(ydatalag_sd == 0) = NaN;
        sig1(sig1 == 0) = NaN;
        xdatalag(xdatalag == 0) = NaN;
        

        faceColor=color_choice(i,:);

        hold on;

        pxx1(i)=plot(xdatalag,ydatalag,'-o','LineWidth',1.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor,'MarkerSize',3);
        plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',1,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w','MarkerSize',3);
        %plot(xdatalag(sig1 > 0.1),ydatalag(sig1 > 0.1), 'o','LineWidth',1,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w','MarkerSize',3);


        %%% plot the uncertainty of y
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
        
        if i == 1 % add the text
            ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}K^{-1})'),'FontSize',10);
            xlabel('Year','FontSize',10);
            set(gca,'box','off');

            xlim([1968 2010]);
            ylim([-1 6]);

            set(gca,'xtick', 1970:10:2010);

            text(1969,5.9,strcat('(',char(96+(p-1)*2 + 1),')'),'Fontsize',12); 
        end
     
                
    end
    
    
      
    % for senW panel
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    subaxis(4,2,p*2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 
    
    ylabel_txt = strcat(char(947),'_{CGR}','^{W} (PgC yr^{-1}100 mm^{-1})');
    
    for i = 1:2 % each panel has two lines
        
        
        if i== 1 && (p == 2 || p == 3 || p == 4) % when use TWS, get a second axis      
            yyaxis right
            set(gca,'YColor','k');
            ylim([-1 0.5]);
            ylabel_txt = strcat(char(947),'_{CGR}','^{TWS} (PgC yr^{-1}TWS^{-1})');        
        elseif i== 2 && (p == 2 || p == 3 || p == 4)
            yyaxis left
            set(gca,'YColor','k');
            ylim([-2 0.5]);
            ylabel_txt = strcat(char(947),'_{CGR}','^{W} (PgC yr^{-1}100 mm^{-1})');
%         elseif i == 1 && p ==1
%             yyaxis left
%             set(gca,'YColor','k');
%             ylim([-2 0.5]);
%             ylabel_txt = strcat(char(947),'_{CGR}','^{W} (PgC yr^{-1}100 mm^{-1})');
        end
        
        tmp = 2*(p-1) + i;
        
        if tmp == 7 % for the last panel, should be 8 and 9 cols
            i = 2;
        elseif tmp == 8
            i = 3;
        end
        
        
        ydatalag=sensCGR.Wsen(:,2*(p-1) + i);
        ydatalag_sd=sensCGR.Wsen_std(:,2*(p-1) + i);
        xdatalag=sensCGR.yearlag(:,2*(p-1) + i);
        sig1 = sensCGR.Wsig(:,2*(p-1) + i);
        
        % remove 0 values
        ydatalag(ydatalag == 0) = NaN;
        ydatalag_sd(ydatalag_sd == 0) = NaN;
        sig1(sig1 == 0) = NaN;
        xdatalag(xdatalag == 0) = NaN;

        faceColor=color_choice(i,:);

        hold on;

        pxx1(i)=plot(xdatalag,ydatalag,'-o','LineWidth',1.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor,'MarkerSize',3);
        plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',1,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w','MarkerSize',3);


        %%% plot the uncertainty of y
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
        
         ylabel(ylabel_txt,'FontSize',10);
         xlabel('Year','FontSize',10);
         set(gca,'box','off');
         
         xlim([1968 2010]);
         
         
         set(gca, 'YDir','reverse');

         set(gca,'xtick', 1970:10:2010);

         

                
    end
    
    if p == 1
       text(1969,-1.2,strcat('(',char(96+(p-1)*2 + 2),')'),'Fontsize',12); 
    else
         text(1969,-2,strcat('(',char(96+(p-1)*2 + 2),')'),'Fontsize',12); 
    end
    
    if p < 4  
        legend([pxx1(1),pxx1(2)],{char(ldtext{2*(p-1)+1}),char(ldtext{2*(p-1)+2})},'FontSize',10,'box','off','Orientation','vertical','Location',[0.25,0.05 + (4-p)*0.245,0.5,0.05]);
    else
        legend([pxx1(1),pxx1(2),pxx1(3)],{char(ldtext{7}),char(ldtext{8}),char(ldtext{9})},'FontSize',10,'box','off','Orientation','vertical','Location',[0.25,0.03,0.5,0.05]);
    end
end


set(f1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure1sense','-djpeg','-r600');     




sensCGR.Trsq(sensCGR.Trsq==0) = NaN;
disp(nanmean(sensCGR.Trsq));

sensCGR.Taic(sensCGR.Taic==0) = NaN;
disp(nanmean(sensCGR.Taic));

sensCGR.Tvif(sensCGR.Tvif==0) = NaN;
disp(nanmean(sensCGR.Tvif));



%%%%%%%%%%%%%%%%%%% End: Figure S1. T and W sensitivity, but different methods %%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%% Figure 1. T and W sensitivity %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% panel (a) and (b) T and W sensitivity and modelled results.
%%%%% panel (c) shows normalised IAVs
%%%%% panel (d) relationship between sensitivity and IAV_CGR.
f1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 20], ...
    'OuterPosition', [2, 2, 25, 20]);

% plot out T_sensitivity

%color_choice=[1,0.4,0.4;0.4,0.4,1;0,0,0;];

color_choice2=[255,100,100;100,100,255;5,113,176; 150,150,150;0,0,0]./255;

color_choice=[0,0,0;200,100,100;150,150,150; 194,165,207;146,197,222]./255;
lag_choice=[15,20,25];

%consider the effect of volcaneo
%vol_years=[1963; 1964; 1982; 1985; 1991; 1992; 1993]; %according to P.Cox
vol_years=[1991; 1992; 1993];

% sub_panel figure
SpacingVert = 0.10;
SpacingHoriz = 0.11;
MR = 0.08;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.15;


cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

% plot observed gamma_T
for ttt=1:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPDT(end-56:end) aData.TPDT(end-56:end)];


        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end) aData.TPDT(end-36:end)];


    end

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
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);
    
    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    
    if ttt == 1
        gamma_T = ydatalag;
    end
end



%%%% plot modelled gamma_T
%calculate TRENDY a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;

% ydata_sen = detrend(squeeze(Trendy_v6.NBP(end-57:end,3,:))).*(-1); % results from S2 simulation, 10 models
% ydata_sen = ydata_sen(end-56:end,:);

ydata_sen = detrend(squeeze(GCPdata.GCPmodels(end-57:end,4:18))).*(-1); % results from S2 simulation, 10 models 4:18
ydata_sen = ydata_sen(end-56:end,:);

yearGCB = 1960:2016;

num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT(end-56:end) aData.TPPDT(end-56:end) aData.TPDT(end-56:end)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

td_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%calculate FluxCom a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=FC.NEEDT.*(1); %NBP to NEE

%ydata_sen = ydata_sen(end-33:end,:);

yearFC=1980:2013;
num_models=size(ydata_sen,2);


%change lag to 20 to get more values, originally it was 10
lag=20;
for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT(end-36:end-3) aData.TPPDT(end-36:end-3) aData.TPDT(end-36:end-3)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearFC(i)+round(lag*0.5);
        
    end

end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

fc_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];
 
%plot a2
hold on;

%trendy
facecolor = color_choice(4,:);
px2=plot(td_a2(:,1),td_a2(:,2),'-','LineWidth',2.5,'color', facecolor);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(td_a2(:,1),td_a2(:,2),td_a2(:,3),{'-'},0);

set(H.patch,'FaceColor', facecolor,'EdgeColor', facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color', facecolor);
set(H.edge,'color','none');
 
%fluxcom
facecolor = color_choice(5,:);
px4=plot(fc_a2(:,1),fc_a2(:,2),'-','LineWidth',2.5,'color',facecolor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(fc_a2(:,1),fc_a2(:,2),fc_a2(:,3),{'-'},0);

set(H.patch,'FaceColor',facecolor,'EdgeColor',facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color',facecolor);
set(H.edge,'color','none');





%%%%% panels axis elements    
 ylabel(strcat(char(947),'_{CGR}','^T',...
     {' '},'or',{' '}, char(947),'_{NEE}','^T',...
     {' '},'(PgC yr^{-1} K^{-1})'),'FontSize',12);
 
 xlabel('Year','FontSize',12);
 set(gca,'box','off');

 %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
 xlim([1968 2010]);
 ylim([0 6.5]);

 set(gca,'xtick', 1970:10:2010);

 text(1969,6.4,'(a)','Fontsize',12); 
     
  
 
 %%%%%%% plot W_sensitivity
 
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

 
for ttt=1:1
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100
            
            
        case 2 % use MAT and TWS
            yearGCB = 1979:2016;
            ydata_sen=aData.CO2DT(end-37:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE aData.TTDT(end-36:end) aData.TPDT(end-36:end)];
            
    end

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
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(1,:);

    hold on;

    pxx1=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %gamma_W = vertcat(nan(1,1),ydatalag);
    gamma_W = ydatalag;
 end

 
%%%%%%%%%% plot modelled gamma_W
%calculate TRENDY a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;

%ydata_sen = detrend(squeeze(Trendy_v6.NBP(end-57:end,3,:))).*(-1); % results from S2 simulation, 10 models
ydata_sen = detrend(squeeze(GCPdata.GCPmodels(end-57:end,4:18))).*(-1);
ydata_sen = ydata_sen(end-56:end,:);

yearGCB = 1960:2016;

num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

td_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%calculate FluxCom a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=FC.NEEDT.*(1); %NBP to NEE

%ydata_sen = ydata_sen(end-33:end,:);

yearFC=1980:2013;
num_models=size(ydata_sen,2);


%change lag to 20 to get more values, originally it was 10
lag=20;
for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TPPDT(end-36:end-3)./100 aData.TTDT(end-36:end-3) aData.TPDT(end-36:end-3)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearFC(i)+round(lag*0.5);
        
    end

end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

fc_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%plot a2
hold on;

%color_choice=[1 ,0.2 ,0.2; 0,0,0; 0.2,0.2,1; 0,0.8,0];

%trendy
facecolor = color_choice(4,:);
pxx2=plot(td_a2(:,1),td_a2(:,2),'-','LineWidth',2.5,'color', facecolor);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(td_a2(:,1),td_a2(:,2),td_a2(:,3),{'-'},0);

set(H.patch,'FaceColor', facecolor,'EdgeColor', facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color', facecolor);
set(H.edge,'color','none');

%fluxcom
facecolor = color_choice(5,:);
pxx3=plot(fc_a2(:,1),fc_a2(:,2),'-','LineWidth',2.5,'color',facecolor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(fc_a2(:,1),fc_a2(:,2),fc_a2(:,3),{'-'},0);

set(H.patch,'FaceColor',facecolor,'EdgeColor',facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color',facecolor);
set(H.edge,'color','none');
 
 
%%%% axis elements of the figure
ylabel(strcat(char(947),'_{CGR}','^W',...
    {' '}, 'or', {' '},char(947),'_{NEE}','^W',...
    {' '},'(PgC yr^{-1}100 mm^{-1})'),'FontSize',12);

xlabel('Year','FontSize',12);
set(gca,'box','off');

set(gca, 'YDir','reverse');
xlim([1968 2010]);
ylim([-2 1.5]);

set(gca,'xtick', 1970:10:2010);
text(1969,-1.9,'(b)','Fontsize',12); 

 
 
 
 
 %%%%%%%%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 for ttt=2:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end);  
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100
            
        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) aData.TPDT(end-36:end)];


    end

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
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    faceColor=color_choice(2,:);
    
    ydatalag(1,1) = NaN;

    hold on;

    pxx1b=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');


    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
 end
 ylabel(strcat(char(947),'_{CGR}','^{TWS} (PgC yr^{-1}TWS^{-1})'),'FontSize',12);
 %xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-0.5 0]); 
 
 

 legend([pxx1,pxx1b,pxx2,pxx3],{'obs.(M1)','obs. (based on TWS; M3)','DGVMs','FLUXCOM'},'FontSize',12,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);
 
 
 
 
%%%%% panel C show normalised IAV
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,3,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

% compile all fluxes
annual_fluxes = [];
tmp1 = squeeze(GCPdata.GCPmodels(end-57:end,4:18));
%tmp1 = squeeze(Trendy_v6.NBP(end-57:end,3,:)); % results from S2 simulation, 10 models
tmp2 = vertcat(nan(21,3),FC.annualNEE,nan(3,3)); % refill fc data to 1959-2016

annual_fluxes = horzcat(GCPdata.growthRateGCP, GCPdata.oceanUptake, GCPdata.lucEmissions,tmp1,tmp2);

% remove volcanoe years
[vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');
annual_fluxe(vol_ind,1) = NaN;


% add GPP to the explored fluxes
tmp1 = squeeze(Trendy_v6.GPP(end-57:end,3,:)); 
tmp2 = vertcat(nan(21,3),FC.annualGPP,nan(3,3)); % refill fc data to 1959-2016
tmp3 = vertcat(nan(23,1),aData.lai3g); % 1982-2016;
annual_fluxes = horzcat(annual_fluxes,tmp1,tmp2,tmp3);

% add Reco to the explored fluxes
tmp1 = squeeze(Trendy_v6.GPP(end-57:end,3,:) - Trendy_v6.NBP(end-57:end,3,:)); %tmp1 = tmp1(end-57:end);
tmp2 = vertcat(nan(21,3),FC.annualGPP + FC.annualNEE,nan(3,3));% 1980:2013
annual_fluxes = horzcat(annual_fluxes,tmp1,tmp2);

% % remove vol years
% annual_fluxes(33:35,1) = NaN;

%%%% add global T, PPlag and TWS
tmp1 = aData.globalT;
tmp2 = aData.globalPP;
tmp3 = vertcat(nan(20,1),TWS.y_GRACE); % 1979-2016

annual_fluxes = horzcat(annual_fluxes,tmp1,tmp2,tmp3);



%%%% add tropical T, PPlag and TWS
tmp1 = aData.tropicalT;
tmp2 = aData.tropicalPP;
tmp3 = vertcat(nan(20,1),TWS.y_tGRACE);
tmp4 = vertcat(nan(1,1),aData.tropicalPPlag);

annual_fluxes = horzcat(annual_fluxes,tmp1,tmp2,tmp3,tmp4);

IAV_fluxes = [];
nIAV_fluxes = [];


lag = 20;

for i = 1:  size(annual_fluxes,2) % get IAV (per 20 yrs) of each flux
    
    clear tmp;
    
    ydata = annual_fluxes(:,i);
    
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
%     
    %nIAV_fluxes(:,i) = IAV_fluxes(:,i)./nanmean(IAV_fluxes(:,i));
end

hold on;

years = 1970:2007;
%color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
%color_choice=[1 ,0.2 ,0.2; 0,0,0; 0.2,0.2,1; 0,0.8,0];


% line 1: IAV of CGR
xdata = years;
ydata = nIAV_fluxes(:,1);
    
faceColor = color_choice(1,:);
ppx1 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

% line 2: IAV of TRENDY
xdata = years;
% ydata = nanmean(nIAV_fluxes(:,4:13),2);
% ydata_sd = nanstd(nIAV_fluxes(:,4:13),1,2);

ydata = nanmean(nIAV_fluxes(:,4:18),2);
ydata_sd = nanstd(nIAV_fluxes(:,4:18),1,2);

% ydata = nanmean(nIAV_fluxes(:,14),2);
% ydata_sd = nanstd(nIAV_fluxes(:,14),1,2);

faceColor = color_choice(4,:);
ppx2 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(10),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% line 3: IAV of FLUXCOM
xdata = years;
% ydata = nanmean(nIAV_fluxes(:,14:16),2);
% ydata_sd = nanstd(nIAV_fluxes(:,14:16),1,2);

ydata = nanmean(nIAV_fluxes(:,19:21),2);
ydata_sd = nanstd(nIAV_fluxes(:,19:21),1,2);

faceColor = color_choice(5,:);
ppx3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(3),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% % line 4: IAV of climate
% xdata = years;
% ydata = nIAV_fluxes(:,47); % tropical T
% faceColor = color_choice(2,:);
% ppx1 = plot(xdata,ydata,'--','LineWidth',2.5,'color',faceColor);
% 
% 
% xdata = years;
% ydata = nIAV_fluxes(:,48); % tropical W
% faceColor = color_choice(3,:);
% ppx1 = plot(xdata,ydata,'--','LineWidth',2.5,'color',faceColor);


set(gca,'box','off');
ylabel('normalized STD_{CGR} or STD_{NEE}','FontSize',12);
xlabel('Year','FontSize',12);
text(1971,1.48,'(c)','fontsize',12);


%legend([ppx1,ppx2,ppx3],{'CGR','DGVMs','FLUXCOM'},'Orientation','vertical','Position',[0.6    0.6    0.1    0.1],'Box','off');


%%% show correlation between r_T and IAV_CGR
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,4,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

clear sxx;
%color_choice=[0.8,0,0;0,0,0.8];
xdata = IAV_fluxes(2:end,1);
xdata_all=[ones(size(xdata)) xdata];

faceColor = color_choice2(1,:);
for i = 1:1 % plot against r_T 
    
    switch i
        case 1
            ydata = gamma_T;
        case 2
            ydata = gamma_W;
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


ylim([0 6.5]);
xlim([0.8 1.4]);
text(1.18,5.8,strcat('r^2=',num2str(round(mdl.Rsquared.Ordinary,2)),', p<0.01'),'FontSize',12,'color',faceColor);


ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1} K^{-1})'),'FontSize',12);
xlabel({'STD_{CGR} (PgC yr^-^1)'},'FontSize',12); 
text(0.82,6.4,'(d)','fontsize',12);



 %%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 faceColor = color_choice2(2,:);
 for i = 2:2 % plot against  r_W
    
    switch i
        case 1
            ydata = gamma_T;
        case 2
            ydata = gamma_W;
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
%xlim([0.8 1.4]);
ylim([-2 1.5]); 
text(1.18,0.5,strcat('r^2=',num2str(round(mdl.Rsquared.Ordinary,2)),', p>0.05'),'FontSize',12,'color',faceColor);

ylabel(strcat(char(947),'_{CGR}','^W (PgC yr^{-1}100 mm^{-1})'),'FontSize',12);
xlabel({'STD_{CGR} (PgC yr^-^1)'},'FontSize',12); 
 


legend(sxx,{strcat(char(947),'_{CGR}','^T'); strcat(char(947),'_{CGR}','^W')},'location', 'southwest' ,'box','off','FontSize',8);




set(f1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure1','-djpeg','-r600');     
     
  

%%%%%%%%%%%%%%%%%%%%%% End: Figure 1. T and W sensitivity %%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%% Figure R1. T and W sensitivity, with P = 0.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% panel (a) and (b) T and W sensitivity and modelled results.
%%%%% panel (c) shows normalised IAVs
%%%%% panel (d) relationship between sensitivity and IAV_CGR.
f1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 12], ...
    'OuterPosition', [2, 2, 25, 12]);

% plot out T_sensitivity

%color_choice=[1,0.4,0.4;0.4,0.4,1;0,0,0;];

color_choice2=[255,100,100;100,100,255;5,113,176; 150,150,150;0,0,0]./255;

color_choice=[0,0,0;200,100,100;150,150,150; 194,165,207;146,197,222]./255;
lag_choice=[15,20,25];

%consider the effect of volcaneo
%vol_years=[1963; 1964; 1982; 1985; 1991; 1992; 1993]; %according to P.Cox
vol_years=[1991; 1992; 1993];

% sub_panel figure
SpacingVert = 0.10;
SpacingHoriz = 0.11;
MR = 0.08;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.25;


cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

% plot observed gamma_T
for ttt=1:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPDT(end-56:end) aData.TPDT(end-56:end)];


        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end) aData.TPDT(end-36:end)];


    end

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
    %ydatalag_sd=r1datalag(:,2);
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);
    
    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.1),ydatalag(sig1 > 0.1), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    
    if ttt == 1
        gamma_T = ydatalag;
    end
end



%%%% plot modelled gamma_T
%calculate TRENDY a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;

% ydata_sen = detrend(squeeze(Trendy_v6.NBP(end-57:end,3,:))).*(-1); % results from S2 simulation, 10 models
% ydata_sen = ydata_sen(end-56:end,:);

ydata_sen = detrend(squeeze(GCPdata.GCPmodels(end-57:end,4:18))).*(-1); % results from S2 simulation, 10 models 4:18
ydata_sen = ydata_sen(end-56:end,:);

yearGCB = 1960:2016;

num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT(end-56:end) aData.TPPDT(end-56:end) aData.TPDT(end-56:end)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

td_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%td_a2=[yearlag',en_mean_a2,en_sd_a2];

%calculate FluxCom a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=FC.NEEDT.*(1); %NBP to NEE

%ydata_sen = ydata_sen(end-33:end,:);

yearFC=1980:2013;
num_models=size(ydata_sen,2);


%change lag to 20 to get more values, originally it was 10
lag=20;
for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT(end-36:end-3) aData.TPPDT(end-36:end-3) aData.TPDT(end-36:end-3)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearFC(i)+round(lag*0.5);
        
    end

end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

fc_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%fc_a2=[yearlag',en_mean_a2,en_sd_a2];
 
%plot a2
hold on;

%trendy
facecolor = color_choice(4,:);
px2=plot(td_a2(:,1),td_a2(:,2),'-','LineWidth',2.5,'color', facecolor);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(td_a2(:,1),td_a2(:,2),td_a2(:,3),{'-'},0);

set(H.patch,'FaceColor', facecolor,'EdgeColor', facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color', facecolor);
set(H.edge,'color','none');
 
%fluxcom
facecolor = color_choice(5,:);
px4=plot(fc_a2(:,1),fc_a2(:,2),'-','LineWidth',2.5,'color',facecolor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(fc_a2(:,1),fc_a2(:,2),fc_a2(:,3),{'-'},0);

set(H.patch,'FaceColor',facecolor,'EdgeColor',facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color',facecolor);
set(H.edge,'color','none');





%%%%% panels axis elements    
 ylabel(strcat(char(947),'_{CGR}','^T',...
     {' '},'or',{' '}, char(947),'_{NEE}','^T',...
     {' '},'(PgC yr^{-1} K^{-1})'),'FontSize',12);
 
 xlabel('Year','FontSize',12);
 set(gca,'box','off');

 %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
 xlim([1968 2010]);
 ylim([0 6.5]);

 set(gca,'xtick', 1970:10:2010);

 text(1969,6.4,'(a)','Fontsize',12); 
     
  
 
 %%%%%%% plot W_sensitivity
 
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

 
for ttt=1:1
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100
            
            
        case 2 % use MAT and TWS
            yearGCB = 1979:2016;
            ydata_sen=aData.CO2DT(end-37:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE aData.TTDT(end-36:end) aData.TPDT(end-36:end)];
            
    end

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
    %ydatalag_sd=r1datalag(:,2);
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(1,:);

    hold on;

    pxx1=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.1),ydatalag(sig1 > 0.1), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %gamma_W = vertcat(nan(1,1),ydatalag);
    gamma_W = ydatalag;
 end

 
%%%%%%%%%% plot modelled gamma_W
%calculate TRENDY a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;

%ydata_sen = detrend(squeeze(Trendy_v6.NBP(end-57:end,3,:))).*(-1); % results from S2 simulation, 10 models
ydata_sen = detrend(squeeze(GCPdata.GCPmodels(end-57:end,4:18))).*(-1);
ydata_sen = ydata_sen(end-56:end,:);

yearGCB = 1960:2016;

num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
    
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

td_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

%td_a2=[yearlag',en_mean_a2,en_sd_a2];

%calculate FluxCom a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=FC.NEEDT.*(1); %NBP to NEE

%ydata_sen = ydata_sen(end-33:end,:);

yearFC=1980:2013;
num_models=size(ydata_sen,2);


%change lag to 20 to get more values, originally it was 10
lag=20;
for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TPPDT(end-36:end-3)./100 aData.TTDT(end-36:end-3) aData.TPDT(end-36:end-3)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearFC(i)+round(lag*0.5);
        
    end

end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

fc_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];
%fc_a2=[yearlag',en_mean_a2,en_sd_a2];

%plot a2
hold on;

%color_choice=[1 ,0.2 ,0.2; 0,0,0; 0.2,0.2,1; 0,0.8,0];

%trendy
facecolor = color_choice(4,:);
pxx2=plot(td_a2(:,1),td_a2(:,2),'-','LineWidth',2.5,'color', facecolor);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(td_a2(:,1),td_a2(:,2),td_a2(:,3),{'-'},0);

set(H.patch,'FaceColor', facecolor,'EdgeColor', facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color', facecolor);
set(H.edge,'color','none');

%fluxcom
facecolor = color_choice(5,:);
pxx3=plot(fc_a2(:,1),fc_a2(:,2),'-','LineWidth',2.5,'color',facecolor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(fc_a2(:,1),fc_a2(:,2),fc_a2(:,3),{'-'},0);

set(H.patch,'FaceColor',facecolor,'EdgeColor',facecolor,'FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color',facecolor);
set(H.edge,'color','none');
 
 
%%%% axis elements of the figure
ylabel(strcat(char(947),'_{CGR}','^W',...
    {' '}, 'or', {' '},char(947),'_{NEE}','^W',...
    {' '},'(PgC yr^{-1}100 mm^{-1})'),'FontSize',12);

xlabel('Year','FontSize',12);
set(gca,'box','off');

set(gca, 'YDir','reverse');
xlim([1968 2010]);
ylim([-2 1.5]);

set(gca,'xtick', 1970:10:2010);
text(1969,-1.9,'(b)','Fontsize',12); 

 
 
 
 
 %%%%%%%%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 for ttt=2:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end);  
            xdata_sen=[ones(size(ydata_sen)) aData.TPPDT(end-56:end)./100 aData.TTDT(end-56:end) aData.TPDT(end-56:end)]; %for pp-based sensitivity, times 100
            
        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) aData.TPDT(end-36:end)];


    end

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
    %ydatalag_sd=r1datalag(:,2);
    ydatalag_sd=r1datalag(:,2)./sqrt(100);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    faceColor=color_choice(2,:);
    
    ydatalag(1,1) = NaN;

    hold on;

    pxx1b=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');


    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
 end
 ylabel(strcat(char(947),'_{CGR}','^{TWS} (PgC yr^{-1}TWS^{-1})'),'FontSize',12);
 %xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-0.5 0]); 
 
 

legend([pxx1,pxx1b,pxx2,pxx3],{'obs.(M1)','obs. (based on TWS; M3)','DGVMs','FLUXCOM'},'FontSize',12,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);
 
set(f1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureR1','-djpeg','-r600');  
%%%%%%%%%%%%%% End, Figure R1 %%%%%%%%%%%





%%%% specify a group of variables that link 2D and 1D dataset

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
filename=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/data4Remi/Pmodel_data/LUE_out_m.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename);
lonlatland=lonlat;
indXClimate=ismember(lonlat2,lonlatland,'rows');

%load the 0.5 degree LC
filename=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/data4remi/LCMODIS.mat');
load(filename)

ppp_names={'ENF','MF','DF','EBF','SAV','GRA','SH','CRO'};
ppp_nums=[1,5,4,2,9,10,7,12];

LC05D=landCoverPadded;
LC05D(LC05D==6)=7; %merge CSH and OSH
LC05D(LC05D==8)=9; %merge WSA and SAV
LC05D(LC05D==3)=4; %merge DNF and DBF

% load in area grid
filename_area=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/areaGrid.mat'); % the file 'LUE_out_m2' only has 31 years
load(filename_area);
areaGrid2 = areaGrid;

% tmp = reshape(areaGrid2(areaGrid2>0),[],1);
% area1D = tmp;

filename_area=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/area.mat'); % the file 'LUE_out_m2' only has 31 years

load(filename_area);
area1D = area; % 1 dimensional area

LC1D = reshape(LC05D(areaGrid2>0),[],1); % change LC map to 1 dimensional map


LC_Tropical = nan(size(LC1D));
LC_Tropical(LC1D == 2) = 1; % EBF
LC_Tropical(LC1D == 7 | LC1D == 9) = 2; % SH and SAV, (arid)
LC_Tropical(LC1D == 12 | LC1D == 14) = 3; % CRO

%%%%% end: specify a group of variables that link 2D and 1D dataset











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


for  i = 1: size(tmp2,2)
    
    tmp_r1 = corrcoef(xdata(1:end-1,i),xdata(2:end,i),'rows','pairwise');
    tmp_r2(i,1) = tmp_r1(1,2);
    
end

non_auto_xdata = xdata(1:end-1,:) - xdata(2:end,:).*tmp_r2';
non_auto_ydata = ydata(1:end-1,:) - 0.895*ydata(2:end,:);


for i = 1:size(tmp2,2)
    mdl = fitlm(1:37,non_auto_xdata(:,i)); 
    [p1,tmp_DW1] = dwtest(mdl,'exact','both');
    DW1(i,1)=round(tmp_DW1.*100)./100;
end


mdl = fitlm(1:37,non_auto_ydata(:,1)); 
[p1,tmp_DW1] = dwtest(mdl,'exact','both');
disp(tmp_DW1);

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
text([1:3]-0.17,reog_bar(1:3,1,1)+0.1,'*','FontSize',12);
text(1+0.12,reog_bar(1,2,1)+0.1,'***','FontSize',12);
text(3+0.12,reog_bar(3,2,1)+0.1,'*','FontSize',12);

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

text(0.7,0.85,'(a)','fontsize',12);

 
 
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

    
  







%%%%%%%%%%%%%%%%% removed later Figure 3. Drought on IAV %%%%%%%%%%%%%%%%%%%%%%%

% how does water influence IAV of CGR?
% two reasons:
% 1. water influence CGR through stomatal control and drought-induced mortality.
% stomatal control is short term control and its sensitivity is largely known.

% 2. implied non-linear response to climate
% not regular climate (as fluxcom considered it, with extremes)
% extreme conditions, where non-linear response is strongest.

% water works through drought.
% intensity of drought defined by percentage of rainfall 1%, 5%, 10%, 20%, 30%, 40%, 50%

% correlation between drought area (experience at least one month of extreme drought, 10%), against drought intensity.
% GPP (LAI) and Reco 


%%%%%% prepare data for convert 2d to 1d, just need to run once
% % concantenate monthly data
% env.temp_month = [];
% env.pre_month = [];
% years = 1959:2016;
% 
% for year=years
%     
%     yearstr=num2str(year);
%     
%     % load temp (is in degrees C)
%     filename=strcat(inputdirClim,'tmp2/cru_tmp_',yearstr,'.mat');
%     tmp=load(filename);
%     names=fieldnames(tmp);
%     env.temp=tmp.(names{:}); 
%     
%     env.temp_month = horzcat(env.temp_month,env.temp(:,3:end));
%     
%     % load PP
%     filename=strcat(inputdirClim,'pre2/cru_pre_',yearstr,'.mat');
%     tmp=load(filename);
%     names=fieldnames(tmp);
%     env.pre=tmp.(names{:}); 
%     env.pre_month = horzcat(env.pre_month,env.pre(:,3:end));
% 
% end
% 
% 
% % deseason and detrend TT
% env.dd_temp = [];
% for i = 1: size(env.temp_month,1) % for every vegetated pixel
%     
%     raw_var = env.temp_month(i,1:end)';
%     
%     seasonal_mean = repmat(mean(reshape(raw_var,[12,58]),2),58,1); % get the monthly mean over the period, assign it to every year
%     deseason_val = raw_var-(seasonal_mean-mean(seasonal_mean)); % remove the annual baseline and seasonal baseline, get anomalies
%     
%     dd_val = detrend(deseason_val); % 
%     
%     env.dd_temp(i,:) = dd_val;
% 
% end
% 
% 
% % deseason and detrend precipitation
% env.dd_pre = [];
% for i = 1: size(env.pre_month,1) % for every vegetated pixel
%     
%     raw_var = env.pre_month(i,1:end)';
%     
%     seasonal_mean = repmat(mean(reshape(raw_var,[12,58]),2),58,1); % get the monthly mean over the period
%     deseason_val = raw_var-(seasonal_mean-mean(seasonal_mean));
%     
%     dd_val = detrend(deseason_val); % 
%     
%     env.dd_pre(i,:) = dd_val;
% 
% end
% 

% 
% 
% %%%% specify a group of variables that link 2D and 1D dataset
% 
% % read in those index that change 2D to 1D
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
% % find location of land in whole grid
% filename=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/data4Remi/Pmodel_data/LUE_out_m.mat'); % the file 'LUE_out_m2' only has 31 years
% load(filename);
% lonlatland=lonlat;
% indXClimate=ismember(lonlat2,lonlatland,'rows');
% 
% %load the 0.5 degree LC
% filename=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/data4remi/LCMODIS.mat');
% load(filename)
% 
% ppp_names={'ENF','MF','DF','EBF','SAV','GRA','SH','CRO'};
% ppp_nums=[1,5,4,2,9,10,7,12];
% 
% LC05D=landCoverPadded;
% LC05D(LC05D==6)=7; %merge CSH and OSH
% LC05D(LC05D==8)=9; %merge WSA and SAV
% LC05D(LC05D==3)=4; %merge DNF and DBF
% 
% % load in area grid
% filename_area=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/areaGrid.mat'); % the file 'LUE_out_m2' only has 31 years
% load(filename_area);
% areaGrid2 = areaGrid;
% 
% % tmp = reshape(areaGrid2(areaGrid2>0),[],1);
% % area1D = tmp;
% 
% filename_area=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/area.mat'); % the file 'LUE_out_m2' only has 31 years
% 
% load(filename_area);
% area1D = area;
% 
% LC1D = reshape(LC05D(areaGrid2>0),[],1); % change LC map to 1 dimensional map
% 
% %%%%% end: specify a group of variables that link 2D and 1D dataset
% 
% 
%%%%%%%%%%% data process script, just need to run once %%%%%%%%%%%%%%%%
% quantify drought area, duration (frequency) over the area, intensity over the area for each year. (multi-year just sum it up)
% use precip only

% env.drought_duration = [];
% env.drought_area = [];
% env.drought_intensity = [];
% 
% cutoff_v = [1,5,10,20,50];
% 
% for i = 1: size(env.pre_month,1) % for every pixel;
%     
%     
%     dd_val = reshape(env.pre_month(i,:),[12,58]);
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
%         drought_area = drought_duration > 0; % if there is a drought that happens, at least 1.
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


% %%%%% calculate the r2 between drought area and nIAV CGR
% 
% r_IAV = []; % iniate the matrix to store results, cutoff, model, r.
% 
% for j = 1:size(cutoff_v,2)
%     
%     indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0; % index of the tropical land
% 
%     aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,j)),1);
%     aData.drought_area = nansum(squeeze(env.drought_area(indX,:,j)),1);
%     aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,j)),1);
%     
%     %M = movmean(aData.drought_area,20,'omitnan');
%     total_area = nansum(squeeze(area1D(indX,:)),1);
%     
%     
%     tmp = nansum(squeeze(env.drought_duration(indX,:,j)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
%     aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
%     aData.drought_area2 = tmp./total_area/12; % the drought affected area per month
%     M = movmean(aData.drought_area2,20,'omitnan');
% %     
%     
%     % get the correlation coefficient from CGR, TRENDY NBP and Fluxcom NBP
%     % add TRENDY GPP, fluxcom GPP, LAI3g
%     % add TRENDY Reco, fluxcom Reco.
%     for mm = 1:size(IAV_fluxes,2)
%        
%         ydata = IAV_fluxes(:,mm);
% 
%         [r1,r2] = corrcoef(M(11:48),ydata,'rows','pairwise');
%         sig=round(r2(2)*1000)./1000;    
%         rsq=round(r1(2)^2*100)./100; 
%         
%         r_IAV(mm,j) = rsq;
%     end
%     
% 
% end
% 
% 
% % plot the lines
% 
% f3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 12, 18], ...
%     'OuterPosition', [2, 2, 12, 18]);
% 
% % sub_panel figure
% SpacingVert = 0.1;
% SpacingHoriz = 0.05;
% MR = 0.05;
% ML = 0.15;
% MarginTop = 0.03;
% MarginBottom = 0.1;
% 
% 
% color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
% 
% % panel a
% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
% subaxis(3,1,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 
% 
% 
% for mm = 1:3
%     
%     switch mm % CGR, TRENDY, FLUXCOM
%         case 1
%             ydata = r_IAV(1,:);
%         case 2
%             ydata = nanmean(r_IAV(4:13,:),1);
%             ydata_sd = nanstd(r_IAV(4:13,:),1,1);
%         case 3
%             ydata = nanmean(r_IAV(14:16,:),1);
%             ydata_sd = nanstd(r_IAV(14:16,:),1,1);
%     end
%     
%     faceColor = color_choice(mm,:);
%     
%     xdata = 1:5;
%     ppx(mm) = plot(1:5, flip(ydata) ,'-','LineWidth',2.5,'color',faceColor);
%     
%     hold on;
%     
%     if mm == 2 | mm ==3
%         cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
%         H=shadedErrorBar(xdata,flip(ydata),flip(ydata_sd),{'-'},0);
%         set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
%         set(H.mainLine,'color','none');
%         set(H.edge,'color','none');
%     end
%     
% end
% 
% factor_names = {'50%';'20%';'10%';'5%';'1%'};
% 
% ylabel('IAV_{CGR} or IAV_{NEE} explained (r2)','FontSize',9);
% h = gca;
% h.XTick = 1:length(factor_names);
% h.XTickLabel = factor_names;
% %h.XTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% xlabel('Drought Intensity (ppt deficit extremes)','FontSize',10);
% 
% ylim([0 1]);
% 
% legend([ppx(1),ppx(2),ppx(3)],{'CGR','TRENDY','FLUXCOM'},'Orientation','vertical','Location','best','Box','off');
% 
% 
% % panel b
% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
% subaxis(3,1,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 
% 
% 
% for mm = 1:3
%     
%     switch mm
%         
%         case 1
%             ydata = nanmean(r_IAV(17:26,:),1);
%             ydata_sd = nanstd(r_IAV(17:26,:),1,1);
%         case 2
%             ydata = nanmean(r_IAV(27:29,:),1);
%             ydata_sd = nanstd(r_IAV(27:29,:),1,1);
%         case 3
%             ydata = r_IAV(30,:);
%     end
%     
%     faceColor = color_choice(mm+1,:);
%     
%     xdata = 1:5;
%     ppx(mm) = plot(1:5, flip(ydata) ,'-','LineWidth',2.5,'color',faceColor);
%     
%     hold on;
%     
%     if mm == 1 | mm ==2
%         cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
%         H=shadedErrorBar(xdata,flip(ydata),flip(ydata_sd),{'-'},0);
%         set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
%         set(H.mainLine,'color','none');
%         set(H.edge,'color','none');
%     end
%     
% end
% 
% factor_names = {'50%';'20%';'10%';'5%';'1%'};
% 
% ylabel('IAV_{GPP} or IAV_{LAI} explained','FontSize',9);
% 
% ylim([0 1]);
% 
% h = gca;
% h.XTick = 1:length(factor_names);
% h.XTickLabel = factor_names;
% %h.XTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% xlabel('Drought Intensity','FontSize',10);
% legend([ppx(1),ppx(2),ppx(3)],{'TRENDY','FLUXCOM','LAI3g'},'Orientation','vertical','Location','best','Box','off');
% 
% 
% % panel 3
% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
% subaxis(3,1,3,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 
% 
% 
% for mm = 1:2
%     
%     switch mm
%         
%         case 1
%             ydata = nanmean(r_IAV(30:39,:),1);
%             ydata_sd = nanstd(r_IAV(30:39,:),1,1);
%         case 2
%             ydata = nanmean(r_IAV(40:42,:),1);
%             ydata_sd = nanstd(r_IAV(40:42,:),1,1);
% 
%     end
%     
%     faceColor = color_choice(mm+1,:);
%     
%     xdata = 1:5;
%     ppx(mm) = plot(1:5, flip(ydata) ,'-','LineWidth',2.5,'color',faceColor);
%     
%     hold on;
%     
%     if mm == 1 | mm ==2
%         cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
%         H=shadedErrorBar(xdata,flip(ydata),flip(ydata_sd),{'-'},0);
%         set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
%         set(H.mainLine,'color','none');
%         set(H.edge,'color','none');
%     end
%     
% end
% 
% factor_names = {'50%';'20%';'10%';'5%';'1%'};
% 
% ylabel('IAV_{Reco} explained','FontSize',9);
% 
% ylim([0 1]);
% 
% h = gca;
% h.XTick = 1:length(factor_names);
% h.XTickLabel = factor_names;
% %h.XTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% xlabel('Drought Intensity','FontSize',10);
% 
% legend([ppx(1),ppx(2)],{'TRENDY','FLUXCOM'},'Orientation','vertical','Location','best','Box','off');
% 
% 
% set(f3,'PaperPositionMode','auto');
% print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure3','-djpeg','-r600'); 
% 
% 
% %%%%%%%%%%%%%%%%% removed later End: Figure 3. Drought on IAV %%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%% Figure 3. Evolution of drough occurenance %%%%%%%%%

% drought area distribution (south america at the beginning, then southeast asia, then move back to south america?)
% contribution of tropical forest and semi-arid region (bar shows the area)

% display drought duration of 10% precip.
% 10% precip indicate drought intensity, map indicate area, then we can only show duration.



f4=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 15], ...
    'OuterPosition', [2, 2, 15, 15]);

% sub_panel figure
SpacingVert = 0.05;
SpacingHoriz = 0.11;
MR = 0.08;
ML = 0.10;
MarginTop = 0.03;
MarginBottom = 0.13;

colorscheme = [254,229,217;...
252,174,145;...
251,106,74;...
222,45,38;...
165,15,21]./255;
% colorscheme = colorscheme./255;
% colorscheme = [1,0,0; 0,1,0];
% colorscheme = autumn(4);

tropical_LC = LC05D;
tropical_LC(tropical_LC ~=2) = 1; % non EBF all semi-arid and arid system
tropical_LC([1:133,227:360],:) = 0; % focus on tropical region

for i = 1:3
    
    switch i
        case 1 %1960-1979
            tmp = env.drought_duration(:,2:21,3); % choose 10% cutoff
            time_title = '1960-1979';
        case 2 %1980-1999
            tmp = env.drought_duration(:,22:41,3);
            time_title = '1980-1999';
        case 3 %1997-2016
            tmp = env.drought_duration(:,39:58,3);
            time_title = '1997-2016';
    end
    
    tmp2 = nanmean(squeeze(tmp),2);
    
    tmp3 = nan(360*720,1);
    
    tmp3(indXClimate) = tmp2;
    
    duration_dis = reshape(tmp3,720,360); 
    
    duration_dis = duration_dis./12; % weighted by time 
    
    % change desert to NaN;
    tmp_LC = rot90(LC05D,3);
    duration_dis(tmp_LC == 16) = NaN;
    
    %duration_dis(duration_dis < 1) = NaN;
    
    %duration_dis([1:133,227:360],:) = NaN;
    
     mapdata=flip(rot90(duration_dis,1));
     mapdata([1:133,227:360],:) = NaN;
     
     %mapdata([1:120,240:360],:) = NaN;
%     
%     duration_dis2 = rot90(duration_dis,1);
%     
%     tropical_LC(duration_dis2 < 1) = NaN;
%     
%     mapdata = flip(tropical_LC);
    
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    subaxis(4,2,(i-1)*2+1:i*2,'SpacingVert',SpacingVert-0.03,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    %mapdata= smooth2a(mapdata,1,1); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;
    
    mapdata_T = mapdata(120:240,:);
    lat_T = lat(120:240);

    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions/m_map');
    m_proj('robinson','lon',[-180 180],'lat',[-23 23]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat_T,mapdata_T); % draw your map here
    %m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colormap(gca,colorscheme);

     caxis([0 0.25]); 
% 
    h=colorbar;%
    set(h,'YTick',0:0.05:0.25);
    set(h,'TickLabels',{'0','5','10','15','20','25'});
    
    if i == 2
        ylabel(h,{'Probability of'; 'extreme droughts (%)'},'FontSize',9);
    end
    
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');cbfreeze(h);

    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions/m_map');
    
    if i < 3
            m_grid('box','on','linewi',0.5,'linest','none', 'tickdir','in','xticklabels',[],'fontsize',7);
    else
            m_grid('box','on','linewi',0.5,'linest','none', 'tickdir','in','fontsize',7);

    end
%     
%     m_grid('box','on','linewi',0.5,'linest','none', 'tickdir','in','fontsize',7);
    
%     if i < 3
%             m_grid('box','on','linewi',0.5,'linest','none', 'tickdir','in','xticklabels',[],'fontsize',7);
%     end
%    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);

    m_text(-170,19,char(96+i),'fontsize',10);
    m_text(-170,10,time_title,'fontsize',10);
    hold off;
    
    
    
end

%%% dynamic of drought affected area in rainforest and semi-arid regions
area_drought = [];
r_IAV = [];

indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0;
total_area = nansum(squeeze(area1D(indX,:)),1);

for rr = 1:5 % for different regions
    
    switch rr
        case 1 % EBF
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 1;
        case 2 % semi-arid system
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 2; %(9 or 7)
        case 3 % tropical America
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<-45 & env.temp(:,1)>-120 & LC_Tropical > 0;
        case 4 % tropical Africa
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<60 & env.temp(:,1)>-45 & LC_Tropical > 0;
        case 5 % tropical Asia
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<180 & env.temp(:,1)>60 & LC_Tropical > 0;
    end
    
    % get all the indicators under 10%
    aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,3)),1); 
    aData.drought_area = nansum(squeeze(env.drought_area(indX,:,3)),1);
    aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,3)),1);
    
%     M = movmean(aData.drought_duration/12,20,'omitnan');M(11:48); display(M(1));display(M(21));display(M(38));
%     M = movstd(aData.drought_duration/12,20,'omitnan');M(11:48); display(M(1));display(M(21));display(M(38));

    
    tmp = nansum(squeeze(env.drought_duration(indX,:,3)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
    aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
    aData.drought_area2 = tmp./total_area/12; % the drought affected area per month
    
    %M = movmean(aData.drought_duration/12,20,'omitnan');M(11:48);
    %M = movstd(aData.drought_duration/12,20,'omitnan');M(11:48);

    
    M = movmean(aData.drought_area2,20,'omitnan');
    M2 = movstd(aData.drought_area2,20,'omitnan');
    
    area_drought(rr,:,1) = M(11:48);
    area_drought(rr,:,2) = M2(11:48);
    
    for mm = 1:size(IAV_fluxes,2)

        ydata = IAV_fluxes(:,mm);

        [r1,r2] = corrcoef(M(11:48),ydata,'rows','pairwise');
        sig=round(r2(2)*1000)./1000;    
        rsq=round(r1(2)^2*100)./100; 
        
        r_IAV(mm,rr) = rsq;
    end
    
    
end

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(4,2,7,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML+0.01,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

color_choice = [90,180,172;216,179,101]./255;

for rr = 1:2
    
    xdata = 1969:2006;
    ydata = squeeze(area_drought(rr,:,1)).*10^2; % to 100%
    ydata_sd = squeeze(area_drought(rr,:,2)).*10^2; % 
    
    hold on;
    faceColor = color_choice(rr,:);
    pxx1(rr)=plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %title_region2{rr} = [title_region{rr},' r = ',num2str(r_IAV(1,rr))];
end

xlim([1969 2010]);
ylim([0 5]);
set(gca,'box','off');
text(1970,5,'d');

ylabel({'Drought'; 'affected area (%)'},'FontSize',10);
xlabel('Year','FontSize',10); 

legend([pxx1],{'T.forests','T.semi-arid'},'Orientation','vertical','Location','southeast','Box','off');

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(4,2,8,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom+0.01); 

bar(r_IAV(1,:),'FaceColor',[150 150 150]/255,'EdgeColor','none');
ylim([0 0.9]);
h = gca;
h.XTick = 1:5;
h.XTickLabel = {'T.forests','T.semi-arid','T.America','T.Africa','T.Asia'};
h.XTickLabelRotation = 40;
set(gca,'XTickLabel',{'T.forests','T.semi-arid','T.America','T.Africa','T.Asia'},'fontsize',9)

ylabel('r^2 with STD_{CGR}','FontSize',10);
text(0.2,0.85,'e');

set(f4,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure4','-djpeg','-r600'); 

%%%%%%%%%% End: Figure 4. Hotspots of drought during the three periods %%%%%%%%%








% Sptial for temporal replacement constraint: the relationship between drought-area and STD
% the extreme drought-induced STD


% Sptial for temporal replacement constraint: the relationship between drought-area and STD
% the extreme drought-induced STD


%%%% get the EC (fixed the intercept at zero?)

    % std_drought_induced = 0.2*area(*100)

%%%% get the drought affected area (tropic, regional, arid, forests)
    LC_Tropical = nan(size(LC1D));
    LC_Tropical(LC1D == 2) = 1; % EBF
    LC_Tropical(LC1D == 7 | LC1D == 9) = 2; % SH and SAV, (arid)
    LC_Tropical(LC1D == 12 | LC1D == 14) = 3; % CRO
    
    
    area_drought = [];
    r_IAV = [];

    indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0;
    total_area = nansum(squeeze(area1D(indX,:)),1);
    
    % add the weight of flux per pixel
    flux_weight2D = transpose(squeeze(nanmean(nanmean(FC.NEE2,3),4)));
    tmp = reshape(flux_weight2D,[],1);
    flux_weight = tmp(indXClimate);
    aData.F_NEE = [];
    

    for rr = 1:6 % for different regions

        switch rr
            case 1 % EBF
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 1;
            case 2 % semi-arid system
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 2; %(9 or 7)
            case 3 % tropical America
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<-45 & env.temp(:,1)>-120 & LC_Tropical > 0;
            case 4 % tropical Africa
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<60 & env.temp(:,1)>-45 & LC_Tropical > 0;
            case 5 % tropical Asia
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<180 & env.temp(:,1)>60 & LC_Tropical > 0;
            case 6 % pan tropic
                indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0;
        end

        % get all the indicators under 10%
        aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,3)),1); 
        aData.drought_area = nansum(squeeze(env.drought_area(indX,:,3)),1);
        aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,3)),1);


        tmp = nansum(squeeze(env.drought_duration(indX,:,3)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
        aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
        aData.drought_area2 = tmp./total_area/12; % the drought affected area per month

        M = movmean(aData.drought_area2,20,'omitnan');
        M2 = movstd(aData.drought_area2,20,'omitnan');

        area_drought(rr,:,1) = M(11:48);
        area_drought(rr,:,2) = M2(11:48);

        %%%% get the drought-induced STD (from emergent contraint)

        std_drought_induced(rr,:,1) = 0.142 * M(11:48) * 100;
        std_drought_induced(rr,:,2) = 0.142 * M2(11:48) * 100;
        
        
        %%%% get the weight by Flux...the mean NEE for each region, also
        %%%% need to consider area and duration
        %aData.F_NEE(rr) = nanmean(flux_weight(indX)); 
        
%         tmp1 = env.drought_duration(indX,:,3) > 0;
%         tmp = nanmean(squeeze(tmp1).*flux_weight(indX),1);

        tmp = nansum(squeeze(env.drought_duration(indX,:,3)).*squeeze(area1D(indX,:)).*flux_weight(indX),1);
        %tmp = nanmean(squeeze(env.drought_duration(indX,:,3)).*flux_weight(indX),1);
        M3 = movmean(tmp,20,'omitnan');
        aData.F_NEE(rr,:) = M3(11:48);
    end
    
    %%%%% adjusted by weighted
    aData.Fweight = aData.F_NEE./aData.F_NEE(6,:);
    
    for rr = 1: 6
        
        std_drought_induced(rr,:,1) = std_drought_induced(6,:,1).* aData.Fweight(rr,:);
        std_drought_induced(rr,:,2) = std_drought_induced(6,:,2).* aData.Fweight(rr,:);
        
    end




    
    
    
%%%%%% ATTENTION: read this one only for the first time run the algorithm    
%%%%% get the NBP of drought affected area (TRENDY)

% %%%%%%%%%%% load in TRENDY v6 (1901-2016) %%%%%
dataHome = '/Users/xiangzhongluo/Documents/Work_folder/project13_ModelGPPsen/data/TRENDYv6/';

Trendy_v6.Models={'CLM4.5','DLEM','ISAM','JULES','LPJ-wsl','ORCHIDEE','ORCHIDEE-MICT',...
    'SDGVM','VEGAS','VISIT'};

scenarios = {'S0','S1','S2','S3'};


areaGrid2 = rot90(flip(areaGrid),3);


%%%%% load in NBP;

    
    for iii = 1: length(Trendy_v6.Models)
        
        filename = strcat(dataHome,Trendy_v6.Models{iii},'_',scenarios{4},'_nbp.mat');
        
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
        
        % NBP spatial
        Trendy_v6.NBP2(:,:,:,iii) = tmp2; % g/m2/yr      
        
    end
%%%%%% ATTENTION: read this one only for the first time run the algorithm    
    



%%%%% get the fire-induced STD

filename = '/Users/xiangzhongluo/Documents/Work_folder/Data_hub/GFED4/GFED4.1s_C.txt';

tmp = readtable(filename,'Delimiter', '\t');
tmp1 = tmp.Var1(18,1);

% from 1997 to 2020
fire_em = [ 3031   2892   2226   1894   1962   2343   2199   2213   2252   2207   2202   1862   1862   2150   1872   2047   1773   2044   2289   1869   1927   1832   2325   1611];
fire_em = fire_em./1000;



%%%%% get the NBP of drought affected area (Fluxcom & TRENDY)
    
    % get the drought affect area every year.
    for yy = 1: 58
        % convert 1d to 2d data
        % if there is extreme drought, weighted by duration
        tmp2 = nan(size(indXClimate));
        tmp2(indXClimate) = squeeze(env.drought_duration(:,yy,3)); % but these give you exact area
        tmp3 = reshape(tmp2,[720,360]);

        tmp4 = tmp3./12; % this is the two D map that used to times NBP
        
        drought_area_2d(:,:,yy) = tmp4;
    end
    
    % get the NBP of every drought affected area (FLUXCOM)
    
    % check direction
    imagesc(drought_area_2d(:,:,1));
    imagesc(FC.NEE2(:,:,1,1));
    
    dr_nbp = nan(34,3);
    
    for i = 1:3 % for every model
        
        for j = 1:34 % for every year  1980 - 2013
            
            tmp = squeeze(FC.NEE2(:,:,j,i))./(10^15);
            
            
            fl_data = flip(tmp,2);
            dr_data = squeeze(drought_area_2d(:,:,21+j)); % drought affected area
            
            tmp2 = reshape(fl_data.*dr_data,[],1);
            
            dr_nbp(j,i) = nansum(tmp2);
                         
        end
    end
    
    
    
    % get the NBP of every drought affected area (TRENDY v6)
    
    % check direction
    imagesc(drought_area_2d(:,:,1))
    imagesc(Trendy_v6.NBP2(:,:,1,1))
    
    dr_tr_nbp = nan(58,10);
    
    for i = 1:10 % for every model
        
        for j = 1:58 % for every year  1959 - 2016
            
            tmp = squeeze(Trendy_v6.NBP2(:,:,58+j,i)).*areaGrid/(10^15);%  trendy v6 1901-2016
            
            
            fl_data = rot90(tmp,3);
            dr_tr_data = squeeze(drought_area_2d(:,:,j)); % drought affected area
            
            tmp2 = reshape(fl_data.*dr_data,[],1);
            
            dr_tr_nbp(j,i) = nansum(tmp2);
                         
        end
    end
    
    
    
    
    
    
    %%%% put all std data together, bar figure
    
    % rainforest, semi-arid, T.AM, T.AF, T.AS, Pan T. CGR
    
    bar_data = [];
    bar_data_std = [];
    
    bar_data(:,1) = vertcat(squeeze(std_drought_induced(:,1,1)),std(aData.CO2DT(1:20)),nanmean(nanstd(dr_tr_nbp(1:20,:))),NaN,NaN);
    bar_data_std(:,1) = vertcat(squeeze(std_drought_induced(:,1,2)),NaN,nanstd(nanstd(dr_tr_nbp(1:20,:))),NaN,NaN);
    
    bar_data(:,2) = vertcat(squeeze(std_drought_induced(:,21,1)),std(aData.CO2DT(21:40)),nanmean(nanstd(dr_tr_nbp(21:40,:))),nanmean(nanstd(dr_nbp(1:20,:))),NaN);
    bar_data_std(:,2) = vertcat(squeeze(std_drought_induced(:,21,2)),NaN,nanstd(nanstd(dr_tr_nbp(21:40,:))),nanstd(nanstd(dr_nbp(1:20,:))),NaN);
    
    bar_data(:,3) = vertcat(squeeze(std_drought_induced(:,38,1)),std(aData.CO2DT(39:58)),nanmean(nanstd(dr_tr_nbp(39:58,:))),nanmean(nanstd(dr_nbp(21:end,:))),std(fire_em));
    bar_data_std(:,3) = vertcat(squeeze(std_drought_induced(:,38,2)),NaN,nanstd(nanstd(dr_tr_nbp(39:58,:))),nanstd(nanstd(dr_nbp(21:end,:))),NaN);
    
    
    
    fs1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 20], ...
    'OuterPosition', [2, 2, 20, 20]);

    % sub_panel figure
    SpacingVert = 0.05;
    SpacingHoriz = 0.06;
    MR = 0.01;
    ML = 0.15;
    MarginTop = 0.03;
    MarginBottom = 0.15;
    
    for tt = 1:3
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        subaxis(3,2,tt*2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

        
        contri_all = bar_data(:,tt);
        contri_all_std = bar_data_std(:,tt);
        
        hold on;
        
        tmp = contri_all_std;
        tmp(isnan(tmp)) = 0;
        tmp2 = contri_all + tmp;
        
        text(0.6:9.6,tmp2 + 0.12,num2str(round(contri_all,2)));
        
        
        b = bar(contri_all,'FaceColor',[166 206 227]/255,'EdgeColor','none');

        b.FaceColor = 'flat';
        for i = 1:6 
            b.CData(i,:) = [166 206 227]/255;
        end
        
        for i = 7:7 
            b.CData(i,:) = [150 150 150]/255;
        end
        
        for i = 8:9 
            b.CData(i,:) = [206 227 166]/255;
        end
        
        for i = 10:10 
            b.CData(i,:) = [227 166 206]/255;
        end

        hold on;

        er = errorbar(contri_all,contri_all_std);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.CapSize = 10;

        set(gca, 'box', 'off');

        ylim([0,2]);

        factor_names = {'T.forests','T.semi-arid','T.America','T.Africa','T.Asia', 'Pan Tropic', 'CGR','DGVMs','FLUXCOM','Fire'};

        ylabel({'Tropical extreme droughts'; 'induced STD_{NEE} (PgC yr^{-1})'},'FontSize',10);

        if tt == 3
            h = gca;
            %h.XTickLabel = Mdl.PredictorNames;
            h.XTick = 1:length(factor_names);
            h.XTickLabel = factor_names;
            h.XTickLabelRotation = 45;
            h.TickLabelInterpreter = 'none';
            %xlabel('Variables','FontSize',10);
        else
            h = gca;
            h.XTick = 1:length(factor_names);
            h.XTickLabel = {'','','','','','','','','',''};
        end

        
        switch tt
            case 1
                text(0.2,1.8,'(b) 1960 - 1979','fontsize',10);
            case 2
                 text(0.2,1.8,'(d) 1980 - 1999','fontsize',10); 
            case 3
                text(0.2,1.8,'(f) 1997 - 2016','fontsize',10);
        end
        
        
    end
    
    
    
% the change of biomes, rainforest and semi-arid    
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(3,2,1,'SpacingVert',SpacingVert +  0.05,'SpacingHoriz',SpacingHoriz + 0.1,'MR',MR,'ML',ML-0.05,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

color_choice = [90,180,172;216,179,101]./255;

for rr = 1:2
    
    xdata = 1969:2006;
    ydata = squeeze(area_drought(rr,:,1)).*10^2; % to 100%
    ydata_sd = squeeze(area_drought(rr,:,2)).*10^2; % 
    
    hold on;
    faceColor = color_choice(rr,:);
    pxx1(rr)=plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %title_region2{rr} = [title_region{rr},' r = ',num2str(r_IAV(1,rr))];
end

xlim([1969 2010]);
ylim([0 5.5]);
set(gca,'box','off');
text(1970,5,'(a)');

ylabel({'Extreme drought'; 'affected area (%)'},'FontSize',10);
xlabel('Year','FontSize',10); 

legend([pxx1(1),pxx1(2)],{'T.forests','T.semi-arid'},'Orientation','vertical','Location','southeast','Box','off');



% the changes of drought-affected area...three continents
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(3,2,3,'SpacingVert',SpacingVert +  0.05,'SpacingHoriz',SpacingHoriz + 0.1,'MR',MR,'ML',ML-0.05,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

color_choice = [127,201,127; 190,174,212; 253,192,134]./255;

for rr = 3:5
    
    xdata = 1969:2006;
    ydata = squeeze(area_drought(rr,:,1)).*10^2; % to 100%
    ydata_sd = squeeze(area_drought(rr,:,2)).*10^2; % 
    
    hold on;
    faceColor = color_choice(rr-2,:);
    pxx2(rr-2)=plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %title_region2{rr} = [title_region{rr},' r = ',num2str(r_IAV(1,rr))];
end

xlim([1969 2010]);
ylim([0 5.5]);
set(gca,'box','off');
text(1970,5,'(c)');

ylabel({'Extreme drought'; 'affected area (%)'},'FontSize',10);
xlabel('Year','FontSize',10); 

legend([pxx2(1),pxx2(2),pxx2(3)],{'T.America','T.Africa','T.Asia'},'Orientation','vertical','Location','southeast','Box','off');

set(fs1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure_s1_part1','-djpeg','-r600'); 




% panel show the progress of droughts

colorscheme = [1,1,1; 1,0.5,0.5];


tropical_LC = LC05D;
tropical_LC(tropical_LC ~=2) = 1; % non EBF all semi-arid and arid system
tropical_LC([1:133,227:360],:) = 0; % focus on tropical region


% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
% subaxis(3,2,5,'SpacingVert',SpacingVert - 0.05,'SpacingHoriz',SpacingHoriz,'MR',MR + 0.1,'ML',ML-0.05,'MarginTop',MarginTop,'MarginBottom',MarginBottom - 0.1); 
fs2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 20], ...
    'OuterPosition', [2, 2, 15, 10]);

hold on;

for i = 1:3
    
    switch i
        case 1 %1960-1979
            tmp = env.drought_duration(:,2:21,3); % choose 10% cutoff
            time_title = '1960-1979';
            %colorscheme = [1,1,1; 1,0,0];
        case 2 %1980-1999
            tmp = env.drought_duration(:,22:41,3);
            time_title = '1980-1999';
            %colorscheme = [1,1,1; 0,1,0];
        case 3 %1997-2016
            tmp = env.drought_duration(:,39:58,3);
            time_title = '1997-2016';
            %colorscheme = [1,1,1; 0,0,1];
    end
    
    tmp2 = nanmean(squeeze(tmp),2);
    
    tmp3 = nan(360*720,1);
    
    tmp3(indXClimate) = tmp2;
    
    
    duration_dis = reshape(tmp3,720,360); 
    
    duration_dis = duration_dis./12; % weighted by time 
    
    % change desert to NaN;
    tmp_LC = rot90(LC05D,3);
    duration_dis(tmp_LC == 16) = NaN;
    
    
    duration_dis(duration_dis<0.12) = NaN; %only consider those has 10% extreme droughts
    
     mapdata=flip(rot90(duration_dis,1));
     mapdata([1:133,227:360],:) = NaN;
     
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    subaxis(3,1,i,'SpacingVert',SpacingVert-0.03,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML-0.1,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    %mapdata= smooth2a(mapdata,1,1); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;
    
    mapdata_T = mapdata(120:240,:);
    lat_T = lat(120:240);

    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions/m_map');
    m_proj('robinson','lon',[-180 180],'lat',[-23 23]); % map projections, range of lon and lat of your data
    
    %m_proj('hammer-aitoff');
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat_T,mapdata_T); % draw your map here
    %m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colormap(gca,colorscheme);

%      caxis([0 0.25]); 
% % 
%     h=colorbar;%
%     set(h,'YTick',0:0.05:0.25);
%     set(h,'TickLabels',{'0','5','10','15','20','25'});
%     
%     if i == 2
%         ylabel(h,{'Probability of'; 'extreme droughts (%)'},'FontSize',9);
%     end
%     
%     cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');cbfreeze(h);

    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions/m_map');
    

    m_grid('box','on','linewi',0.5,'linest','none', 'tickdir','in','xticklabels',[],'yticklabels',[],'fontsize',7);
    
    %m_grid('xaxis','middle');


%     m_text(-170,19,char(96+i),'fontsize',10);
    m_text(-170,10,time_title,'fontsize',10);
    hold off;
    
    
    
end

set(fs2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/Figure_s1_part2','-djpeg','-r600');  
    

    
 
    

    













%%%%%%%%%%%%%% Figure S1: test the sig of multi-regression, using MAT and MAPlag %%%%%%%%%%%%%

%%% several regression:
% apparent T
% apparent W
% T and W
% T and W and T*W

% for each, a figure of overall P, P for T and W if there is.
% for each, a figure of R2, overall adjust R2.
% for each, a figure of T
% for each, a figure of W
% for each, a figure of T*W

stats_sen = []; % statistics for sensitivity analysis
mdl = [];

for i = 1:4 % regressions tested
    
    yearGCB = 1960:2016;
    
    switch i
        case 1
            ydata=aData.CO2DT(end-56:end); 
            xdata=[aData.TTDT(end-56:end)];
        case 2
            ydata=aData.CO2DT(end-56:end); 
            xdata=aData.TPPlagDT(end-56:end)./100;
        case 3
            ydata=aData.CO2DT(end-56:end); 
            xdata=[aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100];
        case 4
            ydata=aData.CO2DT(end-56:end); 
            xdata=[aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100 normalize(aData.TTDT(end-56:end)).*normalize(aData.TPPlagDT(end-56:end)./100) ];

            %xdata=[ones(size(ydata)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end).*aData.TPPlagDT(end-56:end)./100 ];
    end
    
    [vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');
    
    ydata(vol_ind) = NaN;
    xdata(vol_ind,:) = NaN;
    
    lag = 20;
    for j = 1:length(ydata)-lag
        
        mdl = fitlm(xdata(j:j+lag,:),ydata(j:j+lag));
        
        % overall R2 (adjust)
        stats_sen(i,j,1) = mdl.Rsquared.Ordinary;
       
        % overall p value
        stats_sen(i,j,2) = coefTest(mdl);
        
        % p value for each predictors
        
        for vv = 1:size(mdl.Coefficients.pValue)-1
            stats_sen(i,j,vv+2) = mdl.Coefficients.pValue(vv+1);
        end
        
        stats_sen(i,j,6) = mdl.ModelCriterion.AIC;
        
    end


end

stats_sen(2,:,4) = stats_sen(2,:,3); % for W sens calculation, only get gamma_W.
stats_sen(2,:,3) = 0;


fs1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 20], ...
    'OuterPosition', [2, 2, 25, 20]);

% sub_panel figure
SpacingVert = 0.08;
SpacingHoriz = 0.06;
MR = 0.03;
ML = 0.06;
MarginTop = 0.03;
MarginBottom = 0.1;

indicator_types = {'R^2','overall p','p for T','p for W','p for TxW','AIC'};

for ii = 1:4
    
    for jj = 1:6
        
        nn = (ii-1)*6 + jj;
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        subaxis(4,6,nn,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

        plot(1970:2006, squeeze(stats_sen(ii,:,jj)));
        
        ylabel(indicator_types{jj},'FontSize',9);
        
%         if jj < 2
%             ylim([0 0.8]);
%         elseif jj == 6
%             ylim([0 50]);
%         else
%             ylim([0 0.1]);
%             %break_axis('position', [0.5 0.95]);
%             
%         end
        
        xlim([1968 2010]);
        set(gca,'xtick', 1970:20:2010);
        
    end
end


    
set(fs1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS1','-djpeg','-r600'); 

%%%%%%%%%%%% End: Figure S1, statistics of univariate and multivariate regression %%%%%%%%%%%%%%%





%%%%%%%%%%%%%% Figure S2: test the sig of multi-regression, using MAT and TWS %%%%%%%%%%%%%

%%% several regression:
% apparent T
% apparent W
% T and W
% T and W and T*W

% for each, a figure of overall P, P for T and W if there is.
% for each, a figure of R2, overall adjust R2.
% for each, a figure of T
% for each, a figure of W
% for each, a figure of T*W

stats_sen = []; % statistics for sensitivity analysis
mdl = [];

for i = 1:4 % regressions tested
    
    yearGCB = 1979:2016;
    
    switch i
        case 1
            ydata=aData.CO2DT(end-37:end); 
            xdata=aData.TTDT(end-37:end);
        case 2
            ydata=aData.CO2DT(end-37:end); 
            xdata=TWS.dt_tGRACE;
        case 3
            ydata=aData.CO2DT(end-37:end); 
            xdata=[aData.TTDT(end-37:end) TWS.dt_tGRACE];
        case 4
            ydata=aData.CO2DT(end-37:end); 
            xdata=[aData.TTDT(end-37:end) TWS.dt_tGRACE normalize(aData.TTDT(end-37:end)).*normalize(TWS.dt_tGRACE)];

            %xdata=[ones(size(ydata)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end).*aData.TPPlagDT(end-56:end)./100 ];
    end
    
    [vol_ind,Locb] = ismember(yearGCB',vol_years,'rows');
    
    ydata(vol_ind) = NaN;
    xdata(vol_ind,:) = NaN;
    
    lag = 20;
    for j = 1:length(ydata)-lag
        
        mdl = fitlm(xdata(j:j+lag,:),ydata(j:j+lag));
        
        % overall R2 (adjust)
        stats_sen(i,j,1) = mdl.Rsquared.Ordinary;
       
        % overall p value
        stats_sen(i,j,2) = coefTest(mdl);
        
        % p value for each predictors
        
        for vv = 1:size(mdl.Coefficients.pValue)-1
            stats_sen(i,j,vv+2) = mdl.Coefficients.pValue(vv+1);
        end
        
        stats_sen(i,j,6) = mdl.ModelCriterion.AIC;
        
    end


end

stats_sen(2,:,4) = stats_sen(2,:,3); % for W sens calculation, only get gamma_W.
stats_sen(2,:,3) = 0;


fs2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 20], ...
    'OuterPosition', [2, 2, 25, 20]);

% sub_panel figure
SpacingVert = 0.08;
SpacingHoriz = 0.06;
MR = 0.03;
ML = 0.06;
MarginTop = 0.03;
MarginBottom = 0.1;

indicator_types = {'R^2','overall p','p for T','p for W','p for TxW','AIC'};

for ii = 1:4
    
    for jj = 1:6
        
        nn = (ii-1)*6 + jj;
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        subaxis(4,6,nn,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

        plot(1989:2006, squeeze(stats_sen(ii,:,jj)));
        
        ylabel(indicator_types{jj},'FontSize',9);
        
%         if jj < 2
%             ylim([0 0.8]);
%         elseif jj == 6
%             ylim([0 50]);
%         else
%             ylim([0 0.1]);
%             %break_axis('position', [0.5 0.95]);
%             
%         end
        
        xlim([1988 2010]);
        set(gca,'xtick', 1988:10:2010);
        
    end
end


    
set(fs2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS2','-djpeg','-r600'); 

%%%%%%%%%%%% End: Figure S2, statistics of univariate and multivariate regression, using MAT and TWS %%%%%%%%%%%%%%%




%%%%%%%%%%%% Figure S3: percentage of changes in area for different regions %%%%%%%%%%%%%%
%%% EBF, Arid, Tropical America, Tropical Africa and Tropical Asia.

LC_Tropical = nan(size(LC1D));
LC_Tropical(LC1D == 2) = 1; % EBF
LC_Tropical(LC1D == 7 | LC1D == 9) = 2; % SH and SAV, (arid)
LC_Tropical(LC1D == 12 | LC1D == 14) = 3; % CRO

area_drought = [];
r_IAV = [];

for rr = 1:5 % for different regions
    
    switch rr
        case 1 % EBF
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 1;
        case 2 % semi-arid system
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical == 2; %(9 or 7)
        case 3 % tropical America
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<-45 & env.temp(:,1)>-120 & LC_Tropical > 0;
        case 4 % tropical Africa
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<60 & env.temp(:,1)>-45 & LC_Tropical > 0;
        case 5 % tropical Asia
            indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & env.temp(:,1)<180 & env.temp(:,1)>60 & LC_Tropical > 0;
    end
    
    % get all the indicators under 10%
    aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,3)),1); 
    aData.drought_area = nansum(squeeze(env.drought_area(indX,:,3)),1);
    aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,3)),1);
    
    total_area = nansum(squeeze(area1D(indX,:)),1);

    
    tmp = nansum(squeeze(env.drought_duration(indX,:,3)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
    aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
    aData.drought_area2 = tmp./total_area/12; % the drought affected area per month

    
    M = movmean(aData.drought_area2,20,'omitnan');
    M2 = movstd(aData.drought_area2,20,'omitnan');
    
    area_drought(rr,:,1) = M(11:48);
    area_drought(rr,:,2) = M2(11:48);
    
    for mm = 1:size(IAV_fluxes,2)

        ydata = IAV_fluxes(:,mm);

        [r1,r2] = corrcoef(M(11:48),ydata,'rows','pairwise');
        sig=round(r2(2)*1000)./1000;    
        rsq=round(r1(2)^2*100)./100; 
        
        r_IAV(mm,rr) = rsq;
    end
    
    
end

fs3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 15], ...
    'OuterPosition', [2, 2, 15, 15]);

% sub_panel figure
SpacingVert = 0.08;
SpacingHoriz = 0.06;
MR = 0.03;
ML = 0.06;
MarginTop = 0.03;
MarginBottom = 0.1;

title_region = {'tropical forests','arid and semi-arid','tropical America','tropical Africa','tropical Asia'};

color_choice=[0.8,0,0;0,0.8,0;...
    0,0,0;0,0,0.8;0,0.7,0.8];



for rr = 1:5
    
    xdata = 1969:2006;
    ydata = squeeze(area_drought(rr,:,1)).*10^2; % to 100%
    ydata_sd = squeeze(area_drought(rr,:,2)).*10^2; % 
    
    hold on;
    faceColor = color_choice(rr,:);
    pxx1(rr)=plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);

    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdata,ydata,ydata_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    title_region2{rr} = [title_region{rr},' r = ',num2str(r_IAV(1,rr))];
end

ylabel('drought affected area (%)','FontSize',10);
xlabel('year','FontSize',10); 


legend([pxx1],title_region2,'Orientation','vertical','Position',[0.7    0.85    0.1    0.1],'Box','off');
 

set(fs3,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS3','-djpeg','-r600'); 



% %%%%%%%%%%%% End: Figure S3: percentage of changes in area for different regions %%%%%%%%%%%%%%



% %%%%%%%%%%%% Figure S4: examplary relationship between drought area/duration and IAV_CGR %%%%%%%%%%%%%%
% %%% 
% 
% fs4=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 15], ...
%     'OuterPosition', [2, 2, 20, 15]);
% 
% SpacingVert = 0.1;
% SpacingHoriz = 0.1;
% MR = 0.03;
% ML = 0.1;
% MarginTop = 0.03;
% MarginBottom = 0.1;
% 
% title_drought = {'1%','5%','10%','20%','50%'};
% 
% 
% color_choice=[0.8,0,0;0,0,0;0,0,0.8;0,0.8,0];
% 
% 
% %%%%% remove volcanoe years %%%%%%%
% % env.drought_duration(:,33:35,:) = NaN;
% % env.drought_area(:,33:35,:) = NaN;
% % env.drought_intensity(:,33:35,:) = NaN;
% 
% for j = 1:size(cutoff_v,2)
%     
%     indX=env.temp(:,2)<23 & env.temp(:,2)>-23 & LC_Tropical > 0; % & LC1D == 9; % index of the tropical land
% 
%     
%     total_area = nansum(squeeze(area1D(indX,:)),1);
%     
%     aData.drought_duration = nanmean(squeeze(env.drought_duration(indX,:,j)),1);
%     aData.drought_area = nansum(squeeze(env.drought_area(indX,:,j)),1)./total_area;
%     aData.drought_intensity = nanmean(squeeze(env.drought_intensity(indX,:,j)),1);
%     
%     tmp = nansum(squeeze(env.drought_duration(indX,:,j)).*squeeze(area1D(indX,:)),1); % drought affected area * drought duration
%     aData.drought_duration2 = tmp./total_area; % the average drough duration (area-weighted);
%     aData.drought_area2 = tmp./total_area/12; % the drought affected area per month
% 
%     
%     
% %     tmp2= env.drought_area;%.*env.drought_duration./nansum(env.drought_duration,1);
% %     tmp2(env.drought_duration<2 & env.drought_duration > 0)= NaN;    
% %     aData.drought_area3 = nansum(squeeze(tmp2(indX,:,j)),1)./total_area; % only have at least 1 month of drought
% %     
% 
%     % remove vol years
%     aData.drought_area2(1,33:35) = NaN;
%     M = movmean(aData.drought_area2,20,'omitnan');
%     ydata = IAV_fluxes(:,1);
%     xdata = M(11:48)'.*100; % to 100%
%     
% 
%     [r1,r2] = corrcoef(xdata,ydata,'rows','pairwise');
%     sig=round(r2(2)*1000)./1000;    
%     rsq=round(r1(2)^2*100)./100; 
%         
%    
%     cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
%     subaxis(2,3,j,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 
% 
%     
%     faceColor = color_choice(2,:);
%     %ppx3 = plot(xdata,ydata,'.','LineWidth',2.5,'color',faceColor);
%     
%     %scatter(xdata, ydata, 10,faceColor,'filled');
%     xdata_all=[ones(size(xdata)) xdata];
%     faceColor = color_choice(2,:);
%     ssx2 = scatter(xdata, ydata, 10,faceColor,'filled');
% 
%     [b,bint,~,~,~] = regress(ydata,xdata_all,0.05);
% 
%     mdl = fitlm(xdata,ydata);
%     [ypred,yci] = predict(mdl);
%     ydata_low=yci(:,1);
%     ydata_up=yci(:,2);
% 
%     pxx=plot(xdata, xdata.*b(2)+b(1),'-','color','k');
% 
%     hold on;
% 
%     [xdata_s,I] = sort(xdata);
%     h=fill([xdata_s; flip(xdata_s)], [ydata_low(I); flip(ydata_up(I))], 'k');
%     set(h,'facealpha',.2);
%     set(h,'LineStyle','none');
% 
% 
%     text(nanmin(xdata)+0.01*(nanmax(xdata)-nanmin(xdata)),1.35,title_drought{j});
%     text(nanmin(xdata)+0.2*(nanmax(xdata)-nanmin(xdata)),1.35,strcat('r2 =', num2str(rsq)));
% 
%     ylabel('STD_{CGR}','FontSize',9);
%     xlabel('drought affected area (%)','FontSize',9);
% 
% end
% 
% set(fs4,'PaperPositionMode','auto');
% print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureS4','-djpeg','-r600'); 
% 
% 
% %%%%%%%%%%%% End: Figure S4: percentage of changes in area for different regions %%%%%%%%%%%%%%






%%%%%%%%%%%%% Extended Data Fig. 2, STD of many fluxes %%%%%%%%%%%%%%%%%%%%%%%%%

fE2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 12], ...
    'OuterPosition', [2, 2, 15, 12]);

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,1,1,'SpacingVert',0,'SpacingHoriz',0,'MR',0.15,'ML',0.15,'MarginTop',0.05,'MarginBottom',0.2); 


color_choice=[0,0,0; 141,160,203; 252,141,98;228,26,28 ;102,194,165]./255;

annual_fluxes2 = [];

tmp1 = GCPdata.data.Land0x2DUseChangeEmissions(:,4:5); % LUC emissions from two book keeping model 
tmp2 = GCPdata.data.OceanSink(:,4:11); % refill fc data to 1959-2016
tmp3 = vertcat(nan(23,1),aData.Tlai3g); % 1982-2016;
tmp4 = vertcat(nan(38,1),fire_em(1:20)'); % 1997-2020 (update to 2016)

annual_fluxes2 = horzcat(GCPdata.growthRateGCP, GCPdata.oceanUptake, GCPdata.lucEmissions,tmp1,tmp2,tmp4,tmp3);

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
    
    
%     if i == 14 % for fire emission as it is short use lag of 10 years
%         
%         lag2 = 10;
%         
%         for j = 1:length(ydata)-lag
%             
%             %if sum(isnan(ydata(j:j+lag2)))>2 %only if there are more than 3 records
%                 tmp(j) = nanstd(ydata(j:j+lag2)); 
%             %end
% 
%         end
%     
%     end
        
      
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
ydata = nanmean(IAV_fluxes2(:,6:13),2);
ydata_sd = nanstd(IAV_fluxes2(:,6:13),1,2);

faceColor = color_choice(2,:);
ppx2 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(8),{'-'},0);
set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
set(H.mainLine,'color','none');
set(H.edge,'color','none');


% line 3: IAV of land use emission
% xdata = years;
% ydata = nanmean(nIAV_fluxes2(:,14:16),2);
% ydata_sd = nanstd(nIAV_fluxes2(:,14:16),1,2);
% 
% faceColor = color_choice(5,:);
% ppx3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
% cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
% H=shadedErrorBar(xdata,ydata,ydata_sd./sqrt(3),{'-'},0);
% set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
% set(H.mainLine,'color','none');
% set(H.edge,'color','none');

xdata = years;
ydata = IAV_fluxes2(:,3);
    
faceColor = color_choice(3,:);
ppx3 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);


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


%%%%%%%%% add another axis
yyaxis right
set(gca,'YColor','k');

% line 5: IAV of LAI
xdata = years;
ydata = IAV_fluxes2(:,15);
    
faceColor = color_choice(5,:);
ppx5 = plot(xdata,ydata,'-','LineWidth',2.5,'color',faceColor);
ylim([0 0.1]);

ylabel('STD_{LAI}','FontSize',12);





legend([ppx1,ppx2,ppx3,ppx4,ppx5],{'CGR','Ocean','LUC','Fire','LAI'},'FontSize',11,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);


set(fE2,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureE2','-djpeg','-r600');     
     



%%%%%%%%%%%%% End: Extended Data Fig. 2 %%%%%%%%%%%%%%%%%%%%%%%%% 



%%%%%%%%%%%%% Extended Data Fig. 1 %%%%%%%%%%%%%%%%%%%%%%%%%



% panel (a)apparent temporal sensitivity
fE1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 10], ...
    'OuterPosition', [2, 2, 25, 10]);

% panel (a) 
% BEST, CRUTEM, GISS, UDEL

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 


color_choice=lines(4);

lag_choice=[15,20,25];

yearGCB = 1959:2016;

for ttt=1:4
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(2);
        

    ydata_sen=aData.CO2DT; 
        
    switch ttt
        case 1
            xdata_sen=[ones(size(ydata_sen)) apData.BESTDT_T'];
        case 2
            xdata_sen=[ones(size(ydata_sen)) apData.CRUTEMDT_T'];
        case 3
            xdata_sen=[ones(size(ydata_sen)) apData.GISSDT_T'];
        case 4
            xdata_sen=[ones(size(ydata_sen)) apData.UDELDT_T'];
    end
    
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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(i,1)=round(tmp2(2)*1000)./1000;

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    %sort x and y, in order to put in uncertainty
    [xdatalag,sortInd]=sort(xdatalag);
    ydatalag=ydatalag(sortInd);
    ydatalag_sd=ydatalag_sd(sortInd);

    %plot the change of y to x
    %faceColor=[0.8,0.8,0.8];
    faceColor=color_choice(ttt,:);

    pxx1(ttt)=plot(xdatalag,ydatalag,'-','LineWidth',2.5,'color',faceColor);
    hold on;

    %%plot the uncertainty of y
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.2,'EdgeAlpha',0.2)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


 ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1} K^{-1})'),'FontSize',10);
 xlabel('Year','FontSize',12);
 set(gca,'box','off');

 %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
 xlim([1965.5 2010.5]);
 ylim([0 8]);

 set(gca,'xtick', 1970:10:2010);

 text(1967,7.8,'(a)','Fontsize',12);
     
 legend([pxx1(1),pxx1(2),pxx1(3),pxx1(4)],{'BEST','CRUTEM','GISS','UDEL'},'Orientation','vertical','Position',[0.08    0.7    0.1    0.1],'Box','off');


% panel (b)
% GPCC, PRECL, UDEL

cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 

for ttt=1:3
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(2);
        

    ydata_sen=aData.CO2DT; 
        
    switch ttt
        case 1
            xdata_sen=[ones(size(ydata_sen)) apData.GPCClagDT_P'./100];
        case 2
            xdata_sen=[ones(size(ydata_sen)) apData.PRECLlagDT_P'./100];
        case 3
            xdata_sen=[ones(size(ydata_sen)) apData.UDELlagDT_P'./100];
    end
    
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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(i,1)=round(tmp2(2)*1000)./1000;

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    %sort x and y, in order to put in uncertainty
    [xdatalag,sortInd]=sort(xdatalag);
    ydatalag=ydatalag(sortInd);
    ydatalag_sd=ydatalag_sd(sortInd);

    %plot the change of y to x
    %faceColor=[0.8,0.8,0.8];
    faceColor=color_choice(ttt,:);

    pxx2(ttt)=plot(xdatalag,ydatalag,'-','LineWidth',2.5,'color',faceColor);
    hold on;

    %%plot the uncertainty of y
        cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.2,'EdgeAlpha',0.2)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


ylabel(strcat(char(947),'_{CGR}','^W (PgC yr^{-1}100 mm^{-1})'),'FontSize',10);
xlabel('Year','FontSize',12);
set(gca,'box','off');

%xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
xlim([1965.5 2010.5]);
ylim([-2.5 0]);

set(gca,'xtick', 1970:10:2010);

set(gca, 'YDir','reverse');


text(1967,-2.4,'(b)','Fontsize',12);
legend([pxx2(1),pxx2(2),pxx2(3)],{'GPCC','PRECL','UDEL'},'Orientation','vertical','Position',[0.56    0.22    0.1    0.1],'Box','off');

hold off;

set(fE1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureE1','-djpeg','-r600');  

%%%%%%%%%%%%% End: Extended Data Fig. 1 %%%%%%%%%%%%%%%%%%%%%%%%%




%%%% Extended Data Fig. S3, gamma_T and W in multivariate regression %%%%%%%


fe3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 20], ...
    'OuterPosition', [2, 2, 25, 20]);


oceanDT = detrend(GCPdata.oceanUptake);

% plot out T_sensitivity

color_choice=[1,0.4,0.4;0,0,0;0.4,0.4,1];
lag_choice=[15,20,25];

%consider the effect of volcaneo
%vol_years=[1963; 1964; 1982; 1985; 1991; 1992; 1993]; %according to P.Cox
vol_years=[1991; 1992; 1993];

% sub_panel figure
SpacingVert = 0.08;
SpacingHoriz = 0.11;
MR = 0.08;
ML = 0.08;
MarginTop = 0.03;
MarginBottom = 0.1;


cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


for ttt=1:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 

            %xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)];
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)];


        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
 
            %xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end)];
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end)];


    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);
              

%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'Rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);
    
    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    
    if ttt == 1
        gamma_T = ydatalag;
    end
end

    
 ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}?C^{-1})'),'FontSize',12);
 xlabel('Year','FontSize',12);
 set(gca,'box','off');

 %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
 xlim([1968 2010]);
 ylim([-1 6]);

 set(gca,'xtick', 1970:10:2010);

 text(1969,5.9,'(a)','Fontsize',12); 
     

 legend([pxx1(2),pxx1(1)],...
     {strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS'),...
     strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}')},...
     'Orientation','vertical','Position',[0.1    0.59    0.4    0.1],'Box','off');
 
%  
 
 

 
 
 %%%% plot W_sensitivity
 
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,2,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

 for ttt=2:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) ]; %for pp-based sensitivity, times 100
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end)];
            
        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end)];
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end)];


    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;
        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    faceColor=color_choice(ttt,:);
    
    ydatalag(1,1) = NaN;

    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');


    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
 end
 ylabel(strcat(char(947),'_{CGR}','^{TWS} (PgC yr^{-1}TWS^{-1})'),'FontSize',12);
 %xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-1 0.5]); 


 
 %%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 for ttt=1:1
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) ]; %for pp-based sensitivity, times 100
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) normalize(aData.TPPlagDT(end-56:end)./100).*normalize(aData.TTDT(end-56:end)) ];
            
            
        case 2 % use MAT and TWS
            yearGCB = 1979:2016;
            ydata_sen=aData.CO2DT(end-37:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE aData.TTDT(end-37:end)];
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) normalize(TWS.dt_tGRACE(end-36:end)).* normalize(aData.TTDT(end-36:end))];

    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);
% 
%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'Rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;
        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable


        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);

    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %gamma_W = vertcat(nan(1,1),ydatalag);
    gamma_W = ydatalag;
end

    
 ylabel(strcat(char(947),'_{CGR}','^W (PgC yr^{-1}100 mm^{-1})'),'FontSize',12);
 xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-2 0.5]);

 set(gca,'xtick', 1970:10:2010);
 text(1969,-1.9,'(b)','Fontsize',12); 
 
 
 
 
 cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,3,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


for ttt=1:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 

            %xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end)];
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-56:end) aData.TPPlagDT(end-56:end) normalize(aData.TTDT(end-56:end)).*normalize(aData.TPPlagDT(end-56:end))];


        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
 
            %xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end)];
            xdata_sen=[ones(size(ydata_sen)) aData.TTDT(end-36:end) TWS.dt_tGRACE(end-36:end) normalize(aData.TTDT(end-36:end)).*normalize(TWS.dt_tGRACE(end-36:end))];


    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);
              

%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'Rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;

        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);
    
    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    
    if ttt == 1
        gamma_T = ydatalag;
    end
end

    
 ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}?C^{-1})'),'FontSize',12);
 xlabel('Year','FontSize',12);
 set(gca,'box','off');

 %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
 xlim([1968 2010]);
 ylim([-1 6]);

 set(gca,'xtick', 1970:10:2010);

 text(1969,5.9,'(c)','Fontsize',12); 
     

 legend([pxx1(2),pxx1(1)],...
     {strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^{TWS}',char(916),'TWS',' + ', char(947),'_inorm(MATxTWS)' ),...
     strcat(char(916),'CGR =', char(947),'_{CGR}','^T',char(916),'MAT + ',char(947),'_{CGR}','^W',char(916),'MAP_{lag}',' + ', char(947),'_inorm(MATxMAP_{lag})')},...
     'Orientation','vertical','Position',[0.1    0.1    0.4    0.1],'Box','off');
 
%  
 
 

 
 
 %%%% plot W_sensitivity
 
cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
subaxis(2,2,4,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

 for ttt=2:2
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) ]; %for pp-based sensitivity, times 100
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) normalize(aData.TPPlagDT(end-56:end)./100).*normalize(aData.TTDT(end-56:end)) ];
            
        case 2 % use MAT and TWS
            yearGCB = 1980:2016;
            ydata_sen=aData.CO2DT(end-36:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end)];
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) normalize(TWS.dt_tGRACE(end-36:end)).* normalize(aData.TTDT(end-36:end))];


    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);

%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;
        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable

        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);

    faceColor=color_choice(ttt,:);
    
    ydatalag(1,1) = NaN;

    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');


    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
 end
 ylabel(strcat(char(947),'_{CGR}','^{TWS} (PgC yr^{-1}TWS^{-1})'),'FontSize',12);
 %xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-1 0.5]); 


 
 %%% add another axis
 yyaxis right
 set(gca,'YColor','k');
 
 for ttt=1:1
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd sig1;
    lag=lag_choice(2);
    
    switch ttt
        case 1 % use MAT and lagged PP
            yearGCB = 1960:2016;
            ydata_sen=aData.CO2DT(end-56:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) ]; %for pp-based sensitivity, times 100
            xdata_sen=[ones(size(ydata_sen)) aData.TPPlagDT(end-56:end)./100 aData.TTDT(end-56:end) normalize(aData.TPPlagDT(end-56:end)./100).*normalize(aData.TTDT(end-56:end)) ];
            
            
        case 2 % use MAT and TWS
            yearGCB = 1979:2016;
            ydata_sen=aData.CO2DT(end-37:end); 
            
            %xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE aData.TTDT(end-37:end)];
            xdata_sen=[ones(size(ydata_sen)) TWS.dt_tGRACE(end-36:end) aData.TTDT(end-36:end) normalize(TWS.dt_tGRACE(end-36:end)).* normalize(aData.TTDT(end-36:end))];

    end

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
        
        cd('/Users/xiangzhongluo/Documents/Work_folder/project7/code');
        r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,500);
% 
%         [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag,2),ydata_sen(i:i+lag),'Rows','pairwise');
%         sig1(i,1)=round(tmp2(2)*1000)./1000;
        mdl = fitlm(xdata_sen(i:i+lag,:),ydata_sen(i:i+lag));
        sig1(i,1)= round(mdl.Coefficients.pValue(3)*1000)./1000; % not the overall p, but the p for the first variable


        yearlag(i)=yearGCB(i)+round(lag*0.5);

    end
         
    ydatalag=r1datalag(:,1);
    ydatalag_sd=r1datalag(:,2);

    xdatalag=yearlag';


    if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
        xlim_value(1)=min(xdatalag);
    end
    if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
        xlim_value(2)=max(xdatalag);
    end

    xlim_value(3)=xlim_value(2)-xlim_value(1);


    faceColor=color_choice(ttt,:);

    hold on;

    pxx1(ttt)=plot(xdatalag,ydatalag,'-o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor',faceColor);
    %plot(xdatalag(sig1 > 0.01),ydatalag(sig1 > 0.01), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');
    plot(xdatalag(sig1 > 0.05),ydatalag(sig1 > 0.05), 'o','LineWidth',2.5,'color',faceColor,'MarkerEdgeColor',faceColor,'MarkerFaceColor','w');

    
    %%% plot the uncertainty of y
    cd('/Users/xiangzhongluo/Documents/Work_folder/project6/code4Remi/functions');
    H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

    set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
    set(H.mainLine,'color','none');
    set(H.edge,'color','none');
    
    %gamma_W = vertcat(nan(1,1),ydatalag);
    gamma_W = ydatalag;
end

    
 ylabel(strcat(char(947),'_{CGR}','^W (PgC yr^{-1}100 mm^{-1})'),'FontSize',12);
 xlabel('Year','FontSize',12);
 set(gca,'box','off');
 
 set(gca, 'YDir','reverse');
 xlim([1968 2010]);
 ylim([-2 0.5]);

 set(gca,'xtick', 1970:10:2010);
 text(1969,-1.9,'(d)','Fontsize',12); 
 
set(fe3,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureE3','-djpeg','-r600');     
     
  
 

%%%% End Extended Data Fig. S3, gamma_T and W in multivariate regression %%%%%%%




%%%% Res Fig. 4 comparison with TRMM and CRU-NCEP %%%%%%%

% load in CRU-NCEP, from 1980-2016
filename = '/Users/xiangzhongluo/Documents/Work_folder/project7/Data/CRU_NCEP_V7_1980_2016_Prec.nc';
rawdata = ncread(filename,'Prec');

tmp = permute(rawdata,[3,2,1]);
CRUNCEP_P = [];

% calculate the annual values
for i = 1:37   
    CRUNCEP_P(:,:,i) = mean(tmp(:,:,(i-1)*12+1: i*12),3).*365.*24.*3600;    
end

% calculate tropical land mean value
CRUNCEP_tropicalP = [];
for i = 1:37   
    tmp =  squeeze(CRUNCEP_P(:,:,i));
    tmp(LC05D < 1) = NaN;
    CRUNCEP_tropicalP(i,1) = nanmean(reshape(tmp(134:226,:),[],1));
end


% load in CRU-JRA, from 1959-2018
filename = '/Users/xiangzhongluo/Documents/Work_folder/project7/Data/cru_jra_yy.mat'; % variable name is CRU_JRA_yy
%
load(filename);

% calculate tropical land mean value
CRUJRA_tropicalP = [];
for i = 1:58 % only get 1959-2016
    tmp =  squeeze(CRU_JRA_yy(:,:,i));
    tmp(LC05D < 1) = NaN;
    CRUJRA_tropicalP(i,1) = nanmean(reshape(tmp(134:226,:),[],1));
end


% load in TRMM, from 1998-2016
filename=strcat('/Users/xiangzhongluo/Documents/Work_folder/project6/datasets_pro/pp_trmm.mat');    
tmp = load(filename);
TRMM_P = permute(tmp.yearData,[2,3,1]);

% calculate tropical land mean value
TRMM_tropicalP = [];
for i = 1:19   
    tmp =  squeeze(TRMM_P(:,:,i));
    tmp(LC05D < 1) = NaN;
    TRMM_tropicalP(i,1) = nanmean(reshape(tmp(134:226,:),[],1));
end

%%% 

fR1=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 12], ...
    'OuterPosition', [2, 2, 15, 12]);
color_choice=[0,0,0; 141,160,203; 252,141,98;102,194,165]./255;

hold on;
%plot(1980:2016, CRUNCEP_tropicalP);
bias = nanmean(CRUJRA_tropicalP(end-36:end,1) - CRUNCEP_tropicalP);
ppx1 = plot(1998:2016, TRMM_tropicalP,'-','LineWidth',2.5,'color',color_choice(1,:));
ppx2 = plot(1959:2016, CRUJRA_tropicalP - bias,'-','LineWidth',2.5,'color',color_choice(2,:));
ppx3 = plot(1959:2016, aData.tropicalPP,'-','LineWidth',2.5,'color',color_choice(3,:));


ylim([1000 1500]);
xlim([1955 2020]);
ylabel('Tropical Precip (mm/yr)','FontSize',12);
xlabel('Year','FontSize',12);


legend([ppx1,ppx2,ppx3],{'TRMM','CRU-NCEP','CRU'},'FontSize',12,'box','off','Orientation','horizontal','Location',[0.25,0.15,0.5,0.1]);

set(fR1,'PaperPositionMode','auto');
print('/Users/xiangzhongluo/Documents/Work_folder/project7/Figures/FigureR1','-djpeg','-r600');

%%%% Res Fig. 4 comparison with TRMM and CRU-NCEP %%%%%%%








