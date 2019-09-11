



%%%Figure 1, Atmospheric CO2 growth (ACG) verus time

f1=figure('Name','Atmospheric CO2 growth (ACG) verus time','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 10], ...
    'OuterPosition', [2, 2, 15, 10]);
ydata=aData.CO2gr;
%ydata=sum(aData.CO2EMD(1:3,:),1)';
xdata=years';

%consider the effect of volcaneo
% vol_years=[1963; 1964; 1982; 1985; 1992; 1993]; %according to P.Cox
% for i=1:length(vol_years)
%     index=xdata==vol_years(i);
%     xdata(index)=-9999;
% end
% index=xdata>0;
% ydata=ydata(index);
% xdata=xdata(index);

pxx1=plot(xdata,ydata,'-','LineWidth',1,'color',[0.3,0.3,0.3]);

hold on;

%add regression line
p=polyfit(xdata,ydata,1);
[r1,r2]=corrcoef(xdata,ydata);
sig=round(r2(2)*1000)./1000;
h=refline(p(1),p(2));
h.Color = [1,0.3,0.3];
h.LineWidth=1;
h.LineStyle='--';

p(1)=round(p(1).*100)./100;
p(2)=round(p(2).*100)./100;
sig=round(r2(2)*1000)./1000; 
xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  p=',num2str(sig)); 
%text(1960,5,xfun);


%add EMD lines
hold on;
h2=plot(years,aData.CO2EMD(4,:));

xlim([1959,2016]);
ylim([0.5 6.5]);
ylabel({'CGR (PgC yr^{-1})'},'FontSize',10);
xlabel({'Year'},'FontSize',10);
legend([h,h2],{'linear trend','adaptive trend'},'Orientation','vertical','Location', 'northwest','Box','off')

set(f1,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/FigureS2','-djpeg','-r600');   


%




%create figure s2 which exhibit the correlation between EMD sensitivity and
%climate variables
f6=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 9, 18], ...
'OuterPosition', [2, 2, 9, 18]);

hold on;

%month_ft=12;
yearGCB=years;

%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

for fign=1:3
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 1 % x t-tmp
        ydata=sum(aData.TTEMD(1:end-1,:),1);
        ytitle='$\overline\Delta\sf\overline{MAT} (^\circ C)$';
         
%         case 2 % x t-par
%         ydata=aData.tropicalPAR;
%         ytitle='Trop Par (W m^{-2})';
        
        case 2 % x t-pp
        ydata=sum(aData.TPPEMD(1:end-1,:),1);
        ytitle='$\overline\Delta\sf\overline{MAP} (mm\ yr^{-1})$';
        
%         case 3 % x t-vpd
%         ydata=sum(aData.TDEMD(1:end-1,:),1);
%         ytitle='$\overline\Delta\sf\overline{MAD} (kPa)$';
        
                
        case 3 % x t-swc
        ydata=sum(aData.TSWCEMD(1:end-1,:),1);
        ytitle='$\overline\Delta\sf\overline{MASW} (\%)$';


    end
    
%     ylim_value(1)=min(ydata);
%     ylim_value(2)=max(ydata);
%     ylim_value(3)=ylim_value(2)-ylim_value(1);
%     
%     ylim([ylim_value(1)-0.01*ylim_value(3) ylim_value(2)+0.01*ylim_value(3)]);
        
   
          
    %plotting response curve 
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(3,1,fign,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.05,'ML',0.15,'MarginTop',.03,'MarginBottom',.1);   
    hold on;
    xlim_value(1)=min(ydata);
    xlim_value(2)=max(ydata);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        %ydata_sen=mean(aData.CO2EMD(1:3,:),1)'; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        %ydata_sen=aData.CO2gr-aData.CO2EMD(4,:)';
        ydata_sen=sum(aData.CO2EMD(1:end-1,:),1)';
        xdata_sen=[ones(size(ydata_sen)) sum(aData.TTEMD(1:3,:),1)' sum(aData.TPEMD(1:3,:),1)' sum(aData.TPPEMD(1:3,:),1)'];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            %r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
       
%         if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
%             xlim_value(1)=min(xdatalag);
%         end
%         if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
%             xlim_value(2)=max(xdatalag);
%         end
%         
%         xlim_value(3)=xlim_value(2)-xlim_value(1);
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);

        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);

            
            %%plot the uncertainty of y
%             cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
%             H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);
% 
%             set(H.patch,'FaceColor',faceColor,'EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2)
%             set(H.mainLine,'color','none');
%             set(H.edge,'color','none');
%         
        %refline(r1datalag(:,1),r1datalag(:,3));
        
            p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig<0.01
             ptext='p<0.01';
         elseif sig<0.05
             ptext='p<0.05';
         else
             ptext='p>0.1';
         end
         
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end
        %xtitle={xtitle;xfun};
        %end of ploting
        %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel(ytitle,'FontSize',10,'Interpreter','Latex');
     text(xlim_value(1),1.7,xfun);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([1 7]);
     
     text(xlim_value(1)-0.05*xlim_value(3),6.5,strcat('(',char(96+fign),')'));
     clear xlim_value
     
     
end

hold off;

set(f6,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/FigureS3','-djpeg','-r600');  









%%% Figure S2, apparent temporal sensitivity
f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 11, 10], ...
    'OuterPosition', [2, 2, 11, 10]);


% cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
% subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.08,'MR',0.03,'ML',0.08,'MarginTop',.03,'MarginBottom',.18); 
%color_choice=[0.6,0,0;0.3,0.3,0.3;0,0,1];
color_choice=[1,0.4,0.4;0,0,0;0.4,0.4,1];
lag_choice=[15,20,25];

for ttt=1:3
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(ttt);
    %lag=20;
    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    %ydata_sen=sum(aData.CO2EMD(1:3,:),1)';
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT];
    xlim_value=nan(3,1);

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Volumes/RemiLBNL/project7/code');
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
        %cd('/Volumes/RemiLBNL\xl\code\functions');
        cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


    
     ylabel(strcat('dCGR/','dMAT (PgC yr^{-1}°C^{-1})'),'FontSize',12);
     xlabel('Year','FontSize',12);
     set(gca,'box','off');
     
     %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     xlim([1965.5 2010.5]);
     ylim([1 8]);
     
     set(gca,'xtick', 1970:10:2010);
     
     %text(1967,7.3,'(a)','Fontsize',12);
     
 legend([pxx1(1),pxx1(2),pxx1(3)],{'15-yr window','20-yr window','25-yr window'},'Orientation','vertical','Position',[0.26    0.8    0.1    0.1],'Box','off');
set(f2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/FigureS1','-djpeg','-r600'); 
 
 
 




%%temperature sensitivity
%%% Figure 2, three periods and dynamic of a2
% temperature sensitivity

f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 10], ...
    'OuterPosition', [2, 2, 25, 10]);

 
%creat subplot 2
%subplot(1,2,2);
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 


% f3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 11, 10], ...
%     'OuterPosition', [2, 2, 11, 10]);

%color_choice=[0.6,0,0;0.3,0.3,0.3;0,0,1];
color_choice=[1,0.4,0.4;0,0,0;0.4,0.4,1];
lag_choice=[15,20,25];

for ttt=1:3
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(ttt);
    %lag=20;
    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    %ydata_sen=sum(aData.CO2EMD(1:3,:),1)';
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPDT aData.TPPDT];
    xlim_value=nan(3,1);

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Volumes/RemiLBNL/project7/code');
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
        %cd('/Volumes/RemiLBNL\xl\code\functions');
        cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',12);
     xlabel('Year','FontSize',12);
     set(gca,'box','off');
     
     %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     xlim([1965.5 2010.5]);
     ylim([1 8]);
     
     set(gca,'xtick', 1970:10:2010);
     
     text(1967,7.9,'(a)','Fontsize',12); 
     
     %text(1967,7.3,'(b)','Fontsize',12);
     
 legend([pxx1(1),pxx1(2),pxx1(3)],{'15-yr window','20-yr window','25-yr window'},'Orientation','vertical','Position',[0.09    0.8    0.1    0.1],'Box','off');
 
 
 %creat subplot 1
%subplot(1,2,1);
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 

ydata=aData.CO2DT;
%ydata=sum(aData.CO2EMD(1:3,:),1)';
xdata=aData.TTDT;

yearGCB=years;

% pxx1=scatter(xdata,ydata,60,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],...
%                       'LineWidth',1);
hold on;


% pxx1=scatter(xdata(1:20),ydata(1:20),50,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],...
%                       'LineWidth',1);
% pxx2=scatter(xdata(21:40),ydata(21:40),40,'MarkerEdgeColor',[0.4,0.4,0.4],'MarkerFaceColor','none',...
%   'LineWidth',2);
% pxx3=scatter(xdata(41:end),ydata(41:end),50,'MarkerEdgeColor',[0.4,0.4,0.4],'MarkerFaceColor',[0.4,0.4,0.4],...
%   'LineWidth',1);

for nn=1:length(xdata)
    tmp=num2str(yearGCB(nn));
    ratio=1-(yearGCB(nn)-min(yearGCB))./(max(yearGCB)-min(yearGCB));
    %text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0.7,0.7,0.7].*ratio,'FontSize',7.2);
    
    %to avoid overlap, special cases
    if yearGCB(nn)==1967 || yearGCB(nn)==1981 || yearGCB(nn)==1990  || yearGCB(nn)==1997 || yearGCB(nn)==2002 || yearGCB(nn)==2006
        %text(xdata(nn)+0.03,ydata(nn)+0.1,tmp(3:4),'Color',[0.7,0.7,0.7].*ratio,'FontSize',7.2);
            if nn <=22
                text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0.8,0.8,0.8],'FontSize',7.2);
            elseif nn<=38 && nn>22 %n<=42
                text(xdata(nn),ydata(nn),strcat('\it',tmp(3:4)),'Color',[0.5,0.5,0.5],'FontSize',7.2);
            else
                %text(xdata(nn),ydata(nn),strcat('\bf',tmp(3:4)),'Color',[0,0,0],'FontSize',10);
                text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0,0,0],'FontSize',7.2);
            end
    else
        %text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0.7,0.7,0.7].*ratio,'FontSize',7.2);
        
            if nn <=22
                text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0.8,0.8,0.8],'FontSize',7.2);
            elseif nn<=42 && nn>22
                text(xdata(nn),ydata(nn),strcat('\it',tmp(3:4)),'Color',[0.5,0.5,0.5],'FontSize',7.2);
            else
                %text(xdata(nn),ydata(nn),strcat('\bf',tmp(3:4)),'Color',[0,0,0],'FontSize',10);
                text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0,0,0],'FontSize',7.2);
            end
    end
%     if nn <=22
%         text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0.8,0.8,0.8],'FontSize',7.2);
%     elseif nn<=42 && nn>22
%         text(xdata(nn),ydata(nn),strcat('\it',tmp(3:4)),'Color',[0.5,0.5,0.5],'FontSize',7.2);
%     else
%         %text(xdata(nn),ydata(nn),strcat('\bf',tmp(3:4)),'Color',[0,0,0],'FontSize',10);
%         text(xdata(nn),ydata(nn),tmp(3:4),'Color',[0,0,0],'FontSize',7.2);
%     end

end


%pick the 1st period
ydata=aData.CO2DT(2:21);
xdata=aData.TTDT(2:21);

p=polyfit(xdata,ydata,1);
[r1,r2]=corrcoef(xdata,ydata);
p(1)=round(p(1).*100)./100;
p(2)=round(p(2).*100)./100;
sig=round(r2(2)*1000)./1000; 
% if p(2)>0 && sig<0.05
%     xfun1=strcat(char(949),'CGR_{1959-1978}=',num2str(p(1)),char(949),'MAT+',num2str(p(2)), ';  p<0.05'); 
% elseif p(2)<0 && sig>=0.05
%     xfun1=strcat(char(949),'CGR_{1959-1978}=',num2str(p(1)),char(949),'MAT',num2str(p(2)), ';  p>0.05'); 
% elseif p(2)>0 && sig>0.05
%     xfun1=strcat(char(949),'CGR_{1959-1978}=',num2str(p(1)),char(949),'MAT+',num2str(p(2)), ';  p>0.05'); 
% else
%     xfun1=strcat(char(949),'CGR_{1959-1978}=',num2str(p(1)),char(949),'MAT',num2str(p(2)), ';  p<0.05'); 
% end

xdata_sim=min(xdata):0.01:max(xdata);
ydata_sim=p(1)*xdata_sim+p(2);
%pxx2r=plot(xdata_sim,ydata_sim,'-','LineWidth',3,'color',[1,0.3,0.3]);
pxx2r=plot(xdata_sim,ydata_sim,'-','LineWidth',1.5,'color','k');


%pick the 2st period
ydata=aData.CO2DT(22:41);
xdata=aData.TTDT(22:41);

p=polyfit(xdata,ydata,1);
[r1,r2]=corrcoef(xdata,ydata);
p(1)=round(p(1).*100)./100;
p(2)=round(p(2).*100)./100;
sig=round(r2(2)*1000)./1000; 
% if p(2)>0 && sig<0.05
%     xfun2=strcat(char(949),'AGR_{1979-1998}=',num2str(p(1)),char(949),'MATT+',num2str(p(2)), ';  p<0.05'); 
% elseif p(2)<0 && sig>=0.05
%     xfun2=strcat(char(949),'AGR_{1979-1998}=',num2str(p(1)),char(949),'MATT',num2str(p(2)), ';  p>0.05'); 
% elseif p(2)>0 && sig>0.05
%     xfun2=strcat(char(949),'AGR_{1979-1998}=',num2str(p(1)),char(949),'MATT+',num2str(p(2)), ';  p>0.05'); 
% else
%     xfun2=strcat(char(949),'AGR_{1979-1998}=',num2str(p(1)),char(949),'MATT',num2str(p(2)), ';  p<0.05'); 
% end
    
xdata_sim=min(xdata):0.01:max(xdata);
ydata_sim=p(1)*xdata_sim+p(2);
%pxx3r=plot(xdata_sim,ydata_sim,'-','LineWidth',3,'color',[0.3,1,0.3]);
pxx3r=plot(xdata_sim,ydata_sim,'-','LineWidth',1.5,'color','k');

    
    
%pick the 3rd period
ydata=aData.CO2DT(end-19:end);
xdata=aData.TTDT(end-19:end);

p=polyfit(xdata,ydata,1);
[r1,r2]=corrcoef(xdata,ydata);
p(1)=round(p(1).*100)./100;
p(2)=round(p(2).*100)./100;
sig=round(r2(2)*1000)./1000; 
% if p(2)>0 && sig<0.05
%     xfun3=strcat('\epsilonAGR_{1999-2016}=',num2str(p(1)),'\epsilonMATT+',num2str(p(2)), ';  p<0.05'); 
% elseif p(2)<0 && sig>=0.05
%     xfun3=strcat('\epsilonAGR_{1999-2016}=',num2str(p(1)),'\epsilonMATT',num2str(p(2)), ';  p>0.05'); 
% elseif p(2)>0 && sig>0.05
%     xfun3=strcat('\epsilonAGR_{1999-2016}=',num2str(p(1)),'\epsilonMATT+',num2str(p(2)), ';  p>0.05'); 
% else
%     xfun3=strcat('\epsilonAGR_{1999-2016}=',num2str(p(1)),'\epsilonMATT',num2str(p(2)), ';  p<0.05'); 
% end

xdata_sim=min(xdata):0.01:max(xdata);
ydata_sim=p(1)*xdata_sim+p(2);
%pxx4r=plot(xdata_sim,ydata_sim,'-','LineWidth',3,'color',[0.3,0.3,1]);
pxx4r=plot(xdata_sim,ydata_sim,'-','LineWidth',1.5,'color','k');

    
% text(-0.45,2.8,xfun1);
% text(-0.45,2.0,xfun2);
% text(-0.45,1.2,xfun3);


text(0.44,2.4,'1980-1999');
text(0.44,1.7,'1997-2016');
text(0.38,0.49,'1960-1979');

xlim([-0.5 0.7]);
ylim([-2.5 3]);
ylabel(strcat(char(916),'CGR (PgC yr^{-1})'),'FontSize',12);
xlabel(strcat(char(916),'MAT (°C)'),'FontSize',12); 
 
text(-0.47,2.9,'(b)','Fontsize',12); 
 
set(gca,'box','off');
     
set(f2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure1','-djpeg','-r600');     
     
     

     
     
     
     
     

%%% Figure 3. plot climatic variables and its influence on sensitivity
%%% T (linear), PP, swc, VPD
f6=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 20], ...
'OuterPosition', [2, 2, 18, 20]);

hold on;

%month_ft=12;
yearGCB=years;

%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

for fign=1:3
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 1 % x t-tmp
        ydata=aData.tropicalT;
        ytitle='$\sf{MAT}\ and\ \overline{MAT} (^\circ C)$'; 
        
        ydata2=aData.TTDT;
        ytitle2='$\overline\Delta\sf\overline{MAT} (^\circ C)$';
         
        
        case 2 % x t-pp
        ydata=aData.tropicalPP;
        ytitle='$\sf{MAP}\ and\ \overline{MAP}(mm\ yr^{-1})$';
        
        ydata2=aData.TPPDT;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';
        
%         case 3 % x t-vpd
%         ydata=aData.tropicalVPD;
%         ytitle='$\sf{MAD}\ and\ \overline{MAD}(kPa)$';
%         
%         ydata2=aData.TVPDDT;
%         ytitle2='$\overline\Delta \sf\overline{MAD} (kPa)$';
        
                
        case 3 % x t-swc
        ydata=aData.tropicalSWC;
        ytitle='$\sf{MASW}\ and\ \overline{MASW} (\%)$';
        
        ydata2=aData.TSWCDT;
        ytitle2='$\overline\Delta \sf\overline{MASW} (\%)$';


    end
        
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    %%%%start ploting, long term meteo
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    %cd('/Volumes/RemiLBNL\xl\code\functions');
    subaxis(3,2,(fign-1)*2+1,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08); 
    
    pxx2=plot(yearGCB,ydata,'-','LineWidth',2,'color',[0.8 0.8 0.8]);
    hold on;
    
    ylim_value(1)=min(ydata);
    ylim_value(2)=max(ydata);
    ylim_value(3)=ylim_value(2)-ylim_value(1);
    
    
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT];
        clear yearlag ydatalag temp_ydatalag;

        for i=1:length(ydata_sen)-lag
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_ydatalag(i,1)=mean(ydata(i:i+lag));
            
        end

        ydatalag=temp_ydatalag;
        
        if min(ydatalag)<ylim_value(1) || isnan(ylim_value(1))
            ylim_value(1)=min(ydatalag);
        end
        if max(ydatalag)>ylim_value(2) || isnan(ylim_value(2))
            ylim_value(2)=max(ydatalag);
        end
        
        ylim_value(3)=ylim_value(2)-ylim_value(1);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        pxx1=plot(yearlag',ydatalag,'-','LineWidth',2,'color',faceColor);
        hold on;
            
        %end of ploting
        clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
    end   
    
     xlabel({'Year'},'FontSize',10);
     ylabel({ytitle},'FontSize',10,'Interpreter','Latex');
     
     ylim([ylim_value(1)-0.01*ylim_value(3) ylim_value(2)+0.01*ylim_value(3)]);
     xlim([1959 2016]);
     text(1961,ylim_value(2)-0.05*ylim_value(3),strcat('(',char(96+(fign-1)*2+1),')'));
     
     
    %plotting response curve 
    subaxis(3,2,(fign-1)*2+2,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08);   
    hold on;
    xlim_value(1)=min(ydata2);
    xlim_value(2)=max(ydata2);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata2(i:i+lag)); %meteo was assigned to ydata
            
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
       
%         if min(xdatalag)<xlim_value(1) || isnan(xlim_value(1))
%             xlim_value(1)=min(xdatalag);
%         end
%         if max(xdatalag)>xlim_value(2) || isnan(xlim_value(2))
%             xlim_value(2)=max(xdatalag);
%         end
%         
%         xlim_value(3)=xlim_value(2)-xlim_value(1);
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);
        
        
        
        
        
        
%          % added part for auto-correlation
%              mdl = fitlm(yearlag(1:end-1),ydatalag(1:end-1)-ydatalag(2:end)); 
%             [p1,DW] = dwtest(mdl,'exact','both');
%             disp(p1);disp(DW);
%             
%             mdl = fitlm(yearlag(1:end-1),xdatalag(1:end-1)-xdatalag(2:end));
%             [p1,DW] = dwtest(mdl,'exact','both');
%             disp(p1);disp(DW);
%         
%         
%         no_auto_y=ydatalag(1:end-1)-ydatalag(2:end);
%         no_auto_x=xdatalag(1:end-1)-xdatalag(2:end);
%         no_auto_ysd=ydatalag_sd(1:end-1);
%         
%         
%         xdatalag=[];ydatalag=[];ydatalag_sd=[];
%         [xdatalag,sortInd]=sort(no_auto_x);
%         ydatalag=no_auto_y(sortInd);
%         ydatalag_sd=no_auto_ysd(sortInd);
        
        
        
        

        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;
            
            %%plot the uncertainty of y
%             cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
%             H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);
% 
%             set(H.patch,'FaceColor',faceColor,'EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2)
%             set(H.mainLine,'color','none');
%             set(H.edge,'color','none');

        
            p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig < 0.01
             ptext='p<0.01';
         elseif sig < 0.05
             ptext='p<0.05';
         else
             ptext='p>0.1';
         end
            
  
            
        %xfun=strcat('Y_{',num2str(lag),'}=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p=',num2str(sig)); 
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end
        %xtitle={xtitle;xfun};
        %end of ploting
        %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel(ytitle2,'FontSize',10,'Interpreter','Latex');
     text(xlim_value(1),0.7,xfun);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([0 7.5]);
     
     text(xlim_value(1)-0.05*xlim_value(3),7,strcat('(',char(96+(fign-1)*2+2),')'));
     clear xlim_value
end

hold off;

set(f6,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure2','-djpeg','-r600');  







 
 
 
 
 

%Figure 8, after select the best model, simulate a2 and ACG
%y=ydatalag;
 
%get x and y values first
%prepare the matrix
lag_choice=[15,20,25];
for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
    lag=lag_choice(tt);
    %1, in the time window, calculate sensitivity

    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
    clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=b(3,1); %coefficient for TPP
        r1datalag(i,5)=b(4,1);% for TP

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(i,1)=round(tmp2(2)*1000)./1000;

        yearlag(i)=yearGCB(i)+round(lag*0.5);
        xdatalag(i,1)=mean(aData.TTDT(i:i+lag)); %get average meteo
        xdatalag(i,2)=mean(aData.TPPDT(i:i+lag));
        xdatalag(i,3)=mean(aData.TSWCDT(i:i+lag));
        xdatalag(i,4)=mean(aData.TPDT(i:i+lag));
        %xdatalag(i,5)=yearlag(i)-1968;

         %xdatalag(i,5)=mean(aData.tropicalT(i:i+lag)); %get average meteo
         xdatalag(i,5)=i;

    end

     ydatalag=r1datalag(:,1);

    %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
end
 x1=xdatalag(:,2);
 x2=xdatalag(:,5);
 x3=xdatalag(:,1);
 y=ydatalag;
 
 x=xdatalag(:,[2,5]);
 %x=horzcat(xdatalag(:,[1,i+1]),xdatalag(:,1).*xdatalag(:,i+1));
 [b,~,stats] = glmfit(x,ydatalag,'normal');

 
 %y=b1+b2*PP+b3*i;
 %set baseline in a nominal year with 24C and 1300 mm
 f7=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 20], ...
'OuterPosition', [2, 2, 20, 20]);
 hold on;
 
 %plot a2
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 subaxis(3,1,1,'SpacingVert',0.1,'SpacingHoriz',0.08,'MR',0.03,'ML',0.1,'MarginTop',.03,'MarginBottom',.1); 
 hold on;
 %observation
 px1=plot(yearlag,ydatalag,'-','LineWidth',2,'color',[0 0 0]);
 
%  %simulation
  %y=b(2)*x1+b(3)*x2+b(1) - (b(2)*x1(1)+b(3)*x2(1)+b(1));
  y=b(2)*x1+b(3)*x2+b(1);
  px2=plot(yearlag,y,'--','LineWidth',2,'color',[0.5 0.5 0.5]);
 
 %contribution from P
 %y1=b(2)*x1+b(1) - (b(2)*x1(1)+b(1)); %only from PP
 y1=b(2)*x1+b(1) + b(3)*x2(1); %only from PP
 px3=plot(yearlag,y1,'--','LineWidth',2,'color',[0.5,0.5,1]);
 
 %y2=b(3)*x2 - (b(3)*x2(1)); %adjustment put by MAT
 y2=b(3)*x2 +b(1) + (b(2)*x1(1)); %adjustment put by MAT
 px4=plot(yearlag,y2,'--','LineWidth',2,'color',[1, 0.6, 0.2]);
 y12=b(2)*x1+b(3)*x2+b(1);
 
%   y3=b(4)*x3+b(1); %only from time
%  px5=plot(yearlag,y3,'--','LineWidth',2,'color',[0.3 0.8 0.3]);
 
 h=refline(0,2);
 h.Color = 'k';


 legend([px1,px2,px3,px4],{'$\sf{obs.}$', '$\sf{pred.}$','$f\big(\overline\Delta\sf\overline{MAP}\big)$','$f\big(\sf\textit{t}\big)$'},'Orientation','horizontal','Location', 'south','Box','off','Interpreter','latex')
 
 %xlim([1959 2016]);
 xlim([1968 2008]);
 ylim([-2 8]);
 ylabel({strcat(char(947),'_{CGR}','^T'); '(PgC yr^{-1}°C^{-1})'},'FontSize',10);
 %set(gca,'XTickLabel',[]);
 xlabel('The center year of each 20-year window','FontSize',10);
 text(1969,7.2,'(a)','FontSize',12);
 
 
 %plot natural anomaly of ACG
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 subaxis(3,1,2:3,'SpacingVert',0.1,'SpacingHoriz',0.08,'MR',0.03,'ML',0.1,'MarginTop',.03,'MarginBottom',.1); 
 
%   %add volcane and fires
%  indX = years==1963 | years==1982 | years==1991;  % ?Mount Agung (in 1963), El Chichon (in 1982) and Mount Pinatubo (in 1991)
%  px_v=plot([years(indX)' years(indX)'],[-3 2.1],...
%     '-','LineWidth',3,'color',[1 0.8 0.8],'DisplayName','Vocalnoes');
%  %px_h=plot(years,detrend(GCPdata.ffEmissions-GCPdata.oceanUptake));
 
 
 
 hold on;
 
 
 epsilon_acg=nan(58,38); %58 years in total, but use 20 years window
 epsilon_acg_T=nan(58,38); 
 epsilon_acg_P=nan(58,38); 
 epsilon_acg_R=nan(58,38);
 epsilon_acg_TP=nan(58,38);
 epsilon_acg_TT=nan(58,38);
 epsilon_acg_un=nan(58,38);
 %epsilon_acg_base=nan(58,38);
 
 lag=20;
 
 for i=1:length(ydatalag)

    epsilon_acg(i:i+20,i)=(y(i))*aData.TTDT(i:i+20,1)+r1datalag(i,3)+r1datalag(i,4)*aData.TPPDT(i:i+20,1)...
        +r1datalag(i,5)*aData.TPDT(i:i+20,1);

    epsilon_acg_TP(i:i+20,i)=(b(2)*x1(i) )*(aData.TTDT(i:i+20,1));%+r1datalag(i,4)*aData.TPPDT(i:i+20,1)...
        %+r1datalag(i,5)*aData.TPDT(i:i+20,1);
    epsilon_acg_TT(i:i+20,i)=(b(3)*x2(i) - b(3)*x2(1))*(aData.TTDT(i:i+20,1));%+r1datalag(i,4)*aData.TPPDT(i:i+20,1)...
        %+r1datalag(i,5)*aData.TPDT(i:i+20,1);
    epsilon_acg_T(i:i+20,i)=(b(1) + b(3)*x2(1))*(aData.TTDT(i:i+20,1));%+r1datalag(i,4)*aData.TPPDT(i:i+20,1)...
        %+r1datalag(i,5)*aData.TPDT(i:i+20,1);
    epsilon_acg_P(i:i+20,i)=r1datalag(i,4)*aData.TPPDT(i:i+20,1);
    epsilon_acg_R(i:i+20,i)=r1datalag(i,5)*aData.TPDT(i:i+20,1);
    epsilon_acg_un(i:i+20,i)=r1datalag(i,3); %unexplained residue
    
 end

%  px22=plot(years,nanmean(epsilon_acg,2),'-','LineWidth',1,'color',[0.5 0.5 0.5]);
%  px23=plot(years,nanmean(epsilon_acg_T,2),'-','LineWidth',1,'color',[0.5 0 0]);
%  px24=plot(years,nanmean(epsilon_acg_P,2),'-','LineWidth',1,'color',[0 0 0.5]);
%  plot(years,nanmean(epsilon_acg_T,2),'-','LineWidth',2,'color',[1 0 0]);
%  plot(years,nanmean(epsilon_acg_P,2),'-','LineWidth',2,'color',[0 0 1]);
%  plot(years,nanmean(epsilon_acg_TP,2),'-','LineWidth',2,'color',[0 1 0]);

hold on;
epsilon_total=[nanmean(epsilon_acg_T,2),nanmean(epsilon_acg_R,2),nanmean(epsilon_acg_P,2),nanmean(epsilon_acg_un,2), nanmean(epsilon_acg_TP,2),nanmean(epsilon_acg_TT,2)];

epsilon_total1=epsilon_total;
epsilon_total1(epsilon_total1<0)=NaN;
epsilon_total2=epsilon_total;
epsilon_total2(epsilon_total2>0)=NaN;

bar_color={[1,0,0.1]; [0.9,0.9,0]; [0,0,1]; [0.9,0.9,0.9]; [0.5,0.5,1]; [1, 0.6, 0.2]};

b1=bar(epsilon_total1,'stacked','BarWidth',0.7);
set(b1,{'FaceColor'},bar_color);
set(b1,{'EdgeColor'},bar_color);

b2=bar(epsilon_total2,'stacked','BarWidth',0.7);
set(b2,{'FaceColor'},bar_color);
set(b2,{'EdgeColor'},bar_color);

px21=plot(years-1958,aData.CO2DT,'-','LineWidth',1,'color',[0 0 0]);
px22=plot(years-1958,nanmean(epsilon_acg,2),'--','LineWidth',1,'color',[0 0 0]);
 

 

legend([b2(1),b2(2),b2(3),b2(5),b2(6),b2(4)],{'$\Delta\sf{MAT}$','$\Delta\sf{MAR}$','$\Delta\sf{MAP}$','$\Delta\sf{MAT}\times \rm\overline\Delta \sf\overline{MAP}$','$\Delta\sf{MAT}\times\sf\textit{t}$','$\sf{others}$'},'Orientation','horizontal','Location', 'south','Box','off','Interpreter','latex');
% \Delta can only display in \rm form
%legend([px21, px22],{'obs.','pred.'},'Orientation','vertical','Location',
%'northeast','Box','off');

 
%xlim([1959 2016]);
set(gca,'XTick',3:10:60);
set(gca,'xticklabels',{'1960','1970','1980','1990','2000','2010'});

ylim([-3 2.5]);
ylabel(strcat(char(916),'CGR (PgC yr^{-1})'),'FontSize',10);
xlabel('Year','FontSize',10); 
text(1.5,2.3,'(b)','FontSize',12);

ax2=axes('position',get(gca,'position'),'visible','off');
legend(ax2,[px21, px22],{'obs.','pred.'},'Orientation','horizontal','Position', [0.14    0.15    0.2    0.08],'Box','off');

set(f7,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure3','-djpeg','-r600');  









%Figure 9. plot the obs a2, TRENDY a2, CMIP5 a2, Fluxcom a2, just for NEE.
%plot the distribution of their a,b,c values (a2=a*T+b*PP+c)

%get the a2
lag_choice=[10,20,30];

lag=lag_choice(2);
%1, in the time window, calculate sensitivity

%calculate obs a2
ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

for i=1:length(ydata_sen)-lag

    %sensitivity calculation
    [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
    r1datalag(i,1)=b(2,1);
    %r1datalag(i,2)=bint(2,1)-b(2,1);
    cd('/Volumes/RemiLBNL/project7/code');
    r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
    r1datalag(i,3)=b(1,1);
    
    yearlag(i)=yearGCB(i)+round(lag*0.5);
    xdatalag(i,1)=mean(aData.TPPDT(i:i+lag)); %get average meteo
    %xdatalag(i,2)=mean(aData.tropicalT(i:i+lag));
    xdatalag(i,2)=yearlag(i);

end

obs_a2=[yearlag',r1datalag(:,1),r1datalag(:,2)./2./sqrt(1000)]; %change sd to se

tmp_y=[];

clear obs_b123 td_b123 cmip_b123 fc_b123

for j=1:500
    
    tmp_y=normrnd(r1datalag(:,1),r1datalag(:,2));
    %[b,~,stats] = glmfit(xdatalag,r1datalag(:,1),'normal');
    [b,~,stats] = glmfit(xdatalag,tmp_y,'normal');

    obs_b123(:,j)=[b;stats.se];
end

%calculate the coefficient for simulating a2


%calculate TRENDY a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=aData.modelzDT.*(-1); %convert NEP to NEE

num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT aData.TPPDT aData.TPDT];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearGCB(i)+round(lag*0.5);
        xdatalag(i,1)=mean(aData.TPPDT(i:i+lag)); %get average meteo
        %xdatalag(i,2)=mean(aData.tropicalT(i:i+lag));
        xdatalag(i,2)=yearlag(i);

    end
    
    [b,~,stats] = glmfit(xdatalag,r1datalag(:,1,j),'normal');

    td_b123(:,j)=[b; stats.se];
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

td_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];



%calculate CMIP5 a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=CMIP5.NEPDT.*(-1); %NBP to NEE

yearCMIP=1959:2005;
num_models=size(ydata_sen,2);

for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) bData.TTDT(1:size(ydata_sen(:,j)),j) bData.TPPDT(1:size(ydata_sen(:,j)),j) aData.TPDT(1:47)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearCMIP(i)+round(lag*0.5);
        
        xdatalag(i,1)=mean(bData.TPPDT(i:i+lag,j))-273.15; %get average meteo, change from K to C
        %xdatalag(i,2)=mean(bData.tropicalT(i:i+lag,j));
        xdatalag(i,2)=yearlag(i);

    end
    
    [b,~,stats] = glmfit(xdatalag,r1datalag(:,1,j),'normal');

    cmip_b123(:,j)=[b;stats.se];
end

cmip_b123(:,8)=NaN; %one model is not coupled, inmcm4

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

cmip_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];




%calculate FluxCom a2
clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
ydata_sen=FC.NEEDT.*(1); %NBP to NEE

yearFC=1980:2013;
num_models=size(ydata_sen,2);

% 1959:2017
% 1980:2013

%change lag to 10 to get more values
lag=10;
for j=1:num_models
    xdata_sen=[ones(size(ydata_sen(:,j))) aData.TTDT(22:55) aData.TPPDT(22:55) aData.TPDT(22:55)];
    %xdata_sen=[ones(size(ydata_sen(:,j))) cData.TTDT(2:35) cData.TPPDT(2:35) aData.TPDT(22:55)];

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag,j),xdata_sen(i:i+lag,:));
        r1datalag(i,1,j)=b(2,1);
        r1datalag(i,2,j)=bint(2,1)-b(2,1);
        r1datalag(i,3,j)=b(1,1);

        yearlag(i)=yearFC(i)+round(lag*0.5);
        
        xdatalag(i,1)=mean(aData.TPPDT(1+i:1+i+lag)); %get average meteo
        %xdatalag(i,2)=mean(aData.tropicalT(1+i:1+i+lag));
        xdatalag(i,2)=yearlag(i);
        
%         xdatalag(i,1)=mean(aData.tropicalT(21+i:21+i+lag)); %get average meteo
%         xdatalag(i,2)=mean(aData.tropicalPP(21+i:21+i+lag));

    end
    
    [b,~,stats] = glmfit(xdatalag,r1datalag(:,1,j),'normal');
    fc_b123(:,j)=[b;stats.se];
end

en_mean_a2=squeeze(nanmean(r1datalag(:,1,:),3));
en_sd_a2=squeeze(nanstd(r1datalag(:,1,:),1,3));

fc_a2=[yearlag',en_mean_a2,en_sd_a2./sqrt(num_models)];

 
f9=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 22, 12], ...
'OuterPosition', [2, 2, 22, 12]);
 hold on;
 
 %plot a2
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.05,'ML',0.1,'MarginTop',.03,'MarginBottom',.12); 
 hold on;
 %observation
 px1=plot(obs_a2(:,1),obs_a2(:,2),'-','LineWidth',2,'color',[0.2 0.2 0.2]);
 
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(obs_a2(:,1),obs_a2(:,2),obs_a2(:,3),{'-'},0);

        set(H.patch,'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2],'FaceAlpha',0.2,'EdgeAlpha',0.2)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
 
 %trendy
 px2=plot(td_a2(:,1),td_a2(:,2),'-','LineWidth',1,'color',[1 0.2 0.2]);
 
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(td_a2(:,1),td_a2(:,2),td_a2(:,3),{'-'},0);

        set(H.patch,'FaceColor',[1 0.2 0.2],'EdgeColor',[1 0.2 0.2],'FaceAlpha',0.5,'EdgeAlpha',0.2)
        set(H.mainLine,'color',[1 0.2 0.2]);
        set(H.edge,'color','none');
 
 %cmip5
 px3=plot(cmip_a2(:,1),cmip_a2(:,2),'-','LineWidth',1,'color',[1 0.8 0.2]);
  cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(cmip_a2(:,1),cmip_a2(:,2),cmip_a2(:,3),{'-'},0);

        set(H.patch,'FaceColor',[1 0.8 0.2],'EdgeColor',[1 0.8 0.2],'FaceAlpha',0.5,'EdgeAlpha',0.2)
        set(H.mainLine,'color',[1 0.8 0.2]);
        set(H.edge,'color','none');
 
 %fluxcom
 px4=plot(fc_a2(:,1),fc_a2(:,2),'-','LineWidth',1,'color',[0.2 0.2 1]);
   cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(fc_a2(:,1),fc_a2(:,2),fc_a2(:,3),{'-'},0);

        set(H.patch,'FaceColor',[0.2 0.2 1],'EdgeColor',[0.2 0.2 1],'FaceAlpha',0.5,'EdgeAlpha',0.2)
        set(H.mainLine,'color',[0.2 0.2 1]);
        set(H.edge,'color','none');
        
        
 xlim([1969 2009]);
 ylim([-0.5 7]);
 tmp_yl=strcat(char(947),'_{CGR}','^T',{' '}, 'or',{' '},char(947),'_{NEE}','^T');
 ylabel({tmp_yl{1,1}; '(PgC yr^{-1}°C^{-1})'},'FontSize',12);
 xlabel('Year','FontSize',12);
 text(1971.3,6.6,'(a)','FontSize',12);
 
 legend([px1,px2,px3,px4],{'obs.','DGVMs','CMIP5','FLUXCOM'},'Orientation','vertical','Position',[0.13    0.17    0.1    0.1],'Box','off');
 
 %plot b(2)
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 subaxis(1,2,2,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.05,'ML',0.1,'MarginTop',.05,'MarginBottom',.05); 
 hold on;
 
 s1=scatter(td_b123(2,:), td_b123(3,:), 50 ,[1 0.2 0.2]);
 s1.LineWidth = 1;
 s1.MarkerEdgeColor = [1 0.2 0.2];
 s1.MarkerFaceColor = [1 0.2 0.2];
 s1.MarkerFaceAlpha = 0.8;
 s1.MarkerEdgeAlpha = 0.8;
 
 % add error bar
 e1=errorbar(td_b123(2,:), td_b123(3,:),td_b123(6,:),'o');
 e2=herrorbar(td_b123(2,:), td_b123(3,:),td_b123(5,:),'o');
 
 e1.Marker = 'o';
 e1.MarkerSize = 1;
 e1.MarkerEdgeColor=[1 0.2 0.2];
 e1.MarkerFaceColor='none';
 e1.Color = [1 0.2 0.2];

 e2(2,1).Marker = 'o';
 e2(2,1).MarkerEdgeColor=[1 0.2 0.2];
 e2(2,1).MarkerFaceColor='none';
 e2(2,1).MarkerSize = 1;
 e2(1,1).Color = [1 0.2 0.2];
 
 
 
 s2=scatter(cmip_b123(2,:), cmip_b123(3,:), 50 ,[0.5 0.5 0]);
 s2.LineWidth = 1;
 s2.MarkerEdgeColor = [1 0.8 0.2];
 s2.MarkerFaceColor = [1 0.8 0.2];
 s2.MarkerFaceAlpha = 0.8;
 s2.MarkerEdgeAlpha = 0.8;
 
 
 % add error bar
 e1=errorbar(cmip_b123(2,:), cmip_b123(3,:),cmip_b123(6,:),'o');
 e2=herrorbar(cmip_b123(2,:), cmip_b123(3,:),cmip_b123(5,:),'o');
 
 e1.Marker = 'o';
 e1.MarkerSize = 1;
 e1.MarkerEdgeColor=[1 0.8 0.2];
 e1.MarkerFaceColor='none';
 e1.Color = [1 0.8 0.2];

 e2(2,1).Marker = 'o';
 e2(2,1).MarkerEdgeColor=[1 0.8 0.2];
 e2(2,1).MarkerFaceColor='none';
 e2(2,1).MarkerSize = 1;
 e2(1,1).Color = [1 0.8 0.2];
 
 
 
 s3=scatter(fc_b123(2,:), fc_b123(3,:), 50,[0 0 0.5]);
 s3.LineWidth = 1;
 s3.MarkerEdgeColor = [0.2 0.2 1];
 s3.MarkerFaceColor = [0.2 0.2 1];
 s3.MarkerFaceAlpha = 0.8;
 s3.MarkerEdgeAlpha = 0.8;
 
 % add error bar
 e1=errorbar(fc_b123(2,:), fc_b123(3,:),fc_b123(6,:),'o');
 e2=herrorbar(fc_b123(2,:), fc_b123(3,:),fc_b123(5,:),'o');
 
 e1.Marker = 'o';
 e1.MarkerSize = 1;
 e1.MarkerEdgeColor=[0.2 0.2 1];
 e1.MarkerFaceColor='none';
 e1.Color = [0.2 0.2 1];

 e2(2,1).Marker = 'o';
 e2(2,1).MarkerEdgeColor=[0.2 0.2 1];
 e2(2,1).MarkerFaceColor='none';
 e2(2,1).MarkerSize = 1;
 e2(1,1).Color = [0.2 0.2 1];
 
 
 s4=scatter(mean(obs_b123(2,:)), mean(obs_b123(3,:)), 70,[0 0 0],'p');
 s4.LineWidth = 1;
 s4.MarkerEdgeColor = [0.2 0.2 0.2];
 s4.MarkerFaceColor = [0.2 0.2 0.2];
 s4.MarkerFaceAlpha = 0.8;
 s4.MarkerEdgeAlpha = 0.8;
 
 % add error bar
 e1=errorbar(mean(obs_b123(2,:)), mean(obs_b123(3,:)),mean(obs_b123(6,:)),'o');
 e2=herrorbar(mean(obs_b123(2,:)), mean(obs_b123(3,:)),mean(obs_b123(5,:)),'o');
 
 e1.Marker = 'o';
 e1.MarkerSize = 1;
 e1.MarkerEdgeColor=[0.2 0.2 0.2];
 e1.MarkerFaceColor='none';
 e1.Color = [0.2 0.2 0.2];

 e2(2,1).Marker = 'o';
 e2(2,1).MarkerEdgeColor=[0.2 0.2 0.2];
 e2(2,1).MarkerFaceColor='none';
 e2(2,1).MarkerSize = 1;
 e2(1,1).Color = [0.2 0.2 0.2];
 
 ylim([-0.12 0.15]);
 
 PlotAxisAtOrigin(mean(obs_b123(2,:)), mean(obs_b123(3,:)));
 %errorbar(mean(obs_b123(2,:)), mean(obs_b123(3,:)),std(obs_b123(3,:),1));
 %herrorbar(mean(obs_b123(2,:)), mean(obs_b123(3,:)),std(obs_b123(2,:),1),std(obs_b123(2,:),1)); %horizontal option in errorbar is available after 2016b

%  xlim([-0.1 0.06]);

  text(-0.1,0.142,'(b)','Fontsize',12);
  
  text(-0.1,0.04,'parameter \bf\ita','Fontsize',10);
  text(-0.1,0.02,'\rm(PgC °C^{-1} mm^{-1})','Fontsize',10);
  
  text(-0.017,0.088,'parameter \bf\itb','Fontsize',10,'Rotation',90);
  text(-0.007,0.088,'\rm(PgC °C^{-1} yr^{-1})','Fontsize',10,'Rotation',90);
 %ylim([-6 5]);
  % ylim([-0.15 0.2]);
  
 ylabel({'parameter a (PgC °C^{-1} °C^{-1})'},'FontSize',12);
 legend([s4,s1,s2,s3],{'obs.','DGVMs','CMIP5','FLUXCOM'},'Orientation','vertical','Position',[0.57    0.13    0.1    0.1]);

 
 set(f9,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure4','-djpeg','-r600'); 
 
 
 
 
 
 
 
 
 
 
 
 
% GRACE vs CGR
f10=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 9], ...
'OuterPosition', [2, 2, 15, 9]);

hold on;

yearGCB=years;

color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

%plot GRACE against CGR
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.18);   
hold on;

ydata=[nan(20,1);y_GRACE(:,1)];

%ydata=aData.tropicalPP;

xlim_value(1)=min(ydata);
xlim_value(2)=max(ydata);
xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;
        
        
%         for i=1:length(ydata_sen)-lag
%             
%             yearlag(i)=yearGCB(i)+round(lag*0.5);
%             xdatalag(i,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata, if one year is NAN, then the decadal average is NAN;
%             ydatalag(i,1)=mean(ydata_sen(i:i+lag));
%         end
        
        
        xdatalag=ydata;
        ydatalag=ydata_sen;
        %ydatalag(ydatalag<-1.5)=NaN;
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;
            
        
            p=polyfit(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999),1);
            [r1,r2]=corrcoef(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999));
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p>0.05'); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),', p>0.05'); 
        end
        %xtitle={xtitle;xfun};
        %end of ploting
        %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
        
        
        %plot the change of y2 to x, exclude volcano years.
        faceColor='r';%color_choice(2,:);
        ydatalag(33:35,1)=NaN;
        
        pxx3=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',faceColor,...
    'MarkerFaceColor',faceColor);
        hold on;
            
        
            p=polyfit(xdatalag(xdatalag>-9999 & ydatalag>-9999),ydatalag(xdatalag>-9999 & ydatalag>-9999),1);
            [r1,r2]=corrcoef(xdatalag(xdatalag>-9999 & ydatalag>-9999),ydatalag(xdatalag>-9999 & ydatalag>-9999));
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
        if p(2)>0
        xfun2=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.05'); 
        else
            xfun2=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.05'); 
        end
        
        
        
        
    end   
    
     ylabel(strcat('\DeltaCGR (PgC yr^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel('$\sf{GRACE-rec}(cm)$','FontSize',10,'Interpreter','latex');
     text(xlim_value(1),-3.5,xfun);
     text(xlim_value(1),-2.5,xfun2,'Color','r');
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([-4 3]);
     
     text(xlim_value(1)-0.05*xlim_value(3),2.5,'(a)');
     clear xlim_value




%plot GRACE against sensitivity
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.18);   
hold on;

ydata=[nan(20,1);y_GRACE(:,1)];
%ydata=aData.tropicalPP;

xlim_value(1)=min(ydata);
xlim_value(2)=max(ydata);
xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata, if one year is NAN, then the decadal average is NAN;
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag(:,1);
         ydatalag_sd=r1datalag(:,2);
       
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);

        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;
            
            %%plot the uncertainty of y
            cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
%             H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);
% 
%             set(H.patch,'FaceColor',faceColor,'EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2)
%             set(H.mainLine,'color','none');
%             set(H.edge,'color','none');

        
            p=polyfit(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999),1);
            [r1,r2]=corrcoef(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999));
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.01'); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.01'); 
        end
        %xtitle={xtitle;xfun};
        %end of ploting
        %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel('$\sf\overline{GRACE-rec}(cm)$','FontSize',10,'Interpreter','latex');
     text(xlim_value(1),1.7,xfun);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([1 7]);
     
     %text(xlim_value(1)-0.05*xlim_value(3),6.5,strcat('(',char(96+(fign-1)*2+2),')'));
     text(xlim_value(1)-0.05*xlim_value(3),6.6,'(b)');
     clear xlim_value
 
     
     set(f10,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure10','-djpeg','-r600');  








%select the best model to simulate rT

%select the best model to simulate rT

load('/Volumes/RemiLBNL/project7/code/aic_matrix.mat');

% model_list={...
% '$\Delta \sf\overline{MAT}$', '$\itf\rm(\DeltaMAP)$','$\itf\rm(\DeltaMASW)$','$\itf\rm(Year)$', ...
% '$\itf\rm(\DeltaMAT,\DeltaMAP)$', '$\itf\rm(\DeltaMAT,\DeltaMASW)$',...
% '$\itf\rm(\DeltaMAP,\DeltaMASW)$',...
% '$\itf\rm(\DeltaMAT,Year)$',...
% '$\bf\itf\rm\bf(\DeltaMAP,Year)',...
% '\itf\rm(\DeltaMASW,Year)$', ...
% '$\itf\rm(\DeltaMAT,\DeltaMAP,\DeltaMASW)$', ...
% '$\itf\rm(\DeltaMAT,\DeltaMAP,Year)$',...
% '$\itf\rm(\DeltaMAT,\DeltaMASW,Year)$',...
% '$\itf\rm(\DeltaMAP,\DeltaMASW,Year)$',...
% '$\itf\rm(\DeltaMAT,\DeltaMAP,\DeltaMASW, Year)$'
% };

model_list={...
'\itf\rm(\DeltaMAT)', '\itf\rm(\DeltaMAP)','\itf\rm(\DeltaMASW)','\itf\rm(\itt\rm)', ...
'\itf\rm(\DeltaMAT,\DeltaMAP)', '\itf\rm(\DeltaMAT,\DeltaMASW)','\itf\rm(\DeltaMAP,\DeltaMASW)',...
'\itf\rm(\DeltaMAT,\itt\rm)','\bf\itf\rm\bf(\DeltaMAP,\itt\rm\bf)','\itf\rm(\DeltaMASW,\itt\rm)', ...
'\itf\rm(\DeltaMAT,\DeltaMAP,\DeltaMASW)', '\itf\rm(\DeltaMAT,\DeltaMAP,\itt\rm)','\itf\rm(\DeltaMAT,\DeltaMASW,\itt\rm)','\itf\rm(\DeltaMAP,\DeltaMASW,\itt\rm)',...
'\itf\rm(\DeltaMAT,\DeltaMAP,\DeltaMASW,\itt\rm)'
};


f11=figure('Name','AIC','Units', 'centimeters','Color','white', 'Position', [2, 2, 22, 15], ...
'OuterPosition', [2, 2, 22, 15]);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,9,1:5,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.28,'MarginTop',.05,'MarginBottom',.22);

zvalue=flip(glm_aic(:,6:end));
zvalue = [zvalue; zvalue(end,:)]; %pcolor will flip and remove the first row and col
zvalue = [zvalue, zvalue(:,end)];

h=pcolor(zvalue);
%flip(zvalues,1)
set(h, 'EdgeColor', 'none');

%color_s=jet(50);
color_s=gray(50);

colormap(gca,color_s(1:50,:)); %or flip(gray(50)) jet(50)
caxis([0 150]);

set(gca,'XTick',1.5:5:16.5);
%set(gca,'xticklabels',{'15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'});
set(gca,'xticklabels',{'15','20','25','30'});
xlabel('Length of the moving window','FontSize',10);
text(1,16.5,'(a)','FontSize',12);

set(gca,'YTick',1.5:15.5);

%set(0,'defaulttextinterpreter','tex');
set(gca,'yticklabels',flip(model_list));


cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,9,6:9,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.22);

zvalue2=flip(gam_aic(:,6:end));
zvalue2 = [zvalue2; zvalue2(end,:)]; %pcolor will flip and remove the first row and col
zvalue2 = [zvalue2, zvalue2(:,end)];

h2=pcolor(zvalue2);
set(h2, 'EdgeColor', 'none');

%color_s=jet(50);
color_s=gray(50);

colormap(gca,color_s(1:50,:)); %or flip(gray(50)) jet(50)
caxis([0 150]);

set(gca,'XTick',1.5:5:16.5);
%set(gca,'xticklabels',{'15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'});
set(gca,'xticklabels',{'15','20','25','30'});

% set(gca,'YTick',1.5:15.5);
set(gca,'yticklabels',[]);
xlabel('Length of the moving window','FontSize',10);
text(1,16.5,'(b)','FontSize',12);

h=colorbar('Position', [0.3  0.08  0.65  0.03], 'Orientation','Horizontal');%('Ticks',[-200,-150,-100,-50,0,50,100,150]); 
set(h,'YTick',0:30:150);
h.Label.String = 'Akaike information criterion (AIC)';
hold off;

print('/Volumes/RemiLBNL/project7/Figures/Figure11','-djpeg','-r600'); 








% AIC figure
f12=figure('Name','AIC','Units', 'centimeters','Color','white', 'Position', [2, 2, 22, 15], ...
'OuterPosition', [2, 2, 22, 15]);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,9,1:5,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.28,'MarginTop',.05,'MarginBottom',.22);

zvalue=flip(glm_r2(:,6:end));
zvalue = [zvalue; zvalue(end,:)]; %pcolor will flip and remove the first row and col
zvalue = [zvalue, zvalue(:,end)];

h=pcolor(zvalue);
%flip(zvalues,1)
set(h, 'EdgeColor', 'none');

%color_s=jet(50);
%color_s=flip(gray(50));
color_s=gray(50);

colormap(gca,color_s(1:50,:)); %or flip(gray(50)) jet(50)
caxis([0 1]);

set(gca,'XTick',1.5:5:16.5);
%set(gca,'xticklabels',{'15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'});
set(gca,'xticklabels',{'15','20','25','30'});
xlabel('Length of the moving window','FontSize',10);
text(1,16.5,'(a)','FontSize',12);

set(gca,'YTick',1.5:15.5);

%set(0,'defaulttextinterpreter','tex');
set(gca,'yticklabels',flip(model_list));


cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,9,6:9,'SpacingVert',0.06,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.22);

zvalue2=flip(gam_r2(:,6:end));
zvalue2 = [zvalue2; zvalue2(end,:)]; %pcolor will flip and remove the first row and col
zvalue2 = [zvalue2, zvalue2(:,end)];

h2=pcolor(zvalue2);
set(h2, 'EdgeColor', 'none');

%color_s=jet(50);
color_s=gray(50);

colormap(gca,color_s(1:50,:)); %or flip(gray(50)) jet(50)
caxis([0 1]);

set(gca,'XTick',1.5:5:16.5);
%set(gca,'xticklabels',{'15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'});
set(gca,'xticklabels',{'15','20','25','30'});

% set(gca,'YTick',1.5:15.5);
set(gca,'yticklabels',[]);
xlabel('Length of the moving window','FontSize',10);
text(1,16.5,'(b)','FontSize',12);

h=colorbar('Position', [0.3  0.08  0.65  0.03], 'Orientation','Horizontal');%('Ticks',[-200,-150,-100,-50,0,50,100,150]); 
set(h,'YTick',0:0.2:1);
h.Label.String = 'r^2';

print('/Volumes/RemiLBNL/project7/Figures/Figure12','-djpeg','-r600'); 










% only VPD
%create figure s2 which exhibit the correlation between EMD sensitivity and
%climate variables
f10=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 9], ...
'OuterPosition', [2, 2, 18, 9]);

hold on;

%month_ft=12;
yearGCB=years;

%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

for fign=1:2
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 2 % x t-vpd
        ydata=sum(aData.TVPDEMD(1:end-1,:),1);
        ydata_sen=sum(aData.CO2EMD(1:end-1,:),1)';
        xdata_sen=[ones(size(ydata_sen)) sum(aData.TTEMD(1:3,:),1)' sum(aData.TPEMD(1:3,:),1)' sum(aData.TPPEMD(1:3,:),1)'];
        ytitle='$\overline\Delta\sf\overline{MAD} (kPa)$';
         
        
        case 1 % x t-vpd
        ydata=aData.TVPDDT;
        ydata_sen=aData.CO2DT;
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPDT aData.TPPDT];
        ytitle='$\overline\Delta\sf\overline{MAD} (kPa)$';


    end
   
          
    %plotting response curve 
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(1,2,fign,'SpacingVert',0.03,'SpacingHoriz',0.1,'MR',0.05,'ML',0.1,'MarginTop',.03,'MarginBottom',.2);   
    hold on;
    xlim_value(1)=min(ydata);
    xlim_value(2)=max(ydata);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity
        
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            %r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
       
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);

        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;
            
            %%plot the uncertainty of y
%             cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
%             H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);
% 
%             set(H.patch,'FaceColor',faceColor,'EdgeColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2)
%             set(H.mainLine,'color','none');
%             set(H.edge,'color','none');

        
            p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig<0.01
             ptext='p<0.01';
         elseif sig<0.05
             ptext='p<0.05';
         else
             ptext='p>0.1';
         end
         
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end
        %xtitle={xtitle;xfun};
        %end of ploting
        %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel(ytitle,'FontSize',10,'Interpreter','Latex');
     text(xlim_value(1),1.7,xfun);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([1 7]);
     
     text(xlim_value(1)-0.05*xlim_value(3),6.5,strcat('(',char(96+fign),')'));
     clear xlim_value
     
     
end

hold off;

set(f10,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/FigureS3B','-djpeg','-r600');  















% round 1 revision added figure
% discrete regression relationships - to deal with auto-correlation






% test regression with other Tair and PP datasets
% panel (a) 
% BEST, CRUTEM, GISS, UDEL

% panel (a)apparent temporal sensitivity
f1r=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 25, 10], ...
    'OuterPosition', [2, 2, 25, 10]);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 


%color_choice=[1,0.4,0.4;0,0,0;0.4,0.4,1;0.4,1,1];

color_choice=lines(4);

lag_choice=[15,20,25];

for ttt=1:4
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(2);
        

    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    %ydata_sen=sum(aData.CO2EMD(1:3,:),1)';
    %xdata_sen=[ones(size(ydata_sen)) aData.TTDT];
        
    switch ttt
        case 1
            xdata_sen=[ones(size(ydata_sen)) apData.BESTDT_T' aData.TPPDT aData.TPDT];
        case 2
            xdata_sen=[ones(size(ydata_sen)) apData.CRUTEMDT_T' aData.TPPDT aData.TPDT];
        case 3
            xdata_sen=[ones(size(ydata_sen)) apData.GISSDT_T' aData.TPPDT aData.TPDT];
        case 4
            xdata_sen=[ones(size(ydata_sen)) apData.UDELDT_T' aData.TPPDT aData.TPDT];
    end
    
    xlim_value=nan(3,1);

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Volumes/RemiLBNL/project7/code');
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
        cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.2,'EdgeAlpha',0.2)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


    
     %ylabel(strcat('dCGR/','dMAT (PgC yr^{-1}°C^{-1})'),'FontSize',12);
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     xlabel('Year','FontSize',12);
     set(gca,'box','off');
     
     %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     xlim([1965.5 2010.5]);
     ylim([0 8]);
     
     set(gca,'xtick', 1970:10:2010);
     
     text(1967,7.8,'(a)','Fontsize',12);
     
 legend([pxx1(1),pxx1(2),pxx1(3),pxx1(4)],{'BEST','CRUTEM','GISS','UDEL'},'Orientation','vertical','Position',[0.08    0.7    0.1    0.1],'Box','off');


% panel (b)
% GPCC, CPC, PRECL, UDEL

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.02,'ML',0.06,'MarginTop',.03,'MarginBottom',.15); 

for ttt=1:4
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (ttt)
        case 1 % x GPCC      
        ydata2=apData.GPCCDT_P;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';
         
        
        case 2 % x CPC
        ydata2=apData.CPCDT_P;
        %ydata2=apData.GHCNDT_P;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';

                
        case 3 % x PRECL       
        ydata2=apData.PRECLDT_P;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';
        
        
        case 4 % x UDEL       
        ydata2=apData.UDELDT_P;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';


    end
        
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    %%%%start ploting, long term meteo
    
     
    %plotting response curve  
    hold on;
    xlim_value(1)=min(ydata2);
    xlim_value(2)=max(ydata2);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata2(i:i+lag)); %meteo was assigned to ydata
            
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
      
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);

        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        %plot the change of y to x
        faceColor=color_choice(ttt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.1,'MarkerFaceColor',faceColor,'MarkerEdgeColor','none');
        
        pxx2(ttt)=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',faceColor,...
    'MarkerFaceColor',faceColor);
        hold on;
            
            %%plot the uncertainty of y
%             cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
%             H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);
% 
%             set(H.patch,'FaceColor',faceColor,'EdgeColor',faceColor,'FaceAlpha',0.2,'EdgeAlpha',0.2)
%             set(H.mainLine,'color','none');
%             set(H.edge,'color','none');

           %errorbar(xdatalag,ydatalag,ydatalag_sd, 'LineStyle','none','LineWidth',0.1,'MarkerFaceColor',faceColor,'MarkerEdgeColor','none');
        
            tmp_ind=isnan(xdatalag);
            
            p=polyfit(xdatalag(~tmp_ind),ydatalag(~tmp_ind),1);
            [r1,r2]=corrcoef(xdatalag(~tmp_ind),ydatalag(~tmp_ind));
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '-';
            h.LineWidth=2;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig < 0.01
             ptext='p<0.01';
         elseif sig < 0.05
             ptext='p<0.05';
         else
             ptext='p>0.1';
         end
            
  
            
        %xfun=strcat('Y_{',num2str(lag),'}=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p=',num2str(sig)); 
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end

        
    end   
    
     
     clear xlim_value
end

xlim([-30 20]);
ylim([0 8]);
text(-28,7.8,'(b)','Fontsize',12);
set(gca,'box','off');

ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);

xlabel('$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$','FontSize',10,'Interpreter','Latex');


legend([pxx2(1),pxx2(2),pxx2(3),pxx2(4)],{'GPCC','CPC','PRECL','UDEL'},'Orientation','vertical','Position',[0.56    0.25    0.1    0.1],'Box','off');

hold off;

set(f1r,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure1r','-djpeg','-r600');  







% revision figure 
% figure used to detect the discrete correlation between MAP and
% sensitivity
fr2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 20], ...
'OuterPosition', [2, 2, 18, 20]);

hold on;

%month_ft=12;
yearGCB=years;

%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

for fign=2:2
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 1 % x t-tmp
        ydata=aData.tropicalT;
        ytitle='$\sf{MAT}\ and\ \overline{MAT} (^\circ C)$'; 
        
        ydata2=aData.TTDT;
        ytitle2='$\overline\Delta\sf\overline{MAT} (^\circ C)$';
         
        
        case 2 % x t-pp
        ydata=aData.TPPDT;
        ytitle='$\sf{MAP}\ and\ \overline{MAP}(mm\ yr^{-1})$';
        
        ydata2=aData.TPPDT;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';
            
                
        case 3 % x t-swc
        ydata=aData.tropicalSWC;
        ytitle='$\sf{MASW}\ and\ \overline{MASW} (\%)$';
        
        ydata2=aData.TSWCDT;
        ytitle2='$\overline\Delta \sf\overline{MASW} (\%)$';


    end
        
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;

    
    xlim_value(1)=min(ydata2);
    xlim_value(2)=max(ydata2);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,r,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            %r1datalag(i,4)=r(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata2(i:i+lag)); %meteo was assigned to ydata
            
        end
        
        
        
        % newly added content for auto-correlation
        auto_res=[];
        count=1;
        
        for auto_lag=1:20
            
            [tmp1,tmp2] = corrcoef (r1datalag(1:end-auto_lag,1),r1datalag(1+auto_lag:end,1));
            auto_res(count,1)=tmp1(2); % R2
            auto_res(count,2)=tmp2(2); % sig
            
             ydatalag=r1datalag(1:auto_lag:end,1);
             xdatalag=temp_xdatalag(1:auto_lag:end);
             p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            auto_res(count,3)=p(1); % slope
            auto_res(count,4)=rsq; % rsq
            auto_res(count,5)=sig; % sig
            
            count=count+1;
        end
        
               
         
        
    end   
    

end


% reorganise p values
tmp=auto_res(:,2);
tmp(tmp<0.05)=0;
tmp(tmp>0.1)=1;
tmp(tmp>0.05 & tmp<0.1)=0.1;
auto_res(:,2)=tmp;

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(3,1,1,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.12,'MarginTop',.01,'MarginBottom',.08);   

hold on;
zvalue2=[auto_res(:,2)';auto_res(:,2)'];
zvalue2 = [zvalue2; zvalue2(end,:)]; %pcolor will flip and remove the first row and col
zvalue2 = [zvalue2, zvalue2(:,end)];

h=pcolor(0.5:1:20.5,-1:1:1,zvalue2);%if interval too large, then the color cannot show
set(h, 'EdgeColor', 'none');
set(gca,'YTick',-1:1:1);
%set(gca,'XTick',1:1:20);

% color_s=[0.6 0.6 0.6;
%     0.8 0.8 0.8;
%     1 1 1];

color_s=[0.6:0.05:1;
    0.6:0.05:1;
    0.6:0.05:1]';


colormap(gca,color_s); 
%caxis([0;0.05;0.1;1]);

h=refline(0,0);
h.Color = 'k';
h.LineStyle = '--';
h.LineWidth=1.5;

plot(1:20, auto_res(:,1),'-','LineWidth',2,'color',[1 0 0]);
xlim([1 20]);

ylabel({'Degree of Autocorrelation'},'FontSize',11);
text(1.5,0.85,'(a)','FontSize',12);
xlabel('The time lag for autocorrelation test','FontSize',10);


% panel (b)
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(3,1,2,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.12,'MarginTop',.01,'MarginBottom',.08); 

hold on;
h=pcolor(0.5:1:20.5,-0.2:0.1:0,zvalue2);%if interval too large, then the color cannot show
set(h, 'EdgeColor', 'none');
set(gca,'YTick',-0.2:0.1:0);

color_s=[0.6:0.05:1;
    0.6:0.05:1;
    0.6:0.05:1]';

colormap(gca,color_s); 


plot(1:20, auto_res(:,3),'-','LineWidth',2,'color',[1 0 0]);
xlim([1 20]);
ylim([-0.2 0]);        
ylabel(strcat('d',char(947),'_{CGR}','^T/d\DeltaMAP'),'FontSize',10);
xlabel(['The time interval used to resample ',char(947),'_{CGR}','^T'],'FontSize',10);
text(1.5,-0.02,'(b)','FontSize',12);
        
       

hold on;






lag_choice=[15,20,25];
for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
    lag=lag_choice(tt);
    %1, in the time window, calculate sensitivity

    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
    clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
    
    count=1;
    
    for i=1:9:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(count,1)=b(2,1);
        r1datalag(count,2)=bint(2,1)-b(2,1);
        r1datalag(count,3)=b(1,1);
        r1datalag(count,4)=b(3,1); %coefficient for TPP
        r1datalag(count,5)=b(4,1);% for TP

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(count,1)=round(tmp2(2)*1000)./1000;

        yearlag(count)=yearGCB(i)+round(lag*0.5);
        xdatalag(count,1)=mean(aData.TTDT(i:i+lag)); %get average meteo
        xdatalag(count,2)=mean(aData.TPPDT(i:i+lag));
        xdatalag(count,3)=mean(aData.TSWCDT(i:i+lag));
        xdatalag(count,4)=mean(aData.TPDT(i:i+lag));
        %xdatalag(count,5)=yearlag(i)-1968;

         %xdatalag(count,5)=mean(aData.tropicalT(i:i+lag)); %get average meteo
        xdatalag(count,5)=i;
         count=count+1;
    end
    
    
        i=38; % add a separate year
        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(count,1)=b(2,1);
        r1datalag(count,2)=bint(2,1)-b(2,1);
        r1datalag(count,3)=b(1,1);
        r1datalag(count,4)=b(3,1); %coefficient for TPP
        r1datalag(count,5)=b(4,1);% for TP

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(count,1)=round(tmp2(2)*1000)./1000;

        yearlag(count)=yearGCB(i)+round(lag*0.5);
        xdatalag(count,1)=mean(aData.TTDT(i:i+lag)); %get average meteo
        xdatalag(count,2)=mean(aData.TPPDT(i:i+lag));
        xdatalag(count,3)=mean(aData.TSWCDT(i:i+lag));
        xdatalag(count,4)=mean(aData.TPDT(i:i+lag));
        %xdatalag(count,5)=yearlag(i)-1968;

        xdatalag(count,5)=mean(aData.tropicalT(i:i+lag)); %get average meteo

         count=count+1;
    
    

    
     ydatalag=r1datalag(:,1);

    %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
end
 x1=xdatalag(:,2);
 x2=xdatalag(:,5);
 x3=xdatalag(:,1);
 y=ydatalag;
 
 x=xdatalag(:,[2,5]);
 %x=horzcat(xdatalag(:,[1,i+1]),xdatalag(:,1).*xdatalag(:,i+1));
 [b,~,stats] = glmfit(x,ydatalag,'normal');

 
 %y=b1+b2*PP+b3*Year;
 %set baseline in a nominal year with 24C and 1300 mm
%  f7=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 20], ...
% 'OuterPosition', [2, 2, 20, 20]);
 hold on;
 
 %plot a2
 cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 subaxis(3,1,3,'SpacingVert',0.1,'SpacingHoriz',0.08,'MR',0.03,'ML',0.12,'MarginTop',.03,'MarginBottom',.1); 
 hold on;
 %observation
 px1=plot(yearlag,ydatalag,'-','LineWidth',2,'color',[0 0 0]);
 
%  %simulation
  %y=b(2)*x1+b(3)*x2+b(1) - (b(2)*x1(1)+b(3)*x2(1)+b(1));
  y=b(2)*x1+b(3)*x2+b(1);
  px2=plot(yearlag,y,'--','LineWidth',2,'color',[0.5 0.5 0.5]);
 
 %contribution from P
 %y1=b(2)*x1+b(1) - (b(2)*x1(1)+b(1)); %only from PP
 y1=b(2)*x1+b(1) + b(3)*x2(1); %only from PP
 px3=plot(yearlag,y1,'--','LineWidth',2,'color',[0.5,0.5,1]);
 
 %y2=b(3)*x2 - (b(3)*x2(1)); %adjustment put by MAT
 y2=b(3)*x2 +b(1) + (b(2)*x1(1)); %adjustment put by MAT
 px4=plot(yearlag,y2,'--','LineWidth',2,'color',[1, 0.6, 0.2]);
 y12=b(2)*x1+b(3)*x2+b(1);
 
%   y3=b(4)*x3+b(1); %only from time
%  px5=plot(yearlag,y3,'--','LineWidth',2,'color',[0.3 0.8 0.3]);
 
 h=refline(0,2);
 h.Color = 'k';


 legend([px1,px2,px3,px4],{'$\sf{obs.}$', '$\sf{pred.}$','$f\big(\overline\Delta\sf\overline{MAP}\big)$','$f\big(\sf\overline{MAT}\big)$'},'Orientation','horizontal','Location', 'south','Box','off','Interpreter','latex')
 
 %xlim([1959 2016]);
 xlim([1968 2008]);
 ylim([-2 8]);
 ylabel({strcat(char(947),'_{CGR}','^T'); '(PgC yr^{-1}°C^{-1})'},'FontSize',10);
 %set(gca,'XTickLabel',[]);
 xlabel('The center year of each 20-year window','FontSize',10);
 text(1969,7.2,'(c)','FontSize',12);





set(fr2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure2r','-djpeg','-r600');  











% revision add a new figure
% add a figure to display future climate and future sensitivity

% revision add a new figure
% add a figure to display future climate and future sensitivity

fr3=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 10], ...
'OuterPosition', [2, 2, 15, 10]);

hold on;

yearGCB=1959:2016;
yearRCP=2006:2100;
yearALL=1959:2100;

CMIP5.Models={'BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
%CMIP5.Models={'BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','HadGEM2-ES','inmcm4','IPSL-CM5A','IPSL-CM5B','MIROC-ESM','MIROC-ESM-CHEM','MPI-ESM','NorESM1'};
    

% panel (a) % plot future P under RCP 4.5
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,1,1,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.3,'ML',0.12,'MarginTop',.03,'MarginBottom',.15); 

clear yearlag ydatalag

lag = 20;
for i = 1: length(yearGCB)-lag
    
    yearlag(i)=yearGCB(i)+round(lag*0.5);
    ydatalag(i,1)=mean(aData.TPPDT(i:i+lag));
    
end

pxx1=plot(yearlag,ydatalag,'-','LineWidth',2,'color','k');
hold on;

% rcp 45 temperature
ydata=[bData.TPPDT; rcpData.TPPDT45];

clear yearlag ydatalag

lag = 20;
for i = 1: length(yearALL)-lag
    
    yearlag(i)=yearALL(i)+round(lag*0.5);
    %ydatalag(i,:)=mean(rcpData.TPPDT45(i:i+lag,:),1);
    ydatalag(i,:)=mean(ydata(i:i+lag,:),1);
    
end

ydatalag(:,[5,8,10,12])=[];
color_s = lines(10);

for j=1:10
    
    if j<8
        pxx2b=plot(yearlag,ydatalag(:,j),'-','LineWidth',0.5,'Color',color_s(j,:));
    else
        pxx2b=plot(yearlag,ydatalag(:,j),'--','LineWidth',0.5,'Color',color_s(j-7,:));
    end
    
end



xlim([1959 2100]);
ylim([-40 40]);
ylabel('$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$','FontSize',10,'Interpreter','Latex');
set(gca,'xticklabels','','box','off');
text(1961,35,'(a)','FontSize',12);


% panel (b) 
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,1,2,'SpacingVert',0.08,'SpacingHoriz',0.05,'MR',0.3,'ML',0.12,'MarginTop',.03,'MarginBottom',.15); 

clear yearlag ydatalag

lag = 20;
for i = 1: length(yearGCB)-lag
    
    yearlag(i)=yearGCB(i)+round(lag*0.5);
    ydatalag(i,1)=mean(aData.TPPDT(i:i+lag));
    
end

pxx2=plot(yearlag,ydatalag,'-','LineWidth',2,'color','k');
hold on;

% rcp 85 temperature
ydata=[bData.TPPDT; rcpData.TPPDT85];

clear yearlag ydatalag

lag = 20;
for i = 1: length(yearALL)-lag
    
    yearlag(i)=yearALL(i)+round(lag*0.5);
    %ydatalag(i,:)=mean(rcpData.TPPDT85(i:i+lag,:),1);
    ydatalag(i,:)=mean(ydata(i:i+lag,:),1);
    
end

ydatalag(:,[5,8,10,12])=[];
color_s = lines(10);

for j=1:10
    
    if j<8
        pxx2b=plot(yearlag,ydatalag(:,j),'-','LineWidth',0.5,'Color',color_s(j,:));
    else
        pxx2b=plot(yearlag,ydatalag(:,j),'--','LineWidth',0.5,'Color',color_s(j-7,:));
    end
    
end


xlim([1959 2100]);
ylim([-40 40]);
ylabel('$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$','FontSize',10,'Interpreter','Latex');
%set(gca,'xticklabels','','box','off');
set(gca,'box','off');
text(1961,35,'(b)','FontSize',12);

legend({'obs.','BNU-ESM','CanESM2','CCSM4','GFDL-ESM2G','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A','MIROC-ESM','MPI-ESM','NorESM1'},'Orientation','vertical','Position', [0.8    0.2    0.1    0.7],'Box','off');
%legend([px1, pxx2b(1), pxx2b(2)],{'obs.','BNU-ESM','CanESM2'});
    
set(fr3,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure3r','-djpeg','-r600');  


















% correct auto, random codes
% two big sections

%%% Figure 3. plot climatic variables and its influence on sensitivity
%%% T (linear), PP, swc, VPD
f6=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 24], ...
'OuterPosition', [2, 2, 18, 24]);

hold on;

%month_ft=12;
yearGCB=years;

%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,20,30];

for fign=1:3
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 1 % x t-tmp        
        ydata2=aData.TTDT;
        ytitle2='$\overline\Delta\sf\overline{MAT_{adj}} (^\circ C)$';
         
        
        case 2 % x t-pp       
        ydata2=aData.TPPDT;
        ytitle2='$\overline\Delta\sf\overline{MAP_{adj}}(mm\ yr^{-1})$';
        
                
        case 3 % x t-swc       
        ydata2=aData.TSWCDT;
        ytitle2='$\overline\Delta \sf\overline{MASW_{adj}} (\%)$';


    end
        
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    %%%%start ploting, long term meteo
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');     
     
    %plotting response curve 
    subaxis(4,2,(fign-1)*2+1,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08);   
    hold on;
    xlim_value(1)=min(ydata2);
    xlim_value(2)=max(ydata2);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata2(i:i+lag)); %meteo was assigned to ydata
            
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
         
         
         auto_pcy = polyfit(ydatalag(1:end-1),ydatalag(2:end),1);
         auto_pcx = polyfit(xdatalag(1:end-1),xdatalag(2:end),1);
        
        
        % added part for auto-correlation
        tmp_r=corrcoef(xdatalag(1:end-1),xdatalag(2:end));
        
             mdl = fitlm(yearlag(1:end-1),ydatalag(1:end-1)-auto_pcy(1).*ydatalag(2:end)); %0.95
            [p1,DW1] = dwtest(mdl,'exact','both');
            DW1=round(DW1.*100)./100;
            %disp(p1);disp(DW1);
            
            mdl = fitlm(yearlag(1:end-1),xdatalag(1:end-1)-auto_pcx(1).*xdatalag(2:end)); %1
            [p2,DW2] = dwtest(mdl,'exact','both');
            DW2=round(DW2.*100)./100;
            %disp(p2);disp(DW2);
        
        
        
        no_auto_y=ydatalag(1:end-1)-0.95.*ydatalag(2:end);
        no_auto_x=xdatalag(1:end-1)-xdatalag(2:end);
        no_auto_ysd=ydatalag_sd(1:end-1);
        
        
        xdatalag=[];ydatalag=[];ydatalag_sd=[];
        [xdatalag,sortInd]=sort(no_auto_x);
        ydatalag=no_auto_y(sortInd);
        ydatalag_sd=no_auto_ysd(sortInd);
        

        xlim_value(1)=min(no_auto_x);
        xlim_value(2)=max(no_auto_x);
        xlim_value(3)=xlim_value(2)-xlim_value(1);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;

        
            p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig < 0.01
             ptext='p<0.01';
         elseif sig < 0.05
             ptext='p<0.05';
         elseif sig < 0.1
             ptext='p<0.1';
         else
             ptext='p>0.1';
         end
            

        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end
        
        
        xfun2=strcat('DW_Y = ',num2str(DW1),{' '},'(p>0.1)'); 
        xfun3=strcat('DW_X = ',num2str(DW2),{' '},'(p>0.1)'); 
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T_{adj} (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     xlabel(ytitle2,'FontSize',10,'Interpreter','Latex');
     text(xlim_value(1),-1.5,xfun);
     text(xlim_value(1),-2.2,xfun2);
     text(xlim_value(1),-2.8,xfun3);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([-3.1 2.2]);
     
     text(xlim_value(1)-0.05*xlim_value(3),1.9,strcat('(',char(96+(fign-1)*2+1),')'));
     clear xlim_value
end



% check the stats of the figure above
lag_choice=[15,20,25];
for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
    lag=lag_choice(tt);
    %1, in the time window, calculate sensitivity

    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
    clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=b(3,1); %coefficient for TPP
        r1datalag(i,5)=b(4,1);% for TP

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(i,1)=round(tmp2(2)*1000)./1000;

        yearlag(i)=yearGCB(i)+round(lag*0.5);
        xdatalag(i,1)=mean(aData.TTDT(i:i+lag)); %get average meteo
        xdatalag(i,2)=mean(aData.TPPDT(i:i+lag));
        xdatalag(i,3)=mean(aData.TSWCDT(i:i+lag));
        xdatalag(i,4)=mean(aData.TPDT(i:i+lag));
        %xdatalag(i,5)=yearlag(i)-1968;

         xdatalag(i,5)=mean(aData.tropicalT(i:i+lag)); %get average meteo

    end

     ydatalag=r1datalag(:,1);

    %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
end
 x1=xdatalag(:,2);
 x2=xdatalag(:,5);
 x3=xdatalag(:,1);
 y=ydatalag;
 
 x=xdatalag(:,[2,5]);
 %x=horzcat(xdatalag(:,[1,i+1]),xdatalag(:,1).*xdatalag(:,i+1));
 [b,~,stats1] = glmfit(x,ydatalag,'normal');
 
no_auto_y=y(2:end)-0.96.*y(1:end-1);
no_auto_x(:,1)=x(2:end,1)-x(1:end-1,1);
%no_auto_x(:,2)=x(2:end,2)-x(1:end-1,2);

no_auto_x(:,2)=x2(1:end-1);


[b2,~,stats2] = glmfit(no_auto_x, no_auto_y,'normal');




cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(4,2,7,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08);   
hold on;

ydata=[nan(20,1);y_GRACE(:,1)];
%ydata=aData.tropicalPP;

xlim_value(1)=min(ydata);
xlim_value(2)=max(ydata);
xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        %lag=lag_choice(tt);
        lag=20;
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        for i=1:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(i,1)=b(2,1);
            r1datalag(i,2)=bint(2,1)-b(2,1);
            r1datalag(i,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(i,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(i,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(i)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(i,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata, if one year is NAN, then the decadal average is NAN;
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag(:,1);
         ydatalag_sd=r1datalag(:,2);
       
         % coefficient of autocorrelation
         
         tmp_x=ydatalag(1:end-1); tmp_y=ydatalag(2:end);
         ind_nan = isnan(tmp_x) | isnan(tmp_y);
          auto_pcy = polyfit(tmp_x(~ind_nan),tmp_y(~ind_nan),1);
          
         tmp_x=xdatalag(1:end-1); tmp_y=xdatalag(2:end);
         ind_nan = isnan(tmp_x) | isnan(tmp_y);
          auto_pcx = polyfit(tmp_x(~ind_nan),tmp_y(~ind_nan),1);          

         
%         auto_pcy = polyfit(ydatalag(1:end-1),ydatalag(2:end),1);
%         auto_pcx = polyfit(xdatalag(1:end-1),xdatalag(2:end),1);
%          
         
         % added part for auto-correlation
        tmp_r=corrcoef(xdatalag(1:end-1),xdatalag(2:end));
        
             mdl = fitlm(yearlag(1:end-1),ydatalag(1:end-1)-auto_pcy(1).*ydatalag(2:end)); 
            [p1,DW1] = dwtest(mdl,'exact','both');
            DW1=round(DW1.*100)./100;
            %disp(p1);disp(DW1);
            
            mdl = fitlm(yearlag(1:end-1),xdatalag(1:end-1)-auto_pcx(1)*xdatalag(2:end));
            [p2,DW2] = dwtest(mdl,'exact','both');
            DW2=round(DW2.*100)./100;
            %disp(p2);disp(DW2);
        
        
        
            no_auto_y=ydatalag(1:end-1)-auto_pcy(1).*ydatalag(2:end);
            no_auto_x=xdatalag(1:end-1)-auto_pcx(1).*xdatalag(2:end);
            no_auto_ysd=ydatalag_sd(1:end-1);


            xdatalag=[];ydatalag=[];ydatalag_sd=[];
            [xdatalag,sortInd]=sort(no_auto_x);
            ydatalag=no_auto_y(sortInd);
            ydatalag_sd=no_auto_ysd(sortInd);


            xlim_value(1)=min(no_auto_x);
            xlim_value(2)=max(no_auto_x);
            xlim_value(3)=xlim_value(2)-xlim_value(1);
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;
            
            %%plot the uncertainty of y
            cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
       
            p=polyfit(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999),1);
            [r1,r2]=corrcoef(xdatalag(xdatalag>-9999),ydatalag(xdatalag>-9999));
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
           h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig < 0.01
             ptext='p<0.01';
         elseif sig < 0.05
             ptext='p<0.05';
         elseif sig < 0.1
             ptext='p<0.1';
         else
             ptext='p>0.1';
         end
            

        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end
        
        
        xfun2=strcat('DW_Y = ',num2str(DW1),{' '},'(p>0.1)'); 
        xfun3=strcat('DW_X = ',num2str(DW2),{' '},'(p=0.07)'); 
        
        
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T_{adj} (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     xlabel('$\sf\overline{GRACE-rec_{adj}}(cm)$','FontSize',10,'Interpreter','latex');
     text(xlim_value(1),-1.5,xfun);
     text(xlim_value(1),-2.2,xfun2);
     text(xlim_value(1),-2.8,xfun3);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([-3.1 2.2]);
     
     text(xlim_value(1)-0.05*xlim_value(3),1.9,strcat('(',char(96+(4-1)*2+1),')'));
     clear xlim_value
 






        
% the second col, interval is 8 years
%color_choice=[0.8,0,0;0,0.8,0;0,0,0.8];
color_choice=[0.8,0,0;0,0,0;0,0,0.8];
lag_choice=[10,8,30];

for fign=1:3
    
   %define the values of x
   ylim_value=nan(3,1);%xmin, xmax, xlengeth
    
   %define the variables used in axises     
    switch (fign)
        case 1 % x t-tmp      
        ydata2=aData.TTDT;
        ytitle2='$\overline\Delta\sf\overline{MAT} (^\circ C)$';
         
        
        case 2 % x t-pp       
        ydata2=aData.TPPDT;
        ytitle2='$\overline\Delta\sf\overline{MAP}(mm\ yr^{-1})$';
        
                
        case 3 % x t-swc        
        ydata2=aData.TSWCDT;
        ytitle2='$\overline\Delta \sf\overline{MASW} (\%)$';


    end
        
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    %%%%start ploting, long term meteo
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
 
    %plotting response curve 
    subaxis(4,2,(fign-1)*2+2,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08);   
    hold on;
    xlim_value(1)=min(ydata2);
    xlim_value(2)=max(ydata2);
    xlim_value(3)=xlim_value(2)-xlim_value(1);
 
    for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
        lag=lag_choice(tt);
        %1, in the time window, calculate sensitivity

        ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
        xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
        clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

        
        count=1;
        
        for i=1:lag:length(ydata_sen)-lag
            
            %sensitivity calculation
            [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
            r1datalag(count,1)=b(2,1);
            r1datalag(count,2)=bint(2,1)-b(2,1);
            r1datalag(count,3)=b(1,1);
            
            cd('/Volumes/RemiLBNL/project7/code');
            r1datalag(count,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
            
            [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
            sig1(count,1)=round(tmp2(2)*1000)./1000;
            
            yearlag(count)=yearGCB(i)+round(lag*0.5);
            temp_xdatalag(count,1)=mean(ydata2(i:i+lag)); %meteo was assigned to ydata
            
            
            if sig1>0.05
                
                r1datalag(count,1)=NaN;
            end
            
            count=count+1;
            
        end
        
         ydatalag=r1datalag(:,1);
         xdatalag=temp_xdatalag;
         ydatalag_sd=r1datalag(:,2);
        
        xlim_value(1)=min(xdatalag);
        xlim_value(2)=max(xdatalag);
        xlim_value(3)=xlim_value(2)-xlim_value(1);
        
 
        %sort x and y, in order to put in uncertainty
        [xdatalag,sortInd]=sort(xdatalag);
        ydatalag=ydatalag(sortInd);
        ydatalag_sd=ydatalag_sd(sortInd);
        
        
        mdl = fitlm(yearlag(1:end-1),ydatalag(1:end-1)-0.95.*ydatalag(2:end)); 
        [p1,DW1] = dwtest(mdl,'exact','both');
        %disp(p1);disp(DW1);

        mdl = fitlm(yearlag(1:end-1),xdatalag(1:end-1)-xdatalag(2:end));
        [p2,DW2] = dwtest(mdl,'exact','both');
        
        
        %plot the change of y to x
        faceColor=color_choice(tt,:);
        
        errorbar(xdatalag,ydatalag,ydatalag_sd, 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
        
        hold on;
        
        pxx2=plot(xdatalag,ydatalag,'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
    'MarkerEdgeColor',[0.2,0.2,0.2],...
    'MarkerFaceColor',[0.2,0.2,0.2]);
        hold on;


        
            p=polyfit(xdatalag,ydatalag,1);
            [r1,r2]=corrcoef(xdatalag,ydatalag);
            sig=round(r2(2)*1000)./1000;    
            rsq=round(r1(2)^2*100)./100; 
            
            
            h=refline(p(1),p(2));
            h.Color = faceColor;
            h.LineStyle = '--';
            h.LineWidth=1.5;
            
            p(1)=round(p(1).*100)./100;
            p(2)=round(p(2).*100)./100;
            
         if sig < 0.05
             ptext='p<0.05';
         elseif sig < 0.1
             ptext='p<0.1';
         else
             ptext='p>0.1';
         end
            
         disp(sig);
  
            
        %xfun=strcat('Y_{',num2str(lag),'}=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p=',num2str(sig)); 
        if p(2)>0
        xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        else
            xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),{', '},ptext); 
        end

%         xfun2=strcat('DW_X = ',num2str(DW1),'(p>0.1)'); 
%         xfun3=strcat('DW_Y = ',num2str(DW2),'(p>0.1)'); 
        
    end   
    
     ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
     %xlabel({ytitle;xfun},'FontSize',10);
     xlabel(ytitle2,'FontSize',10,'Interpreter','Latex');
     text(xlim_value(1)-0.05*xlim_value(3),0.5,xfun);
%      text(xlim_value(1),-0.3,xfun2);
%      text(xlim_value(1),-1.7,xfun3);
     
     xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     ylim([0 8]);
     
     text(xlim_value(1)-0.05*xlim_value(3),7.5,strcat('(',char(96+(fign-1)*2+2),')'));
     clear xlim_value
end




% cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
% subaxis(4,2,8,'SpacingVert',0.08,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.01,'MarginBottom',.08);   
% hold on;
% 
% ydata=[nan(20,1);y_GRACE(:,1)];
% %ydata=aData.tropicalPP;
% 
% xlim_value(1)=min(ydata);
% xlim_value(2)=max(ydata);
% xlim_value(3)=xlim_value(2)-xlim_value(1);
%  
%     for tt=2:2 % the window used for sensitivity analysis, 10, 20, 30 years
%         lag=8;
%         %1, in the time window, calculate sensitivity
% 
%         ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
%         xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
%         clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag r1datalag;
% 
%         for i=1:lag:length(ydata_sen)-lag
%             
%             %sensitivity calculation
%             [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
%             r1datalag(count,1)=b(2,1);
%             r1datalag(count,2)=bint(2,1)-b(2,1);
%             r1datalag(count,3)=b(1,1);
%             
%             cd('/Volumes/RemiLBNL/project7/code');
%             r1datalag(count,2)=2.*bootstrap_slope(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:),0.2,1000);
%             
%             [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
%             sig1(count,1)=round(tmp2(2)*1000)./1000;
%             
%             yearlag(count)=yearGCB(i)+round(lag*0.5);
%             temp_xdatalag(count,1)=mean(ydata(i:i+lag)); %meteo was assigned to ydata
%             
%             
%             if sig1>0.05
%                 
%                 r1datalag(count,1)=NaN;
%             end
%             
%             count=count+1;
%         end
%         
%          ydatalag=r1datalag(:,1);
%          xdatalag=temp_xdatalag(:,1);
%          ydatalag_sd=r1datalag(:,2);
%        
%         
%         xlim_value(1)=min(xdatalag);
%         xlim_value(2)=max(xdatalag);
%         xlim_value(3)=xlim_value(2)-xlim_value(1);
% 
%         %sort x and y, in order to put in uncertainty
%         [xdatalag,sortInd]=sort(xdatalag);
%         ydatalag=ydatalag(sortInd);
%         ydatalag_sd=ydatalag_sd(sortInd);
%         
%         zero_ind = ydatalag == 0 | isnan(xdatalag);
%         
%         
%         %plot the change of y to x
%         faceColor=color_choice(tt,:);
%         
%         errorbar(xdatalag(~zero_ind),ydatalag(~zero_ind),ydatalag_sd(~zero_ind), 'Color',faceColor,'LineStyle','none','LineWidth',0.5,'MarkerFaceColor',faceColor,'MarkerEdgeColor','k');
%         
%         hold on;
%         
%         pxx2=plot(xdatalag(~zero_ind),ydatalag(~zero_ind),'o','LineWidth',1,'color',faceColor,'MarkerSize',3,...
%     'MarkerEdgeColor',[0.2,0.2,0.2],...
%     'MarkerFaceColor',[0.2,0.2,0.2]);
%         hold on;
%             
%             %%plot the uncertainty of y
%             cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
% 
%         
%             
%             p=polyfit(xdatalag(~zero_ind),ydatalag(~zero_ind),1);
%             [r1,r2]=corrcoef(xdatalag(~zero_ind),ydatalag(~zero_ind));
%             sig=round(r2(2)*1000)./1000;    
%             rsq=round(r1(2)^2*100)./100; 
%             
%             
%             h=refline(p(1),p(2));
%             h.Color = faceColor;
%             h.LineStyle = '--';
%             h.LineWidth=1.5;
%             
%             p(1)=round(p(1).*100)./100;
%             p(2)=round(p(2).*100)./100;
%             
%         if p(2)>0
%         xfun=strcat('Y=',num2str(p(1)),'X+',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.01'); 
%         else
%             xfun=strcat('Y=',num2str(p(1)),'X',num2str(p(2)), ';  r^2=',num2str(rsq),', p<0.01'); 
%         end
% 
%         
%     end   
%     
%      ylabel(strcat(char(947),'_{CGR}','^T (PgC yr^{-1}°C^{-1})'),'FontSize',10);
%      %xlabel({ytitle;xfun},'FontSize',10);
%      xlabel('$\sf\overline{GRACE-rec}(cm)$','FontSize',10,'Interpreter','latex');
%      text(xlim_value(1)-0.05*xlim_value(3),0.5,xfun);
%      
%      xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
%      ylim([0 8]);
%      
%      text(xlim_value(1)-0.05*xlim_value(3),7.5,strcat('(',char(96+(4-1)*2+2),')'));
%      clear xlim_value

hold off;

set(f6,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project7/Figures/Figure4r','-djpeg','-r600');  













% not used
lag_choice=[8,20,25];
for tt=1:1 % the window used for sensitivity analysis, 10, 20, 30 years
    lag=lag_choice(tt);
    %1, in the time window, calculate sensitivity

    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    xdata_sen=[ones(size(ydata_sen)) aData.TTDT aData.TPPDT aData.TPDT];
    clear yearlag ydatalag xdatalag temp_ydatalag temp_xdatalag;

    count=1;
    
    for i=1:lag:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(count,1)=b(2,1);
        r1datalag(count,2)=bint(2,1)-b(2,1);
        r1datalag(count,3)=b(1,1);
        r1datalag(count,4)=b(3,1); %coefficient for TPP
        r1datalag(count,5)=b(4,1);% for TP

        [tmp1,tmp2]=corrcoef(xdata_sen(i:i+lag),ydata_sen(i:i+lag));
        sig1(count,1)=round(tmp2(2)*1000)./1000;
        
       if sig1>0.05
                
                r1datalag(count,1)=NaN;
       end

        yearlag(count)=yearGCB(i)+round(lag*0.5);
        xdatalag(count,1)=mean(aData.TTDT(i:i+lag)); %get average meteo
        xdatalag(count,2)=mean(aData.TPPDT(i:i+lag));
        xdatalag(count,3)=mean(aData.TSWCDT(i:i+lag));
        xdatalag(count,4)=mean(aData.TPDT(i:i+lag));
        %xdatalag(i,5)=yearlag(i)-1968;

         xdatalag(count,5)=mean(aData.tropicalT(i:i+lag)); %get average meteo
%         xdatalag(i,2)=mean(aData.tropicalPP(i:i+lag));
%         xdatalag(i,3)=mean(aData.tropicalSWC(i:i+lag));
%         xdatalag(i,4)=mean(aData.tropicalVPD(i:i+lag));

        count=count+1;
    end

     ydatalag=r1datalag(:,1);

    %clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
        
end
 x1=xdatalag(:,2);
 x2=xdatalag(:,5);
 x3=xdatalag(:,1);
 y=ydatalag;
 
 x=xdatalag(:,[2,5]);
 %x=horzcat(xdatalag(:,[1,i+1]),xdatalag(:,1).*xdatalag(:,i+1));
 [b,~,stats] = glmfit(x,ydatalag,'normal');        
       
        
 
 
 
 
 
 % to get water sensitivity
 
 for i = 1: 57
     
     aData.lagtropicalPP(i+1,1)=sum(mData.tropicalPP((i-1)*12+1+4:i*12+4,1));
     
 end
 
 aData.lagtropicalPP(1,1) =  aData.lagtropicalPP(2,1);
 
 f2=figure('Name','global average','Units', 'centimeters','Color','white', 'Position', [2, 2, 11, 10], ...
    'OuterPosition', [2, 2, 11, 10]);


% cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
% subaxis(1,2,1,'SpacingVert',0.08,'SpacingHoriz',0.08,'MR',0.03,'ML',0.08,'MarginTop',.03,'MarginBottom',.18); 
%color_choice=[0.6,0,0;0.3,0.3,0.3;0,0,1];
color_choice=[1,0.4,0.4;0,0,0;0.4,0.4,1];
lag_choice=[15,20,25];

for ttt=1:3
    clear r1datalag xdatalag ydatalag ydatalag_sd temp_xdatalag yearlag sortInd;
    lag=lag_choice(ttt);
    %lag=20;
    ydata_sen=aData.CO2DT; %if not detrend, what you see is actually due to fossil fuel and ocean sink
    %ydata_sen=sum(aData.CO2EMD(1:3,:),1)';
    xdata_sen=[ones(size(ydata_sen)) aData.lagtropicalPP];
    xlim_value=nan(3,1);

    for i=1:length(ydata_sen)-lag

        %sensitivity calculation
        [b,bint,~,~,~] = regress(ydata_sen(i:i+lag),xdata_sen(i:i+lag,:));
        r1datalag(i,1)=b(2,1);
        r1datalag(i,2)=bint(2,1)-b(2,1);
        r1datalag(i,3)=b(1,1);
        r1datalag(i,4)=-(b(2,1)).*nanmean(xdata_sen(i:i+lag,2),1);
        
        cd('/Volumes/RemiLBNL/project7/code');
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
        %cd('/Volumes/RemiLBNL\xl\code\functions');
        cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        H=shadedErrorBar(xdatalag,ydatalag,ydatalag_sd,{'-'},0);

        set(H.patch,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.3,'EdgeAlpha',0.3)
        set(H.mainLine,'color','none');
        set(H.edge,'color','none');
end


    
     ylabel(strcat('dCGR/','dMAT (PgC yr^{-1}°C^{-1})'),'FontSize',12);
     xlabel('Year','FontSize',12);
     set(gca,'box','off');
     
     %xlim([xlim_value(1)-0.1*xlim_value(3) xlim_value(2)+0.1*xlim_value(3)]);
     xlim([1965.5 2010.5]);
     ylim([1 8]);
     
     set(gca,'xtick', 1970:10:2010);
     
     %text(1967,7.3,'(a)','Fontsize',12);
     
 legend([pxx1(1),pxx1(2),pxx1(3)],{'15-yr window','20-yr window','25-yr window'},'Orientation','vertical','Position',[0.26    0.8    0.1    0.1],'Box','off');
set(f2,'PaperPositionMode','auto');
%print('/Volumes/RemiLBNL/project7/Figures/FigureS1','-djpeg','-r600'); 
 