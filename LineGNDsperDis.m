function [LN] = LineGNDsperDis(fo)
% this function calculates the distrubation of dislocations in any 
% selcted fatures (preferably slip band) and then proceed to cauclate 
% the contribuation of each dislcaotion type before proceeding to
% calculate the weighted average angle between the feature and the 
% the grain mean orienation
% the input for this folder should be the directory of an xEBSD_v3 file (fo)

% fo is directory of an xEBSD_v3 file
clc;close all
set(0,'defaultAxesFontSize',25);       set(0,'DefaultLineMarkerSize',14)
% fo = 'A:\OneDrive - Nexus365\Work\EBSD Data\20-12-11 DSS Slip\20_12_11_Slips_1_15nA_20kV_XEBSD';
load(fo,'GND','Data_InputMap','Map_RefID','GrainData',...
    'Grain_Map_A0_sample','Maps','Map_EBSD_MTEX');
%%
if ~exist('Spec','var') && length(GrainData.RefPoint.x)~=1
    s3=subplot(1,1,1);
    imagesc(Data_InputMap.X_axis,Data_InputMap.Y_axis,Map_RefID); hold on
    axis off;           axis image;         axis xy;        s3.YDir='reverse';
    colormap jet;       s3.XDir='reverse';
    title('Respond in the Command Line')
    if isempty(GrainData.RefPoint)
        [GrainData.RefPoint] = to_label(Data_InputMap, Map_RefID);
    else
        for i=1:length(GrainData.RefPoint.x)
            GrainData.RefPoint.prop.labels{i} = num2str(i);
        end
    end
    scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'k','filled');
    scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'w');
    labelpoints(GrainData.RefPoint.x,GrainData.RefPoint.y,...
        GrainData.RefPoint.prop.labels,'FontSize',20);
    hold off;      set(gcf,'position',[30 50 1300 950])
    Spec    = input('Which Grain you want to explore?   '); close
elseif ~exist('Spec','var')
    Spec = 1;
end
% fo = erase(fo,'XEBSD');

%% cut
LN.onI = num2str(randi(100)); LN.fo=fileparts(fo);
LN.fo = fullfile(LN.fo,LN.onI); mkdir(LN.fo)
A = squeeze(Grain_Map_A0_sample(:,:,Spec,1,1)); A(A==0)=NaN;
GND.total(isnan(A))=NaN;
f= figure; s3=subplot(1,1,1);
imagesc(Data_InputMap.X_axis,Data_InputMap.Y_axis,real(log10(GND.total)))%,'LineStyle','none');
brighten(0.2);
c = colorbar;        c.Label.String = ['GNDs [log_{10} m^{-2}]'];
axis image;         axis xy;        s3.YDir='reverse';
colormap jet;       s3.XDir='reverse';
title('Select Feature'); %caxis([11 13.8])
set(gcf,'position',[30 50 1300 950])
[LN.x,LN.y,LN.c] = improfile; axis off;
addScale([1 1 1],[Data_InputMap.XSample(:) Data_InputMap.YSample(:)]);
hold on;
plot(LN.x,LN.y,'-.k','LineWidth',3,'DisplayName','Data line')
hold off; title ''
DirSave = fullfile(LN.fo, [LN.onI '_Lined Data.fig']);      saveas(gcf,DirSave);
DirSave = fullfile(LN.fo, [LN.onI '_Lined Data.tif']);      saveas(gcf,DirSave);     close

%%
LN.MapsGND_screw = zeros(size(GND.total));
LN.MapsGND_edge  = zeros(size(GND.total));
for iV=1:length(GND.Sliplabels{Spec})
    LN.MapsGNDper(:,:,iV)  = abs(reshape(GND.rho(:,iV),length(Data_InputMap.Y_axis),...
        length(Data_InputMap.X_axis)))./GND.total;
    LN.MapsGNDreal(:,:,iV) = reshape(GND.rho(:,iV),length(Data_InputMap.Y_axis),...
        length(Data_InputMap.X_axis));
    
    if strfind(GND.Sliplabels{Spec}{iV}, 'screw')
        LN.MapsGND_screw = LN.MapsGND_screw+reshape(GND.rho(:,iV),length(Data_InputMap.Y_axis),...
            length(Data_InputMap.X_axis));
    elseif strfind(GND.Sliplabels{Spec}{iV}, 'edge')
        LN.MapsGND_edge = LN.MapsGND_edge+reshape(GND.rho(:,iV),length(Data_InputMap.Y_axis),...
            length(Data_InputMap.X_axis));
    else
        dips('What the f**k is wrong with you!')
    end
end
save([LN.fo '\' LN.onI '_Line_Data.mat'],'LN');

%% get data on the line
for i=1:length(LN.x)
    [~,ix(i)]= min(abs(Data_InputMap.X_axis-LN.x(i)));
    [~,iy(i)]= min(abs(Data_InputMap.Y_axis-LN.y(i)));
    LN.Tot(i) = norm([[LN.x(1),LN.y(1)]-[LN.x(i),LN.y(i)]]);
    
    LN.GNDs(i)      = GND.total(iy(i),ix(i));
    LN.GND_screw(i) = LN.MapsGND_screw(iy(i),ix(i));
    LN.GND_edge(i)  = LN.MapsGND_edge(iy(i),ix(i));
    for iV=1:length(GND.Sliplabels{Spec})
        LN.GNDreal(i,iV) = LN.MapsGNDreal(iy(i),ix(i),iV);
    end
    LN.PH(i)  = Maps.PH_2(iy(i),ix(i));
    LN.MAE(i) = Maps.MAE_2(iy(i),ix(i));
end

for iV=1:length(GND.Sliplabels{Spec})
    LN.GND_per(iV)  = sum(abs(LN.GNDreal(:,iV)),'all','omitnan')./...
        sum(abs(LN.GNDreal),'all')*100;
    LN.GND_real(iV) = sum(LN.GNDreal(:,iV),'all','omitnan');
end
LN.DisName = GND.Sliplabels{Spec};
save([LN.fo '\' LN.onI '_Line_Data.mat'],'LN');

%%
bar(1:iV,LN.GND_per,'k');
xticks([1:iV]);      xlabel('Dislocation System'); xlim([0 iV])
ylabel('% of \rho_{GND}');
ylim([0 max([0:5:max(LN.GND_per)+5])]); yticks([0:5:max(LN.GND_per)+5])
set(gcf,'position',[41 195 1800 600]); box off;grid on
set(gca,'GridLineStyle','-.')
saveas(gcf, [LN.fo '\' LN.onI '_GND_Per.fig']);
saveas(gcf, [LN.fo '\' LN.onI '_GND_Per.tiff']);close

%%
bar(1:iV,LN.GND_real,'k');
xticks([1:iV]);      xlabel('Dislocation System'); xlim([0 iV])
ylabel('\rho_{GNDs} [m^{-2}]');
set(gcf,'position',[41 195 1800 600]); box off;grid on
set(gca,'GridLineStyle','-.')
saveas(gcf, [LN.fo '\' LN.onI '_GND_Mag.fig']);
saveas(gcf, [LN.fo '\' LN.onI '_GND_Mag.tiff']);close

%%
for iO = 1:length(LN.DisName)
    if strfind(LN.DisName{iO}, 'screw')
        W = erase(LN.DisName{iO},'screw');
    elseif strfind(LN.DisName{iO}, 'edge')
        W = erase(LN.DisName{iO},'edge');
    end
    try
        O = erase(W,'b');
        O = strrep(O,'t',',');          O = erase(O,'[');
        O = erase(O,']');               O = {erase(O,' ')};
        Od  = cell2mat(cellfun(@str2num,O,'UniformOutput',false));
        b(iO,1:3) = Od(1:3);
    catch
        O = erase(W,'b');
        O = strrep(O,'t',' ');          O = erase(O,'[');
        O = erase(O,']');               O = {erase(O,'t')};
        Od  = cell2mat(cellfun(@str2num,O,'UniformOutput',false));
        b(iO,1:3) = Od(1:3);
    end
    t(iO,1:3) = Od(4:end);
    clear Od O
end
LN.b = b; LN.t = t; 
save([LN.fo '\' LN.onI '_Line_Data.mat'],'LN');

%%
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

ebsd2 = Map_EBSD_MTEX('indexed'); % removing non indexed points
%constructing grains
[grains,ebsd2.grainId,ebsd2.mis2mean]=calcGrains(ebsd2,'angle',5*degree);
grains2=grains(grains.grainSize > 20); % removing noise
% cs = grains2(Spec).CS; % asign crystal system
LN.ori = grains2(Spec).meanOrientation; % grain mean oreination
% r=vector3d.Z; % normal to surface
% h=inv(ori)*r;
% u = h.round;
% u = [u.h; u.k; u.l];
dir = [0 0 1];
for iO=1:1:length(LN.DisName)
    u = LN.ori.matrix*[1 0 0;0 1 0; 0 0 1]*LN.b(iO,1:3)'; 
    
    % two diretcion .. the angle between tow directions
    %     LN.ang(iO) = round(rad2deg(atan2(norm(cross(u,b(iO,1:3))),...
    %         dot(u,b(iO,1:3)))),0);
    LN.ang(iO) = round(rad2deg(atan2(norm(cross(u',dir)),...
        dot(u',dir))),1);
    if LN.ang(iO) > 90
        LN.ang(iO) = LN.ang(iO)-90;
    end
end
LN.Table = table(LN.b,LN.t,LN.GND_real',LN.GND_per',LN.ang',...
    'VariableNames',{'b','t','m-2','%','deg'});
save([LN.fo '\' LN.onI '_Line_Data.mat'],'LN');

%%
A = unique(LN.ang);
Cat = zeros(size(A));
GP = zeros(length(LN.DisName),length(A));
for iO=1:length(LN.DisName)
    [~,B] = ismember(LN.ang(iO),A);
    Cat(B) = Cat(B)+LN.GND_per(iO);
    GP(iO,B) = LN.GND_per(iO);
end

C = find(sum(GP,2)==0);
B = find(sum(GP,2)~=0);
GP(C,:)=[];

%
% fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 1]]);yyaxis left;
% plot(LN.Tot,real(log10(sum(abs(LN.GNDreal),2))),'-k','LineWidth',3,'DisplayName','Total')
subplot(1,2,1)
plot(LN.Tot,LN.GNDs,'-k','LineWidth',3,'DisplayName','Mag.')
hold on;
plot(LN.Tot,LN.GND_screw,'-r','LineWidth',3,'DisplayName','Screw')
plot(LN.Tot,LN.GND_edge,'-g','LineWidth',3,'DisplayName','Edge')
% plot(LN.Tot,real(log10(LN.GND_screw)),'-r','LineWidth',3,'DisplayName','Screw')
% plot(LN.Tot,real(log10(LN.GND_edge)),'-g','LineWidth',3,'DisplayName','Edge')
hold off; ylabel('GNDs [m^{-2}]');
% yyaxis right;
% plot(LN.Tot,LN.PH,'-g','LineWidth',3,'DisplayName','PH');ylabel('PH');
% plot(LN.Tot,LN.MAE,'-b','LineWidth',3,'DisplayName','MAE');ylabel('MAE');

xlabel(['Total Distance [\mum]']);axis tight
[~,objh]=legend('location','northoutside','box','off','fontsize',28,...
    'Orientation','horizontal'); box off
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 12); %// set marker size as desired

subplot(1,2,2)
H=bar(1:size(GP,2),GP,'stacked');
cc = jet(size(GP,1));
for i = 1:size(GP,1)
    H(i).FaceColor = 'flat';
    H(i).CData = cc(i,:);
end
xlabel('\Phi^{\circ}'); xlim([0 size(GP,2)+1]);
for iO = 1:size(GP,2)
    OO{iO}=num2str(A(iO));
end
xticklabels(OO)
ylabel('% of \rho_{GND}');
legend(string(B),'location','best','NumColumns',2)
box off;grid on; set(gca,'GridLineStyle','-.')
title(['$\bar{\Phi} = ' num2str(round(sum(A.*Cat)/sum(Cat),1)) ...
    '^{\circ}$'],'Interpreter','Latex')
% saveas(gcf, [LN.fo '\' LN.onI '_GND_Bar.fig']);
% saveas(gcf, [LN.fo '\' LN.onI '_GND_Bar.tif']);close
set(gcf,'WindowStyle','normal'); set(gcf,'position',[10 45 1700 700]);  box off
saveas(gcf, [LN.fo '\' LN.onI '_GND_line.fig']);
saveas(gcf, [LN.fo '\' LN.onI '_GND_line.tiff']);close

LN.AvgAng = round(sum(A.*Cat)/sum(Cat),0);
save([LN.fo '\' LN.onI '_Line_Data.mat'],'LN');

end
