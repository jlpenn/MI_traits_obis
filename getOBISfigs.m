clear all;close all;

% switches
do.Eoeff = 1; % if dEdT ref Eo to median T
do.F1min = 0.5; % minimum F1 for analyses 
do.qc=1; % filter lab trait?

% paths to cbrewer
addpath(genpath('/Users/jpenn/Desktop/Research/FRES/Cambridge_prisms/Cambridge'))

% laboratory traits from Deutsch et al., (Nature) 2020
load('Phy.mat')

% constants
phic_mn = 3.0; % average phi-crit
dEdT=0.025; % dEo/dT (lab)
Tref = 15; % reference T
kb = 8.6173324*10^-5; % Boltzmann's constant (eV/K)
ss.temp=-3:1:35;% temp
ss.texp = (1./(ss.temp+273.15)-1./(Tref+273.15)).*(1/kb); %Arrhenius temp

% extact lab data
for ii = 1:366
    Par.tmedP(ii) = median(Phy.P.temp(ii,:),'omitnan');% median experimental temperature
end
Par.tmedP=Par.tmedP';
Phy.Eoref=Phy.Eo-dEdT*(Par.tmedP-Tref);%  Eo @ Tref
Phy.Eo(~Phy.Ieo)=nan; % Eo (resting)
Phy.Emet(~Phy.Iem)=nan;
Phy.alphaD(~Phy.Iem)=nan;
aD = Phy.alphaS./Phy.Ao;% species w/alphaD and Pcrit 
idx = ~isnan(aD) & isnan(Phy.alphaD);
Phy.alphaD(idx)=aD(idx);% fill missing alphaD
Phy.Ao(~Phy.Ieo)=nan; % resting hypoxia tolerance (1/atm)
Phy.Ac(~Phy.Ieo)=nan; % active hypoxia tolerance (1/atm)

% Remove non-marine species
% 'Carassius carassius' Crucian carp, freshwater fish
% 'Oreochromis niloticus' Nile fish, largely freshwater
% 'Palaemon pugio' grass shrimp, brackish/fresh 

if do.qc
    %quality control
    Sp_cases={'Carassius carassius','Oreochromis niloticus','Palaemonetes pugio'};
    for ii=1:length(Sp_cases)
        Phy.Eo(strcmp(Phy.Species_unq,Sp_cases(ii)))=NaN;
        Phy.Ao(strcmp(Phy.Species_unq,Sp_cases(ii)))=NaN;
        Phy.Ac(strcmp(Phy.Species_unq,Sp_cases(ii)))=NaN;
        Phy.phic(strcmp(Phy.Species_unq,Sp_cases(ii)))=NaN;
        Phy.alphaD(strcmp(Phy.Species_unq,Sp_cases(ii)))=NaN;

    end

end

% update lab data
Phy.Eo(strcmp(Phy.Species_unq,"Cancer irroratus"))=0.88;% includes all lifestages
Phy.Ao(strcmp(Phy.Species_unq,"Cancer irroratus"))=30.8;% includes all lifestages


% add geographic data for lab species 
load('obis_AE_og.mat')
Phy.Zmed = obis_AE.Zmed;
Phy.Ymed = obis_AE.Ymed;

% load obis traits
load('obis_traits_fulldata.mat')

% load grid
load Goc.mat
prct=Goc.prct;
plim=5;% percentile

% filter traits based on F1
obis.Eo(obis.F1<do.F1min)=NaN; % Eeco (eV)
obis.Ac(obis.F1<do.F1min)=NaN; % Aeco (1/atm)

% Median obis temperature
obis.Tmed = obis.Tprct(:,find(Goc.prct==50));% median

% effective Eo at median T
obis.Eref=obis.Eo; % Eeco(Tref)
if do.Eoeff
    obis.dEdT(isnan(obis.dEdT))=0;% no T-dependence of Eeco for these species
    obis.Eoeff = obis.Eo + obis.dEdT.*(obis.Tmed-15);% Eeco(Tmed)
    obis.Eo=obis.Eoeff;% change obis.Eo(Tref) to obis.Eo(Tmed)
    Phy.Eoref=Phy.Eo; % Eo(Tmed)
end

% occurrences in T-pO2
obis.TObin=double(obis.TObin);
obis.TObin(obis.TObin==0)=nan;

% occurrences in XYZ
obis.XYZbin = double(obis.XYZbin);

% land mask for mapping
lat1 = -87.5:5:87.5;
lon1 = 2.5:5:357.5;

% land mask
load('etopo.mat')

% rearrange longitude
topo2 = zeros(1081,540);
topo2(1:61,:) = topo(1021:end,:);
topo2(62:end,:) = topo(1:1020,:);
land=topo2*nan;
land(topo2>0)=1;
lon3 = zeros(1081,1);
lon3(1:61)=lon(1021:end)-360;
lon3(62:end)=lon(1:1020);

%% Fig. 1A: Species richness (with trait data)
figure

% lat-lon occurrences
obis.XYbin = (max(obis.XYZbin(:,:,:,:),[],4))>0;
Rxy = sum(obis.XYbin,1,'omitnan'); 
Rxy(Rxy==0)=nan;

%%Create two axes
f = figure;
f.Position = [250 250 640*1.5 400*1.25];
ax1 = axes;
pcolor(ax1,lon1,lat1,squeeze(log10(Rxy))')
shading flat
view(2)
ax2 = axes;
pcolor(ax2,lon3,lat,land')
shading flat
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'jet')
colormap(ax2,'gray')
%ylim([-80 87])
ylim(ax1,[-80 87])
ylim(ax2,[-80 87])
xlim(ax1,[1 359])
xlim(ax2,[2 358])
box on
caxis(ax1,[0 3.7])
set(ax1,'fontsize',18)
xticks(ax1,0:60:360)
xticks(ax2,0:60:360)
xticklabels(ax1,{'0^oE','60^oE','120^oE','180^oE','120^oW','60^oW'})
yticks(ax1,[-80 -60 -30 0 30 60 87])
yticks(ax2,[-80 -60 -30 0 30 60 87])
yticklabels(ax1,{'80^oS','60^oS','30^oS','0^o','30^oN','60^oN','90^oN'})
xticklabels(ax1,{'0^oE','60^oE','120^oE','180^oE','120^oW','60^oW'})
cb1 = colorbar(ax1,'Position',[.92 .11 .03 .815]);
hold on

% % load Fig 1C,D species (data downloaded from https://obis.org/)
load('c_striata/obis_c_striata.mat')
dat.lon(dat.lon<0)= dat.lon(dat.lon<0)+360;
% remove 3 land occurrences for plotting
dat.lon(dat.lon<150)=nan;
dat.lon(dat.lon<277 & dat.lat>35)=nan;
scatter(dat.lon,dat.lat,15,[.1 .1 .6]*0,'marker','x')
load('p_antarctica/obis_p_antarctica.mat')
dat.lon(dat.lon<0)= dat.lon(dat.lon<0)+360;
scatter(dat.lon,dat.lat,15,[.1 .6 .1]*0,'marker','o')
load('p_diacanthus/obis_p_diacanthus.mat')
% remove 3 land occurrences for plotting
dat.lon(dat.lon<0)= dat.lon(dat.lon<0)+360;
dat.lon(dat.lon<50 & dat.lat>-.5 & dat.lat<10.5)=nan;
scatter(dat.lon,dat.lat,15,[.6 .1 .1]*0,'marker','^')
cd figs
saveas(gcf,'Fig1A','pdf')
cd ..

%% Fig S2: LAT-LAT 
Yedge = Goc.Yedge;
Ycent = Goc.Ycent;
Ieo = ~isnan(Phy.Eo);

%obis
% southern range
YminO=obis.Yprct(:,find(prct==plim));
% northern range
YmaxO=obis.Yprct(:,find(prct==(100-plim)));
%lab
% southern range
YminL=obis_AE.Yprct(Ieo,find(prct==plim));
%northern range
YmaxL=obis_AE.Yprct(Ieo,find(prct==(100-plim)));
psym='-:::::::';
dy=[0 10 20 40 60 80];
set(groot,'DefaultAxesFontSize', 20)

figure; clf reset;
N=histcounts2(YminO,YmaxO,Yedge,Yedge); N(N==0)=nan;
pcolor(Ycent,Ycent,log10(N)'); shading flat;colorbar
hold on
plot(YminL,YmaxL,'ro','MarkerSize',10,'LineWidth',3)
idx=find(strcmp(obis.species,'Centropristis striata'))
scatter(YminO(idx),YmaxO(idx),200,[0 .6 1],'marker','x','linewidth',4)
idx=find(strcmp(obis.species,'Paraeuchaeta antarctica'))
scatter(YminO(idx),YmaxO(idx),200,[0 .8 0],'marker','o','linewidth',3)
idx=find(strcmp(obis.species,'Pygoplites diacanthus'))
scatter(YminO(idx),YmaxO(idx),200,[.9 .5 0],'marker','^','linewidth',3)

hold on; for i=1:length(dy); plot(Ycent,Ycent+dy(i),psym(i)); end
plot([0 0],[-80 80],'k:',[-80 80],[0 0],'k:')
xlabel('Southern Range Limit (^oN)'); ylabel('Northern Range Limit (^oN)')
colormap(cbrewer('seq','Blues',20))
legend('Biogeographic','Experimental','C. striata','P. antarctica','P. diacanthus','location','southeast')
xlim([-77 87.5])
ylim([-77 87.5])
caxis([0 2.5])

cd figs
saveas(gcf,'FigS2','pdf')
cd ..
%% Fig. 1B,C Metabolic rates and Pcrit vs. T

% species traits
Es = 0.3;
Eo = 0.4;
Ephicrit = -0.35;
Ermr=Eo+Es;
Emmr=Eo+Es+Ephicrit;
Eeco = Eo+Ephicrit;
Ao = 20;
aD = 150/Ao;
phicrit = 3;

% Pcrit
pcrit_rest = (1./Ao).*exp(-Eo.*ss.texp);
pcrit_act = (phicrit./Ao).*exp(-Eeco.*ss.texp);
pcrit_act_const = (phicrit./Ao).*exp(-Eo.*ss.texp);

% Metabolism
rmr = aD.*exp(-Ermr.*ss.texp);
mmr = phicrit*aD.*exp(-Emmr.*ss.texp);
mmr_const = phicrit*aD.*exp(-Ermr.*ss.texp);

figure
plot(ss.temp,pcrit_rest,'r','linewidth',4)
hold on
plot(ss.temp,pcrit_act_const,'color',[.7 .7 1],'linewidth',4)
plot(ss.temp,pcrit_act,'b','linewidth',4)
set(gca,'fontsize',18)
xlabel('Temperature (^oC)')
ylabel('pO_2 (atm)')
ylim([0 0.25])
xlim([-3 30])
cd figs
saveas(gcf,'Fig1C','pdf')
cd ..


figure
plot(ss.temp,rmr./rmr(1),'r','linewidth',4)
hold on
plot(ss.temp,mmr_const./rmr(1),'color',[.7 .7 1],'linewidth',4)
plot(ss.temp,mmr./rmr(1),'b','linewidth',4)
set(gca,'fontsize',18)
xlabel('Temperature (^oC)')
ylabel('Metabolic rate (normalized)')
xlim([-3 30])
cd figs
saveas(gcf,'Fig1B','pdf')
cd ..

%% Fig. 1D,E, Tropical and polar endemics
f = figure;
%%Create two axes
% polar species
idx0 = find(strcmp(obis.species,'Paraeuchaeta antarctica'));
for ii =1:length(idx0)
    idx=idx0(ii);

    %figure
    ax1 = axes;
    pcolor(ax1,Goc.Tcent,Goc.Ocent,log10(squeeze(obis.TObin(idx,:,:)))')
    hold on
    shading flat;
    colormap jet

    % best fit Pcrit (OBIS)
    this_dEdT=obis.dEdT(idx);
    if isnan(this_dEdT)
        this_dEdT=0;
    end
    Eo = obis.Eref(idx)+this_dEdT*(ss.temp-Tref);
    ss.pcrit = exp(-Eo.*ss.texp)./obis.Ac(idx);
    hold on
    plot(ss.temp,ss.pcrit,'color',[0.2 0.6 0.2],'linewidth',3)
    ylim([0 0.25])

end
view(2)
ax2 = axes;

%tropical indo-pacific species
idx0 = find(strcmp(obis.species,'Pygoplites diacanthus'));
for ii =1:length(idx0)
    idx=idx0(ii);

    %figure
    pcolor(ax2,Goc.Tcent,Goc.Ocent,log10(squeeze(obis.TObin(idx,:,:)))')
    hold on
    shading flat;
    colormap jet

    % best fit Pcrit (OBIS)
    this_dEdT=obis.dEdT(idx);
    if isnan(this_dEdT)
        this_dEdT=0;
    end
    Eo = obis.Eref(idx)+this_dEdT*(ss.temp-Tref);
    ss.pcrit = exp(-Eo.*ss.texp)./obis.Ac(idx);
    hold on
    plot(ss.temp,ss.pcrit,'color',[0.6 0.2 0.2],'linewidth',3)
    ylim([0 0.25])

end

linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
%%Give each one its own colormap
colormap(ax1,cbrewer('seq','BuGn',10))
colormap(ax2,cbrewer('seq','OrRd',10))
caxis(ax1,[0 3.5])
caxis(ax2,[0 3.5])
set(ax1,'fontsize',18)
%cb1 = colorbar(ax1,'Position',[.92 .11 .03 .815]);
%cb1 = colorbar(ax2,'Position',[.914 .11 .03 .815]);
xlabel(ax1,'Temperature (^oC)')
ylabel(ax1,'pO_2 (atm)')

cd figs
saveas(gcf,'Fig1E','pdf')
cd ..

%% Fig: 2A,B: Map of Ac
% Set up trait map 
% all depth 
obis.XYbin = max(obis.XYZbin(:,:,:,:),[],4);

% mask
obis.XYbin(obis.XYbin<1)=nan;% no occurrence
obis.XYbin(obis.XYbin>0)=1;% occurrence

% open space
obis.XYeo = zeros(size(obis.XYZbin,1),72,36);
obis.XYac = zeros(size(obis.XYZbin,1),72,36);

% populate species trait wherever occurrence
for ii = 1:size(obis.XYZbin,1)
    obis.XYeo(ii,:,:) = obis.XYbin(ii,:,:)*obis.Eo(ii);
    obis.XYac(ii,:,:) = obis.XYbin(ii,:,:)*obis.Ac(ii);
end

% inter-species median
obis.Eomap = squeeze(median(obis.XYeo,1,'omitnan'));
obis.Eomap2 = squeeze(nanmedian(obis.XYeo,1));
obis.Acmap = squeeze(median(obis.XYac,1,'omitnan'));

% Trait map - Ac
%%Create two axes
f = figure;
f.Position = [250 250 640*1.25 400*1.25];
ax1 = axes;
pcolor(ax1,lon1,lat1,(obis.Acmap)')
shading flat

view(2)
ax2 = axes;
pcolor(ax2,lon3,lat,land')
shading flat
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'jet')
colormap(ax2,'gray')
ylim(ax1,[-80 87])
ylim(ax2,[-80 87])
box on
caxis(ax1,[0 25])
set(ax1,'fontsize',18)
xticks(ax1,0:60:360)
xticks(ax2,0:60:360)
xticklabels(ax1,{'0^oE','60^oE','120^oE','180^oE','120^oW','60^oW'})
xlabel(ax1,'Longitude')
ylabel(ax1,'Latitude (^oN)')
cb1 = colorbar(ax1,'Position',[.92 .11 .03 .815]);
xlim([2 357])

cd figs
saveas(gcf,'Fig2A','eps')
cd ..

% Trait map - Eo
%%Create two axes
f = figure;
f.Position = [250 250 640*1.25 400*1.25];
ax1 = axes;
pcolor(ax1,lon1,lat1,(obis.Eomap)')
shading flat
view(2)
ax2 = axes;
pcolor(ax2,lon3,lat,land')
shading flat
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'jet')
colormap(ax2,'gray')
%ylim([-80 87])
ylim(ax1,[-80 87])
ylim(ax2,[-80 87])
box on
caxis(ax1,[-1 1])
set(ax1,'fontsize',18)
xticks(ax1,0:60:360)
xticks(ax2,0:60:360)
xticklabels(ax1,{'0^oE','60^oE','120^oE','180^oE','120^oW','60^oW'})
xlabel(ax1,'Longitude')
ylabel(ax1,'Latitude (^oN)')
cb1 = colorbar(ax1,'Position',[.92 .11 .03 .815]);
xlim([2 357])

cd figs
saveas(gcf,'Fig2B','eps')
cd ..

%% Histograms of Ac-Eo vs. lat 

% Traits vs. latitude
Ylat = -90:5:90;
Ycent = -87.5:5:87.5;
E = -2:0.2:2;
A = 0:1:35;
Ecent = -2+0.2/2:0.2:2-0.2/2;
Acent = 0+0.5:1:35-0.5;
[Ney]=histcounts2(obis.Ymed,obis.Eo,Ylat,E);
[Nay]=histcounts2(obis.Ymed,obis.Ac,Ylat,A);

% Eo vs lat
for yy =1:37-1
    lidx = obis.Ymed > Ylat(yy) & obis.Ymed<= Ylat(yy+1);
    this_E = obis.Eo(lidx);
    Emy(yy) =median(this_E(:),'omitnan');  
end

% Ac vs lat
for yy = 1:37-1
    lidx = obis.Ymed > Ylat(yy) & obis.Ymed<= Ylat(yy+1);
    this_E = obis.Ac(lidx);
    Amy(yy) =median(this_E(:),'omitnan');
end

%% Fig 2C: Ac vs. latitude
figure
Nmsk = Nay;
Nmsk(Nay==0)=nan;
pcolor(Ycent,Acent,log10(Nmsk'))
hold on
scatter(Phy.Ymed,Phy.Ao/phic_mn,75,'r','linewidth',3)

% contributions from resting and active metabolism
Arest = median(Phy.alphaS,'omitnan')./(Phy.alphaD.*median(Phy.phic,'omitnan'));
Aact = median(Phy.alphaS,'omitnan')./(Phy.alphaD.*Phy.phic);

% variations from resting metabolism
scatter(Phy.Ymed,Arest,75,[.6 .6 .6]*0,'linewidth',1)
% variations from active metabolism
scatter(Phy.Ymed,Aact,75,[.6 .6 .6]*0,'linewidth',1,'marker','^')

% measured Ao/mean phicrit
scatter(Phy.Ymed,Phy.Ao/phic_mn,75,'r','linewidth',3)

% central tendencies
plot(Ycent(2:end),Amy(2:end),'color',[0 1 0],'linewidth',2)

colormap(cbrewer('seq','Blues',20))
set(gca,'fontsize',16)
ylabel('Hypoxia Tolerance, A (atm^-^1)')
legend('Biogeographic (A_e_c_o)','Experimental (A_o)','Resting demand','Active demand','location','northwest')
xlabel('Latitude (^oN)')
caxis([0 2.5])
shading flat
xlim([-78 87])
set(gca,'fontsize',16)
box on
colorbar 
ylim([0 37])

cd figs
saveas(gcf,'Fig2C','pdf')
cd ..

% tropical vs. extra-tropical resting traits
nanmedian(Phy.Ao(abs(Phy.Ymed)<20))
nanmedian(Phy.Ao(abs(Phy.Ymed)>20))


%% Fig 2D Eo vs. latitude
figure
Nmsk = Ney;
Nmsk(Ney==0)=nan;
pcolor(Ycent,Ecent,log10(Nmsk'))
hold on
scatter(Phy.Ymed,Phy.Eoref,75,'r','linewidth',3)
scatter(Phy.Ymed,Phy.Emet,75,[0 0 0],'linewidth',1)
scatter(Phy.Ymed2,Phy.Emet2,75,[0 0 0],'linewidth',1)% additional Emet data
scatter(Phy.Ymed,Phy.Eoref,75,'r','linewidth',3)
plot(Ycent(2:end),Emy(2:end),'color',[0 1 0],'linewidth',2)
colormap(cbrewer('seq','Blues',20))
set(gca,'fontsize',16)
ylabel('Temperature Sensitivity, E (eV)')
legend('Biogeographic (E_e_c_o)','Experimental (E_o)','location','southwest')
xlabel('Latitude (^oN)')
caxis([0 2.5])
shading flat
xlim([-78 87])
ylim([-1.95 1.95])
set(gca,'fontsize',16)
box on
colorbar 

cd figs
saveas(gcf,'Fig2D','pdf')
cd ..

% tropical vs. extra-tropical resting traits
nanmedian(Phy.Eo(abs(Phy.Ymed)<20))
nanmedian(Phy.Eo(abs(Phy.Ymed)>20))

%%  Fig3A:Trait PDF Eo-Ac 
% Eo vs Ac
[N,xe,ye]=histcounts2(obis.Eo,obis.Ac,E,A);

figure
Nmsk = N;
Nmsk(N==0)=nan;
pcolor(Ecent,Acent,log10(Nmsk'))
hold on
shading flat
scatter(Phy.Eoref,Phy.Ao/phic_mn,100,'r','linewidth',2)
colormap(cbrewer('seq','Blues',20))
set(gca,'fontsize',16)
xlabel('Temperature Sensitivity, E (eV)')
ylabel('Hypoxia Tolerance, A (atm^-^1)')
caxis([0 3])
xlim([-1.9 1.9])
ylim([0 35])
set(gca,'fontsize',18)
box on
colorbar 
text(-1.8,26,'Num. Species (log_1_0)','fontsize',16)
legend('Biogeographic (A_e_c_o, E_e_c_o)','Experimental (A_o, E_o)','location','northwest')

cd figs
saveas(gcf,'Fig3A','pdf')
cd ..

%%  Fig3B: Traits and aerobic habitat volume 
if ~do.Eoeff % !!!!!!!!! WARNING: only run with do.Eoeff = 0 (i.e., with Eo(Tref))
figure
% upper 100 m habitat 
load('modern_Hvol_phimax_100m_seasonal.mat')
E = -2+0.1/2:0.1:2-0.1/2;
A = 0+0.5:1:35-0.5;
contourf(E(1:end),A(1:end),V')
hold on
Eo2Ec = 0.07; % median difference between lab Eo(Tref) and obis Eeco(Tref)
scatter(Phy.Eoref(Phy.Zmed<100)-Eo2Ec,Phy.Ac(Phy.Zmed<100),100,'r','linewidth',2)
colormap(cbrewer('seq','Blues',20))
set(gca,'fontsize',16)
xlabel('Temperature Sensitivity, E_e_c_o (eV)')
ylabel('Hypoxia Tolerance, A_e_c_o (atm^-^1)')
xlim([-1.9 1.9])
ylim([0 30])
set(gca,'fontsize',18)
box on
colorbar 
text(-1.8,28,'Habitat Volume (m^3)','fontsize',16)
text(-1.8,26,'(0-100m)','fontsize',16)

cd figs
saveas(gcf,'Fig3B','pdf')
cd ..

% Fig 3C
% upper 1000 m habitat 
load('modern_Hvol_phimax_100_1000m_seasonal.mat')
figure
contourf(E(1:end),A(1:end),V')
hold on
scatter(Phy.Eoref(Phy.Zmed>100 & Phy.Zmed<10000)-Eo2Ec,Phy.Ac(Phy.Zmed>100 & Phy.Zmed<10000),100,'r','linewidth',2)
colormap(cbrewer('seq','Blues',20))
set(gca,'fontsize',16)
xlabel('Temperature Sensitivity, E_e_c_o (eV)')
ylabel('Hypoxia Tolerance, A_e_c_o (atm^-^1)')
xlim([-1.9 1.9])
ylim([0 30])
set(gca,'fontsize',18)
box on
colorbar 
text(-1.8,28,'Habitat Volume (m^3)','fontsize',16)
text(-1.8,26,'(100-1000m)','fontsize',16)

cd figs
saveas(gcf,'Fig3C','pdf')
cd ..

end

%% Fig. 4A: T-warm vs. cold
% T percentiles
TmaxO=obis.Tprct(:,prct==(100-plim));
TmedO=obis.Tprct(:,prct==50);
TminO=obis.Tprct(:,prct==(plim));

% Eo vs T (only for species with dEdT)
dEdT = obis.dEdT;dEdT(obis.dEdT==0)=NaN;
EoO_tmaxO=obis.Eref + dEdT.*(TmaxO-Tref);
EoO_tminO=obis.Eref + dEdT.*(TminO-Tref);
% indexes
Ikp=1:length(obis.F1);
Ikp=obis.F1>0.5;
phimax=obis.phi1e_prct(:,prct==(100-plim));% 95th percentile 
m2m=phimax(Ikp);% relative to Phi-crit =1 
median(m2m,'omitnan')
figure
subplot(121);
dEo = EoO_tmaxO(Ikp)-EoO_tminO(Ikp);
histogram(dEo,0:.1:2); hold on;
xlabel('∆E_e_c_o (eV)');ylabel('# Species')
xlim([0 2])
subplot(122);
histogram(m2m,1:.1:5); xlim([1 3])
xlabel('\Phi_{max}/\Phi_{crit}')
ylim([0 4500])

cd figs
saveas(gcf,'Fig4A','pdf')
cd ..

%% Fig S1A: bivariate F1 vs nobs
if do.F1min == 0
figure
nobs = 1:.25:5;
f1 = 0:0.025:1;
nobsm = 1.125:.25:4.875;
f1m = 0.0125:0.025:.9875;
[N,xe,ye] = histcounts2(obis.F1,log10(double(obis.nobs)),f1,nobs);
N(N==0)=nan;
pcolor(nobsm,f1m,log10(N))
xlabel('log_1_0(Num. of observations)')
ylabel('F_1 score')
colormap jet
colorbar

colormap jet
shading flat

cd figs
saveas(gcf,'FigS1A','eps')
cd ..

%% Fig S1B,C: Eo, Ac vs. F1

fmin = [0 .5:.1:.8];
for ff = 1:length(fmin)
figure(14)
binr = -2:.1:2;
binm=-1.95:.1:1.95;
binc=histcounts(obis.Eo(obis.F1(:)>fmin(ff)),binr);
plot(binm(:),binc(:),'linewidth',2)
hold on
xlabel('Temperature sensitivity (E_o, eV)')
ylabel('Number of species')
legend('F1>0','F1>0.5','F1>0.6','F1>0.7','F1>0.8')
end

cd figs
saveas(gcf,'FigS1C','eps')
cd ..

for ff = 1:length(fmin)

figure(15)
binr = 0:1:35;
binm=.5:1:34.5;
binc=histcounts(obis.Ac(obis.F1(:)>fmin(ff)),binr);
plot(binm(:),binc(:),'linewidth',2)
hold on
xlabel('Active hypoxia tolerance (Ac, atm^-^1)')
ylabel('Number of species')
xlim([0 30])
legend('F1>0','F1>0.5','F1>0.6','F1>0.7','F1>0.8')
end

cd figs
saveas(gcf,'FigS1B','eps')
cd ..

end
%% Fig. S3
% Eo
% low lat
ll.Eo = obis.Eo(abs(obis.Ymed)<30);
ll.Eo=ll.Eo(~isnan(ll.Eo));
% high lat
hl.Eo = obis.Eo(abs(obis.Ymed)>60);
hl.Eo=hl.Eo(~isnan(hl.Eo));

% Ac
% low lat
ll.Ac = obis.Ac(abs(obis.Ymed)<30);
ll.Ac=ll.Ac(~isnan(ll.Ac));
% high lat
hl.Ac = obis.Ac(abs(obis.Ymed)>60);
hl.Ac=hl.Ac(~isnan(hl.Ac));

% ks test
clear E A
[E.h,E.p,E.kstat] = kstest2(ll.Eo,hl.Eo);
[A.h,A.p,A.kstat] = kstest2(ll.Ac,hl.Ac);

% Eo
figure
binr = -2:.1:2;
binm = -1.95:.1:1.95;
binc=histcounts(ll.Eo,binr);
plot(binm,binc./max(binc),'r','linewidth',3)
hold on
binc=histcounts(hl.Eo,binr);
plot(binm,binc./max(binc),'b','linewidth',3)
xlabel('E_o')
ylabel('normalized')
legend('<30 deg.')
legend('>60 deg.')
legend('Low latitude (<30 deg.)','High latitude (>60 deg.)','location','northwest')
ylim([0 1.2])
xlim([-2 2])

cd figs
saveas(gcf,'FigS3A','pdf')
cd ..

% Ac
figure
binr = 0:2:35;
binm = 1:2:34;
binc=histcounts(ll.Ac,binr);
plot(binm,binc./max(binc),'r','linewidth',3)
hold on
binc=histcounts(hl.Ac,binr);
plot(binm,binc./max(binc),'b','linewidth',3)
xlabel('A_c')
ylabel('normalized')
legend('Low latitude (<30 deg.)','High latitude (>60 deg.)','location','northeast')
ylim([0 1.2])

cd figs
saveas(gcf,'FigS3B','pdf')
cd ..

% cdfs 
ebin = -2:.1:2;
Ecdf =zeros(length(ebin),1);
abin = 0:1:35;
Acdf =zeros(length(abin),1);
for aa = 1:length(ebin)
ll.Ecdf(aa)=sum(ll.Eo<ebin(aa));
hl.Ecdf(aa)=sum(hl.Eo<ebin(aa));
end
for aa = 1:length(abin)
ll.Acdf(aa)=sum(ll.Ac<abin(aa));
hl.Acdf(aa)=sum(hl.Ac<abin(aa));
end

figure
plot(abin,ll.Acdf./max(ll.Acdf),'r','linewidth',3)
hold on
plot(abin,hl.Acdf./max(hl.Acdf),'b','linewidth',3)
xlabel('A_c')
ylabel('CDF')
legend('<30 deg.')
legend('>60 deg.')
legend('Low latitude (<30 deg.)','High latitude (>60 deg.)','location','southeast')
ylim([0 1])

cd figs
saveas(gcf,'FigS3C','pdf')
cd ..


figure
plot(ebin,ll.Ecdf./max(ll.Ecdf),'r','linewidth',3)
hold on
plot(ebin,hl.Ecdf./max(hl.Ecdf),'b','linewidth',3)
xlabel('E_o')
ylabel('CDF')
legend('<30 deg.')
legend('>60 deg.')
legend('Low latitude (<30 deg.)','High latitude (>60 deg.)','location','northwest')
ylim([0 1])
xlim([-2 2])


cd figs
saveas(gcf,'FigS3D','pdf')
cd ..


%% Fig. S4: aD vs. lat
figure
scatter(Phy.Ymed,log10(Phy.alphaD),75,'k','filled')
hold on
scatter(Phy.Ymed,log10(Phy.alphaD.*Phy.phic),75,'b','filled')
ylim([-4 4])
set(gca,'fontsize',18)
box on
xlabel('Latitude')
ylabel('Metabolic rate (log10)')
legend('Resting','Sustained','Maximum')
xlim([-80 80])

idx = ~isnan(Phy.Ymed) & ~isnan(Phy.alphaD);
[R,p] = corrcoef(abs(Phy.Ymed(idx)),log10(Phy.alphaD(idx)))

idx = ~isnan(Phy.Ymed) & ~isnan(Phy.alphaD.*Phy.phic);
[R,p] = corrcoef(abs(Phy.Ymed(idx)),log10(Phy.alphaD(idx).*Phy.phic(idx)))


cd figs
saveas(gcf,'FigS4','pdf')
cd ..

%% Fig S7 obis vs lab histograms
figure
%obis
subplot(221)
binr = 0:2:35;
binm = 1:2:34;
idx1 = obis.F1>0.5;
[N,edges] = histcounts(obis.Ac(:),binr);
bar(binm,N,'r','facealpha',0.5)
hold on
xlabel('Active hypoxia tolerance (A_e_c_o, atm^{-1})')
ylabel('Number of species ')
box on
set(gca,'fontsize',18)
ylim([0 8000])

subplot(222)
binr = -2:0.2:2;
binm = -1.9:.2:1.9;
[N,edges] = histcounts(obis.Eo(:),binr);
bar(binm,N,'b','facealpha',0.5)
hold on
ylabel('Number of species')
xlabel('Active Temperature Sensitivity (E_e_c_o, eV)')
set(gca,'fontsize',18)
box on
ylim([0 8000])

% lab
subplot(224)
binr = -2:0.2:2;
binm = -1.9:.2:1.9;
[N,edges] = histcounts(Phy.Eoref(:),binr);
bar(binm,N,'b','facealpha',0.5)
ylabel('Number of species')
xlabel('Resting temperature Sensitivity (E_o, eV)')
set(gca,'fontsize',18)
box on

subplot(223)
binr = 0:5:150;
binm = 2.5:5:147.5;
[N,edges] = histcounts(Phy.Ao(:),binr);
bar(binm,N,'r','facealpha',0.5)
hold on
xlabel('Resting hypoxia tolerance (A_o, atm^{-1})')
ylabel('Number of species ')
box on
set(gca,'fontsize',18)
xlim([0 100])

cd figs
saveas(gcf,'FigS7','pdf')
cd ..


