% Code for Fig. 6 (Ephicrit analysis)
clear all; close all

% obis and lab traits fit over same temperature range from laboratory experiments 
load obis_AE_Tlab.mat

% Remove non-marine species
Sp_cases={'Carassius carassius','Oreochromis niloticus','Palaemon pugio'};
% 'Carassius carassius' Crucian carp, freshwater fish
% 'Oreochromis niloticus' Nile fish, largely freshwater
% 'Palaemon pugio' grass shrimp, brackish/fresh 
% exclude freshwater species
for ii=1:length(Sp_cases)
    obis.Eolab(strcmp(obis.species,Sp_cases(ii)))=NaN;
    obis.Aolab(strcmp(obis.species,Sp_cases(ii)))=NaN;

end

% update lab data
obis.Eolab(strcmp(obis.species,"Cancer irroratus"))=0.88;% includes all lifestages
obis.Aolab(strcmp(obis.species,"Cancer irroratus"))=30.8;% includes all lifestages 

% traits based on fitting Arrhenius over lab T range with dE/dT = 0
obis.Eeco=obis.p1s(:,2); % Eeco
obis.F1=1-obis.f1;% F1 score

% solve for Ephicrit 
obis.Ephicrit = obis.Eeco-obis.Eolab; % eV

% Experimental MMR:RMR data
load mmr_temp_data.mat
input=mmr;clear mmr
input.species_unq = unique(input.species);
for ii = 1:length(input.species_unq)

    idx_mmr = find(strcmp(input.species,input.species_unq(ii)));
    if strcmp(input.species_unq(ii),'Centropristis striata 1') || strcmp(input.species_unq(ii),'Centropristis striata 2') 
        input.species_unq(ii) = 'Centropristis striata';
    end
    
    % obis index
    idxo = find(strcmp(obis.species,input.species_unq(ii)));

   
    %extract data
    temp = input.temp(idx_mmr);
    mmr = input.mmr(idx_mmr);
    rmr = input.rmr(idx_mmr);
    fas = mmr./rmr;

    % constants
    tref=15;% reference T (oC)
    kb = 8.6173324*10^-5; % Boltzmann's constant (eV/K)
    texp = ((1./(temp+273.15))-(1./(tref+273.15)))./kb; % Arrhenius temperature

    % FAS = MMR./RMR 
    % FAS = fas_act*exp(-Efas,T)
    % fit arrhenius to FAS
    y=log(fas);
    x = -texp;
    p = polyfit(x,y,1);
    Efas=p(1);
    fas_ref = exp(p(2));
    fas_fit = fas_ref.*exp(-Efas.*texp);

    %save
    out.fas_ref(ii)=fas_ref; % fas at Tref
    out.Efas(ii)=Efas; % T-dep, eV
    
    % Traits for species with FAS vs T
    out.Eo(ii)=obis.Eolab(idxo); % lab Eo (Eo_rest) 
    out.Ephicrit(ii)=obis.Ephicrit(idxo); % Ephicrit from obis + lab
    out.Eeco(ii) = obis.Eeco(idxo);% E_eco (obis) 
    out.F1(ii)=obis.F1(idxo); % F1-score 
    out.Eeco_lab(ii) = out.Efas(ii)+out.Eo(ii);% Eeco_lab  = Efas + Eo (lab data)
    
end

% Fig. 6
figure
subplot(211)
dE=0.2;
binranges = -2:dE:2;
binm=-2+dE/2:dE:2-dE/2;
filt=obis.F1>0.5;% filtering criteria
binc=histcounts(obis.Eeco(filt),binranges);
plot(binm,binc,'b','linewidth',3);
hold on
binc=histcounts(obis.Ephicrit(filt),binranges);
bar(binm,binc,1,'facecolor',[0 0 0],'facealpha',0.25,'edgealpha',1);
binc=histcounts(obis.Eolab(filt),binranges);
bar(binm,binc,1,'facecolor',[1 0 0],'facealpha',0.25,'edgealpha',1);
hold on
plot([0,0],[0,50],'k')
xlabel('Temperature sensitivity, E (eV)')
ylabel('Number of species')
text(.2,18,'Active (E_e_c_o)','color','b','fontsize',18)
text(.5,13,'Resting (E_o)','color',[1 .7 .7],'fontsize',18)
text(-1.6,11,'Active/Resting (E_s_m_s)','color',[.6 .6 .6],'fontsize',18)
xlim([-1.5 1.5])
ylim([0 22])
title('Biogeographic')
set(gca,'fontsize',16);

subplot(212)
x =[1 2 3 4 5 6];% species number
del=0.25;% offsets for visualization
x2 = x+del;% offsets for visualization
x3 = x2+del; % offsets for visualization

% reorder species by decreasing Efas for plotting
y = [out.Eo(3);out.Eo(5);out.Eo(2) ;out.Eo(1) ;out.Eo(4) ; out.Eo(6)]; % lab Eo
y2 = [out.Eeco_lab(3);out.Eeco_lab(5);out.Eeco_lab(2) ;out.Eeco_lab(1) ;out.Eeco_lab(4) ; out.Eeco_lab(6)]; % lab Eeco (Efas+Eo)
y3 = [ out.Efas(3); out.Efas(5) ;out.Efas(2) ;out.Efas(1) ;out.Efas(4) ; out.Efas(6)]; % lab Ephicrit
y4=[ out.Ephicrit(3); out.Ephicrit(5) ;out.Ephicrit(2) ;out.Ephicrit(1) ;out.Ephicrit(4) ; out.Ephicrit(6)]; % obis Ephicrit

barh(x3,y3,del,'k','facealpha',.25) % lab Ephicrit
hold on
title('Experimental')
scatter(y4,x3,50,[.4 .4 .4]*0,'filled') % obis Ephicrit
barh(x2,y2,del,'b') % lab Eeco
hold on
barh(x,y,del,'r','facealpha',.25)% lab Eo
set(gca,'yticklabels',[])
ylabel('Species')
ylim([0.5 7])
xlabel('Temperature sensitivity, E (eV)')
xlim([-1.5 1.5])
text(-1.4,6.25,'S. ocellatus','fontsize',14)
text(-1.4,5.25,'C. lumpus','fontsize',14)
text(-1.4,4.25,'C. striata (flume)','fontsize',14)
text(-1.4,3.25,'C. striata (chase)','fontsize',14)
text(-1.4,2.5,'P. maximus','fontsize',14)
text(-1.4,1,'C. atripectoralis','fontsize',14)
set(gca,'fontsize',16);

% correlation statistics
[r,p]=corrcoef(out.Efas,out.Ephicrit)
r.^2
p
