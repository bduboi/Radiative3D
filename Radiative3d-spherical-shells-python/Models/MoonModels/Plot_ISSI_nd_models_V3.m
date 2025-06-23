% plot of lunar internal structure models provided for ISSI team

Vs_dens=load('MoonKhan2002_Vs_probability_density.txt');
Vp_dens=load('MoonKhan2002_Vp_probability_density.txt');
dkhan = 0.0 - [1:1:1736];
vpkhan = reshape(Vp_dens(:,3),11,1736);
vpdens = reshape(Vp_dens(:,1),11,1736);
vskhan = reshape(Vs_dens(:,3),10,1736);
vsdens = reshape(Vs_dens(:,1),10,1736);

filemod{1}='MoonToksoz_RGSP_1974.nd'
filemod{2}='MoonNakamura_JGR_1983.nd'
filemod{3}='MoonLognonneC_EPSL_2003.nd'
filemod{4}='MoonBeyneixC_PEPI_2006.nd'
filemod{5}='MoonWeber_Science_2011.nd'
filemod{6}='MoonGarcia_PEPI_2011.nd'
filemod{7}='MoonKhan_JGR_2014_mean.nd'
filemod{8}='MoonMatsumoto_GRL_2015.nd'
filemod{9}='MoonKhan_JGR_2014_std.nd'
filemod{10}='MoonMatsumoto_GRL_2015_std.nd'
filemod{11}='GilletMargerin_attenuation_2016.nd'

textlegend{1}='Toksoz et al. 1974'
textlegend{2}='Nakamura et al., 1983'
textlegend{3}='Lognonne et al., 2003'
textlegend{4}='Beyneix et al., 2006'
textlegend{5}='Weber et al., 2011'
textlegend{6}='Garcia et al., 2011'
textlegend{7}='Khan et al., 2014'
textlegend{8}='Matsumoto et al., 2015'

nkhan=7;
nvpremoon=6;
nmatsu=8;

% load error estimates associated to VPREMOON model
load('MoonGarcia_PEPI_2011_std.mat')

nbfile=length(filemod)

% read model files
for i=1:nbfile 
    [smodel{i},nl(i)]=read_nd_models(filemod{i});
end

% define model colors
clor = lines(nbfile) ;
% clor = hsv(nbfile) ;
%clor = flag(nbfile) ;
clor(8,1)=0.45

% remove the standard deviation file from Amir + scattering model by Gillet
nbfile=nbfile-3

% for i=1:nbfile
%     textdum=filemod{i};
%     nl=length(textdum)
%     textlegend{i}=textdum(5:(nl-3))
% end



% plot density Vp Vs
figure
for i=1:nbfile
model=smodel{i};

ax(1)=subplot(1,3,1)
hold on;
plot(model(:,2),0.0-model(:,1),'-','color',clor(i,:),'LineWidth',2)

ax(2)=subplot(1,3,2)
hold on;
plot(model(:,3),0.0-model(:,1),'-','color',clor(i,:),'LineWidth',2)

ax(3)=subplot(1,3,3)
hold on;
plot(model(:,4),0.0-model(:,1),'-','color',clor(i,:),'LineWidth',2)

end
% include Khan 2002 contour of probability density
subplot(1,3,1)
%contour(vpdens(:,1),dkhan(1:1200),vpkhan(:,1:1200)',4,'k--','LineWidth',1)
contour(vpdens(:,1),dkhan(1:530),vpkhan(:,1:530)',1,'k--','LineWidth',1)

subplot(1,3,2)
%contour(vsdens(:,1),dkhan(1:1200),vskhan(:,1:1200)',4,'k--','LineWidth',1)
contour(vsdens(:,1),dkhan(1:530),vskhan(:,1:530)',1,'k--','LineWidth',1)


% include VPREMOON errors
model=smodel{nvpremoon};
subplot(1,3,1)
plot(model(:,2)-sigvp(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)
plot(model(:,2)+sigvp(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)

subplot(1,3,2)
plot(model(:,3)-sigvs(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)
plot(model(:,3)+sigvs(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)

subplot(1,3,3)
plot(model(:,4)-sigrho(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)
plot(model(:,4)+sigrho(1,:)',0.0-model(:,1),'--','color',clor(nvpremoon,:),'LineWidth',1)


% include Khan errors
model=smodel{nkhan};
modelstd=smodel{nkhan+2};
subplot(1,3,1)
plot(model(:,2)-modelstd(:,2),0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)
plot(model(:,2)+modelstd(:,2),0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)

subplot(1,3,2)
plot(model(:,3)-modelstd(:,3),0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)
plot(model(:,3)+modelstd(:,3),0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)

subplot(1,3,3)
plot(model(:,4)-modelstd(:,4)/1000.0,0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)
plot(model(:,4)+modelstd(:,4)/1000.0,0.0-model(:,1),'--','color',clor(nkhan,:),'LineWidth',1)


% include Matsumoto errors
model=smodel{nmatsu};
modelstd=smodel{nmatsu+2};
subplot(1,3,1)
plot(model(:,2)-modelstd(:,2),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
plot(model(:,2)+modelstd(:,2),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
ylabel('Depth (in km)')
xlabel('Vp (in km/s)')
grid 'on'
xlim([0.9 13])

subplot(1,3,2)
plot(model(:,3)-modelstd(:,3),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
plot(model(:,3)+modelstd(:,3),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
xlabel('Vs (in km/s)')
grid 'on'
xlim([-0.1 7.0])

subplot(1,3,3)
plot(model(:,4)-modelstd(:,4),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
plot(model(:,4)+modelstd(:,4),0.0-model(:,1),'--','color',clor(nmatsu,:),'LineWidth',1)
xlabel('Density (in g/cm^3)')
grid 'on'
xlim([2.4 8.1])

lgd=legend(textlegend{1:nbfile},'Location','SouthEast')
lgd.Interpreter='none'
linkaxes(ax,'y')

ylim([-1738 0])

% % PLot attenuation models
% figure
% for i=1:7:nbfile
% model=smodel{i}
% 
% bx(1)=subplot(1,2,1)
% semilogx(model(:,5),0.0-model(:,1),'-','color',clor(i,:))
% hold on;
% 
% bx(2)=subplot(1,2,2)
% semilogx(model(:,6),0.0-model(:,1),'-','color',clor(i,:))
% hold on;
% 
% 
% end
% 
% % scattering attenuation model by Gillet et al., 2016
% model=smodel{6}
% subplot(1,2,1)
% semilogx(model(:,7),0.0-model(:,1),'-','color',clor(i+1,:))
% ylabel('Depth (in km)')
% xlabel('Qp or Qscatt')
% grid('on')
% 
% subplot(1,2,2)
% semilogx(model(:,6),0.0-model(:,1),'-','color',clor(i+1,:))
% xlabel('Qs or Qintrinsic')
% grid('on')
% legendary={filemod{1:nbfile} filemod{6}};
% %legendary{nbfile+1}=filemod{6};
% legend(legendary,'Location','SouthEast')
% 
% linkaxes(bx,'y')


% Function to read the nd files
function [model,nl]=read_nd_models(filename);
    fid=fopen(filename,'r')

nl=0;
test=0;
while (test<1)
    tline = fgetl(fid);
    if (sum(size(tline))>15)
        A=str2num(tline)
        nl=nl+1;
        model(nl,:)=A(1,:);
    end
    test=feof(fid)
end

end