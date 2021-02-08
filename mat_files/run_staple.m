clear

[vent1,meta] = nrrdread('../../Dalhousie_seg/finished_segs/Utah_ventricles.nrrd');
[vent2,~] = nrrdread('../../Dalhousie_seg/finished_segs/nijmegen_ventricles.nrrd');
[vent3,~] = nrrdread('../../Dalhousie_seg/finished_segs/OHSU_ventricles_masked.nrrd');
[vent4,~] = nrrdread('../../Dalhousie_seg/finished_segs/inria_ventricles.nrrd');
[vent5,~] = nrrdread('../../Dalhousie_seg/finished_segs/bordeaux_ventricles.nrrd');


[torso1,~] = nrrdread('../../Dalhousie_seg/finished_segs/Utah_torso.nrrd');
[torso2,~] = nrrdread('../../Dalhousie_seg/finished_segs/nijmegen_torso.nrrd');
[torso3,~] = nrrdread('../../Dalhousie_seg/finished_segs/OHSU_torso.nrrd');
[torso4,~] = nrrdread('../../Dalhousie_seg/finished_segs/inria_torso.nrrd');
[torso5,~] = nrrdread('../../Dalhousie_seg/finished_segs/bordeaux_torso.nrrd');

[rlung1,~] = nrrdread('../../Dalhousie_seg/finished_segs/Utah_RLung.nrrd');
[rlung2,~] = nrrdread('../../Dalhousie_seg/finished_segs/nijmegen_RLung.nrrd');
[rlung3,~] = nrrdread('../../Dalhousie_seg/finished_segs/bordeaux_RLung.nrrd');


[llung1,~] = nrrdread('../../Dalhousie_seg/finished_segs/Utah_LLung.nrrd');
[llung2,~] = nrrdread('../../Dalhousie_seg/finished_segs/nijmegen_LLung.nrrd');

disp('Data loaded')

[Wv,pv,qv]=STAPLE([vent1(:) vent2(:) vent3(:) vent4(:) vent5(:)]);
%[Wv,pv,qv]=STAPLE([vent1(:) vent2(:) vent4(:)]);
Wv_im = reshape(Wv*255,size(vent1));
nrrdwrite(Wv_im,meta,'../../Dalhousie_seg/finished_segs/Wv.nrrd')
vent_staple = reshape(Wv>=.9,size(vent1));
nrrdwrite(vent_staple,meta,'../../Dalhousie_seg/finished_segs/vent_staple.nrrd')
disp('Ventricle staple done')

% [Wt,pt,qt]=STAPLE([torso1(:) torso2(:) torso3(:) torso4(:) torso5(:)]);
% Wt_im = reshape(Wt*255,size(vent1));
% nrrdwrite(Wt_im,meta,'../../Dalhousie_seg/finished_segs/Wt.nrrd')
% torso_staple=reshape(Wt>=.9,size(vent1));
% nrrdwrite(torso_staple,meta,'../../Dalhousie_seg/finished_segs/toros_staple.nrrd')
% disp('Torso staple done')
% 
% [Wrl,prl,qrl]=STAPLE([rlung1(:) rlung2(:)]);
% Wrl_im = reshape(Wrl*255,size(vent1));
% nrrdwrite(Wrl_im,meta,'../../Dalhousie_seg/finished_segs/Wrl.nrrd')
% rlung_staple=reshape(Wrl>=.9,size(vent1));
% nrrdwrite(rlung_staple,meta,'../../Dalhousie_seg/finished_segs/rlung_staple.nrrd')
% disp('rlung staple done')
% 
% [Wll,pll,qll]=STAPLE([llung1(:) llung2(:)]);
% Wll_im = reshape(Wll*255,size(vent1));
% nrrdwrite(Wll_im,meta,'../../Dalhousie_seg/finished_segs/Wll.nrrd')
% llung_staple=reshape(Wll>=.9,size(vent1));
% nrrdwrite(llung_staple,meta,'../../Dalhousie_seg/finished_segs/llung_staple.nrrd')
% disp('llung staple done')


dv(1)=DICE(vent1,vent_staple);
dv(2)=DICE(vent2,vent_staple);
dv(3)=DICE(vent3,vent_staple);
dv(4)=DICE(vent4,vent_staple);
dv(5)=DICE(vent5,vent_staple);

set(0,'Units','inches')
scrsz_in=get(0,'ScreenSize');
fig_sz=[4 2.5];
compare_1=[ 0 1; scrsz_in(3:4)-fig_sz];
compare_2=[scrsz_in(3:4); fig_sz];
fnt_sz=12;
fnt_nm='Arial';
color1=[40 40 200]/255;
color2=[255 30 0]/255;
color3=[50 255 50]/255;

exp_name = {'seg 1','seg 2','seg 3','seg 4','seg 5',};

figure(2)
set(gcf,'Color',[1 1 1],'PaperPositionMode','auto','Units','inches');
set(gcf,'PaperSize',fig_sz,'Position',...
    [min(compare_1) min(compare_2) ])
% subplot(1,3,1)
% set(gca,'OuterPosition',[0 0 .33 1])
b=bar([dv; pv; qv]');
set(b(1),'FaceColor',color1)
set(b(2),'FaceColor',color2)
set(b(3),'FaceColor',color3)
set(gca,'XTickLabel',exp_name,'Fontname',fnt_nm,'Layer','top','Box','off')
ylabel('Ventricles','Interpreter','LaTex','Fontname',fnt_nm,'FontSize',fnt_sz+2)
% ylabel('$\rho$','Interpreter','LaTex','Fontname',fnt_nm,'FontSize',fnt_sz+2)
axis([0.5 5.5 0 1])
%title('Whole Map Metrics')
% xlabel('Segmentation','Fontname',fnt_nm,'FontSize',fnt_sz+2)
lh=legend('DICE','sensitivity','specificity','Location','best');
% yt=get(gca,'YLim');
% text(0,yt(2),'A','Fontname',fnt_nm,'FontWeight','bold','FontSize',fnt_sz+4)
set(gca,'FontSize',fnt_sz)


% tv(1)=DICE(torso1,torso_staple);
% tv(2)=DICE(torso2,torso_staple);
% tv(3)=DICE(torso3,torso_staple);
% tv(4)=DICE(torso4,torso_staple);
% tv(5)=DICE(torso5,torso_staple);
% 
% drl(1)=DICE(rlung1,rlung_staple);
% drl(2)=DICE(rlung2,rlung_staple);
% 
% dll(1)=DICE(llung1,llung_staple);
% dll(2)=DICE(llung2,llung_staple);