clear

[vent1,meta] = nrrdread('../Dalhousie_seg/finished_segs/Utah_ventricles.nrrd');
[vent2,~] = nrrdread('../Dalhousie_seg/finished_segs/nijmegen_ventricles.nrrd');
[vent3,~] = nrrdread('../Dalhousie_seg/finished_segs/OHSU_ventricles_masked.nrrd');
[vent4,~] = nrrdread('../Dalhousie_seg/finished_segs/inria_ventricles.nrrd');
[vent5,~] = nrrdread('../Dalhousie_seg/finished_segs/bordeaux_ventricles.nrrd');


[torso1,~] = nrrdread('../Dalhousie_seg/finished_segs/Utah_torso.nrrd');
[torso2,~] = nrrdread('../Dalhousie_seg/finished_segs/nijmegen_torso.nrrd');
[torso3,~] = nrrdread('../Dalhousie_seg/finished_segs/OHSU_torso.nrrd');
[torso4,~] = nrrdread('../Dalhousie_seg/finished_segs/inria_torso.nrrd');
[torso5,~] = nrrdread('../Dalhousie_seg/finished_segs/bordeaux_torso.nrrd');

[rlung1,~] = nrrdread('../Dalhousie_seg/finished_segs/Utah_RLung.nrrd');
[rlung2,~] = nrrdread('../Dalhousie_seg/finished_segs/nijmegen_RLung.nrrd');
[rlung3,~] = nrrdread('../Dalhousie_seg/finished_segs/bordeaux_RLung.nrrd');


[llung1,~] = nrrdread('../Dalhousie_seg/finished_segs/Utah_LLung.nrrd');
[llung2,~] = nrrdread('../Dalhousie_seg/finished_segs/nijmegen_LLung.nrrd');

disp('Data loaded')

[Wv,pv,qv]=STAPLE([vent1(:) vent2(:) vent3(:) vent4(:) vent5(:)]);
%[Wv,pv,qv]=STAPLE([vent1(:) vent2(:) vent4(:)]);
Wv_im = reshape(Wv*255,size(vent1));
nrrdwrite(Wv_im,meta,'../Dalhousie_seg/finished_segs/Wv.nrrd')
vent_staple = reshape(Wv>=.9,size(vent1));
nrrdwrite(vent_staple,meta,'../Dalhousie_seg/finished_segs/vent_staple.nrrd')
disp('Ventricle staple done')

[Wt,pt,qt]=STAPLE([torso1(:) torso2(:) torso3(:) torso4(:) torso5(:)]);
Wt_im = reshape(Wt*255,size(vent1));
nrrdwrite(Wt_im,meta,'../Dalhousie_seg/finished_segs/Wt.nrrd')
torso_staple=reshape(Wt>=.9,size(vent1));
nrrdwrite(torso_staple,meta,'../Dalhousie_seg/finished_segs/toros_staple.nrrd')
disp('Torso staple done')

[Wrl,prl,qrl]=STAPLE([rlung1(:) rlung2(:)]);
Wrl_im = reshape(Wrl*255,size(vent1));
nrrdwrite(Wrl_im,meta,'../Dalhousie_seg/finished_segs/Wrl.nrrd')
rlung_staple=reshape(Wrl>=.9,size(vent1));
nrrdwrite(rlung_staple,meta,'../Dalhousie_seg/finished_segs/rlung_staple.nrrd')
disp('rlung staple done')

[Wll,pll,qll]=STAPLE([llung1(:) llung2(:)]);
Wll_im = reshape(Wll*255,size(vent1));
nrrdwrite(Wll_im,meta,'../Dalhousie_seg/finished_segs/Wll.nrrd')
llung_staple=reshape(Wll>=.9,size(vent1));
nrrdwrite(llung_staple,meta,'../Dalhousie_seg/finished_segs/llung_staple.nrrd')
disp('llung staple done')

dv(1)=DICE(vent1,vent_staple);
dv(2)=DICE(vent2,vent_staple);
dv(3)=DICE(vent3,vent_staple);
dv(4)=DICE(vent4,vent_staple);
dv(5)=DICE(vent5,vent_staple);

tv(1)=DICE(torso1,torso_staple);
tv(2)=DICE(torso2,torso_staple);
tv(3)=DICE(torso3,torso_staple);
tv(4)=DICE(torso4,torso_staple);
tv(5)=DICE(torso5,torso_staple);

drl(1)=DICE(rlung1,rlung_staple);
drl(2)=DICE(rlung2,rlung_staple);

dll(1)=DICE(llung1,llung_staple);
dll(2)=DICE(llung2,llung_staple);