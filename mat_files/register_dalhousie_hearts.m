clear

data_dir='/Users/jess/FP/segmentation_error/Dalhousie_seg/meshes/';
data_dir2='/Users/jess/FP/segmentation_error/Dalhousie_seg/finished_segs/';


files = ls('-1', [data_dir '*surface_seg.mat']);

files=regexp(files,'\n','split');
files(end)=[];

t_files = ls('-1', [data_dir2 '*points.pts']);

t_files=regexp(t_files,'\n','split');
t_files(end)=[];



load(['~/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/'...
    'Meshes/torso_surface.mat'])

load(['~/Edgar/database/Unmodified/Human_Pacing_Site/Charles_PSTOV-12-07-27/'...
    'Meshes/2012-07-27Subject_daltorso_registered.mat'])
[D, Z, Trans]= procrustes(torso_surface.node,geom.node(:,geom.channels)');
%geom_small.node = (geom.node(:,geom.channels)'*Trans.T + Trans.c)';
geom_small.node = Z;
save(['~/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/'...
    'Meshes/Dalhousie_torso_subset.mat'],'geom_small')
geom_new = geom;
%geom_new.node =(geom.node'*Trans.T + ones(length(geom.node),1)*Trans.c(1,:))';
geom_new.node =(Trans.b*geom.node'*Trans.T + ones(length(geom.node),1)*Trans.c(1,:))';
save(['~/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/'...
    'Meshes/Dalhousie_torso.mat'],'geom_new')


load(['~/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/'...
    'Meshes/torso_surface.mat'])
load(['~/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/'...
    'Meshes/Dalhousie_torso_morphed.mat'])
full_torso.node = scirunfield.node;
full_torso.face = scirunfield.face;


% load([data_dir 'heart_surface_seg.mat'])
% ref_surface = scirunfield;
% data1=ref_surface.node(:,ref_surface.field==2);
% data2=ref_surface.node(:,ref_surface.field==3);
% data3=ref_surface.node(:,ref_surface.field==4);
% set1.data=data1;
% set2.data=data2;
% set3.data=data3;
% 
% opt.delta = 0;
% opt.mode = 'rigid';

% for k = 1:length(files)
%     
%     filename = files{k};
%     load(filename);
%     hs_seg = scirunfield;
%     model1 = hs_seg.node(:,hs_seg.field==2);
%     model2 = hs_seg.node(:,hs_seg.field==3);
%     model3 = hs_seg.node(:,hs_seg.field==4);
%     
%     set1.model=model1;
%     set2.model=model2;
%     set3.model=model3;
%     
%     [R, S, T] = icp_mult(set1,set2,set3,opt);
%     
%     transform = [R*(eye(3)*S) T; 0 0 0 1];
%     
%     dir_split =regexp(filename,'/','split');
%     fname = dir_split(end);
%     froot = fname{1}(1:end-4);
%     
%     save([data_dir froot '_trans.mat'],'transform')
%     
%     
% end
