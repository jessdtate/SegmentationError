clear

data_dir = '/Users/jess/FP/segmentation_error/Dalhousie_seg/Utah_submission/inverse_solutions/';
reference_dir = '/Users/jess/FP/segmentation_error/Dalhousie_seg/staple_submission/inverse_solutions/';

load([data_dir '../geom/mapping_to_staple.mat'])
mapping = scirunmatrix;
files = ls('-1', [data_dir '*.mat']);

files=regexp(files,'\n','split');
files(end)=[];

tot_pot=[];
tot_ref_pot = [];
tot_corr = [];
tot_pe = [];
tot_rrmse = [];
for k =1:length(files)
    
    filename = files{k};
    load(filename);
    pot = scirunmatrix;
    
    dir_split =regexp(filename,'/','split');
    fname = dir_split(end);
    froot = fname{1}(1:end-4);
    
    if ~exist([reference_dir fname{1}],'file')
        continue
    end
    
    load([reference_dir fname{1}])
    ref_pot = scirunmatrix;
   
    %tot_pot= [tot_pot pot];
    %tot_ref_pot= [tot_ref_pot ref_pot];
    corr=Correlation(ref_pot,mapping*pot);
    pe = pError(ref_pot,mapping*pot);
    rrmse = rRMSError(ref_pot,mapping*pot);
    
    tot_corr = [tot_corr corr];
    tot_pe = [tot_pe pe];
    tot_rrmse = [tot_rrmse rrmse];
    
    

end


mean(tot_corr)
mean(tot_rrmse)



