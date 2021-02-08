clear

%  create regularization matrice for inverse

root_dir = '/Users/jess/FP/segmentation_error/';
segmentation_dir = 'Dalhousie_seg/';

path(path,'/Users/jess/software/activationtimes/computeGradientAndHessian')

files = ls('-1',[root_dir segmentation_dir '*_submission/geom/heart_surface_small.mat']);

file=regexp(files,'\n','split');
file(end) =[];

for k = 1:length(file)
    
    
    tmp_filename = file{k};
   
    slashies = regexp(tmp_filename,'/');
   
    out_file = [tmp_filename(1:slashies(end)) 'dist_laplacian_matrix.mat'];
    Wieghts_file = [tmp_filename(1:slashies(end)) 'weights_matrix.mat'];
   
    load(tmp_filename);
    
    heart.node = scirunfield.node';
    heart.face = scirunfield.face';
%    heart.node = heart.node';
%    heart.face = heart.face';
   
    [D] = PairwiseDistance(heart.node');

	% choose sigma
	sigma = (sort(D));
	sigma = max(sigma,[],2);
	sigma = sigma(3);
    
%     sigma = 30;
    
    iE = exp( -D./(2*sigma));
    
    save(Wieghts_file,'iE')
    

	[wghFcn] = invExplonentialWeight(heart, sigma, 1e-5);

	[cDf, cHf] = meshVolDiffHessMatrix(heart,wghFcn);
	
	LHf = LaplacianMatrixFromHessianMatrix(cHf);
    
    save(out_file,'LHf')
    
    disp(tmp_filename)
    
   
end

   

