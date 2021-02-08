clear

% data_dir = '/Users/jess/FP/segmentation_error/Dalhousie_seg/_submission/inverse_solutions/';

root_dir='/Users/jess/FP/segmentation_error/Dalhousie_seg/';
submissions = {'Utah','Nijmegen','OHSU','Inria','Bordeaux'};
reference_dir=[root_dir 'staple_submission/inverse_solutions/'];

m_tot_corr = zeros(1,length(submissions));
m_tot_rrmse = zeros(1,length(submissions));
m_tot_pe = zeros(1,length(submissions));
std_tot_corr = zeros(1,length(submissions));
std_tot_rrmse = zeros(1,length(submissions));
std_tot_pe = zeros(1,length(submissions));


for s = 1:length(submissions)
    
    data_dir = [root_dir submissions{s} '_submission/inverse_solutions/'];
    
    load([data_dir '../geom/mapping_to_staple.mat'])
    mapping = scirunmatrix;
    files = ls('-1', [data_dir '*dist_lap.mat']);

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


    m_tot_corr(s) = mean(tot_corr);
    m_tot_rrmse(s) = mean(tot_rrmse);
    m_tot_pe(s) = mean(tot_pe);
    std_tot_corr(s) = std(tot_corr);
    std_tot_rrmse(s) = std(tot_rrmse);
    std_tot_pe(s) = std(tot_pe);
    
    
end

set(0,'Units','inches')
scrsz_in=get(0,'ScreenSize');
fig_sz=[12 2.5];
compare_1=[ 0 1; scrsz_in(3:4)-fig_sz];
compare_2=[scrsz_in(3:4); fig_sz];
fnt_sz=12;
fnt_nm='Arial';
color1=[40 40 200]/255;
color2=[255 30 0]/255;



figure(3)

fig_sz = [12 2.5];
set(0,'Units','inches')
compare_1=[ fig_sz(1) 3; scrsz_in(3:4)-fig_sz];
compare_2=[scrsz_in(3:4); fig_sz];
set(gcf,'Color',[1 1 1],'PaperPositionMode','auto','Units','inches');
set(gcf,'PaperSize',fig_sz,'Position',...
    [min(compare_1) min(compare_2) ])


subplot(1,3,1)
set(gca,'OuterPosition',[0 0 .33 1])
[b er]=errorbars(m_tot_corr',std_tot_corr',{},{});
set(b(1),'FaceColor',color2)
axis([0.5 5.5 0.5  1])
set(gca,'Box','off')
set(b,'BarWidth',0.5)
ylabel('$\rho$','Interpreter','LaTex','Fontname',fnt_nm,'FontSize',fnt_sz+2)
%title('Whole Map Metrics')
xlabel('Segmentation','Fontname',fnt_nm,'FontSize',fnt_sz+2)
set(gca,'FontSize',fnt_sz)
yt=get(gca,'YLim');
text(0,yt(2),'A','Fontname',fnt_nm,'FontWeight','bold','FontSize',fnt_sz+4)


subplot(1,3,2)
set(gca,'OuterPosition',[.33 0 .33 1])
[b er]=errorbars(m_tot_rrmse',std_tot_rrmse',{},{});
set(b(1),'FaceColor',color2)
axis([0.5 5.5 0  1.5])
set(gca,'Box','off')
set(b,'BarWidth',0.5)
ylabel('rRMSE','Interpreter','LaTex','Fontname',fnt_nm,'FontSize',fnt_sz+2)
xlabel('Segmentation','Fontname',fnt_nm,'FontSize',fnt_sz+2)
set(gca,'FontSize',fnt_sz)
set(gca,'Box','off')
yt=get(gca,'YLim');
text(-0.25,yt(2),'B','Fontname',fnt_nm,'FontWeight','bold','FontSize',fnt_sz+4)

subplot(1,3,3)
set(gca,'OuterPosition',[.66 0 .33 1])
[b er]=errorbars(m_tot_pe',std_tot_pe',{},{});
set(b(1),'FaceColor',color2)
set(b,'BarWidth',0.5)
axis([0.5 5.5 0  1.5])
ylabel('RE','Interpreter','LaTex','Fontname',fnt_nm,'FontSize',fnt_sz+2)
xlabel('Segmentation','Fontname',fnt_nm,'FontSize',fnt_sz+2)
set(gca,'FontSize',fnt_sz)
set(gca,'Box','off')
yt=get(gca,'YLim');
text(-0.25,yt(2),'C','Fontname',fnt_nm,'FontWeight','bold','FontSize',fnt_sz+4)



