clear

root_dir='/Users/jess/FP/segmentation_error/Dalhousie_seg/';
submissions = {'Utah','Nijmegen','OHSU','Inria','Bordeaux'};

ref_dir=[root_dir 'staple_submission/inverse_solutions/'];

file_num = [1,3,12,18];
beat_label={'sinus','paced','paced (phrenic)','RV paced'};
f_it = 1;

scrsz_in=get(0,'ScreenSize');
lw = 2; %line width
ms = 10; %marker size
fs = 12; %font size

files = ls('-1', [ref_dir '*dist_lap.mat']);
files=regexp(files,'\n','split');
files(end)=[];

t_var_tot = [];
t_std_tot = [];

for f = 1:length(files)

    dir_split =regexp(files{f},'/','split');
    fname = dir_split(end);
    froot = fname{1}(1:end-4);

    load(files{f})
    ref_pot = scirunmatrix;

    ref_rms = rmsVoltage(ref_pot);

    tot_pot = zeros([size(ref_pot) length(submissions)]);

    corr_all = zeros(length(submissions),length(ref_rms));
    pe_all = zeros(length(submissions),length(ref_rms));
    rrmse_all = zeros(length(submissions),length(ref_rms));
    rmsv_all = zeros(length(submissions),length(ref_rms));

    for k=1:length(submissions)
        load([root_dir submissions{k} '_submission/geom/mapping_to_staple.mat'])
        mapping = scirunmatrix;

        if ~exist([root_dir submissions{k} '_submission/inverse_solutions/' fname{1}],'file')
            error([files{k} ' does not exist for all submissions']);
        end
        load([root_dir submissions{k} '_submission/inverse_solutions/' fname{1}])
        pot = scirunmatrix;
        mp = mapping*pot;
        tot_pot(:,:,k) = mp;

        corr=Correlation(ref_pot,mp);
        pe = pError(ref_pot,mp);
        rrmse = rRMSError(ref_pot,mp);
        rmsv = rmsVoltage(mp);

        corr_all(k,:) = corr;
        pe_all(k,:) = pe;
        rrmse_all(k,:) = rrmse;
        rmsv_all(k,:) = rmsv;


    end

   
    t_var = var(tot_pot/1000,[],3);
    t_std = std(tot_pot/1000,[],3);
    var_map= mean(t_var,2);
    std_map= mean(t_std,2);
    var_time = mean(t_var);
    std_time = mean(t_std);
    
    t_var_tot = [t_var_tot t_var];
    t_std_tot = [t_std_tot t_std];
    

    save([ref_dir '../stats/' froot '_var_map.mat'],'var_map')
    save([ref_dir '../stats/' froot '_std_map.mat'],'std_map')
    save([ref_dir '../stats/' froot '_var_time.mat'],'t_var')
    save([ref_dir '../stats/' froot '_std_time.mat'],'t_std')

    max(max((var(tot_pot/1000,[],3)./(max(abs(tot_pot/1000),[],3)-...
        min(abs(tot_pot/1000),[],3)))));
    
    t = (0:length(corr)-1)/2000;

    if f_it<=length(file_num) && f == file_num(f_it)
%     if true
        f_it = f_it+1;
        figure(1)

        fig_sz = [12 2.5];
        set(0,'Units','inches')
        compare_1=[ fig_sz(1) 1; scrsz_in(3:4)-fig_sz];
        compare_2=[scrsz_in(3:4); fig_sz];
        set(gcf,'Color',[1 1 1],'PaperPositionMode','auto','Units','inches');
        set(gcf,'PaperSize',fig_sz,'Position',...
            [min(compare_1) min(compare_2) ])

        subplot(1,3,1)
        set(gca,'OuterPosition',[0 0 .33 1])
        plot(t,ref_rms/1000,'-k','LineWidth',lw)
        hold on
        plot(t,rmsv_all/1000','LineWidth',lw)
        axis([0 max(t) 0  5])
        ylabel('RMS Voltage (mV)')
        xlabel('time (s)')
        set(gca,'FontSize',fs)
        set(gca,'Box','off')
        legend('aggregate','segmentation 1','segmentation 2','segmentation 3',...
            'segmentation 4','segmentation 5','Location','best')

        subplot(1,3,2)
        set(gca,'OuterPosition',[.33 0 .33 1])
        plot(t,rrmse_all','LineWidth',lw)
        axis([0 max(t) 0 1.5 ])
        ylabel('rRMSE')
        xlabel('time (ms)')
        set(gca,'FontSize',fs)
        set(gca,'Box','off')

        subplot(1,3,3)
        set(gca,'OuterPosition',[.66 0 .33 1])
        plot(t,corr_all','LineWidth',lw)
        axis([0 max(t) 0.5 1 ])
        ylabel('\rho')
        xlabel('time (ms)')
        set(gca,'FontSize',fs)
        set(gca,'Box','off')
        
        figure(2)
        
        fig_sz = [4 2.5];
        set(0,'Units','inches')
        compare_1=[ 1 1; scrsz_in(3:4)-fig_sz];
        compare_2=[scrsz_in(3:4); fig_sz];
        set(gcf,'Color',[1 1 1],'PaperPositionMode','auto','Units','inches');
        set(gcf,'PaperSize',fig_sz,'Position',...
            [min(compare_1) min(compare_2) ])
        plot(t,var_time,'-k','LineWidth',lw)
        axis([0 max(t) 0  2.5])
        ylabel('Variance (mV^2)')
        xlabel('time (s)')
        set(gca,'Box','off')

    end


    

    
end

var_tot_map=mean(t_var_tot,2);
std_tot_map=mean(t_std_tot,2);



save([ref_dir '../stats/var_tot_map.mat'],'var_tot_map')
save([ref_dir '../stats/std_tot_map.mat'],'std_tot_map')


    


