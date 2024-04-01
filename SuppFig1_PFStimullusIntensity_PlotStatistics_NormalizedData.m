% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Supp Fig1 PFstimulation_stimulus dependence"
% 3. analyze the data within the folder
% 4. Using ctrl+enter within each function to run each section

%% Section 1: Open files
% 1. ctrl+enter to run this section
% 2. Select "yes" while popup screen show "add new genotype",
%    and then got to the WT foler 
% 2. select SpikeAnalist file with the order from 
%    10uA_SpikeAnalist.mat, 20uA_SpikeAnalist.mat,
%    30uA_SpikeAnalist.mat, to 50uA_SpikeAnalist.mat, and then click
%    "cancel" to end the data selection
% 3. select yes for adding new genotype, and then follow the same steps
%    within L7SK2, and CaMKII folders
% 4. select no for "add new genotype"


clear
mn=1;% mouse#
adding=1; %adding new mouse order (default)
while adding==1
    answer = questdlg('Add new genotype?', ...
        'Add new genotype', ...
        'Yes','No','Yes');
    % Handle response
    switch answer
        case 'Yes'
            FN=[];
            FP=[];
            FNlist={};
            FPlist={};
            condlist=[];
            f=1;
            filefolder=pwd;
            while f
                [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeData files', 'Multiselect', 'on',filefolder);
                filefolder=FolderPath;
                if FolderPath==0;f=0;end

                if iscell(FileName)
                    NewAddFile=size(FileName,2);
                elseif FileName~=0
                    NewAddFile=1;
                    if strfind(FileName,'SpikeAnalist')
                        load([FolderPath,FileName])
                        NewAddFile=size(FileName,2);
                    end
                else
                    NewAddFile=0;
                end

                if NewAddFile~=0;
                    for fnumber=1:NewAddFile
                        if iscell(FileName)
                            FNlist=cat(1,FNlist,FileName{fnumber});
                            FPlist=cat(1,FPlist,FolderPath);
                        elseif strfind(FileName,'SpikeAnalist')
                            FNlist=cat(1,FNlist,FN{fnumber});
                            FPlist=cat(1,FPlist,FP{fnumber});
                        else
                            FNlist=cat(1,FNlist,FileName);
                            FPlist=cat(1,FPlist,FolderPath);
                        end
                        condlist=cat(1,condlist,f);
                    end
                    f=f+1;
                end
            end


            % import data
            clear data;
            for i=1:size(FNlist,1)
                data(i)=load([FPlist{i},FNlist{i}]);
            end

            %
            % pull the data out
            pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
            time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
            smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

            % categorize the accumulated signal (during bust) from individual spike
            % signal
            StimInt=[];
            for i=1:length(data)
                StimInt(i)=str2num(FNlist{i}(end-16:end-15));
            end
            ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

            mouse(mn).smoothBC_signal=smoothBC_signal;
            mouse(mn).pksLOCS=pksLOCS;
            mouse(mn).FNlist=FNlist;
            mouse(mn).StimInt=StimInt;
            mouse(mn).condlist=condlist;
            mouse(mn).time=time;
            mn = mn+1;

        case 'No'
            adding=0; % refuse adding new mouse order
    end
end

ana_window=[0.0 0.2];% timewindow for response spike detection
SPana_window=[-4.5 -0.5];%timewindow for spontaneous activity


%% Section 2: plot statistics for normalized probability
% plotting the probability within "window" sec after stim
prob_by_cond=cell(size(mouse));
StimInten=cell(size(mouse));
for mn=1:length(mouse)
    LOCS=mouse(mn).pksLOCS;
    condlist=mouse(mn).condlist;
    counts=cell(size(condlist));
    SPcounts=cell(size(condlist));
    prob=cell(size(condlist));
    SPprob=cell(size(condlist));
    prob_by_cond{mn}=cell(max(condlist),1);
    SPprob_by_cond{mn}=cell(max(condlist),1);
    Norm_prob_by_cond{mn}=cell(max(condlist),1);
    SPNorm_prob_by_cond{mn}=cell(max(condlist),1);

    time=mouse(mn).time;
    y=[];
    for cond=1:max(condlist)
        %% histogram based on time
        %             cellnum=sum(cellfun(@(a) size(a,2),riseLOCS(condlist==cond)));
        condidx=find(condlist==cond)';
        StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
        %%
        for i=1:length(condidx)
            % average clls
            idx=condidx(i);
            ISinterval=mean(diff(time{idx}));
            st_idx=data(1).stim(1);
            T=time{idx}-time{idx}(st_idx+1);
            for roi=1:size(LOCS{idx},2)

                % because the histcounts include left left edge, so
                % convert the sign, and then left-right flip.
                % because the histcounts include the right edge of the
                % last bin, so add one more bin at right
                counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                prob{idx}(roi,1) = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));

                SPcounts{idx}{roi,1} = cellfun(@(a) (histcounts(fliplr(-T(a)),fliplr(-[SPana_window(1) SPana_window]))./diff(SPana_window)).*diff(ana_window),LOCS{idx}(:,roi),'uni',0);
                SPcounts{idx}{roi,1} = cellfun(@(a) a(1),SPcounts{idx}{roi,1});
            end
            prob_by_cond{mn}{cond}=vertcat(prob_by_cond{mn}{cond},prob{idx});
            SPprob_by_cond{mn}{cond}=vertcat(SPprob_by_cond{mn}{cond},SPcounts{idx});
        end
    end
    SPprob_by_cond{mn}=horzcat(SPprob_by_cond{mn}{:});
    SPprob_by_cond{mn}=arrayfun(@(a) vertcat(SPprob_by_cond{mn}{a,:}),(1:size(SPprob_by_cond{mn},1))','UniformOutput',0);
    SPprob_by_cond{mn}=arrayfun(@(a) nanmean(SPprob_by_cond{mn}{a}),(1:size(SPprob_by_cond{mn},1))');
    SPNorm_prob_by_cond{mn}=SPprob_by_cond{mn}./SPprob_by_cond{mn};
    Norm_prob_by_cond{mn}=arrayfun(@(a) prob_by_cond{mn}{a}./SPprob_by_cond{mn},1:length(prob_by_cond{mn}),'UniformOutput',0);
end

% statistics
% N-way ANOVA stat
ydata=[];
g1=[];
g2=[];
for mn=1:length(mouse)
    ydata=[ydata;SPNorm_prob_by_cond{mn}];
    g1=[g1;repmat(mn,length(SPNorm_prob_by_cond{mn}),1)];%genotype
    g2=[g2;repmat(0,length(SPNorm_prob_by_cond{mn}),1)];%condition
    for cond=1:length(prob_by_cond{mn})
        ydata=[ydata;Norm_prob_by_cond{mn}{cond}];
        g1=[g1;repmat(mn,length(Norm_prob_by_cond{mn}{cond}),1)];%genotype
        g2=[g2;repmat(StimInten{mn}(cond),length(Norm_prob_by_cond{mn}{cond}),1)];%condition
    end
end
[p.anovan,~,stats.anovan] = anovan(ydata,{g1,g2},'model','interaction','varnames',{'genotype','condition'});
c.anovan=multcompare(stats.anovan,'Dimension',[1,2],'Display','off');
g2mat=repmat(1:length(prob_by_cond{mn})+1,3,1);
groupindex=[repmat(1:3,1,length(prob_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups

% plot
   
        co=lines;
        co(6,:)=[.5 .5 .5];
        cnum=[6 4 5];
        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        stim=[0 StimInten{1}];
        figure;hold on;
        x=1:5;
        for mn=1:length(mouse)
            y=mean([SPNorm_prob_by_cond{mn} horzcat(Norm_prob_by_cond{mn}{:})]);
            ysem=nanstd([SPNorm_prob_by_cond{mn} horzcat(Norm_prob_by_cond{mn}{:})],1)./sqrt(length(SPNorm_prob_by_cond{mn}));
            errorbar(x,y,ysem,'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end



        %significance
        sigX=0.5;
        sigY=12.5;
        dsigY=1.7;
        formatSpec = '%.3f';
        FS=16;
        if p.anovan(1)<0.001
            text(sigX,sigY+dsigY*2,'{\it p_g_e_n_e} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY+dsigY*2,['{\it p_g_e_n_e} =' num2str(p.anovan(1),formatSpec)],'FontSize',FS)
        end
        if p.anovan(2)<0.001
            text(sigX,sigY+dsigY,'{\it p_s_t_i_m} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY+dsigY,['{\it p_s_t_i_m} =' num2str(p.anovan(2),formatSpec)],'FontSize',FS)
        end
        if p.anovan(2)<0.001
            text(sigX,sigY,'{\it p_g_e_n_e _x _s_t_i_m} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY,['{\it p_g_e_n_e _x _s_t_i_m} =' num2str(p.anovan(2),formatSpec)],'FontSize',FS)
        end

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(x)
        xlim([0 6])
        yticks([1 5 10])
        ylim([0 17])
        ylabel('Normalized probability')
        xticklabels({'SP' '10' '20' '30' '50'})
        xlabel('PF stimulation (\muA)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])


 
%% Section 3: categorize stimulus intensity-dependent cells according to normalized probability
%% plotting the probability within "window" sec after stim
prob_by_cond=cell(size(mouse));
StimInten=cell(size(mouse));
for mn=1:length(mouse)
    LOCS=mouse(mn).pksLOCS;
    condlist=mouse(mn).condlist;
    counts=cell(size(condlist));
    SPcounts=cell(size(condlist));
    prob=cell(size(condlist));
    SPprob=cell(size(condlist));
    prob_by_cond{mn}=cell(max(condlist),1);
    SPprob_by_cond{mn}=cell(max(condlist),1);
    Norm_prob_by_cond{mn}=cell(max(condlist),1);
    SPNorm_prob_by_cond{mn}=cell(max(condlist),1);

    time=mouse(mn).time;
    y=[];
    for cond=1:max(condlist)
        %% histogram based on time
        %             cellnum=sum(cellfun(@(a) size(a,2),riseLOCS(condlist==cond)));
        condidx=find(condlist==cond)';
        StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
        %%
        for i=1:length(condidx)
            % average clls
            idx=condidx(i);
            ISinterval=mean(diff(time{idx}));
            st_idx=data(1).stim(1);
            T=time{idx}-time{idx}(st_idx+1);
            for roi=1:size(LOCS{idx},2)

                % because the histcounts include left left edge, so
                % convert the sign, and then left-right flip.
                % because the histcounts include the right edge of the
                % last bin, so add one more bin at right
                counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                prob{idx}(roi,1) = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));

                SPcounts{idx}{roi,1} = cellfun(@(a) (histcounts(fliplr(-T(a)),fliplr(-[SPana_window(1) SPana_window]))./diff(SPana_window)).*diff(ana_window),LOCS{idx}(:,roi),'uni',0);
                SPcounts{idx}{roi,1} = cellfun(@(a) a(1),SPcounts{idx}{roi,1});
            end
            prob_by_cond{mn}{cond}=vertcat(prob_by_cond{mn}{cond},prob{idx});
            SPprob_by_cond{mn}{cond}=vertcat(SPprob_by_cond{mn}{cond},SPcounts{idx});
        end
    end
    SPprob_by_cond{mn}=horzcat(SPprob_by_cond{mn}{:});
    SPprob_by_cond{mn}=arrayfun(@(a) vertcat(SPprob_by_cond{mn}{a,:}),(1:size(SPprob_by_cond{mn},1))','UniformOutput',0);
    SPprob_by_cond{mn}=arrayfun(@(a) nanmean(SPprob_by_cond{mn}{a}),(1:size(SPprob_by_cond{mn},1))');
    SPNorm_prob_by_cond{mn}=SPprob_by_cond{mn}./SPprob_by_cond{mn};
    Norm_prob_by_cond{mn}=arrayfun(@(a) prob_by_cond{mn}{a}./SPprob_by_cond{mn},1:length(prob_by_cond{mn}),'UniformOutput',0);
    Norm_prob_by_cond{mn}=horzcat(Norm_prob_by_cond{mn}{:});

end


% plot response-stimulus curve
sigY=1.0;
L_L_dist=0.017; %distance between significant lines
L_D_dist=0.005; %distance between significant line and star symbol

x=1:5;
for mn=1:length(mouse)
    y=[SPNorm_prob_by_cond{mn} Norm_prob_by_cond{mn}];
    [R,P] = arrayfun(@(a) corrcoef(x(~isnan(y(a,:))),rmmissing(y(a,:))),(1:size(y,1))','uni',0);
    RP{mn} = arrayfun(@(a) [R{a}(2) P{a}(2)],(1:length(SPNorm_prob_by_cond{mn}))','uni',0);
    RP{mn}=vertcat(RP{mn}{:});
    corrindex=RP{mn}(:,1)>0.6 & RP{mn}(:,2)<0.5;

    %correlated
    figure;hold on;
    plot(x',[SPNorm_prob_by_cond{mn}(corrindex) Norm_prob_by_cond{mn}(corrindex,:)]','-',...
        'color',[co(cnum(mn),:) 0.3],...
        'MarkerSize',5,...
        'LineWidth',1);
    text(0.5,13,[num2str(sum(corrindex)) '/' num2str(length(corrindex))],'FontSize',18)
    pos=get(gcf,'position');
    pos(3)=350;
    set(gcf,'Position',pos);
    xticks(1:5)
    yticks(0:7:20)
    xlim([0 6])
    ylim([0 14])
    ylabel('Normalized probability')
    xticklabels({'SP' '10' '20' '30' '50'})
    xlabel('PF stimulation (\muA)')
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])

    if sum(~corrindex)
        %non-correlated
        figure;hold on;
        plot(x',[SPNorm_prob_by_cond{mn}(~corrindex) Norm_prob_by_cond{mn}(~corrindex,:)]','-',...
            'color',[co(cnum(mn),:) 0.3],...
            'MarkerSize',5,...
            'LineWidth',1);
        text(0.5,13,[num2str(sum(~corrindex)) '/' num2str(length(corrindex))],'FontSize',18)

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(1:5)
        yticks(0:7:20)
        xlim([0 6])
        ylim([0 14])
        ylabel('Normalized probability')
        xticklabels({'SP' '10' '20' '30' '50'})
        xlabel('PF stimulation (\muA)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
    end
end


%% Section 4: plot statistics for normalized amplitude

%%
        %%
        StimInten=cell(size(mouse));
        amp_by_cond=cell(size(mouse));
        for mn=1:length(mouse)
            smoothBC_signal=mouse(mn).smoothBC_signal;
            LOCS=mouse(mn).pksLOCS;
            condlist=mouse(mn).condlist;

            amps=cell(size(condlist));
            SPamps=cell(size(condlist));
            amps_idx=cell(size(condlist));
            SPamps_idx=cell(size(condlist));
            cellmean=cell(size(condlist));
            SPcellmean=cell(size(condlist));
            amp_by_cond{mn}=cell(max(condlist),1);
            SPamp_by_cond{mn}=cell(max(condlist),1);
            Norm_amp_by_cond{mn}=cell(max(condlist),1);
            SPNorm_amp_by_cond{mn}=cell(max(condlist),1);

            time=mouse(mn).time;

            for cond=1:max(condlist)
                %% histogram based on time
                condidx=find(condlist==cond)';
                StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
                %%
                for i=1:length(condidx)
                    % average clls
                    idx=condidx(i);
                    ISinterval=mean(diff(time{idx}));
                    st_idx=data(1).stim(1);
                    T=time{idx}-time{idx}(st_idx+1);

                    for roi=1:size(LOCS{idx},2)
                       




                        activ=smoothBC_signal{idx}(:,roi,:);
                        ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

                        amps{idx}{roi,1} = activ(T>ana_window(1) & T<=ana_window(2),:,:);%filter out the data out of the ana_window
                        amps{idx}{roi,1} = max(amps{idx}{roi,1},[],1);%find the max within ana_window
                        amps{idx}{roi,1} = amps{idx}{roi,1}(:);
                        cellmean{idx}(roi,1) = nanmean(amps{idx}{roi},1);
         
                        %create new random numbers
                        %saved as suffle_idx_calculateAmp in G:\My Drive\UChicago\Hansel lab\databackup\prob amp summary_20200626\all mice stimulus intensity\WT\SpikeAna
                        rand_idx{idx}(:,roi,:)=randsample(find(T>SPana_window(1) & T<=SPana_window(2)-ana_window(2)),1);%randomly select an index before stim
% rand_idx=rand_idx;
                        SPamps{idx}{roi,1} = activ(rand_idx{idx}(:,roi,:):rand_idx{idx}(:,roi,:)+length(amps{idx}{roi,1})-1,:,:);%filter out the data out of the ana_window
                        SPamps{idx}{roi,1} = max(SPamps{idx}{roi,1},[],1);%find the mean within ana_window

                        SPamps{idx}{roi,1} = SPamps{idx}{roi,1}(:);
                        SPcellmean{idx}{roi} =nanmean(SPamps{idx}{roi},1);  
                    end
                    amp_by_cond{mn}{cond}=vertcat(amp_by_cond{mn}{cond},cellmean{idx});
                    SPamp_by_cond{mn}{cond}=vertcat(SPamp_by_cond{mn}{cond},SPcellmean{idx}');

                end
            end
            SPamp_by_cond{mn}=horzcat(SPamp_by_cond{mn}{:});
            SPamp_by_cond{mn}=arrayfun(@(a) vertcat(SPamp_by_cond{mn}{a,:}),(1:size(SPamp_by_cond{mn},1))','UniformOutput',0);
            SPamp_by_cond{mn}=arrayfun(@(a) nanmean(SPamp_by_cond{mn}{a}),(1:size(SPamp_by_cond{mn},1))');
            SPNorm_amp_by_cond{mn}=SPamp_by_cond{mn}./SPamp_by_cond{mn};
            Norm_amp_by_cond{mn}=horzcat(amp_by_cond{mn}{:});
            Norm_amp_by_cond{mn}=Norm_amp_by_cond{mn}./SPamp_by_cond{mn};
        end
        

        % statistics
        % N-way ANOVA stat
        ydata=[];
        g1=[];
        g2=[];
        for mn=1:length(mouse)
            ydata=[ydata;SPNorm_amp_by_cond{mn}];
            g1=[g1;repmat(mn,length(SPNorm_amp_by_cond{mn}),1)];
            g2=[g2;repmat(0,length(SPNorm_amp_by_cond{mn}),1)];
            for cond=1:size(Norm_amp_by_cond{mn},2)
                ydata=[ydata;Norm_amp_by_cond{mn}(:,cond)];
                g1=[g1;repmat(mn,length(Norm_amp_by_cond{mn}(:,cond)),1)];%genotype
                g2=[g2;repmat(StimInten{mn}(cond),length(Norm_amp_by_cond{mn}(:,cond)),1)];%condition
            end
        end
        [p.anovan,~,stats.anovan] = anovan(ydata,{g1,g2},'model','interaction','varnames',{'genotype','stimulus'});
        c.anovan=multcompare(stats.anovan,'Dimension',[1,2],'Display','off');
        g2mat=repmat(1:length(Norm_amp_by_cond{mn})+1,3,1);
        groupindex=[repmat(1:3,1,length(Norm_amp_by_cond{mn})+1);g2mat(1:end)]';%index of combinations and the corresponding groups
        
        % plot all
        %

        % plot response-stimulus curve
        sigY=1.0;
        L_L_dist=0.017; %distance between significant lines
        L_D_dist=0.005; %distance between significant line and star symbol

        cnum=[6 4 5];

        figure;hold on;
        x=1:5;
        for mn=1:length(mouse)
            y=nanmean([SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}],1);
            ysem=nanstd([SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}],1)./sqrt(length(SPNorm_amp_by_cond{mn}));
            errorbar(x,y,ysem,'o-',...
                'color',co(cnum(mn),:),...
                'CapSize',20,...
                'MarkerSize',5,...
                'LineWidth',3);
        end



        %significance
        sigX=0.5;
        sigY=23;
        dsigY=3;
        formatSpec = '%.3f';
        FS=16;
        if p.anovan(1)<0.001
            text(sigX,sigY+dsigY*2,'{\it p_g_e_n_e} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY+dsigY*2,['{\it p_g_e_n_e} =' num2str(p.anovan(1),formatSpec)],'FontSize',FS)
        end
        if p.anovan(2)<0.001
            text(sigX,sigY+dsigY,'{\it p_s_t_i_m} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY+dsigY,['{\it p_s_t_i_m} =' num2str(p.anovan(2),formatSpec)],'FontSize',FS)
        end
        if p.anovan(2)<0.001
            text(sigX,sigY,'{\it p_g_e_n_e _x _s_t_i_m} < 0.001','FontSize',FS)
        else
            formatSpec = '%.3f';
            text(sigX,sigY,['{\it p_g_e_n_e _x _s_t_i_m} =' num2str(p.anovan(2),formatSpec)],'FontSize',FS)
        end
        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(1:5)
        xlim([0 6])
        yticks([1 10 20])
        ylim([0 30])
        ylabel('Normalized amplitude')
        xticklabels({'SP' '10' '20' '30' '50'})
        xlabel('PF stimulation (\muA)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])

%% Section 5: categorize stimulus intensity-dependent cells according to normalized amplitude
StimInten=cell(size(mouse));
amp_by_cond=cell(size(mouse));
for mn=1:length(mouse)
    smoothBC_signal=mouse(mn).smoothBC_signal;
    LOCS=mouse(mn).pksLOCS;
    condlist=mouse(mn).condlist;

    amps=cell(size(condlist));
    SPamps=cell(size(condlist));
    amps_idx=cell(size(condlist));
    SPamps_idx=cell(size(condlist));
    cellmean=cell(size(condlist));
    SPcellmean=cell(size(condlist));
    amp_by_cond{mn}=cell(max(condlist),1);
    SPamp_by_cond{mn}=cell(max(condlist),1);
    Norm_amp_by_cond{mn}=cell(max(condlist),1);
    SPNorm_amp_by_cond{mn}=cell(max(condlist),1);

    time=mouse(mn).time;

    for cond=1:max(condlist)
        %% histogram based on time
        condidx=find(condlist==cond)';
        StimInten{mn}=horzcat(StimInten{mn},mouse(mn).StimInt(condidx(1)));
        %%
        for i=1:length(condidx)
            % average clls
            idx=condidx(i);
            ISinterval=mean(diff(time{idx}));
            st_idx=data(1).stim(1);
            T=time{idx}-time{idx}(st_idx+1);

            for roi=1:size(LOCS{idx},2)
                activ=smoothBC_signal{idx}(:,roi,:);
                ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

                amps{idx}{roi,1} = activ(T>ana_window(1) & T<=ana_window(2),:,:);%filter out the data out of the ana_window
                amps{idx}{roi,1} = max(amps{idx}{roi,1},[],1);%find the max within ana_window
                amps{idx}{roi,1} = amps{idx}{roi,1}(:);
                cellmean{idx}(roi,1) = nanmean(amps{idx}{roi},1);

                %create new random numbers
                %saved as suffle_idx_calculateAmp in G:\My Drive\UChicago\Hansel lab\databackup\prob amp summary_20200626\all mice stimulus intensity\WT\SpikeAna
                rand_idx{idx}(:,roi,:)=randsample(find(T>SPana_window(1) & T<=SPana_window(2)-ana_window(2)),1);%randomly select an index before stim
                % rand_idx=rand_idx;
                SPamps{idx}{roi,1} = activ(rand_idx{idx}(:,roi,:):rand_idx{idx}(:,roi,:)+length(amps{idx}{roi,1})-1,:,:);%filter out the data out of the ana_window
                SPamps{idx}{roi,1} = max(SPamps{idx}{roi,1},[],1);%find the mean within ana_window

                SPamps{idx}{roi,1} = SPamps{idx}{roi,1}(:);
                SPcellmean{idx}{roi} =nanmean(SPamps{idx}{roi},1);
            end
            amp_by_cond{mn}{cond}=vertcat(amp_by_cond{mn}{cond},cellmean{idx});
            SPamp_by_cond{mn}{cond}=vertcat(SPamp_by_cond{mn}{cond},SPcellmean{idx}');

        end
    end
    SPamp_by_cond{mn}=horzcat(SPamp_by_cond{mn}{:});
    SPamp_by_cond{mn}=arrayfun(@(a) vertcat(SPamp_by_cond{mn}{a,:}),(1:size(SPamp_by_cond{mn},1))','UniformOutput',0);
    SPamp_by_cond{mn}=arrayfun(@(a) nanmean(SPamp_by_cond{mn}{a}),(1:size(SPamp_by_cond{mn},1))');
    SPNorm_amp_by_cond{mn}=SPamp_by_cond{mn}./SPamp_by_cond{mn};
    Norm_amp_by_cond{mn}=horzcat(amp_by_cond{mn}{:});
    Norm_amp_by_cond{mn}=Norm_amp_by_cond{mn}./SPamp_by_cond{mn};
end

% statistics
% N-way ANOVA stat
ydata=[];
g1=[];
g2=[];
for mn=1:length(mouse)
    ydata=[ydata;SPNorm_amp_by_cond{mn}];
    g1=[g1;repmat(mn,length(SPNorm_amp_by_cond{mn}),1)];
    g2=[g2;repmat(0,length(SPNorm_amp_by_cond{mn}),1)];
    for cond=1:size(Norm_amp_by_cond{mn},2)
        ydata=[ydata;Norm_amp_by_cond{mn}(:,cond)];
        g1=[g1;repmat(mn,length(Norm_amp_by_cond{mn}(:,cond)),1)];%genotype
        g2=[g2;repmat(StimInten{mn}(cond),length(Norm_amp_by_cond{mn}(:,cond)),1)];%condition
    end
end
[p.anovan,~,stats.anovan] = anovan(ydata,{g1,g2},'model','interaction','varnames',{'time','genotype'},'display','off');
c.anovan=multcompare(stats.anovan,'Dimension',[1,2],'Display','off');

% plot all
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 4 5 1];
% plot response-stimulus curve
sigY=1.0;
L_L_dist=0.017; %distance between significant lines
L_D_dist=0.005; %distance between significant line and star symbol


% plot response-stimulus curve
sigY=1.0;
L_L_dist=0.017; %distance between significant lines
L_D_dist=0.005; %distance between significant line and star symbol

x=1:5;
for mn=1:length(mouse)
    y=[SPNorm_amp_by_cond{mn} Norm_amp_by_cond{mn}];
    [R,P] = arrayfun(@(a) corrcoef(x(~isnan(y(a,:))),rmmissing(y(a,:))),(1:size(y,1))','uni',0);
    RP{mn} = arrayfun(@(a) [R{a}(2) P{a}(2)],(1:length(SPNorm_amp_by_cond{mn}))','uni',0);
    RP{mn}=vertcat(RP{mn}{:});
    corrindex=RP{mn}(:,1)>0.6 & RP{mn}(:,2)<0.5;
    LinW=3;
    %correlated
    figure;hold on;
    plot(x',[SPNorm_amp_by_cond{mn}(corrindex) Norm_amp_by_cond{mn}(corrindex,:)]','-',...
        'color',[co(cnum(mn),:) .3],...
        'MarkerSize',5,...
        'LineWidth',LinW);
    text(0.5,42,[num2str(sum(corrindex)) '/' num2str(length(corrindex))],'FontSize',18)
    pos=get(gcf,'position');
    pos(3)=350;
    set(gcf,'Position',pos);
    xticks(1:5)
    yticks([1 20 40])
    xlim([0 6])
    ylim([0 45])
    ylabel('Amplitude (\DeltaF/F)')
    xticklabels({'SP' '10' '20' '30' '50'})
    ylabel('Normalized amplitude')
    xlabel('PF stimulation (\muA)')
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])

    if sum(~corrindex)
        figure;hold on;
        plot(x',[SPNorm_amp_by_cond{mn}(~corrindex) Norm_amp_by_cond{mn}(~corrindex,:)]','-',...
            'color',[co(cnum(mn),:) .3],...
            'MarkerSize',5,...
            'LineWidth',LinW);
        text(0.5,9,[num2str(sum(~corrindex)) '/' num2str(length(corrindex))],'FontSize',18)

        pos=get(gcf,'position');
        pos(3)=350;
        set(gcf,'Position',pos);
        xticks(1:5)
        yticks([1 5 10])
        xlim([0 6])
        ylim([0 10])
        ylabel('Normalized amplitude')
        xticklabels({'SP' '10' '20' '30' '50'})
        xlabel('PF stimulation (\muA)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
    end
end

