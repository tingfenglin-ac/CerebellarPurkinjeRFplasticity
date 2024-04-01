% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Fig8 PF x tactile interaction\paw"
% 3. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section

clear
FN=[];
FP=[];
FNlist={};

FPlist={};
condlist=[];


f=1;
while f
    [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeData files', 'Multiselect', 'on');

    if FolderPath==0;f=0;end

    if iscell(FileName)
        NewAddFile=size(FileName,2);
    elseif FileName~=0
        NewAddFile=1;
        if strfind(FileName,'SpikeAnalist')
            load([FolderPath,FileName])
            NewAddFile=size(FN,1);
        end
    else
        NewAddFile=0;
    end
    %%
    if NewAddFile~=0;
        %         if ~size(strfind(FileName,'Spikelist'),1)|iscell(FileName)
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

        %         end
    end
end


% import data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

%
% pull the data out
riseLOCS=arrayfun(@(x) data(x).riseLOCS,1:length(data),'uni',0);
upLOCS=arrayfun(@(x) data(x).upLOCS,1:length(data),'uni',0);
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
PKS=arrayfun(@(x) data(x).PKS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

% categorize the accumulated signal (during bust) from individual spike
% signal
ind_rise_LOCS=cell(size(riseLOCS));
accm_rise_LOCS=cell(size(riseLOCS));
ind_pks_LOCS=cell(size(riseLOCS));
accm_pks_LOCS=cell(size(riseLOCS));
ind_PKS=cell(size(riseLOCS));
accm_PKS=cell(size(riseLOCS));
for i=1:length(data)
    ind_rise_LOCS{i}=cellfun(@(x,y) x(~ismember(x,y)),riseLOCS{i},upLOCS{i},'uni',0);
    accm_rise_LOCS{i}=cellfun(@(x,y) x(ismember(x,y)),riseLOCS{i},upLOCS{i},'uni',0);
    ind_pks_LOCS{i}=cellfun(@(x,y,z) z(~ismember(x,y)),riseLOCS{i},upLOCS{i},pksLOCS{i},'uni',0);
    accm_pks_LOCS{i}=cellfun(@(x,y,z) z(ismember(x,y)),riseLOCS{i},upLOCS{i},pksLOCS{i},'uni',0);
    ind_PKS{i}=cellfun(@(x,y,z) z(~ismember(x,y)),riseLOCS{i},upLOCS{i},PKS{i},'uni',0);
    accm_PKS{i}=cellfun(@(x,y,z) z(ismember(x,y)),riseLOCS{i},upLOCS{i},PKS{i},'uni',0);
end
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section 2 Plot all cells or low responsive cells
% ctrl+enter to run this section

for i=1:2
    load('PFcriteria.mat');
    if i==1; Criteria_ID=logical(PFcriteria+PFcriteriaL); %PFcriteria for low responsive and PFcriteriaL for high
    else Criteria_ID=logical(PFcriteria); %PFcriteria for low responsive and PFcriteriaL for high
    end
    increvalue=0.1;
    [co] = return_colorbrewer('Dark2', 8);
    co=distinguishable_colors(50);
    co=lines;
    cnum=[1 2 4 5];
    amps=cell(size(condlist));

    y=[];
    counts=cell(size(condlist));
    for cond=1:max(condlist)
        condidx=find(condlist==cond)';
        ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
        ind_std_cell=arrayfun(@(a) nanstd(smoothBC_signal{a},0,3),condidx,'UniformOutput',0);

        ind_activ{cond}=horzcat(ind_activ_cell{:});
        ind_std{cond}=horzcat(ind_std_cell{:});
    end
    T=time{1}-time{1}(data(1).stim(1));

    figure; hold on
    for cond=1:max(condlist)
        ave_activ{cond}=nanmean(ind_activ{cond}(:,Criteria_ID),2);
        std_activ{cond}=nanstd(ind_activ{cond},0,2)./sqrt(size(ind_activ{cond}(:,Criteria_ID),2));
    end

    co=lines;
    co(6,:)=[.5 .5 .5];
    cnum=[1 2 4 5];
    cnum=[6 6 6];% when ploting WT tetanus group
    delay=1.6; %delay of each spike
    delay=1.9; %delay of each spike

    display_window=[-0.2 1.3]; % display of the time window for each spike
    display_window=[-0.3 1.5]; % display of the time window for each spike

    window_idx=T>display_window(1) & T<display_window(2);
    linestyle={'-' '-' '-'};

    % draw stim based on data
    stimT=[T(data(1).stim)';T(data(1).stim)'];
    stimT=stimT(:);
    Stim_signal=[0 2 2 0]';

    anaT=[0.0 0.2;0.0 0.2];
    anaT=anaT(:);
    ana_signal=[-1 2 2 -1]';

    for cond=1:max(condlist)
        %plotting stimulus
        b=area(stimT+(cond-1).*delay,Stim_signal,...
            'FaceColor',[255 96 0]/255,...
            'EdgeColor','none',...
            'FaceAlpha',0.4);

        %plotting analysis window
        b=area(anaT+(cond-1).*delay,Stim_signal,...
            'FaceColor',[249, 191, 69]./255,...
            'EdgeColor','none',...
            'FaceAlpha',0.4);

        %plotting error
        incre=increvalue.*(1:size(ind_activ{cond},2));
        ErrArea_Smooth(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),...
            std_activ{cond}(window_idx),...
            [co(cnum(cond),:) .5]);
    end

    for cond=1:max(condlist)
        incre=increvalue.*(1:size(ind_activ{cond},2));

        hold on
        plot(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),linestyle{cond},...
            'color',[co(cnum(cond),:)],...
            'LineWidth',1.5);
    end
    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.2 4.5])
    xlim([-0.5 5.5])
    ylim([0.0 0.15]+0.086)
    ylim([0.0 0.15]+0.02)
    ylim([0.0 0.4])

    xlabel('Time (sec)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end
