% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Fig7 and Supp Fig7 8 Airpuff_map plasticity"
% 3. analyze the data within the subfolder of "Fig7 and Supp Fig7 8 Airpuff_map plasticity"
%  which is WT control, WT tetanus, SK2 tetanus, or CaMKII tetanus
% 4. click "Run" to run the whole script

%% Section 1: Open files 
% Select SpikeAnalist file from the pre- (e.g., WTcontrol_pre_SpikeAnalist.mat),
% early- (e.g., WTcontrol_early_SpikeAnalist.mat), to late-tetanus
% (e.g., WTcontrol_late_SpikeAnalist) conditions.
% ctrl+enter to run this section

clear
FN=[];
FP=[];
FNlist={};

FPlist={};
condlist=[];

f=1;
while f
    [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeDatalist files', 'Multiselect', 'on');
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
    %%
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
clear data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

% pull the data out
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section 2: plot average calcium traces of each stimulus conditions
% ctrl+enter to run this section

% From left to right, plot from pre- to late-tetanus conditions 
%%
LOCS=pksLOCS;
ana_window=[0 0.2];
increvalue=0.1;

[co]=grad_co;
cnum=[2 1;2 2;4 1;4 2;3 1;3 2];
cnum=[1 2;1 2;1 2];

amps=cell(size(condlist));

%
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
catmean=[];
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    mean_activ=[];
    amps=[];
    amps_idx=[];
    cellmean=[];
    for i=1:length(condidx)
        idx=condidx(i);
        incre=increvalue.*(1:size(LOCS{idx},2));
        for roi=1:size(LOCS{idx},2)
            % aveF
            % because the histcounts include left left edge, so
            % convert the sign, and then left-right flip.
            % because the histcounts include the right edge of the
            % last bin, so add one more bin at right
            counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
            counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1});

            activ=smoothBC_signal{idx}(:,roi,:);
            mean_activ{idx}(:,roi,:)=nanmean(activ(:,:,counts{idx}{roi,1}),3);


            ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

            amps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>ana_window(1) & T(ana_idx{a})<=ana_window(2)),:,a),...
                1:length(ana_idx),'uni',0);%filter out the data out of the ana_window

            amps_idx{idx}{roi,1} = arrayfun(@(a) find(amps{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%find all the index of amplitude
            amps_idx{idx}{roi,1} = arrayfun(@(a) min(amps_idx{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%find the first index
            amps{idx}{roi,1} = arrayfun(@(a) amps{idx}{roi,1}{a}(amps_idx{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%first amplitude within the ana_window

            cellmean{idx}{roi,1} = nanmean(vertcat(amps{idx}{roi}{:}),1);
        end
    end
    cellmean=vertcat(cellmean{condidx});
    catmean{cond}=horzcat(cellmean{:});
    ind_activ{cond}=horzcat(mean_activ{:});
end

for cond=1:max(condlist)
    act=1:size(ind_activ{cond},2);
%     act([28 8])=[];
    ave_activ{cond}=nanmedian(ind_activ{cond}(:,act)./catmean{1}(act),2);
    std_activ{cond}=nanstd(ind_activ{cond}(:,act)./catmean{1}(act),0,2)./sqrt(size(ind_activ{cond}(:,act),2));

%     std_activ{cond} = mad(ind_activ{cond}(:,act)./catmean{1}(act),1,2);
end

figure;hold on
delay=1.6; %delay of each spike
delay=1.9; %delay of each spike

display_window=[-0.2 1.3]; % display of the time window for each spike
display_window=[-0.25 1.5]; % display of the time window for each spike

window_idx=T>display_window(1) & T<display_window(2);
linestyle={'-' '-' '-'};

% draw stim based on data
stimT=[T(data(1).stim)';T(data(1).stim)'];
stimT=stimT(:);
Stim_signal=[-2 2 2 -2]';

anaT=[ana_window;ana_window];
anaT=anaT(:);
ana_signal=[0 2 2 0]';

hold on
for cond=1:max(condlist)
    %plotting stimulus
    b=area(stimT+(cond-1).*delay,Stim_signal,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4,'ShowBaseLine','off','BaseValue',-2);

    %plotting analysis window
    b=area(anaT+(cond-1).*delay,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4,'ShowBaseLine','off','BaseValue',-2);

    %plotting error
    incre=increvalue.*(1:size(ind_activ{cond},2));
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),...
        std_activ{cond}(window_idx),...
        [co(cnum(cond,1),:,cnum(cond,2)) .5]);
end

for cond=1:max(condlist)
    %     subplot(1,5,cond)
    incre=increvalue.*(1:size(ind_activ{cond},2));

    hold on
    plot(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),linestyle{cond},...
        'color',[co(cnum(cond,1),:,cnum(cond,2))],...
        'LineWidth',1.5);

    %
    box off
    % b.BaseLine.Visible = 'off'
    xlim([-0.2 4.5])
    xlim([-0.5 5.5])
    ylim([0.085 0.235])
    ylim([0.1 1.7])
    ylim([0 1.6]-0.15)

    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
end




