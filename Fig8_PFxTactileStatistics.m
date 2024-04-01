% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Fig8 PF x tactile interaction"
% 3. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section

clear
co=lines;
co(6,:)=[.5 .5 .5];
[gc]=grad_co;
co(3,:)=gc(2,:,2);
cnum=[6 2 4 5];
legendlist=[];
colorlist=[];
markersize=30;
criteria=[0 100];
legendlist={'WT control' 'WT figure' 'WT figure'};
groupnum=3;
clear GOF

load('PF.mat');
con_data=location2;      

load('paw.mat');
tet_data=location1;



%% Section 2: plot  calcium probability
% ctrl+enter to run this section

% pre_c=prctile(tet_data(1,:),75);
pre_c=0.18;
PFcriteria=con_data(2,:)>con_data(1,:)*1 & con_data(3,:)>con_data(1,:)*1 & tet_data(1,:)<pre_c;%normalize
PFcriteriaL=con_data(2,:)>con_data(1,:)*1 & con_data(3,:)>con_data(1,:)*1 & tet_data(1,:)>=pre_c;%normalize

Fcon_data=tet_data(:,PFcriteriaL);
con_mean=nanmean(Fcon_data,2);
con_sem=nanstd(Fcon_data,0,2)./sqrt(size(Fcon_data,2));
matconncat=[Fcon_data];
concatgroup=ones(size(Fcon_data));

Ftet_data=tet_data(:,PFcriteria);
tet_mean=nanmean(Ftet_data,2);
tet_sem=nanstd(Ftet_data,0,2)./sqrt(size(Ftet_data,2));
matconncat=[matconncat Ftet_data];
concatgroup=[concatgroup 2.*ones(size(Ftet_data))];

% repeated ANOVA stat
data=matconncat(1:end);
subject=repmat(1:length(matconncat),3,1);
subject=subject(:)';
g1=repmat(1:3,1,length(matconncat));%time
g2=concatgroup(1:end);%genotype

genotype_cat=categorical(concatgroup(1,:)');
t2 = table(genotype_cat, matconncat(1, :)', matconncat(2, :)', matconncat(3, :)',...
    'VariableNames', {'genotype', 'pre', 'early', 'late'});
Time = [1 2 3]; % 0: pre, 1: during, 2: post
rm2 = fitrm(t2, 'pre-late ~genotype', ...
    'WithinModel', Time, 'WithinModel', 'separatemeans');
ranovatbl =ranova(rm2);
between_sbj=multcompare(rm2,'genotype','By','Time', 'ComparisonType', 'tukey-kramer');
within_sbj=multcompare(rm2,'Time','By','genotype', 'ComparisonType', 'tukey-kramer');

% plot individual cells

x=repmat((1:3)',1,groupnum);
y=[con_mean tet_mean];
ysem=[con_sem tet_sem];

gcomb=[2];% combinations of groups to compare
symbole_pattern={'o';'o';'^--'};
for i=1:size(gcomb,1)
    figure;subplot(1,3,[2 3]);hold on;


         MarkerSize=100;
         space=0.1;
        arrayfun(@(a) errorbar((1:3)'-space,y(:,a),ysem(:,a),symbole_pattern{a},...
        'color',co(cnum(a),:),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),gcomb(i,:),'uni',0);
         arrayfun(@(a) plot((1:3)'+space,matconncat(:,concatgroup(1,:)==a),...
             '-','color',[co(cnum(a),:) 0.3],'LineWidth',1),gcomb(i,:),'Uni',0);
         arrayfun(@(a) scatter((1:3)'+space,matconncat(:,concatgroup(1,:)==a),MarkerSize,'o',...
             'MarkerFaceColor','flat',...
             'MarkerEdgeColor','flat',...
             'MarkerFaceAlpha',0.3,...
             'MarkerEdgeAlpha',0.3,...
             'CData',co(cnum(a),:)),gcomb(i,:),'Uni',0);

    % plot significance of anovan over time
    Ysig=0.55;Yincre=0.04;sigincre=0.005;

    for g=1%plot significance of each group
        idx=within_sbj.genotype==num2str(gcomb(i,g));%find index of combinations between time
        sig(1,:)=[1 2 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==2)];% x1 x2 and p-value
        sig(2,:)=[1 3 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==3)];% x1 x2 and p-value
        sig(3,:)=[2 3 within_sbj.pValue(idx & within_sbj.Time_1==2 & within_sbj.Time_2==3)];% x1 x2 and p-value

        sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
    end

    % config
    xticks(1:4)
    xticklabels({'Pre','Early','Late'});
    ylim([0 0.6])%1stAmp
    yticks(0:0.3:2)
    xlim([0.8 3.2]-0.05)
    ylabel('Amplitude (\DeltaF/F)')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end
%% Section 3: plot correlations among pre, early, and late
% ctrl+enter to run this section

x=repmat((1:3)',1,groupnum);
y=[con_mean tet_mean];
ysem=[con_sem tet_sem];

gcomb=[2];% combinations of groups to compare
symbole_pattern={'o';'o';'^--'};
for i=1:size(gcomb,1)
    figure;subplot(1,3,[2 3]);hold on;

    plot([0 5],[0 5],'k')
    plot([0.18 0.18],[0 5],'k--')
    plot([0 5],[0.18 0.18],'k--')
    MarkerSize=100;
    space=0.1;

    arrayfun(@(a) scatter(matconncat(1,concatgroup(1,:)==a),matconncat(2,concatgroup(1,:)==a),MarkerSize,'o',...
        'MarkerFaceColor','flat',...
        'MarkerEdgeColor','flat',...
        'MarkerFaceAlpha',0.2,...
        'MarkerEdgeAlpha',0.2,...
        'CData',co(cnum(a),:)),gcomb(i,:),'Uni',0);

    arrayfun(@(a) scatter(matconncat(1,concatgroup(1,:)==a),matconncat(3,concatgroup(1,:)==a),MarkerSize,'o',...
        'MarkerFaceColor','flat',...
        'MarkerEdgeColor','flat',...
        'MarkerFaceAlpha',0,...
        'MarkerEdgeAlpha',1,...
        'CData',co(cnum(a),:),...
        'LineWidth',0.9),gcomb(i,:),'Uni',0);

    % config
    lim=[0 0.9];
    lim=[0 0.55];
    ylim(lim)
    xlim([0 0.2])
    ticks=0:0.4:2;
    ticks=0:0.3:2;
    xticks([0:0.09:0.18])
    yticks([0 0.18 0.4])
    axis square
    ylabel('Post (\DeltaF/F)')
    xlabel('Pre (\DeltaF/F)')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% Section 4: plot normalized calcium amplitude
% ctrl+enter to run this section

pre_c=0.18;
PFcriteria=con_data(2,:)>con_data(1,:)*1 & con_data(3,:)>con_data(1,:)*1 & tet_data(1,:)<pre_c;%normalize
PFcriteriaL=con_data(2,:)>con_data(1,:)*1 & con_data(3,:)>con_data(1,:)*1 & tet_data(1,:)>=pre_c;%normalize

Fcon_data=tet_data(:,PFcriteriaL);
con_mean=nanmean(Fcon_data./Fcon_data(1,:),2);
con_sem=nanstd(Fcon_data./Fcon_data(1,:),0,2)./sqrt(size(Fcon_data,2));
matconncat=[Fcon_data./Fcon_data(1,:)];
concatgroup=ones(size(Fcon_data));

Ftet_data=tet_data(:,PFcriteria);
% Ftet_data(:,12)=[];
tet_mean=nanmean(Ftet_data./Ftet_data(1,:),2);
tet_sem=nanstd(Ftet_data./Ftet_data(1,:),0,2)./sqrt(size(Ftet_data,2));
matconncat=[matconncat Ftet_data./Ftet_data(1,:)];
concatgroup=[concatgroup 2.*ones(size(Ftet_data))];


% repeated ANOVA stat
data=matconncat(1:end);
subject=repmat(1:length(matconncat),3,1);
subject=subject(:)';
g1=repmat(1:3,1,length(matconncat));%time
g2=concatgroup(1:end);%genotype

genotype_cat=categorical(concatgroup(1,:)');
t2 = table(genotype_cat, matconncat(1, :)', matconncat(2, :)', matconncat(3, :)',...
    'VariableNames', {'genotype', 'pre', 'early', 'late'});
Time = [1 2 3]; % 0: pre, 1: during, 2: post
rm2 = fitrm(t2, 'pre-late ~genotype', ...
    'WithinModel', Time, 'WithinModel', 'separatemeans');
ranovatbl =ranova(rm2);
between_sbj=multcompare(rm2,'genotype','By','Time', 'ComparisonType', 'tukey-kramer');
within_sbj=multcompare(rm2,'Time','By','genotype', 'ComparisonType', 'tukey-kramer');

% plot normalized individual cells

x=repmat((1:3)',1,groupnum);
y=[con_mean tet_mean];
ysem=[con_sem tet_sem];

gcomb=[2];% combinations of groups to compare
symbole_pattern={'o';'o';'^--'};
for i=1:size(gcomb,1)
    figure;subplot(1,3,[2 3]);hold on;


         MarkerSize=100;
         space=0.1;
        arrayfun(@(a) errorbar((1:3)'-space,y(:,a),ysem(:,a),symbole_pattern{a},...
        'color',co(cnum(a),:),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),gcomb(i,:),'uni',0);
         arrayfun(@(a) plot((1:3)'+space,matconncat(:,concatgroup(1,:)==a),...
             '-','color',[co(cnum(a),:) 0.3],'LineWidth',1),gcomb(i,:),'Uni',0);
         arrayfun(@(a) scatter((1:3)'+space,matconncat(:,concatgroup(1,:)==a),MarkerSize,'o',...
             'MarkerFaceColor','flat',...
             'MarkerEdgeColor','flat',...
             'MarkerFaceAlpha',0.3,...
             'MarkerEdgeAlpha',0.3,...
             'CData',co(cnum(a),:)),gcomb(i,:),'Uni',0);

    % plot significance of anovan over time
    Ysig=7.5;Yincre=1.1;sigincre=0.15;

    for g=1%plot significance of each group
        idx=within_sbj.genotype==num2str(gcomb(i,g));%find index of combinations between time
        sig(1,:)=[1 2 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==2)];% x1 x2 and p-value
        sig(2,:)=[1 3 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==3)];% x1 x2 and p-value
        sig(3,:)=[2 3 within_sbj.pValue(idx & within_sbj.Time_1==2 & within_sbj.Time_2==3)];% x1 x2 and p-value

        sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
    end


    % config
    xticks(1:4)
    xticklabels({'Pre','Early','Late'});
    ylim([0.6 10])%1stAmp
    yticks([0.7 1 5])
    xlim([0.8 3.2]-0.05)
set(gca, 'YScale', 'log')
    ylabel('Normalized amplitude')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

