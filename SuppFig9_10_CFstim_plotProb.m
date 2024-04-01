% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Supp Fig9 CFstimulation"
% 3. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section

clear
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 2 4 5];
legendlist=[];
colorlist=[];
markersize=30;
criteria=[0 100];
legendlist={'WT control' 'WT tetanus' 'SK2 KO' 'CaMKII TT305/6VA'};
groupnum=4;
clear GOF


load('prob.mat');

tet_data=ymat;
tet_data=tet_data(:,tet_data(1,:)>=criteria(1) & tet_data(1,:)<criteria(2));

%% Section 2: plot  calcium probability
% ctrl+enter to run this section


tet_mean=nanmean(tet_data,2);
tet_sem=nanstd(tet_data,0,2)./sqrt(size(tet_data,2));
matconncat=[tet_data];
concatgroup=ones(size(tet_data));

figure;
    hold on;

    MarkerSize=100;
    meanSize=10;
    space=0.17;
    jrange=0.3;
    errorbar(1+space,tet_mean(1),tet_sem(1),'o',...
        'color',co(6,:),...
        'CapSize',30,...
        'MarkerSize',meanSize,...
        'LineWidth',3);
    text(1+space,tet_mean(1)+tet_sem(1)+0.1,num2str(tet_mean(1),'%4.2f'),...
        'color',co(6,:),...
        'HorizontalAlignment','center',...
        'FontSize',20);
    
    scatter(1-space-0.5*jrange+jrange*rand(1,length(tet_data)),tet_data(1,:),MarkerSize,'o',...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(6,:));

        xticks(1:4)
        xticklabels([])
        yticks(0:0.5:10)
    ylim([0 1])
    xlim([0.5 1.5])
    ylabel('Probability')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
%% Section 3 plot normalized mean+sem

NZtet_data=tet_data(:,tet_data(1,:)>0);%remove pre is zero
Ntet_data=NZtet_data./NZtet_data(1,:);%normalize
tet_mean=nanmean(Ntet_data,2);
tet_sem=nanstd(Ntet_data,0,2)./sqrt(size(Ntet_data,2));
matconncat=[Ntet_data];
concatgroup=[ones(size(Ntet_data))];

% 1-way ANOVA stat
[p,tbl,stats] = anova1(matconncat');
c.anovan=multcompare(stats);



x=repmat((1:3)',1,groupnum);
y=[tet_mean];
ysem=[tet_sem];
    figure;subplot(1,3,[2 3]);hold on;
    arrayfun(@(a) errorbar((1:3)',y(:,a),ysem(:,a),'o-',...
        'color',co(2,:),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),1,'uni',0);

    % plot significance of anovan over time
        Ysig=1.1;Yincre=0.045;sigincre=0.007;
    for g=1%plot significance of each group
        idx=[1 2 3];%find index of combinations between time
        sig(1,:)=[1 2 c.anovan(c.anovan(:,1)==idx(1) & c.anovan(:,2)==idx(2),end)];% x1 x2 and p-value
        sig(2,:)=[1 3 c.anovan(c.anovan(:,1)==idx(1) & c.anovan(:,2)==idx(3),end)];% x1 x2 and p-value
        sig(3,:)=[2 3 c.anovan(c.anovan(:,1)==idx(2) & c.anovan(:,2)==idx(3),end)];% x1 x2 and p-value

        sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(2,:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(2,:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(2,:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(2,:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2),[Ysig Ysig],'color',co(2,:),'linewidth',3)
                text(mean(sig(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                    'HorizontalAlignment','center',...
                    'color',co(2,:));
                Ysig=Ysig+Yincre;
            end
        end
    end

   

    % config
    %     legend(legendlist{gcomb(i,:)},'Location','southwest')
    xticks(1:4)
    xticklabels({'Pre','Early','Late'});
    yticks(0.7:0.3:2)
    ylim([0.5 1.3])%1stAmp
    xlim([0.8 3.2])
    ylabel('Normalized probability')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])



