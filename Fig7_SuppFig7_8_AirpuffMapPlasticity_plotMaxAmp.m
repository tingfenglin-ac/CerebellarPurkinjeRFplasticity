% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Fig7 and Supp Fig7 8 Airpuff_map plasticity"
% 3. analyze the data within the folder of "Fig7 and Supp Fig7 8 Airpuff_map plasticity\probability"
% 4. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section

clear
[co]=grad_co;
cnum=[6 2 4 5];
legendlist=[];
colorlist=[];
markersize=30;
criteria=[0 100];
legendlist={'WT finger control' 'WT wrist control' 'WT finger non-tetanized' 'WT wrist non-tetanized' 'WT finger tetanus' 'WT wrist tetanus'};
groupnum=3;
clear GOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load control data
% combine old and new collected data

%0-200ms window
% TRUE control

load('location2_OfforNoPuffTetanus.mat');
loc2_data=location2;
load('location1_OfforNoPuffTetanus.mat');
loc1_data=location1;
con_data=mean(cat(3,loc2_data,loc1_data),3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% untetanized area
load('location2.mat');
offloc2_data=location2;
load('location1.mat');
offloc1_data=location1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load location2 tetanus
induced=[];
load('location2_tet.mat');
tet_loc2=induced;

%%%%%%%%%%%%%%%%
% load location1 tetanus
load('location1_tet.mat');
% induced(:,28)=[];
tet_loc1=induced;

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% CaMKII

% untetanized area
load('camk2_location2_off.mat');
Coffloc2_data=location2;
load('camk2_location1_off.mat');
Coffloc1_data=location1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load location2 tetanus
induced=[];
load('camk2_location2_tet.mat');
Ctet_loc2=induced;

%%%%%%%%%%%%%%%%
% load location1 tetanus
load('camk2_location1_tet.mat');
Ctet_loc1=induced;

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% SK2
% untetanized area
load('sk2_location2_off.mat');
Soffloc2_data=location2;
load('sk2_location1_off.mat');
Soffloc1_data=location1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load location2 tetanus
induced=[];
load('sk2_location2_tet.mat');
Stet_loc2=induced;

%%%%%%%%%%%%%%%%
% load location1 tetanus
load('sk2_location1_tet.mat');
Stet_loc1=induced;

%% Section 2: plot original calcium probability of each genotype before airpuff tetanization
% ctrl+enter to run this section

combine_loc2_data=[loc2_data offloc2_data tet_loc2];
combine_loc2_mean=nanmean(combine_loc2_data,2);
combine_loc2_sem=nanstd(combine_loc2_data,0,2)./sqrt(size(combine_loc2_data,2));
matconncat=[combine_loc2_data];
concatgroup=[1.*ones(size(combine_loc2_data))];

combine_loc1_data=[loc1_data offloc1_data tet_loc1];
combine_loc1_mean=nanmean(combine_loc1_data,2);
combine_loc1_sem=nanstd(combine_loc1_data,0,2)./sqrt(size(combine_loc1_data,2));
matconncat=[matconncat combine_loc1_data];
concatgroup=[concatgroup 2.*ones(size(combine_loc1_data))];

Ccombine_loc2_data=[Coffloc2_data Ctet_loc2];
Ccombine_loc2_mean=nanmean(Ccombine_loc2_data,2);
Ccombine_loc2_sem=nanstd(Ccombine_loc2_data,0,2)./sqrt(size(Ccombine_loc2_data,2));
matconncat=[matconncat Ccombine_loc2_data];
concatgroup=[concatgroup 3.*ones(size(Ccombine_loc2_data))];

Ccombine_loc1_data=[Coffloc1_data Ctet_loc1];
Ccombine_loc1_mean=nanmean(Ccombine_loc1_data,2);
Ccombine_loc1_sem=nanstd(Ccombine_loc1_data,0,2)./sqrt(size(Ccombine_loc1_data,2));
matconncat=[matconncat Ccombine_loc1_data];
concatgroup=[concatgroup 4.*ones(size(Ccombine_loc1_data))];

Scombine_loc2_data=[Soffloc2_data Stet_loc2];
Scombine_loc2_mean=nanmean(Scombine_loc2_data,2);
Scombine_loc2_sem=nanstd(Scombine_loc2_data,0,2)./sqrt(size(Scombine_loc2_data,2));
matconncat=[matconncat Scombine_loc2_data];
concatgroup=[concatgroup 5.*ones(size(Scombine_loc2_data))];

Scombine_loc1_data=[Soffloc1_data Stet_loc1];
Scombine_loc1_mean=nanmean(Scombine_loc1_data,2);
Scombine_loc1_sem=nanstd(Scombine_loc1_data,0,2)./sqrt(size(Scombine_loc1_data,2));
matconncat=[matconncat Scombine_loc1_data];
concatgroup=[concatgroup 6.*ones(size(Scombine_loc1_data))];


% plot
% 1-way ANOVA stat
data=matconncat(1,concatgroup(1,:)==1|concatgroup(1,:)==2|concatgroup(1,:)==3|concatgroup(1,:)==4|concatgroup(1,:)==5|concatgroup(1,:)==6);
% g1=repmat(1:3,1,length(matconncat(concatgroup==6|concatgroup==3|concatgroup==5)));%time
g2=concatgroup(1,concatgroup(1,:)==1|concatgroup(1,:)==2|concatgroup(1,:)==3|concatgroup(1,:)==4|concatgroup(1,:)==5|concatgroup(1,:)==6);%genotype
[p.anova1,~,stats.anova1] = anova1(data,g2);
c.anova1=multcompare(stats.anova1);
    

cnum=[1 2;1 1;4 2;4 1;3 2;3 1];
y=[combine_loc2_mean(1) combine_loc1_mean(1) Ccombine_loc2_mean(1) Ccombine_loc1_mean(1) Scombine_loc2_mean(1) Scombine_loc1_mean(1)];
ysem=[combine_loc2_sem(1) combine_loc1_sem(1) Ccombine_loc2_sem(1) Ccombine_loc1_sem(1) Scombine_loc2_sem(1) Scombine_loc1_sem(1)];
figure;
    hold on;

    MarkerSize=100;
    meanSize=10;
    space=0.2;
    jrange=0.3;
    symbole_pattern={'^--';'o-';'^--';'o-';'^--';'o-'};
    arrayfun(@(a) errorbar(a+space,y(a),ysem(a),symbole_pattern{a},...
        'color',co(cnum(a,1),:,cnum(a,2)),...
        'CapSize',30,...
        'MarkerSize',meanSize,...
        'LineWidth',3),1:6,'uni',0);
    arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.3,num2str(y(a),'%4.2f'),...
        'color',co(cnum(a,1),:,cnum(a,2)),...
        'HorizontalAlignment','center',...
        'FontSize',15),1:6,'uni',0);

    compgroup=[1:6];
    symbole_pattern={'^';'o';'^';'o';'^';'o'};
    arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,symbole_pattern{a},...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a,1),:,cnum(a,2))),1:6,'uni',0);

    %significance
    sigX=1;
    sigY=0.04;
    dsigY=0.04;
    formatSpec = '%.3f';
    FS=16;

        Ysig=1.2;Yincre=0.1;sigincre=0.008;

            sigidx=find(c.anova1(:,end)<0.05 & c.anova1(:,end)>=0.01)';
            if sigidx
                for sn=sigidx
                    plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                    text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
                        'HorizontalAlignment','center',...
                        'color',[0 0 0]);
                    Ysig=Ysig+Yincre;
                end
            end
            sigidx=find(c.anova1(:,end)<0.01 & c.anova1(:,end)>=0.001)';
            if sigidx
                for sn=sigidx
                    plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                    text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
                        'HorizontalAlignment','center',...
                        'color',[0 0 0]);
                    Ysig=Ysig+Yincre;
                end
            end
            sigidx=find(c.anova1(:,end)<0.001)';
            if sigidx
                for sn=sigidx
                    plot(c.anova1(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                    text(mean(c.anova1(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
                        'HorizontalAlignment','center',...
                        'color',[0 0 0]);
                    Ysig=Ysig+Yincre;
                end
            end

    % config
%     legend(legendlist)
    xticks(1:6)
    yticks(0:1:5)
    xticklabels([]);
%     xticklabels({'WT','SK2 KO','CaMKII TT305/6VA'});
%     xtickangle(15)
    ylim([0 2])
%     yticks(0:0.2:0.4)
    xlim([0.5 6.5])
    ylabel('Amplitude (\DeltaF/F)')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])





%% Section 3-1: statistics
% ctrl+enter to run this section
NZcon_data=con_data(:,con_data(1,:)>0);%remove pre is zero
Ncon_data=NZcon_data./NZcon_data(1,:);%normalize
con_mean=nanmean(Ncon_data,2);
con_sem=nanstd(Ncon_data,0,2)./sqrt(size(Ncon_data,2));
matconncat=[Ncon_data];
concatgroup=ones(size(Ncon_data));

NZoffloc2_data=offloc2_data(:,offloc2_data(1,:)>0);%remove pre is zero
Noffloc2_data=NZoffloc2_data./NZoffloc2_data(1,:);%normalize
offloc2_mean=nanmean(Noffloc2_data,2);
offloc2_sem=nanstd(Noffloc2_data,0,2)./sqrt(size(Noffloc2_data,2));
matconncat=[matconncat Noffloc2_data];
concatgroup=[concatgroup 2.*ones(size(Noffloc2_data))];

NZoffloc1_data=offloc1_data(:,offloc1_data(1,:)>0);%remove pre is zero
Noffloc1_data=NZoffloc1_data./NZoffloc1_data(1,:);%normalize
offloc1_mean=nanmean(Noffloc1_data,2);
offloc1_sem=nanstd(Noffloc1_data,0,2)./sqrt(size(Noffloc1_data,2));
matconncat=[matconncat Noffloc1_data];
concatgroup=[concatgroup 3.*ones(size(Noffloc1_data))];

NZtet2_data=tet_loc2(:,tet_loc2(1,:)>0);%remove pre is zero
Ntet2_data=NZtet2_data./NZtet2_data(1,:);%normalize
tet2_mean=nanmean(Ntet2_data,2);
tet2_sem=nanstd(Ntet2_data,0,2)./sqrt(size(Ntet2_data,2));
matconncat=[matconncat Ntet2_data];
concatgroup=[concatgroup 4.*ones(size(Ntet2_data))];

NZtet1_data=tet_loc1(:,tet_loc1(1,:)>0);%remove pre is zero
Ntet1_data=NZtet1_data./NZtet1_data(1,:);%normalize
tet1_mean=nanmean(Ntet1_data,2);
tet1_sem=nanstd(Ntet1_data,0,2)./sqrt(size(Ntet1_data,2));
matconncat=[matconncat Ntet1_data];
concatgroup=[concatgroup 5.*ones(size(Ntet1_data))];
%CaMKII
NZCoffloc2_data=Coffloc2_data(:,Coffloc2_data(1,:)>0);%remove pre is zero
NCoffloc2_data=NZCoffloc2_data./NZCoffloc2_data(1,:);%normalize
Coffloc2_mean=nanmean(NCoffloc2_data,2);
Coffloc2_sem=nanstd(NCoffloc2_data,0,2)./sqrt(size(NCoffloc2_data,2));
matconncat=[matconncat NCoffloc2_data];
concatgroup=[concatgroup 6.*ones(size(NCoffloc2_data))];

NZCoffloc1_data=Coffloc1_data(:,Coffloc1_data(1,:)>0);%remove pre is zero
NCoffloc1_data=NZCoffloc1_data./NZCoffloc1_data(1,:);%normalize
Coffloc1_mean=nanmean(NCoffloc1_data,2);
Coffloc1_sem=nanstd(NCoffloc1_data,0,2)./sqrt(size(NCoffloc1_data,2));
matconncat=[matconncat NCoffloc1_data];
concatgroup=[concatgroup 7.*ones(size(NCoffloc1_data))];

NZCtet2_data=Ctet_loc2(:,Ctet_loc2(1,:)>0);%remove pre is zero
NCtet2_data=NZCtet2_data./NZCtet2_data(1,:);%normalize
Ctet2_mean=nanmean(NCtet2_data,2);
Ctet2_sem=nanstd(NCtet2_data,0,2)./sqrt(size(NCtet2_data,2));
matconncat=[matconncat NCtet2_data];
concatgroup=[concatgroup 8.*ones(size(NCtet2_data))];

NZCtet1_data=Ctet_loc1(:,Ctet_loc1(1,:)>0);%remove pre is zero
NCtet1_data=NZCtet1_data./NZCtet1_data(1,:);%normalize
Ctet1_mean=nanmean(NCtet1_data,2);
Ctet1_sem=nanstd(NCtet1_data,0,2)./sqrt(size(NCtet1_data,2));
matconncat=[matconncat NCtet1_data];
concatgroup=[concatgroup 9.*ones(size(NCtet1_data))];
%SK2
NZSoffloc2_data=Soffloc2_data(:,Soffloc2_data(1,:)>0);%remove pre is zero
NSoffloc2_data=NZSoffloc2_data./NZSoffloc2_data(1,:);%normalize
Soffloc2_mean=nanmean(NSoffloc2_data,2);
Soffloc2_sem=nanstd(NSoffloc2_data,0,2)./sqrt(size(NSoffloc2_data,2));
matconncat=[matconncat NSoffloc2_data];
concatgroup=[concatgroup 10.*ones(size(NSoffloc2_data))];

NZSoffloc1_data=Soffloc1_data(:,Soffloc1_data(1,:)>0);%remove pre is zero
NSoffloc1_data=NZSoffloc1_data./NZSoffloc1_data(1,:);%normalize
Soffloc1_mean=nanmean(NSoffloc1_data,2);
Soffloc1_sem=nanstd(NSoffloc1_data,0,2)./sqrt(size(NSoffloc1_data,2));
matconncat=[matconncat NSoffloc1_data];
concatgroup=[concatgroup 11.*ones(size(NSoffloc1_data))];

NZStet2_data=Stet_loc2(:,Stet_loc2(1,:)>0);%remove pre is zero
NStet2_data=NZStet2_data./NZStet2_data(1,:);%normalize
Stet2_mean=nanmean(NStet2_data,2);
Stet2_sem=nanstd(NStet2_data,0,2)./sqrt(size(NStet2_data,2));
matconncat=[matconncat NStet2_data];
concatgroup=[concatgroup 12.*ones(size(NStet2_data))];

NZStet1_data=Stet_loc1(:,Stet_loc1(1,:)>0);%remove pre is zero
NStet1_data=NZStet1_data./NZStet1_data(1,:);%normalize
Stet1_mean=nanmean(NStet1_data,2);
Stet1_sem=nanstd(NStet1_data,0,2)./sqrt(size(NStet1_data,2));
matconncat=[matconncat NStet1_data];
concatgroup=[concatgroup 13.*ones(size(NStet1_data))];

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


%% Section 3-2; plot all
% ctrl+enter to run this section

for fw=1:2
cnum=[1 2;2 2;2 2;2 1;2 1;4 3;4 3;4 1;4 1;3 3;3 3;3 1;3 1];
symbole_pattern={'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-'};


    x=repmat((1:3)',1,groupnum);
    y=[con_mean offloc2_mean offloc1_mean tet2_mean tet1_mean...
        Coffloc2_mean Coffloc1_mean Ctet2_mean Ctet1_mean...
        Soffloc2_mean Soffloc1_mean Stet2_mean Stet1_mean];
    ysem=[con_sem offloc2_sem offloc1_sem tet2_sem tet1_sem...
        Coffloc2_sem Coffloc1_sem Ctet2_sem Ctet1_sem...
        Soffloc2_sem Soffloc1_sem Stet2_sem Stet1_sem];
    figure;subplot(1,3,[2 3]);hold on;
    if fw==1; Dataid=[1 9 13 5];%wrist
    else Dataid=[1 8 12 4];%finger
    end
    arrayfun(@(a) errorbar((1:3)',y(:,a),ysem(:,a),symbole_pattern{a},...
        'color',co(cnum(a,1),:,cnum(a,2)),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),Dataid,'uni',0);
%1 9 13 5
%1 8 12 4

    %significance
    sigX=1;
    sigY=1.6;
    dsigY=0.05;
    formatSpec = '%.3f';
    FS=16;
    if ranovatbl.pValue(1)<0.001
        text(sigX,sigY+dsigY*2,'{\it p_s_t_i_m} < 0.001','FontSize',FS)
    else
        text(sigX,sigY+dsigY*2,['{\it p_s_t_i_m} = ' num2str(ranovatbl.pValue(1),formatSpec)],'FontSize',FS)
    end
    if ranovatbl.pValue(2)<0.001
        text(sigX,sigY,'{\it p_s_t_i_m _x _g_e_n_e} < 0.001','FontSize',FS)
    else
        text(sigX,sigY,['{\it p_s_t_i_m _x _g_e_n_e} = ' num2str(ranovatbl.pValue(2),formatSpec)],'FontSize',FS)
    end

    % config
    xticks(1:6)
    xticklabels({'Pre','Early','Late'});
        ylim([0.6 1.75])
    yticks(0.6:0.4:6)
    xlim([0.8 3.5])
    ylabel('Normalized amplitude')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% Section 3-3; plot compbinations
% ctrl+enter to run this section

for fw=1:2
    if fw==1;gcomb=[1 3 5;7 9 5;11 13 5];% wrist
    else gcomb=[1 2 4;6 8 4;10 12 4];% finger
    end


for i=1:size(gcomb,1)
    figure;subplot(1,3,[2 3]);hold on;
        symbole_pattern={'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-';'^:';'o-'};
%         symbole_pattern={'o-';'^-';'o-';'^-';'o-';'^-';'o-';'^-';'o-';'^-';'o-';'^-';'o-'};

        arrayfun(@(a) errorbar((1:3)',y(:,a),ysem(:,a),symbole_pattern{a},...
        'color',co(cnum(a,1),:,cnum(a,2)),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3),gcomb(i,:),'uni',0);

        % plot significance of anovan over time
        Ysig=1.55;Yincre=0.053;sigincre=0.006;
        for g=1:3%plot significance of each group
%         for g=1:4 %plot significance of each group %when comparing 4 groups

            idx=within_sbj.genotype==num2str(gcomb(i,g));%find index of combinations between time
            sig(1,:)=[1 2 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==2)];% x1 x2 and p-value
            sig(2,:)=[1 3 within_sbj.pValue(idx & within_sbj.Time_1==1 & within_sbj.Time_2==3)];% x1 x2 and p-value
            sig(3,:)=[2 3 within_sbj.pValue(idx & within_sbj.Time_1==2 & within_sbj.Time_2==3)];% x1 x2 and p-value

            sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
            %          symbole_pattern={'-';'-';':';'-';':';'-';':';'-';':';'-';':';'-';':'};

            if sigidx
                for sn=sigidx
                    plot(sig(sn,1:2),[Ysig Ysig],symbole_pattern{gcomb(i,g)}(2),'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)),'linewidth',3)
                    text(mean(sig(sn,1:2)),Ysig+sigincre,'*','FontSize',24,...
                        'HorizontalAlignment','center',...
                        'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)));
                    Ysig=Ysig+Yincre;
                end
            end
            sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
            if sigidx
                for sn=sigidx
                    plot(sig(sn,1:2),[Ysig Ysig],symbole_pattern{gcomb(i,g)}(2),'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)),'linewidth',3)
                    text(mean(sig(sn,1:2)),Ysig+sigincre,'**','FontSize',24,...
                        'HorizontalAlignment','center',...
                        'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)));
                    Ysig=Ysig+Yincre;
                end
            end
            sigidx=find(sig(:,end)<0.001)';
            if sigidx
                for sn=sigidx
                    plot(sig(sn,1:2),[Ysig Ysig],symbole_pattern{gcomb(i,g)}(2),'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)),'linewidth',3)
                    text(mean(sig(sn,1:2)),Ysig+sigincre,'***','FontSize',24,...
                        'HorizontalAlignment','center',...
                        'color',co(cnum(gcomb(i,g),1),:,cnum(gcomb(i,g),2)));
                    Ysig=Ysig+Yincre;
                end
            end
        end

        % plot significance of anovan between groups
        compa_group=[1 2;1 3;1 4;2 3;2 4;3 4]; %when comparing 4 groups
        compa_group=[1 2;1 3;2 3];

        for t=2:3%plot significance of each time
        xsin=0.16;xincre=0.15;xsigincre=0.125;
            for cp=1:size(compa_group,1)
                cp1=compa_group(cp,1);
                cp2=compa_group(cp,2);
                idx=find(between_sbj.genotype_1==num2str(gcomb(i,cp1)) &...
                    between_sbj.genotype_2==num2str(gcomb(i,cp2)) &...
                    between_sbj.Time==t);%find index of combinations between time
                bw_pvalue=between_sbj.pValue(idx);


                if bw_pvalue<0.05 & bw_pvalue>=0.01
                    if bw_pvalue~=0
                        plot([t t]+xsin,y(t,gcomb(i,[cp1 cp2])),'color',[.3 .3 .3],'linewidth',3)
                        tx=text(t+xsin+xsigincre,mean(y(t,gcomb(i,[cp1 cp2]))),'*','FontSize',24,...
                            'HorizontalAlignment','center','Color',[.3 .3 .3])
                        tx.Rotation=90;
                        xsin=xsin+xincre;
                    end
                end
                if bw_pvalue<0.01 & bw_pvalue>=0.001
                    if bw_pvalue~=0
                        plot([t t]+xsin,y(t,gcomb(i,[cp1 cp2])),'color',[.3 .3 .3],'linewidth',3)
                        tx=text(t+xsin+xsigincre,mean(y(t,gcomb(i,:))),'**','FontSize',24,...
                            'HorizontalAlignment','center','Color',[.3 .3 .3])
                        tx.Rotation=90;
                        xsin=xsin+xincre;

                    end
                end
                if bw_pvalue<0.001
                    if bw_pvalue~=0
                        plot([t t]+xsin,y(t,gcomb(i,[cp1 cp2])),'color',[.3 .3 .3],'linewidth',3)
                        tx=text(t+xsin+xsigincre,mean(y(t,gcomb(i,:))),'***','FontSize',24,...
                            'HorizontalAlignment','center','Color',[.3 .3 .3])
                        tx.Rotation=90;
                        xsin=xsin+xincre;
                    end
                end
            end
        end


    % config
    %         legend(legendlist{gcomb(i,:)})
    xticks(1:6)
    xticklabels({'Pre','Early','Late'});
    ylim([0.8 1.2])%spon
        ylim([0.6 1.75])
    yticks(0.5:1:6)
        yticks(0.6:0.4:6)
    xlim([0.8 3.5])
    ylabel('Normalized amplitude')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

end
