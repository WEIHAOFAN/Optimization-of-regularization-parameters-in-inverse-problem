clear
close all
clc


wl = 850;
nn = [3 6 9 12 15 18 24 30];
path = '/data/eggebrecht/data1/Weihao/CSF/semi_homo/new/voxel_1mm';
d1 = zeros(8,8);
d2 = zeros(8,1);
mdx = [0.43,0.51,0.56,0.63,0.68,16.25,0.174648,0.4844617,0.381169];
mdy = [0.04,0.03,0.06,0.12,0.27,0.03,0.0848,0.056155];

for m=1:8
    if m<5
        frequency = 0;
    else
        frequency = 300;
    end
    load([path,'/',num2str(nn(m)),'/',num2str(frequency),'/lambda',num2str(wl),'_noise_',num2str(nn(m)),'_',num2str(frequency),'MHz.mat'])
    color = [
             0.9 0.1 0.1;
             0.9 0.6 0.1;
             0.9 0.9 0.1;
             0.1 0.9 0.1;
             0.1 0.8 0.9; 
             0.1 0.1 0.9;
             0.7 0.1 0.7;           
             ];
    mark = ['o','s','d','h','v','p','>'];
    nsr = 1./snr;
    if m==1
        ifov = (17*34^2)./fov;
    else
        ifov = (17*44^2)./fov;
    end

    %% Lambda1
    fw = fw(1:7,1:7);
    nsr = nsr(1:7,1:7);
    fw = fw./max(fw(:));
    nsr = nsr./max(nsr(:));
    fig1=figure;
    for i=1:7
        plot(fw(i,1:7),nsr(i,1:7),'--','LineWidth',3.5,'Color',[color(i,1),color(i,2),color(i,3)])
        hold on
    end
    for i=1:7
        for j=1:7
            plot(fw(i,j),nsr(i,j),mark(j),'MarkerSize',14,...
                'MarkerEdgeColor',[color(i,1),color(i,2),color(i,3)],...
                'MarkerFaceColor',[color(i,1),color(i,2),color(i,3)],...
                'LineWidth',2,'Color',[color(i,1),color(i,2),color(i,3)])
            hold on
            d1(i,j) = sqrt((fw(i,j)-mdx(m))^2+(nsr(i,j)-mdy(m))^2);
        end
    end
    % % l=legend('0.3','0.1','0.06','0.03','0.01','0.006','0.003','0.001','Location','West','FontSize',20);
    % % title(l,"{\lambda}_2",'FontSize',25);
    xticks(0.2:0.05:1)
%     xticklabels({})
%     yticks(0:0.2:1)
%     yticklabels({})
    % grid on
    set(gca,'FontName','Arial','fontsize',26,'LineWidth',3.5,'gridlinestyle','--')
        set(gcf,'position',[1000 1000 650 500])
    % ylabel('Normalized NSR','FontSize',15)
    % xlabel('Normalized FWHM','FontSize',15)
    % title([num2str(nn(m)),'mm ',num2str(frequency),'MHz'],'FontSize',15)
%     xlim([0.7 1.0])
%     ylim([0 1])
%     saveas(fig1,[path,'/figures/L_curve/L1_',num2str(nn(m)),'_',num2str(wl),'_',num2str(frequency),'.epsc'])

    % figure 
    % for i=1:8 
    %     plot(1:8,d1(i,:),'.-','MarkerSize',10,'LineWidth',3.5,'color',[color(i,1),color(i,2),color(i,3)])
    %     hold on 
    % end 
    % legend('0.3','0.1','0.06','0.03','0.01','0.006','0.003','0.001','FontSize',20)
    % xlabel('lambda2','FontSize',15) % ylabel('mean distance','FontSize',15)
    % xticks([1 2 3 4 5 6 7 8]) % xtickangle(45) %
    % xticklabels({'0.3','0.1','0.06','0.03','0.01','0.006','0.003','0.001'});
    % ylabel('Distance','FontSize',15)
    % title([num2str(nn(m)),'mm ',num2str(frequency),'MHz'],'FontSize',15)
    % [aa ab] = min(d1(:));
    % [ac ad] = ind2sub([8 8],ab)
end