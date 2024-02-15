clear
close all
clc

wl = 850;
nn = [3 6 9 12 15 18 24 30];
p = [6 6 5 5 3 4 1 1];
path = '/data/eggebrecht/data1/Weihao/CSF/semi_homo/new/voxel_1mm/';

color = [
         0.9 0.9 0.1;
         0.9 0.6 0.1;
         0.1 0.9 0.1;
         0.1 0.9 0.1;
         0.1 0.8 0.9; 
         0.1 0.1 0.9;
         0.7 0.1 0.7;           
         ];
mark = ['<','_','x','+','o','^','s','d','h','v','p','>','|','*'];
for m=3
    if m < 5
        frequency = 0;
    else
        frequency = 300;
    end
    i1 = p(m); 
    load([path,num2str(nn(m)),'/',num2str(frequency),'/lambda',num2str(wl),'_more_noise_',num2str(nn(m)),'_',num2str(frequency),'MHz.mat'])
    datax = zeros(size(noise,2),1);
    datay = zeros(size(noise,2),1);
    datax(:,1) = noise(i1,:);
    datay(:,1) = fov(6,i1,:); 
    
    if m==1
        datay = (41*41*30)./datay;
    else
        datay = (41*41*31)./datay;
    end
    d2 = zeros(size(datax,1),1);

    %% Lambda2
    xx=3;
    datay = datay(1:end-xx,1);
    datax = datax(1:end-xx,1);
    datay = datay./max(datay(1:end));
    datax = datax./max(datax(1:end));
    fig2 = figure;
    plot(datax(1:end,1),100*datay(1:end,1),'--','LineWidth',4,'Color',[color(i1,1),color(i1,2),color(i1,3)])
    hold on
    for j=1:length(datax)
        plot(datax(j,1),100*datay(j,1),mark(j),'MarkerSize',16,...
             'MarkerEdgeColor',[color(i1,1),color(i1,2),color(i1,3)],...
             'MarkerFaceColor',[color(i1,1),color(i1,2),color(i1,3)],...
             'LineWidth',2,'Color',[color(i1,1),color(i1,2),color(i1,3)])
        hold on
%         d2(j,1) = sqrt((datax(j,1)-mdx(m))^2+(datay(j,1)-mdy(m))^2);
    end
%     grid on
    set(gca,'FontName','Arial','fontsize',26,'LineWidth',3.5,'gridlinestyle','--')
    set(gcf,'position',[1000 1000 650 500])
%     ylabel('Normalized FOV/GFFR %','FontSize',15)
%     xlabel('Normalized Image noise','FontSize',15)
%     title([num2str(nn(m)),'mm ',num2str(frequency),'MHz'],'FontSize',15)
    xticks(0:0.2:1)
%     xticklabels({})
%     yticks(84:4:100)
%     yticklabels({})
    xlim([0.2 1])
    ylim([97 100])
    saveas(fig2,[path,'/figures/L_curve/L2_',num2str(nn(m)),'_',num2str(wl),'_',num2str(frequency),'.epsc'])

%     [~,idx] = min(d2(1:end,1))
%     figure
%     plot(1:size(d2,1),d2(:,1),'.-','MarkerSize',10,'LineWidth',3.5,'color',[color(p(m),1),color(p(m),2),color(p(m),3)])
%     % legend('0.3','0.1','0.06','0.03','0.01','0.006','0.003','0.001','Location','West','FontSize',20)
%     xlabel('lambda2','FontSize',15)
%     ylabel('distance','FontSize',15)
%     xticks([1 2 3 4 5 6 7 8])
%     xtickangle(45)
%     xticklabels({'0.3','0.1','0.06','0.03','0.01','0.006','0.003','0.001'});

end