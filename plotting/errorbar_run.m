filename = '/home/xc25/Data_dealing/CASP13/backup_tr884'
FF = load(filename)
means = FF(:,1)
variances = FF(:,2)
x= linspace(1, length(means)+1,length(means))
%e=errorbar(x,means,variances,'-s','MarkerSize',3,'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',0.5);hold on;
figure;
plot(x,variances,'r','LineWidth',2)

%e.LineWidth = 0.5
xname = 'frame'
yname = 'variance'
fsize = 20
    xlabel(xname, 'fontsize', fsize); 
    ylabel(yname, 'fontsize', fsize); 
    set(gca, 'FontSize', fsize);
    xlim([min(x),max(x)]);
    %XTick = [1:5];
    %set(gca,'xtick',XTick)
    %set(gca,'Color','w')
    %caxis([0,6])
    set(gca,'fontsize',fsize)
    set(gca, 'FontName', 'Helvetica')