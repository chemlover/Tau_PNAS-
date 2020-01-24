%pdb_array= {'Abeta42_2' 'Abeta420' 'Abeta40_2'
pdb_array= {'tau'};
%sim_labels = [4 4 4];
sim_labels = [12];

% what data to loadout_file.txt
oligomer ='heptamer'
q_name={['p_total-',oligomer],['e_total-',oligomer]};
T=350
%xsym=' energy '
xsym = 'Q fibril'
ysym='inter para hydrogen bond'
ysym = 'energy'
%ysym = 'ecoul'
%ysym = 'contact number'
%xy_limit=1; % whose range to use for plotting, 1 or 2?

n_contour=30; % # of contour lines
fsize=20; tsize=20; mr=1; mc=length(sim_labels);
cutoff=80.1;scrnsize = get(0,'ScreenSize'); 
figure('position', [1 scrnsize(4) 0.25*mc*scrnsize(3) 0.35*scrnsize(4)]);

for i_label=1:length(sim_labels)
    protein_i = i_label; pdbID_upper = pdb_array{protein_i};
    path = sprintf('/home/xc25/Data_dealing/%s',pdbID_upper);
    %path = sprintf('/home/xun/Downloads/download/fis/specifity-2/small')
    sim_label = sim_labels(i_label);
    qa_name=q_name{1}; qb_name=q_name{2};
    
    filename = sprintf('%s/%s',path, qa_name); qa = load(filename);
    filename = sprintf('%s/%s',path, qb_name); qb = load(filename);
    if strcmp(q_name{1},'dih')
        qa=(mean(qa'))';
    end
    if strcmp(q_name{2},'dih')
        qb=(mean(qb'))';
    end
    if i_label==1
        qa_min=min(qa); qa_max=max(qa);
        qb_min=min(qb); qb_max=max(qb);
    end
    Nsample=size(qa,1);
    assert(Nsample==size(qb,1));
    %Nsample =length(qa)
    %load pmf file and calculate pi_sample
    T_i=T; T=T_i;    
    filename = sprintf('%s/p_total-%s',path,oligomer); q = load(filename);
    filename=sprintf('%s/%s_%d_pmf.dat',path,pdbID_upper, T_i);
    %filename=sprintf('%s/small_300_pmf.dat',path);
    FF=load(filename); qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
    dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
    Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
    pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
    %calculate pi_sample
    for i_bin= 1:nbin
        qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
        ids = find( q >= qi_min & q < qi_max ) ;    
        ni_sample(i_bin) = length(ids);        
        if ni_sample(i_bin) > 0
            pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
        end
    end    
    fprintf('probability = %.3f\n', sum(pi_sample));
        
    binN=20;ObinN=8;
    qa_lin=linspace(min(qa), max(qa),ObinN); qb_lin=linspace(min(qb), max(qb),binN); H=zeros(ObinN,binN);
    [~, bin_index_x] = histc(qa, qa_lin); [~, bin_index_y] = histc(qb, qb_lin);
  % change for -nktln(c1/c0)
    %for i_sample = 1:Nsample
     %   x=bin_index_x(i_sample); y=bin_index_y(i_sample);
      %  if qa_lin(x) ~=  5
       %   qb(i_sample)=qb(i_sample)-qa_lin(x)*0.001987*T*log(0.0002705*(5-qa_lin(x)));
        %  pi_sample(i_sample)=pi_sample(i_sample)*0.0002705*(5-qa_lin(x))*exp(-qa_lin(x)*0.001987*T);
        %end
    %end
    %qb_lin=linspace(min(qb), max(qb),binN);
    %[~, bin_index_x] = histc(qa, qa_lin); [~, bin_index_y] = histc(qb, qb_lin);
    for i_sample = 1:Nsample
        x=bin_index_x(i_sample); y=bin_index_y(i_sample);
        H(x,y) = H(x,y) + pi_sample(i_sample);
      end
    H=H'; fprintf('sum(sum(H))=%.3f\n', sum(sum(H)));
    
    F=-0.001987*T*log(H); 
    ids = (F>= cutoff); F(ids) = -inf; 
    subplot(mr,mc,i_label)
    [~,h] = contourf(qa_lin, qb_lin,F,n_contour,'edgecolor','none'); shading flat,
    colormap(jet), col=colorbar, %set(col,'ylim',[0 10])
    %title([num2str(pdbID_upper), ' T= ', num2str(T)],'fontsize', fsize);
    title([ 'T= ', num2str(T),'-',oligomer],'fontsize', fsize);
    %%%fill the top area white
    ccc = get(h,'children'); max_cdata = -inf; cdata_list=zeros(size(ccc,1), 1);
    
    for k=1:size(ccc,1)
        cd1 = get(ccc(k), 'cdata');
        if cd1 > max_cdata
            max_cdata = cd1 ;
        end
        cdata_list(k) = get(ccc(k),'cdata');
    end
    id = find(cdata_list == max_cdata);
    disp(ccc(id));
    for k=1:size(id,1)
        set(ccc(id(k)), 'facecolor', 'white');
    end
    
    %xlabel(q_name{1}, 'fontsize', fsize), 
    %ylabe'
  
    %ysym ='energy'
    xlabel(xsym, 'fontsize', fsize)
    ylabel(ysym, 'fontsize', fsize)
    %set(gca, 'FontSize', fsize);
    xlim([qa_min, qa_max])
    ylim([qb_min, qb_max])
    %caxis([0,60])
    %xlim([0,180])
    saveas(gcf,['/home/xc25/Data_picture/tau_aggregation_1/',xsym,'-',ysym,'-',oligomer,'.png'])
end