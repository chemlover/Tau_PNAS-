%pdb_array= {'Abeta42_2' 'Abeta420' 'Abeta40_2'
pdb_array= {'phos'};

sim_labels = [12];
% what data to load
%q_name={  'largest_oligomer_total' ,'tc_total'};xname = 'largest oligomer size',yname='contact number'
%q_name={  'largest_oligomer_total' ,'e_total'};xname = 'largest oligomer size',yname='energy'
%q_name={  'largest_oligomer_total' ,'p_total'};xname = 'largest oligomer size',yname='Qw'
q_name={  'largest_oligomer_total' ,'inter_antipara_total'};xname = 'oligomer size',yname='#intrachain antiparalell hydrogen bond'
%q_name={  'largest_oligomer_total' ,'ecoul_total'};xname = 'oligomer size',yname='electrostatic'
q_name={  'largest_oligomer_total' ,'inter_para_registered_total'};xname = 'oligomer size',yname='#interchain parallel registered hydrogen bond'
q_name={  'inter_para_registered_total' ,'inter_antipara_total'};xname = '#interchain parallel registered hydrogen bond',yname='#intrachain antiparallel hydrogen bond'
%q_name = {'charge_ll','dist_protein_dna_bindsite_total'};xname = 'angle';yname = 'dist between helix and dna'
%q_name={  'largest_oligomer_total' ,'inter_para_total'};xname = 'largest oligomer size',yname='interpara HB'
%q_name={  'largest_oligomer_total' ,'intra_para_total'};xname = 'largest oligomer size',yname='intrapara HB'
%q_name={  'largest_oligomer_total' ,'inter_total'};xname = 'largest oligomer size',yname='inter HB'
%q_name={  'largest_oligomer_total' ,'intra_total'};xname = 'largest oligomer size',yname='intra HB'
%q_name={  'inter_para_total-pentamer' ,'inter_antipara_total-pentamer'};xname = '#interchain parallel hydorgen bond',yname='#interchain antiparallel hydorgen bond'
q_name={  'first_oligomer_total' ,'inter_para_total'};xname = 'oligomer size',yname='#interchain parallel hydrogen bond'
%q_name={  'inter_para_registered_total' ,'inter_antipara_total'};xname = '#interchain parallel registered hydrogen bond',yname='#intrachain antiparallel hydrogen bond'
q_name={'48_angle.dat','48_ecoul.dat'};xname = 'DNA bending angle at the binding site',yname='E_{elec} at DNA binding site'

T=300;
%xy_limit=1; % whose range to use for plotting, 1 or 2?
newpoints =15;Onewpoints=12
n_contour=30; % # of contour lines
fsize=20; tsize=20; mr=1; mc=length(sim_labels);
cutoff=4.5;scrnsize = get(0,'ScreenSize'); 
figure('position', [1 scrnsize(4) 0.25*mc*scrnsize(3) 0.35*scrnsize(4)]);



for i_label=1:length(sim_labels)
    protein_i = i_label; pdbID_upper = pdb_array{protein_i};
    path = sprintf('.');
    sim_label = sim_labels(i_label);
    qa_name=q_name{1}; qb_name=q_name{2};contourf
   % path = sprintf('/home/xc25/Data_dealing/tau_R134/tau/phosnew2/phos5/')
    path = sprintf('/home/xc25/Data_dealing/tau_agggregagtion/tau/phosnew2/phos4/')
    path = sprintf('/home/xc25/Data_dealing/PU_1/update5/design/bias/setup2/test5')
    filename = sprintf('%s/%s',path, qa_name); 
    qa =180 - load(filename); 
    filename = sprintf('%s/%s',path, qb_name); 
    %rg = load('rg_total');
    qb = load(filename);
    %qb = qb./rg.^3
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
    %assert(Nsample==size(qb,1));
    %load pmf file and calculate pi_sample
    T_i=T; T=T_i;    
    %path = sprintf('/home/xc25/Data_dealing/tau_agggregagtion/tau_9th_3rd')
    filename = sprintf('%s/p_total',path); q = load(filename);
    filename=sprintf('%s/setup2_%d_pmf.dat',path, T_i);
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
    binN=30;ObinN=20;
    qa_lin=linspace(min(qa), max(qa),ObinN); qb_lin=linspace(min(qb), max(qb),binN); H=zeros(ObinN,binN);
    [~, bin_index_x] = histc(qa, qa_lin); [~, bin_index_y] = histc(qb, qb_lin);
    for i_sample = 1:Nsample
        x=bin_index_x(i_sample); y=bin_index_y(i_sample);
        H(x,y) = H(x,y) + pi_sample(i_sample);
    end
    H=H'; fprintf('sum(sum(H))=%.3f\n', sum(sum(H)));    
    F=-0.001987*T*log(H);%ids = (F>= cutoff); F(ids) = cutoff;
   %[xnew,ynew]=meshgrid(linspace(min(qa),max(qa),newpoints),linspace(min(qb),max(qb),Onewpoints))
   ids = (F>= cutoff); F(ids) = cutoff+1;
   % znew=interp2(qa_lin,qb_lin,F,xnew,ynew,'spline')
   % ids = (znew>= cutoff); znew(ids) = cutoff+1;
    %[~,h] = contourf(xnew, ynew,znew,n_contour,'edgecolor','none');
    [~,h] = contourf(qa_lin, qb_lin,F,n_contour);% shading flat,
    colormap('jet'), colorbar, 
    
    colormap(jet), col=colorbar, %set(col,'ylim',[0 10])name
    cm=colormap;
    cm(256, :) = [1 1 1]
    %cm(63, :) = [1 1 1];
    colormap(cm);
    %title([num2str(pdbID_upper), ' T= ', num2str(T)],'fontsize', fsize);
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
    fsize = 50
    xlabel(xname, 'fontsize', fsize); 
    ylabel(yname, 'fontsize', fsize); 
    set(gca, 'FontSize', fsize);
    xlim([min(qa_lin),max(qa_lin)]);
    %XTick = [1:5];
    %set(gca,'xtick',XTick)
    %set(gca,'Color','w')
    %caxis([0,6])
    set(gca,'fontsize',fsize/2)
    set(gca, 'FontName', 'Helvetica')
    %xlim([1,5])
    %ylim([0,80])
    %title('nophosphoylation','fontsize', fsize)
    %saveas(gcf,[path,'/',xname,'-',yname,'.eps'])
    %untitled.png
end
