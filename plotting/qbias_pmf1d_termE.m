%pdb_array= { 'Abeta42_1' 'Abeta40_1'};
function qbias_pmf1d_termE
figure(1);


energy_array=load ('/home/xc25/Data_dealing/ER_project/ER_3_old/bias/coa/e_test');

%array = [2, 4, 5, 7, 8, 9,10, 17,19,12,18]

%Name = {'Chain','Chi','Rama','DSSP','PAP','Water', 'burial', 'Elec','V_{total}', 'amh-go', 'Q-S'}

for i=1
    %subplot(3,3,i);
    qbias_pmf1d(energy_array,'e_total');
end

function qbias_pmf1d(energy_array,name)


% Burial Chain Chi DSSP Frag_Mem Helix P_AP Rama Water
T=300;
fsize = 15;
for T = 300
    color={'r','m','g','b'};


    for i_label=1
        path=sprintf('.');
        %title([pdbID_upper, ' ', num2str(T), 'K'], 'fontsize', fsize);
        
        %filename = sprintf('%s/%s',path,ea_name); 
        ea = energy_array;        
        path = '/home/xc25/Data_dealing/ER_project/ER_3_old/bias/coa'
        %load q
        filename = sprintf('%s/p_total',path); q = load(filename);
        %if qo_flag == 1
        %    filename = sprintf('%s/qo_%d',path, sim_label); qo = load(filename); 
        %end
        Nsample = length(q);    

        %load pmf file and calculate pi_sample
        filename=sprintf('%s/coa_%d_pmf.dat',path,T);

        FF=load(filename); qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
        dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
        Py=exp(-Fy/(0.001987*T)); P_norm = sum(Py); Py=Py/P_norm;
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

        %calculate the average energy term value in each window.
        h=zeros(nbin,1);
        for i_bin= 1:nbin
            term=0.0;
            prob=0.0;
            qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
            ids = find( q >= qi_min & q < qi_max ) ;   
            
            ni_sample(i_bin) = length(ids);  
            for j=1:ni_sample(i_bin)
                term=term + pi_sample(ids(j))*ea(ids(j));
                prob=prob+pi_sample(ids(j));
            end
            h(i_bin)=term/prob;
        end   
        plot(qx, h+Fy, ':ko', 'linewidth', 2.0); 
        
        set(gca,'FontSize',fsize) ;
        %xlim([-90, 90]);
        ylim([min(h),max(h)])
        xlabel('Twisting Angle', 'fontsize', fsize); ylabel('kcal/mol', 'fontsize', fsize);title(name);
        hold on; %plot(-180*ones(numel(qx),1),h,'b','linewidth',3);

    end
end
