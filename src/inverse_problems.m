clear all;
close all;
clc;

load ../data/data;
A=G;

%% Parameters

corelate_noise = true;     %Correlate Noise
algo = "MNE";   %Sissy, MNE, Gibbs



%lambda = [0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500,1000];
%lambda =[ 0.01, 0.2, 0.3, 0.4,  0.5,0.6, 0.7, 0.8, 0.9, 1:1:10, 10:10:100, 100:100:1000];
%lambda =[ 0.01, 0.05, 0.1,  0.5,1:1:30, 50, 100:100:1000];
%lambda =[  0.05, 0.1,  0.5,1:1:100, 100:100:1000];
%lambda = 1:0.1:30;
%lambda=[10^-10, 10^-9, 10^-8,  10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,0.1:0.1:1, 1:1:100,100:10:500, 1000, 10^4, 10^5, 10^6];  %Complet
%lambda=[  10^-7,10^-6,10^-5,10^-4,10^-3,10^-2, 0.1, 1, 10, 100,1000,10^5,10^6,10^7,10^8,10^9];
%lambda=[  10^-3,10^-2, 0.1:0.01:1,1:0.1:10, 10:1:100,100,200:10:400,400:100:1000];
lambda=10
alpha = 0.1;
%alpha = 0.1:0.1:1;


SNR=0.1;    %Signal to noise ratio


% Heuristics
f=zeros(1,length(lambda));
norm_err = zeros(1,length(lambda));
norm_s = zeros(1,length(lambda));
discrepancy = true;




%% Noise generation

%Generate linear mixture of source signals
Xs=A*S;

%Determine maximum of the signal of interest (here an epileptic spike) to
%apply source localization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
%plot_brain(S, path, mesh);


%Generate Gaussian random noise or correlate one
if corelate_noise==true
    fprintf("Bruit Corrélé\n");
    Sint = (S+1);
    Sint(Sint~=1) =0; 
    Snoize = (S*.0+randn(size(S))).*Sint;
    Noise = A*Snoize;
    bruit = "_corr_";

else 
    fprintf("Bruit non Corrélé\n");
    Noise=randn(size(Xs));
    bruit = "";
end
    
%Normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');


%% Execution Algo

if algo == "Sissy"
    % Construction de T
    T = variation_operator(mesh,'face');
end

if algo=="Gibbs"
    lambda = length(find(S(:,id)~=0))/length(S(:,id));
    sigma_s2 = var(S(:,id));
    sigma_n2 =var(Noise(:,id))/norm(Xs,'fro');  
end

i=1;
  
 for s=SNR
     %Generate noisy data according to given SNR
            X=Xs+1/sqrt(s)*Noise;
    for l=lambda
         for a=alpha
        
                if algo=="Sissy"
                  S= sissy(X(:,id), A, T, l, a);
                  t = ["Sissy | SNR  =  "+num2str(SNR)+" | \lambda = "+l+ " | \alpha = "+ a]; 
                end
                if algo =="Gibbs"
                    S = Gibbs_sampler(X(:,id),A, sigma_s2, sigma_n2, lambda);
                    t = ["Gibbs | SNR  =  "+num2str(SNR)]; 
                end
                if algo=="MNE"
                    S=MNE(X(:,id),A, lambda(i));
                    t = ["MNE | SNR  =  "+num2str(SNR)+" | \lambda = "+l]; 
                end
                
                close all;

                
               path = [algo+bruit+"SNR_"+num2str(SNR)+"lamb_"+l+ "alp_"+ a+'.png'] ;
              plot_brain(S, path, mesh, t);
                close all;
        
               %% heuristic
                        
                norm_s(i)=sqrt(sum(S.^2));
                norm_err(i)=sqrt(sum(((X(:,id)-A*S).^2)));
              
                   if algo=="Sissy"
                        s = S;
                        s(s<0.01*max(s))=0;
                        n2=nnz(s);

                        ts =  T*s;
                        n1 = nnz(ts);

                        f(i) = n1+a*n2;  
                   end
                   
                    i=i+1;
        
    end
end
 end

 

 %% L0
 if algo=="Sissy"
    figure();
    semilogx(lambda,f);
    xlabel("\lambda",'Fontsize',18)
    ylabel("f", 'Fontsize',18);

    %Lambda recovery
    [val, inx]=min(f);
    l_f = lambda(inx);

    line('XData', [l_f l_f], 'YData', [4000 min(f) ], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
     str2 = {"\lambda = "+num2str(l_f)};
     text(5,4100,str2,'Color','red','Fontsize',16);
 end

%% Discrepancy priciple

if discrepancy

    %Lambda recovery
    n_noise = norm(Noise(:,id), 'fro')^2;
    nr =norm_err.^2;
    [ closestValue ind ] = min(abs(nr-n_noise));
    l_disc = lambda(ind);

    figure();
   plot(nr, lambda);
    line('XData', [n_noise n_noise ], 'YData', [0 30 ], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
     line('XData', [min(norm_err.^2) max(norm_err.^2) ], 'YData', [l_disc l_disc ], 'LineStyle', '--', ...
         'LineWidth', 1,'Color','red');
    str2 = {"\lambda = "+num2str(l_disc)};
   text(3*10^4,200,str2,'Color','red','Fontsize',16);


    xlabel("||x-As||_2^2", 'FontSize', 18);
    ylabel("\lambda", 'FontSize', 18);
    title("Discrepancy principle","FontSize",20);
end

%% Lcurve criterion

if algo=="MNE"
     figure();
     loglog(norm_err,norm_s);
     xlabel("||x-As||_2", 'FontSize', 18);
     ylabel("||s||_2", 'FontSize', 18);
     title("L-curve criterion", 'FontSize', 20);
     
     nl =9;
    line('YData', [norm_s(nl) norm_s(nl)], 'XData', [10^-2 10^3 ], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
    
    line('XData', [norm_err(nl) norm_err(nl)], 'YData', [10^-1 10^5], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
    str2 = {"\lambda = "+num2str(lambda(nl))};
    text(405,100,str2,'Color','red','Fontsize',16);
     
end

%% Generalized cross-validation

if algo=="Sissy"
    gcv = zeros (1, length(lambda));
    j=1;
    for l=lambda
        gcv(j)=norm_err(j)^2/trace(eye(91)-A*A'*inv(A*A'+l*eye(91)))^2;
        j=j+1;
    end

    %Lambda recovery
    [val, inx]=min(gcv);
    l_gcv = lambda(inx);

    plot(lambda,gcv)
     xlabel("\lambda", 'FontSize', 18);
     ylabel("GCV(\lambda)", 'FontSize', 18);

     line('YData', [min(gcv) min(gcv) ], 'XData', [0 l_gcv ], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
    line('YData', [0   min(gcv)], 'XData', [l_gcv l_gcv ], 'LineStyle', '--', ...
        'LineWidth', 1,'Color','red');
    str2 = {"\lambda = "+num2str(l_gcv)};
       text(40, 20, str2,'Color','red','Fontsize',16);
end



