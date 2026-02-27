function demBinaryClassification_Marg_fixedhypers_pMALA(rep)


addpath aGrad/code/toolbox/;
addpath aGrad/data/;

%dataName = 'Heart';
for dataName = {'Australian' 'German' 'Heart' 'Pima' 'Ripley'}% 'Caravan'}
     
   disp(['Running with ' dataName ' dataset']);
   %switch(dataName)
    switch(dataName{1})
   case 'Caravan'
    
        %Load and prepare train & test data
        load('caravan.mat');
        
        % End column contains binary response data
        t=X(:,end);
        X(:,end) = [];

   case 'Australian'
    
        %Load and prepare train & test data
        load('australian.mat');
        t=X(:,end);
        X(:,end)=[];  
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernAustralian;
      
   case 'German'

        %Load and prepare train & test data
        load('german.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;
          
        load('initkern_binaryclass.mat'); 
        initkern = initkernGerman; 
              
   case 'Heart'
        
        %Load and prepare train & test data
        load('heart.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;  
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernHeart;
       
   case 'Pima'
       
        %Load and prepare train & test data
        load('pima.mat');
        t=X(:,end);
        X(:,end)=[];  
          
        load('initkern_binaryclass.mat'); 
        initkern = initkernPima;
        
   case 'Ripley'

        %Load and prepare train & test data
        load('ripley.mat');
        t=X(:,end);
        X(:,end)=[];
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernRipley;
   end   
   


   Y = 2*t - 1;
   
   [n D] = size(X);
   


   % normalize the data so that the inputs have unit variance and zero mean  
   meanX = mean(X);
   sqrtvarX = sqrt(var(X)); 
   X = X - repmat(meanX, n, 1);
   X = X./repmat(sqrtvarX, n, 1);
   
   % model options
   options = gpsampOptions('classification'); 
   options.kern = 'se';
   
   model = gpsampCreate(Y, X, options); 
   model.GP.jitter = 1e-8;
   model.GP.logtheta = initkern; 
   
   model.XX = -2*(X*X') + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   
   model.constraints.kernHyper = 'free';
   model.constraints.likHyper = 'free';
   model.delta =2/sqrt(size(X,1));
   mcmcoptions.T = 5000;
   mcmcoptions.Burnin = 5000;
   mcmcoptions.StoreEvery = 1;
   mcmcoptions.Langevin = 1;

   mcmcoptions.opt =0.574;
  
   % precompute
   model.K = kernCompute(model.GP, model.X);
   [model.U, model.Lambda, tmp] = svd(model.K);
   model.Lambda = diag(model.Lambda);



   tic;
   [model samples accRates] = gpsampAuxMarg_fixedhypers_pMALA(model, mcmcoptions);
   elapsedTime = toc;



 % compute statistics 
   summaryMarg = summaryStatistics(samples);
   ess = summaryMarg.eff_F;
   LogL= summaryMarg.LogL;
   summaryMarg.elapsed = elapsedTime;
   summaryMarg.accRates = accRates;
   summaryMarg.delta = model.delta; 
   eff_LogL = mcmc_ess(samples.LogL);
   save(['results/LogisticRegression_GP/' dataName{1} '_repeat' num2str(rep) '_Marg_fixedhypers_pMALA.mat'], ...
       'elapsedTime','accRates','eff_LogL','ess','LogL');
end   

end




