function demBinaryClassification_Marg_fixedhypers_gimala_VR(rep)



addpath aGrad/code/toolbox/; 
addpath aGrad/data/;


%dataName = 'German';
for dataName = {'Australian' 'German' 'Heart' 'Pima' 'Ripley'}% 'Caravan'}
%     
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
   
b2 = zeros(1,n);
means = zeros(1, n);
means_cv = zeros(1, n);
%var_reduction = zeros(1,n);

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
   mcmcoptions.T = 1000;
   mcmcoptions.Burnin = 5000;
   mcmcoptions.StoreEvery = 1;
   mcmcoptions.Langevin = 1;

   mcmcoptions.opt =0.75;
  
   % precompute
   model.K = kernCompute(model.GP, model.X);
   [model.U, model.Lambda, tmp] = svd(model.K);
   model.Lambda = diag(model.Lambda);



   tic;
   [model, samples, accRates] = gpsampAuxMarg_fixedhypers_gimala(model, mcmcoptions);
   elapsedTime = toc;


 myd = size(samples.F,2);
 for jj=1:myd
  
  thetaCV = NaN(2, mcmcoptions.T);


       g1 = samples.accprob(:,1).*( samples.Fy(:,jj)-samples.F(:,jj) )/model.delta;
       g2 = samples.Fy(:,jj)/model.delta - (1.0-model.delta)*samples.F(:,jj)/model.delta -samples.PGhat(:,jj);


       K11 = ergMean(g1.*g1);
       K12 = ergMean(g1.*g2);
       K21 = ergMean(g1.*g2);
       K22 = ergMean(g2.*g2);
       
       th1 = ergMean(samples.F(:,jj).*g1)-ergMean(samples.F(:,jj) ).*ergMean(g1);
       th2 = ergMean(samples.F(:,jj).*g2)-ergMean(samples.F(:,jj) ).*ergMean(g2);
       
  lambda = 1e-8;
tol_rcond = 1e-12;

for i = 100:mcmcoptions.T
    Kmat = [K11(i), K12(i); K21(i), K22(i)];
    Kmat = (Kmat + Kmat')/2;

    if rcond(Kmat) < tol_rcond
        thetaCV(:, i) = thetaCV(:, i-1);
    else
        thetaCV(:, i) = (Kmat + lambda*eye(2)) \ [th1(i); th2(i)];
    end
end
 means(1,jj) = mean(samples.F(:,jj));
 means_cv(1,jj)= mean(samples.F(:,jj)) - thetaCV(1,mcmcoptions.T)*mean(g1)- thetaCV(2,mcmcoptions.T)*mean(g2);

end

   aprob = accRates;
   save(['results/LogisticRegression_GP/' dataName '_repeat' num2str(rep) '_Marg_fixedhypers_VR_' num2str(mcmcoptions.T) '.mat'], ...
       'aprob','means','means_cv');
end









