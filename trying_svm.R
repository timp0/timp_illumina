
Error: could not find function "ksvm"
> library(kernlab)
> 
Error in data.frame(y, x) : 
  arguments imply differing number of rows: 64, 329
> nac
 [1] Adenoma   Carcinoma Normal    Normal    Carcinoma Normal    Normal   
 [8] Adenoma   Normal    Adenoma   Normal    Adenoma   Normal    Adenoma  
[15] Normal    Carcinoma Normal    Carcinoma Carcinoma Adenoma   Adenoma  
[22] Normal    Adenoma   Normal    Adenoma   Normal    Adenoma   Normal   
[29] Carcinoma Normal    Adenoma   Carcinoma Adenoma   Normal    Carcinoma
[36] Carcinoma Carcinoma Normal    Adenoma   Carcinoma Normal    Normal   
[43] Adenoma   Carcinoma Adenoma   Carcinoma Normal    Normal    Normal   
[50] Normal    Adenoma   Carcinoma Carcinoma Normal    Adenoma   Carcinoma
[57] Adenoma   Normal    Adenoma   Carcinoma Carcinoma Normal    Adenoma  
[64] Normal   
Levels: Normal Adenoma Carcinoma
> tdata=data$fqbeta[,th]
> dim(tdata)
[1] 384  64
> tclass=ksvm(x=t(tdata), y=nac, cross=length(nac))
Using automatic sigma estimation (sigest) for RBF or laplace kernel 
> tclass
Support Vector Machine object of class "ksvm" 

SV type: C-svc  (classification) 
 parameter : cost C = 1 

Gaussian Radial Basis kernel function. 
 Hyperparameter : sigma =  0.00149819026920112 

Number of Support Vectors : 58 

Objective Function Value : -13.9556 -9.2787 -20.9742 
Training error : 0 
Cross validation error : 0.25 
> 

tissue_density(data$fqbeta,data$fsamp, data$probes, interest=hector$X[1:20],namey="top20universe.pdf")

