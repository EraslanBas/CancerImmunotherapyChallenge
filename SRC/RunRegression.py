from libraries import *
from parameters import *

from warnings import filterwarnings
filterwarnings('ignore')
from sklearn import linear_model
import statsmodels.formula.api as smf
import multiprocessing
import time

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')



def RunRegression( setIndex_1, setIndex_2, 
                    expressionMatrix, guideMatrix, my_formula):
    
    allRes = pd.DataFrame()


    for i in range(setIndex_1,setIndex_2,1):
        #print(i)
        guideMatrix["y"] = np.array(expressionMatrix)[:,i]

        mod1 = smf.ols(formula=my_formula, data=guideMatrix).fit()

        k = mod1.summary()
        res = pd.DataFrame( k.tables[1])
        res["respGene"] = expressionMatrix.columns[i]

        allRes = pd.concat([allRes,res])
        
    allRes.to_csv("./TmpReg/GuideRegTest_"+str(setIndex_1)+"_"+str(setIndex_2)+".csv", index=False)
    return



if __name__ == "__main__":
    
    
    os.getcwd()
    os.chdir(projectDir)
    
    start_time = time.perf_counter()
    
    
    adata = sc.read_h5ad('./DATA/sc_training.h5ad')
    adata.obs
    
    
    
    guideMatrix = pd.get_dummies(data=adata.obs.condition, drop_first=False)
    guideMatrix = guideMatrix.drop(['Unperturbed'], axis=1)
    guideMatrix = guideMatrix.join(adata.obs[["n_genes", "mt_frac"]])
#     stateMat = pd.get_dummies(data=adata.obs.state, drop_first=False)
#     stateMat.columns=["cycling", "effector", "other", "progenitor", "exhausted"]
#     stateMat = stateMat.drop(['progenitor'], axis=1)
#     guideMatrix = guideMatrix.join(stateMat)
    
    
    expressionMatrix = pd.DataFrame(adata.X.A)
    expressionMatrix.columns = adata.var_names
    expressionMatrix.index = adata.obs.index
    
    
    myFormula = "+".join(guideMatrix.columns)
    my_formula = "y~" + myFormula
    my_formula
    
    
    par_test_guide_interval = 500
    processes = []

    for i in range(0,len(expressionMatrix.columns),par_test_guide_interval):
        setIndex_1 = i

        if setIndex_1 +  par_test_guide_interval < len(expressionMatrix.columns):
            setIndex_2 = setIndex_1 + par_test_guide_interval
        else:
            setIndex_2 = len(expressionMatrix.columns)


        p = multiprocessing.Process(target = RunRegression, args=(setIndex_1, 
                                                                  setIndex_2,
                                                                  expressionMatrix, guideMatrix,
                                                                  my_formula))
        p.start()
        processes.append(p)




    # Joins all the processes 
    for p in processes:
        p.join()
        
    finish_time = time.perf_counter()
    
    print(f"Program finished in {finish_time-start_time} seconds")

    ##Combine all results
    allRes = pd.DataFrame()
    for i in range(0,len(expressionMatrix.columns),par_test_guide_interval):
        setIndex_1 = i

        if setIndex_1 +  par_test_guide_interval < len(expressionMatrix.columns):
            setIndex_2 = setIndex_1 + par_test_guide_interval
        else:
            setIndex_2 = len(expressionMatrix.columns)

        res =pd.read_csv("./TmpReg/GuideRegTest_"+str(setIndex_1)+"_"+str(setIndex_2)+".csv")
        allRes = pd.concat([allRes,res])
        #os.remove("./TmpReg/GuideRegTest_"+str(setIndex_1)+"_"+str(setIndex_2)+".csv")
        
    allRes.columns = ['guides', 'coef', 'stderr', 'z', 'pval', '[0.025', '0.975', 'respGene'] 
    allRes = allRes[~allRes.guides.isin(['Intercept', 'n_genes', 'mt_frac', 'NaN']) & 
                    ~allRes.guides.isnull() &  
                    ~allRes.guides.str.contains("leiden", case=False, na=False) ]
 
    
    allRes.to_csv("AllResults_noCellStates.csv", index=False)
 