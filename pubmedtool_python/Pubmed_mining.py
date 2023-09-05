import metapub
import pandas as pd
from urllib.request import urlretrieve
import textract
from metapub import PubMedFetcher
#from scihub import SciHub
import numpy as np 
import os
#from scihub import SciHub


def download_pubmed_paper(PMIDs, major_path, sub_path):
    count = 0
    log = ["log start"]
    sh = SciHub()
    fetch = PubMedFetcher()
    #"---------------------------------------------"
    exist_ = os.listdir(major_path+sub_path)
    nnn=[]
    for x in exist_: 
        nnn.append(x.replace('.pdf',""))
    #"---------------------------------------------"
    PMIDs_  = pd.DataFrame(set(PMIDs[0].values.tolist())-set(nnn))#扣除掉已經存在資料夾的
    exist_pdf = set(PMIDs[0].values.tolist())&set(nnn) #轉成set才能比對是否有聯集
    print("existed PMIDs:",str(len(exist_pdf)),"; unique PMIDs:",str(len(PMIDs_)))      
    print("start")
    #"---------------------------------------------"
    for i in range(len(PMIDs_)):
        pmid = int(PMIDs_.iloc[i].values) #為了因應有些轉不了str的情況，必須先value再轉int
        #雖然原本單純輸入PMIDs ok 但如果來源是多人種的表單去調出來的，資料型態問題，會一併把欄列名str化，無法單獨取出PMID
        #"---------------------------------------------"
        if str(pmid) in nnn:
            print(pmid,"exist")
            #"---------------------------------------------"
        else:
            try:
                url = metapub.FindIt(str(pmid)).url
                print(pmid)
                
                if not pd.isnull(url):
                    path_pdf  = major_path + sub_path+"/"+str(pmid)+".pdf"
                    #"---------------------------------------------"
                    try:
                        urlretrieve(url, path_pdf)
                        count += 1
                        log = log+[str(pmid)]
                        print(str(count)+" "+str(pmid)+" by Pubmed")
                        log_output = {'log': log}
                        pd.DataFrame(log_output ).to_csv(major_path + sub_path + "_log.csv", sep= ',', mode="a")
                    except:
                        pass
                else:
                    #假如url空的，就啟動sci-hub       
                    try:     
                        count += 1
                        print(str(count)+" "+str(pmid)+" by Sci-hub")
                        log = log+[str(pmid)+"Sci-hub"]
                        path_pdf  = major_path + sub_path+"/"+ str(pmid)+".pdf"
                        #如果拿到PMID，撈不到PDF，則用這個方法下載
                        result = sh.download("https://pubmed.ncbi.nlm.nih.gov/"+str(pmid) +"/", path= path_pdf)
                        log_output = {'log': log}
                        pd.DataFrame(log_output ).to_csv(major_path + sub_path + "_log.csv", sep= ',', mode="a")
                    except:
                        pass
            except:
                pass #先pass 以後再修正成轉去sci-hub 20230112 metapub: 'NoneType' object is not subscriptable
