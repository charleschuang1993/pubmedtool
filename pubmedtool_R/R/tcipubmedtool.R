#' Pubmed文獻搜尋爬蟲: 兩組關鍵子彼此排列組合進行搜尋
#'
#' 可設定不同關鍵字組合，搜尋Pubmed，並把搜尋結果彙整成表格
#' SNPs為主要關鍵字, factors為次要關鍵字。
#' 範例說明: 例如不同SNPs v.s. 不同人種，會有多個排列組合 [SNP] AND [Country]
#' ouput file: Yes (.csv)
#' 
#' @param path 想要輸出到哪個路徑資料夾儲存
#' @param SNPs 主要關鍵字/SNP清單
#' @param factors 次要關鍵字
#' @param filename 依照不同搜尋結果存放在不同資料夾
#' @return 包含PMID的CSV file
#' @export

SNPs_pubmed_for_multiple_factors <- function(path, SNPs, factors, filename){ #先不代入疾病/或種族，先廣再收斂
    count <- 0
    factors_df <- array()
    for(f in factors){
       SNP_df <- data.frame(matrix(ncol=2,nrow=0))
       print(f)
       for(SNP in SNPs){
        count <- count+1 #確認即便有error，迴圈是否繼續進行
        query_ <- paste0(f, " AND ", SNP )
        snp_on_pubmed <- get_pubmed_ids(query_)
        if(snp_on_pubmed$Count >= 1 ){
            print(SNP)
            snp_abstracts_xml <- fetch_pubmed_data(snp_on_pubmed ) #如果沒結果，會卡很久，需要直接跳過不能再繼續執行
            snp_abstracts_list <- articles_to_list(snp_abstracts_xml )
            PMIDs <- paste(fetch_all_pubmed_ids( snp_on_pubmed ), collapse = "; ")
            SNP_df <- rbind(SNP_df, c(SNP,PMIDs))
        }else{
            print("No")
            SNP_df <- rbind(SNP_df,c(SNP,"No"))
        }
        colnames( SNP_df )<- c("rsID",paste0("PMID_",f))     
     }

    factors_df <- cbind(factors_df, SNP_df)
    print(head(factors_df))
   }
    write.csv(factors_df, paste0(path, filename,".csv") )
}

#' Pubmed文獻搜尋爬蟲: 已經預設配對好的關鍵字組合
#'
#' 針對G2 384突變位點所設計，可設定不同關鍵字組合，搜尋Pubmed，並把搜尋結果彙整成表格
#' SNPs為主要關鍵字, factors為次要關鍵字。
#' 範例說明: 例如不同SNPs v.s. 不同人種，會有多個排列組合 [SNP] AND [Country]
#' ouput file: Yes (.csv)
#' 
#' @param path 想要輸出到哪個路徑資料夾儲存
#' @param disease_col 欄位名稱:第一組預設好的關鍵字組
#' @param snp_col 欄位名稱:第二組預設好的關鍵字組
#' @param data 包含兩組預設關鍵字組的完整檔案清單
#' @param filename 輸出的檔案名稱
#' @return CSV file
#' @export


SNPs_pubmed_for_paired_disease_and_SNP <- function(path, disease_col, snp_col, data, filename){ #先不代入疾病/或種族，先廣再收斂
    count <- 0
    result_df <- array()
    data_ <- data[,c(disease_col, snp_col )]

    for(line_ in seq_len(nrow(data_))){
        print(line_ )
        Disease <- data_[line_, disease_col] 
        SNP <- data_[line_, snp_col]        
        print(paste0(c(Disease, SNP), collapse="_"))
        count <- count+1 #確認即便有error，迴圈是否繼續進行
        query_ <- paste0(Disease, " AND ", SNP )
        snp_on_pubmed <- get_pubmed_ids(query_)
        if(snp_on_pubmed$Count >= 1 ){
            print(c(SNP,snp_on_pubmed$Count))
            snp_abstracts_xml <- fetch_pubmed_data(snp_on_pubmed ) #如果沒結果，會卡很久，需要直接跳過不能再繼續執行
            snp_abstracts_list <- articles_to_list(snp_abstracts_xml )
            PMIDs <- paste(fetch_all_pubmed_ids( snp_on_pubmed ), collapse = "; ")
            number_ <- length(unlist(strsplit(split="; ",PMIDs)))
            result_df <- rbind(result_df , c(SNP, Disease, number_,  PMIDs))
            
        }else{
            print("No")
            result_df <- rbind(result_df , c(SNP, Disease, 0 ,"No result"))
        }
        colnames( result_df)<- c("rsID","Disease", "Number", "PMIDs")}
        result_df<-result_df[-1,]
    #print(head(result_df))
    print(c("finished, preparing output results","Total",nrow(result_df)))   
    write.csv(result_df , paste0(path, filename,".csv") )
 
}


#' 輸入已確認過的PMIDs檢查是否自動搜尋和手動搜尋結果一致
#' 
#' 這個function會將ID另外output成一個檔，並且給出驗證訊息
#' @param path 想要輸出到哪個路徑資料夾儲存
#' @param SNP 只輸入單一個SNP rsID/ 第一關鍵字
#' @param race 只輸入單一個人種/第二關鍵字
#' @param disease_type 只輸入一個疾病/第三關鍵字
#' @param PMIDs_manual 將先前手動搜尋得到的PMIDs拿來驗證
#' @return 輸出檔案，其中包含PMID和驗證訊息
#' @export


Return_pubmedid_and_validation <- function(path, disease_type, race, SNP, PMIDs_manual ){
  validation_ <- unique(PMIDs_manual)  
  query_ <- paste0(race, " AND ", SNP," AND (", disease_type,")" )
  snp_on_pubmed <- get_pubmed_ids(query_)
  snp_abstracts_xml <- fetch_pubmed_data(snp_on_pubmed )
  snp_abstracts_list <- articles_to_list(snp_abstracts_xml )
  PMIDs <- data.frame(fetch_all_pubmed_ids( snp_on_pubmed )) #回傳相關的文獻PMIDs
  colnames(PMIDs) <- "PMID"
  variable_ <- gsub(" ", "_", disease_type) #作為輸出檔名，加上底線
  write.csv(PMIDs,paste0(path, variable_, "_PMIDs.txt"), row.names=FALSE)
  "----------------------------------------------------------------"
  match_result <- PMIDs[[1]][match(validation_, PMIDs[[1]])]
  assign(variable_, list("Validation" = validation_ , "PMIDs" = PMIDs[[1]], "Result" = match_result, "NA" = sum(is.na(match_result)) ), envir=.GlobalEnv)
  status_ <- if(sum(is.element(validation_, PMIDs[[1]]))==length(validation_ )){print("OK")}else{validation_[!is.element(validation_ , PMIDs[[1]])]}
  result_ <- data.frame(cbind("auto_PMIDs" = length(PMIDs[["PMID"]]),"manual_PMIDs"=length(validation_ ), "NA" = sum(is.na(match_result)) ))
  result_$"Status" = paste0(status_ , collapse =", ")
  rownames(result_) = variable_ 
  return(result_)
 }

#' 多重搜尋驗證
#' 
#' @param path 想要輸出到哪個路徑資料夾儲存
#' @param disease 第一關鍵字
#' @param SNP 第二組關鍵字
#' @param race 第三組關鍵字
#' @param PMIDs_manual 將先前手動搜尋得到的PMIDs拿來驗證
#' @return txt file
#' @export

Multiple_search_validation_message <- function(path, diseases, SNP, race, PMIDs_manual){
  result_ <- NULL
  for(x in seq_len(length(diseases))){
  result <- Return_pubmedid_and_validation(disease_type = diseases[x], SNP=SNP, race= race, PMIDs_manual[[x]])
  result_ <- rbind(result_,result )}
  write.csv(result_ , paste0(path,"Check_table","_PMIDs.txt"))
  return(result_)}

#' 從句子中萃取rsID
#' 
#' @param only_rsID 輸入句子
#' @return 回傳句子中以rs作為前輟的單字
#' @export

rsID_extraction <- function(only_rsID){
    words <- tokenize_words(only_rsID)
    tab <- table(words[[1]])
    tab_ <- data.frame(word = names(tab), count = as.numeric(tab))
    extracted_rsID <- names(tab)[grepl("^rs",names(tab))] #^代表鎖定起始
    return(extracted_rsID )
  }

#' 輸入句子和關鍵字，回傳確實存在句子中的關鍵字
#' 預設忽略英文字母大小寫
#' 
#' @param sentences 輸入句子
#' @param keywords 輸入關鍵字
#' @return 回傳確實存在於句子中的關鍵字
#' @export

return_value_for_each_sentance = function(sentences, keywords){
  return_=c()
  for(s in sentences){
    result <- paste0(keywords[sapply(tolower(keywords), function(k) grepl(k, s, fixed = FALSE,ignore.case=TRUE),simplify = TRUE,USE.NAMES = FALSE)],collapse=", ")
    return_ <- c(return_, result)
  }
  return(return_)
}

#' 輸入句子和關鍵字，回傳確實存在句子中的關鍵字
#' 預設忽略英文字母大小寫
#' 針對單個pmid，而一組rsID+Disease有多個pmids; 一個pmid有多個句子
#' 輸入單一個SNP
#'  
#' @param path 檔案存放路徑
#' @param pmid_list 用來搜尋資料夾中檔名為PMID的PDF檔
#' @param rsID 第一關鍵字
#' @param SNP 第二關鍵字
#' @param race 第三關鍵字
#' @param disease 第四關鍵字
#' @param odds 第五關鍵字
#' @return 輸出大表格，包含rsID對應疾病
#' @export


text_mining_pubmed_pdf <- function(path_, pmid_list, rsID, SNP, race, disease, odds){ #輸入單一個SNP
 final_merged<-data.frame()
 for(pmid_file in pmid_list){ #針對單個pmid，而一組rsID+Disease有多個pmids; 一個pmid有多個句子
      print(pmid_file )      
      target_pdf <- paste0(path_, pmid_file)
      test_pdf <- pdftools::pdf_text(pdf = target_pdf)
      PDF_ <- tokenize_sentences(test_pdf)
      for(o in seq_len(length(PDF_))){
        i = length(PDF_[[o]])
        for(x in seq_len(i)){
        PDF_[[o]][[x]]= str_squish(PDF_[[o]][[x]])}
      }
      PDF_ <- Filter(length, PDF_ ) # remove empty lists
      PDF_unlist <- unlist(PDF_)
      "先回傳Only的"
      rsID_boolean <- grepl(paste0(rsID,collapse="|"), PDF_unlist, ignore.case=TRUE)
      SNP_boolean <- grepl(paste0(SNP,collapse="|"), PDF_unlist, ignore.case=TRUE)
      race_boolean <- grepl(paste0(race,collapse="|"), PDF_unlist, ignore.case=TRUE)
      disease_boolean <- grepl(paste0(disease,collapse="|"), PDF_unlist, ignore.case=TRUE)
      odds_boolean <- grepl(paste0(odds,collapse="|"), PDF_unlist, ignore.case=FALSE)
      "Caucasians" %in% race
      disease_and_rsID <- PDF_unlist[disease_boolean & rsID_boolean] #key1
      odds_and_rsID <- PDF_unlist[odds_boolean & rsID_boolean]#key2
      disease_and_race <- PDF_unlist[disease_boolean & race_boolean] #key3
      race_and_SNP <- PDF_unlist[race_boolean & SNP_boolean] #key4
      race_and_rsID <- PDF_unlist[race_boolean & rsID_boolean] #key5
      odds_and_SNP <- PDF_unlist[odds_boolean & SNP_boolean] #key6
      element_ <- c("key_rsID","key_Disease","PMID", "Disease","Race", "rsID", "OR","SNP","sentence")
      index_ <- data.frame(matrix(ncol = 1, nrow = length(element_ )))
      rownames(index_) <- element_
      "key_1"
      if(length(disease_and_rsID)!=0){
        print("key1")
        key_1 <- data.frame("sentence" = disease_and_rsID)
        key_1$"Disease" <- return_value_for_each_sentance(disease_and_rsID, disease)
        #key_1$"rsID" <- paste(unique(unlist(sapply(disease_and_rsID, rsID_extraction))),collapse = "; ")
        key_1$"rsID" <- return_value_for_each_sentance(disease_and_rsID, rsID) 
        index_ <- merge(index_,t(key_1), by="row.names", all=TRUE)
        rownames(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_1<-"NA"}
      "key_2"
      if(length(odds_and_rsID)!=0){
        print("key2")
        key_2 <- data.frame("sentence" = odds_and_rsID)
        key_2$"OR" <- return_value_for_each_sentance(odds_and_rsID, odds)
        #key_2$"rsID" <- paste(unlist(rsID_extraction(odds_and_rsID)),collapse = "; ") #萃取所有rsID
        key_2$"rsID" <- return_value_for_each_sentance(odds_and_rsID, rsID)
        index_ <- merge(index_,t(key_2), by="row.names", all=TRUE)
        row.names(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_2<-"NA"}
      "key_3"
      if(length(disease_and_race)!=0){
        print("key3")
        key_3 <- data.frame("sentence" = disease_and_race)
        key_3$"Disease" <- return_value_for_each_sentance(disease_and_race, disease)
        key_3$"Race" <- return_value_for_each_sentance(disease_and_race, race)
        index_ <- merge(index_,t(key_3), by="row.names", all=TRUE)
        rownames(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_3<-"NA"}
      "key_4"
      if(length(race_and_SNP)!=0){
        print("key4")
        key_4 <- data.frame("sentence" = race_and_SNP)
        key_4$"Race" <- return_value_for_each_sentance(race_and_SNP, race)
        #key_4$"SNP" <- paste0(rsID_extraction(race_and_SNP),collapse=";")
        key_4$"SNP" <- "Yes"
        index_ <- merge(index_,t(key_4), by="row.names", all=TRUE)
        rownames(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_4<-"NA"}
      "key_5"
      if(length(race_and_rsID)!=0){
        print("key5")
        key_5 <- data.frame("sentence" = race_and_rsID)
        key_5$"Race" <- return_value_for_each_sentance(race_and_rsID, race)
        #key_5$"rsID" <- paste(unlist(rsID_extraction(race_and_rsID)),collapse = "; ") #萃取所有rsID
        key_5$"rsID" <- return_value_for_each_sentance(race_and_rsID, rsID) #只回傳搜尋的rsID
        index_ <- merge(index_,t(key_5), by="row.names", all=TRUE)
        rownames(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_5<-"NA"}
      "key_6"
      if(length(odds_and_SNP)!=0){
        print("key6")
        key_6 <- data.frame("sentence" = odds_and_SNP)
        key_6$"OR" <- return_value_for_each_sentance(odds_and_SNP, odds)
        #key_6$"SNP" <- paste0(rsID_extraction(odds_and_SNP),collapse = ";")
        key_6$"SNP" <- "Yes"
        index_ <- merge(index_,t(key_6), by="row.names", all=TRUE)
        rownames(index_) <- index_$"Row.names"
        index_<- index_[,-1]
        colnames(index_) <- paste0("x",seq_len(ncol(index_)))
      }else{key_6 <-"NA"}

      index_v2 <- as.data.frame(index_[,-1]) #如果只有一行，扣掉會取消data frame，所以先轉成v2再處理
      rownames(index_v2) <- rownames(index_) 
      final_index <- t(index_v2[c("key_rsID","key_Disease","PMID","Disease","Race", "rsID", "OR","SNP","sentence"),])
      colnames(final_index) <-  c("key_rsID","key_Disease","PMID","Disease","Race", "rsID", "OR","SNP","sentence")
      pmid <- gsub(".pdf","",pmid_file); final_index[,"PMID"] <- pmid
      rownames(final_index ) = NULL
      final_index[,"key_rsID"] <- rsID
      final_index[,"key_Disease"] <- disease
      #新增G2 SNP和disease到最前欄
      final_merged <- rbind(final_merged, final_index)

    }
  final_merged <- final_merged[-1,]
  return(final_merged)
  #write.csv(final_merged ,paste0(path_, disease_type, "_Summary_result_test.csv"), row.names=FALSE)
}

#' 針對已經建立好的rsID對應疾病的G2清單，共384組合進行表格輸出
#' 
#' @param numbers 匯入表格的row numbers
#' @param data 完整的表格，包含疾病和rsID對應組合
#' @param pmid_list 
#' @param file_names 輸出的檔案名稱
#' @return 輸出大表格
#' @export


Processing_text_mining_from_entire_table <- function(numbers, data = G2_rsID_Disease_PMIDs_list_clean, pmid_list, file_names){
  result_ <- NULL
  for(i in numbers){
    rsID <- data[i, "rsID"]
    disease <- data[i, "Disease"]
    race <- list_country_and_race 
    pmid_list <- sapply(strsplit(data[i, "PMIDs"], "; ", fixed = TRUE), function(x) paste0(x,".pdf"))
    pmid_list_existed <- pmid_list[is.element(pmid_list,list.files(path_G2_pdf ))] #比對真的有存在在資料夾的pdf檔案
    each_result <- text_mining_pubmed_pdf(path_G2_pdf, pmid_list_existed, rsID, SNP, race, disease, odds)
    result_ <- rbind(result_, each_result ) #輸入檔案位置，匯入pdf檔
    #pmid 58個
    #write.csv(each_result, file = paste0(path_G2_panel,file_names), append = TRUE)
    print(paste(i, disease, rsID ))
  }
 write.csv(result_ , file = paste0(path_G2_panel, file_names))
 return(result_ )
}



#' pubmed retrieve PMIDs 通用版本，可自定義多種關鍵字
#' 此版本從先前自己寫的SNPs_pubmed_for_multiple_factors啟發和延伸，目的是為變成通用版
#' 有時候get_pubmed_ids會出現warning: error in erl: 無法開啟...連結
#' 20230215 已新增 return retrieve PMIDs 
#' 20230215 已解決 "ERBB2 AND Lung cancer AND drug" >> Note that only 0 PubMed IDs were retrieved (1170 were expected).
#' @param keyword_col_list 檔案表格關鍵字的欄名，也是關鍵字的類別，例如: 生物標誌、癌症、藥物
#' @param data 完整的表格，包含不同類型關鍵字，以欄(column)分別呈現，例如: EGFR, Lung cancer, herceptin
#' @return 輸出不同關鍵字排列組合對應的PMID列表
#' @export
#' 
Pubmed_ID_retrieve <- function(keyword_col_list, data){
  cat("可能輸入錯誤的欄位: ", keyword_col_list[!is.element(keyword_col_list, colnames(data))])
  for(i in keyword_col_list){
    assign(i, data[, i] %>% stri_remove_empty(., na_empty = FALSE), envir=.GlobalEnv)
  } #將個別欄位重新assign給稱作欄名的變數，利用stri_remove_empty把空字串排除
  print("欄位參數化")
  params <- as.list(keyword_col_list) #將欄位名作為參數
  print("提取內容")
  var_list <- lapply(params, get) #將個別變數用get回傳內容物 
  print("組合關鍵字")
  keyword.df <- var_list %>% expand.grid()  #將多個欄位的關鍵字排列組合
  retrieve_df <- data.frame() #建立收集PMID的空集合
  print("準備進入迴圈")
  for(i in 1:nrow(keyword.df)){ #nrow(keyword.df)
    query <- paste0(unlist(keyword.df[i,]), collapse=" AND ")
    print(query)
    pubmed_id <- get_pubmed_ids(query)
    print("連結頁面")
    tryCatch({ if(pubmed_id$Count >= 1 ){
                    abstracts_xml <- fetch_pubmed_data(pubmed_id ) #如果沒結果，會卡很久，需要直接跳過不能再繼續執行
                    print("取得xml")
                    abstracts_list <- articles_to_list(abstracts_xml )
                    print("取得abstract_list")
                    retrieve_PMIDs <- paste(fetch_all_pubmed_ids(pubmed_id), collapse = "; ")
                    print("獲取PMIDs")
                    retrieve_df  <- rbind(retrieve_df, c(query, retrieve_PMIDs))
                }else{
                    print("No")
                    retrieve_df  <- rbind(retrieve_df, c(query,"No"))
                }
        },#some expression
        error=function(e) {
            message('An Error Occurred') #if an error occurs, tell me the error
            print(e)
        },
        warning=function(w) {
            message('A Warning Occurred') #if a warning occurs, tell me the warning
            print(w)
            return(NA)
        })

    cat("\r",round(i/nrow(keyword.df)*100,2), '%     ', query) #進度條
  }
 colnames(retrieve_df) <- c("Keywords", "PMIDs") 
 return(retrieve_df) 
}


#' R search PMIDs
#' @param host_ input number
#' @param port_ input port
#' @return run shiny app
#' @export

run_easypubmed_shiny <- function(host_, port_){
    library(easyPubMed)
    library(tidyverse)
    library(stringi)
    library(shiny)
    ui <- fluidPage(
      titlePanel("easyPubMed 搜尋器"),
      #----------------------------------------------------
      #session 1
      # 創建一個用於上傳CSV檔案的fileInput UI元素
      fileInput("file", "上傳CSV檔案"),
      # 用verbatimTextOutput顯示CSV檔案的欄位名稱
      verbatimTextOutput("columns"),
      #----------------------------------------------------
      #session 2
      # 輸入欄位名
      textInput(inputId = "selected_columns", label = "輸入欲查詢之關鍵字欄位名稱，以[;]或[,]符號分隔"),
      # 點選按鈕
      actionButton(inputId = "selection", label = "Search PMID"),
      # 接受server結果後回傳到頁面
      mainPanel(
          tableOutput("resultTable")
      ),
      #----------------------------------------------------
      #Session 3
      downloadButton("downloadData", "Download table"))

    # Define server logic
    server <- function(input, output) {
        #-----------------------------------------------------
        # 創建一個觀察器，用於處理文件上傳事件
        observeEvent(input$file, {
          # 讀取上傳的CSV檔案，並提取欄位名稱
          df <- read.csv(input$file$datapath)
          columns <- colnames(df)
          # 將欄位名稱傳回到UI中，以更新verbatimTextOutput UI元素的內容
          output$columns <- renderPrint(columns)
        
        })
        #-----------------------------------------------------

        result_table <- eventReactive(input$selection,{ 
            df <- read.csv(input$file$datapath)
            selected_columns <- input$selected_columns
            params <- unlist(strsplit(selected_columns, "[,!?; ]+"))
            PMIDs_table <- Pubmed_ID_retrieve(params, df)
            return(PMIDs_table)})     
        output$resultTable <- renderTable({result_table()})
        # Downloadable csv of selected dataset ----
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("test", ".csv", sep = "")
          },
          content = function(file){
            write.csv(result_table(), file, row.names = FALSE)})


    }
    # Run the application
    for_run <- shinyApp(ui = ui, server = server)
    runApp(for_run, host = host_, port = port_)
 }


#' R download PDF
#' @param host_ input number
#' @param port_ input port
#' @return r shiny app for paper download
#' @export

run_paper_download <- function(host_, port_ ){
    library(easyPubMed)
    library(tidyverse)
    library(stringi)
    library(shiny)
    ui <- fluidPage(
          textInput("pmid", "輸入 PMID"),
          actionButton(inputId = "download_pdf", "下載 PDF"),
          fileInput("file", "上傳關鍵字組合及PMIDs清單 (CSV檔案)"),
          actionButton(inputId = "download_pdf_2", "下載 PDF_2"),
          textInput("path", "欲下載至指定路徑")  
        )

    server <- function(input, output) {
      observeEvent(input$download_pdf ,{ 
          input_pmids <- as.data.frame(strsplit(input$pmid,"; "))
          colnames(input_pmids) <- "PMIDs" #一定要先給欄位名，以便後續可以從data.frame抽出變成series data type
          download_pubmed_paper(r_to_py(input_pmids)["PMIDs"], input$path, "") #一定要把pandas.data.frame轉成pandas.Series
      })

      observeEvent(input$download_pdf_2 ,{ 
        data_file <- read.csv(input$file$datapath)
        PMIDs_clean <- data_file[data_file$"PMIDs"!="No",] #排除掉No result
        for(i in 1:nrow(PMIDs_clean)){
          test_ <- as.data.frame(unlist(strsplit(PMIDs_clean[i,"PMIDs"],"; ")))
          colnames(test_ ) = "PMIDs" #一定要先給欄位名，以便後續可以從data.frame抽出變成series data type
          download_pubmed_paper(r_to_py(test_)["PMIDs"], input$path, "") #一定要把pandas.data.frame轉成pandas.Series
        }

      })
    }

    for_run <- shinyApp(ui = ui, server = server)

    runApp(for_run, host = host_, port = port_)

 }