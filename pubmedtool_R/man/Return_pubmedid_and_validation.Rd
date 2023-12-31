% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcipubmedtool.R
\name{Return_pubmedid_and_validation}
\alias{Return_pubmedid_and_validation}
\title{輸入已確認過的PMIDs檢查是否自動搜尋和手動搜尋結果一致}
\usage{
Return_pubmedid_and_validation(path, disease_type, race, SNP, PMIDs_manual)
}
\arguments{
\item{path}{想要輸出到哪個路徑資料夾儲存}

\item{disease_type}{只輸入一個疾病/第三關鍵字}

\item{race}{只輸入單一個人種/第二關鍵字}

\item{SNP}{只輸入單一個SNP rsID/ 第一關鍵字}

\item{PMIDs_manual}{將先前手動搜尋得到的PMIDs拿來驗證}
}
\value{
輸出檔案，其中包含PMID和驗證訊息
}
\description{
這個function會將ID另外output成一個檔，並且給出驗證訊息
}
