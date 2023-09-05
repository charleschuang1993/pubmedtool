# 此為Charles Chuang 自定義模組

20230327

主要進行關鍵字組合來搜尋相關文獻的PMID，連接的學術網站為Pubmed。

搭配R shiny 可以直接部屬內網，記得填入host ip

需要事先library("shiny")

未來會再發布根據PMID下載文獻PDF之功能，預計6月發布

# 安裝模組請使用下列指令
# tcipubmedtool

```R
devtools::install_gitlab(auth_token = "ZVZBAxqHmRGxX8t_pzFf", host = "http://192.168.2.86:8088/", "charleschuang/tcipubmedtool")
```
