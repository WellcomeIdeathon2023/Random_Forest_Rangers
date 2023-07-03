#Function to look at uploaded cvs files and five information with the data columns
  
  library(httr)
  library(jsonlite)
  
  read_headers <- function(file) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    headers <- colnames(data)
    return(headers)
  }
  
  openai_prompt <- function(prompt) {
    url <- "https://api.openai.com/v1/chat/completions"
    auth_header <- add_headers(`Authorization` = paste("Bearer", "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"),
                               `Content-Type` = "application/json")
    
    body <- toJSON(list(`model` = "gpt-3.5-turbo",
                        `prompt` = prompt,
                        `max_tokens` = 60))
    
    response <- POST(url, headers = auth_header, body = body, encode = "json")
    
    print(content(response, "text"))  # print the response
    
    result <- fromJSON(content(response, "text"))
    return(result$choices[[1]]$text)
  }
  
  interpret_headers <- function(file) {
    headers <- read_headers(file)
    descriptions <- sapply(headers, function(header) {
      description <- openai_prompt(paste("Describe the following data column header:", header))
      return(description)
    })
    names(descriptions) <- headers
    return(descriptions)
  }
  
  file_path <- "/Users/tutg/Documents/DataSciComp/Random_Forest_Rangers/data/sdy180/resultfiles/mbaa_result.csv"
  headers <- read_headers(file_path)
  print(headers)
  
  descriptions <- interpret_headers(file_path)
  print(descriptions)
  
  
  ##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
  