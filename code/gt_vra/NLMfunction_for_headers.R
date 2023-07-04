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
    
    body <- list(
      `model` = "gpt-3.5-turbo",
      `messages` = list(
        list(`role` = "system", `content` = "You are a helpful assistant."),
        list(`role` = "user", `content` = paste("Describe the following data column header:", prompt))
      )
    )
    
    response <- POST(url, 
                     body = body, 
                     encode = "json", 
                     add_headers(c(Authorization = paste0("Bearer ", "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"),
                                   'Content-Type' = 'application/json')))
    
    response_content <- content(response, "parsed")
    assistant_message <- response_content$choices[[1]]$message$content
    
    return(assistant_message)
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
  
  ##Personal API Key###
  ##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
  