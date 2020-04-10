# building files in the right location
#Anthony Davidson
#Oct2019

# put the single output into the githubpages directory
#rmarkdown build
rmarkdown::render(input = "./index.Rmd", 
                  output_dir = "./docs", 
                  output_file = "test_output.docx")

# bookdown build
bookdown::render_book(input = c("index.Rmd", 
                                "00-status.Rmd", 
                                "00-proposal.Rmd", 
                                "03-methods.Rmd"))

rmarkdown::render(input = "./02-data_input.Rmd", 
                  output_dir = "./docs", 
                  output_file = "test_output.docx")
