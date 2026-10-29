rmarkdown::render(
  "exploration/Diversity_Check_Absent_Species.Rmd",
  params = list(model = 1),
  output_file = "Model1.html"
)

rmarkdown::render(
  "exploration/Diversity_Check_Absent_Species.Rmd",
  params = list(model = 2),
  output_file = "Model2.html"
)

rmarkdown::render(
  "exploration/Diversity_Check_Absent_Species.Rmd",
  params = list(model = 3),
  output_file = "Model3.html"
)

rmarkdown::render(
  input = "exploration/Solve Model Convergence Issue.Rmd",
  output_file = "exploration/Solve-Model-Convergence-Issue_r64.html"
)
