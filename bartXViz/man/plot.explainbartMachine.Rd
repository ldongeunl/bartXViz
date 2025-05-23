% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_explainbartMachine.R
\name{plot.explainbartMachine}
\alias{plot.explainbartMachine}
\title{plot.explainbartMachine}
\usage{
\method{plot}{explainbartMachine}(
  object,
  average = NULL,
  adjust = FALSE,
  type = NULL,
  num_post = NULL,
  plot.flag = TRUE,
  probs = 0.95,
  title = NULL,
  ...
)
}
\arguments{
\item{object}{An explainbartMachine class object containing the Shapley values of the BART model.}

\item{average}{Input the reference value for calculating the mean of the object's phi list.
"obs" represents the average based on observations (#post by #variable),
while "post" represents the average based on posterior samples (#obs by #variable).
If "both" is entered, calculations are performed based on both observation and posterior sample criteria.}

\item{adjust}{The default value is FALSE.
Enter TRUE to check the Shapley values adjusted based on the model's average contribution.}

\item{type}{"bar" represents a bar chart that includes the average contribution of each variable,
while "bee" represents a summary plot, allowing you to determine the graph's format.}

\item{num_post}{To check the contribution of variables for a single posterior sample,
enter a value within the number of posterior samples.}

\item{plot.flag}{If average = "obs", the quantile interval of each variable's is provided by default.}

\item{probs}{Enter the probability for the quantile interval. The default value is 0.95.}

\item{title}{The title of the plot, with a default value of NULL.}
}
\value{
Return plot:
If average is "obs" or "post", a bar plot or summary plot is generated based on
the selected averaging criterion. If average is "both", a boxplot is displayed to
show the distribution of Shapley values computed using both criteria.
If adjust is TRUE, the adjusted Shapley values are displayed.
If num_post is specified, a bar plot or summary plot for the selected posterior sample is generated.
}
\description{
This function is implemented to visualize the computed Shapley values in various ways.
The values and format used in the graph are determined based on the input parameters.
}
